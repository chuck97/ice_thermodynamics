import numpy as np
import matplotlib.pyplot as plt

# == общие ==
mu = 0.054
sigma = 5.67e-8
T0 = 273.16
emissivity = 0.99
C_sh = 1e-3
C_lh = 1e-3
mmd_to_ms = 1e-3/3600/24 # перевод осадков из миллиметров в сутки в метры в секунду

# == воздух ==
rho_a = 1.28
c_pa = 1.01e3
P_surf = 1013.25

# == вода ==
c_w = 3.99e3
c_pw = 4.17e3
rho_w = 1023.0
c1_w = 17.27
c2_w = 35.86

# == лёд ==
c0_i = 2.06e3
rho_i = 917.0
L0_i = 3.34e5
kappa_i = 1.5
albedo_i = 0.6
i0_i = 0.15
c1_i = 21.87
c2_i = 7.66

# == снег ==
c0_s = 2.06e3
rho_s = 330.0
#k0_s = 0.31
k0_s = 0.5
L_f0 = 3.33e5
L_s0 = 2.83e6
kappa_s = 10
albedo_s = 0.8
i0_s = 0.08

# == функции для льда ==
Tf_i = lambda S_i: -mu*S_i
k_i = lambda T, S_i: rho_i/917.0 * (2.11 - 0.011*T + (0.09*S_i/T if T != 0 else 0))
#k_i = lambda T, S_i: (2.03 + 0.1172 * (S_i/T if T != 0.0 else 0.0))
c_i = lambda T, T_old, S_i: c0_i - L0_i*Tf_i(S_i)/(T*T_old)
E_i = lambda T, S_i: c0_i*(T - Tf_i(S_i)) - L0_i*(1.0 - Tf_i(S_i)/T) + c_w*Tf_i(S_i)
L_i = lambda T, S_i: c0_i*(T - Tf_i(S_i)) - L0_i*(1.0 - Tf_i(S_i)/T)

# == функции для снега ==
k_s = lambda T, S_s=0: k0_s
c_s = lambda T, T_old, S_s=0: c0_s
E_s = lambda T, S_s=0: c0_s*T - L_f0
L_s = lambda T, S_s=0: c0_s*T - L_f0

class Process:
    def __init__(self,
                 dzi_arr_init=np.array([]), dzs_arr_init=np.array([]),
                 timeline_init=np.array([]),
                 temp_oi_arr_init=np.array([]), temp_ice_arr_init=np.array([]), temp_is_arr_init=np.array([]),
                 temp_snow_arr_init=np.array([]), temp_sa_arr_init=np.array([]),
                 rho_ice_arr_init=np.array([]),
                 snow_filter_init=np.array([])):
        
        assert len(dzi_arr_init) == len(dzs_arr_init) == len(timeline_init) == len(temp_ice_arr_init) \
            == len(temp_ice_arr_init) == len(temp_is_arr_init) == len(temp_snow_arr_init) == len(temp_sa_arr_init) \
            == len(rho_ice_arr_init) == len(snow_filter_init), "Lengths of init arrays are not equal!"
        
        self.ice_dz_history = np.array(dzi_arr_init)
        self.snow_dz_history = np.array(dzs_arr_init)
        self.timeline = np.array(timeline_init)
        self.oi_temp_history = np.array(temp_oi_arr_init)
        self.ice_temp_history = np.array(temp_ice_arr_init)
        self.is_temp_history = np.array(temp_is_arr_init)
        self.snow_temp_history = np.array(temp_snow_arr_init)
        self.sa_temp_history = np.array(temp_sa_arr_init)
        self.ice_density_history = np.array(rho_ice_arr_init)
        self.snow_presence_history = np.array(snow_filter_init)
        
    def __getitem__(self, key):
        if isinstance(key, slice):
            return Process(self.ice_dz_history[key],
                           self.snow_dz_history[key],
                           self.timeline[key],
                           self.oi_temp_history[key],
                           self.ice_temp_history[key],
                           self.is_temp_history[key],
                           self.snow_temp_history[key],
                           self.sa_temp_history[key],
                           self.ice_density_history[key],
                           self.snow_presence_history[key])
        elif isinstance(key, int):
            return Process(self.ice_dz_history[[key]],
                           self.snow_dz_history[[key]],
                           self.timeline[[key]],
                           self.oi_temp_history[[key]],
                           self.ice_temp_history[[key]],
                           self.is_temp_history[[key]],
                           self.snow_temp_history[[key]],
                           self.sa_temp_history[[key]],
                           self.ice_density_history[[key]],
                           self.snow_presence_history[[key]])
    
    def get_zip(self):
        return zip(*self.__dict__.values()) 
# zip(self.ice_dz_history[clip_start:clip_end], self.snow_dz_history[clip_start:clip_end],
#                    self.timeline[clip_start:clip_end],
#                    self.oi_temp_history[clip_start:clip_end], self.ice_temp_history[clip_start:clip_end],
#                    self.is_temp_history[clip_start:clip_end],
#                    self.snow_temp_history[clip_start:clip_end], self.sa_temp_history[clip_start:clip_end],
#                    self.ice_density_history[clip_start:clip_end],
#                    self.snow_presence_history[clip_start:clip_end])
    
    def get_temp(self, mode):
        temp_arrs = [self.oi_temp_history, self.ice_temp_history, self.is_temp_history,
                     self.snow_temp_history, self.sa_temp_history]
        if mode == 'min':
            return min(np.nanmin(temp_arr) for temp_arr in temp_arrs if not np.isnan(temp_arr).all())
        elif mode == 'max':
            return max(np.nanmax(temp_arr) for temp_arr in temp_arrs if not np.isnan(temp_arr).all())
        else:
            raise Exception("'mode' should be either 'min' or 'max'!")
            
    def get_length(self):
        return self.timeline.size
    
    def get_nodes_num(self):
        return len(self.ice_dz_history), len(self.snow_dz_history)
            
def process_from_data(levels,
                      temp_array, Tib, Tis, Tss,
                      h_ib, h_is, h_ss,
                      dsigma_ice, dsigma_snow,
                      timeline, rho_ice,
                      snow_thickness_threshold=1e-3):
    
    assert len(temp_array) == len(Tss) == len(Tis) == len(Tib) == len(h_ib) == len(h_is) == len(h_ss) == len(timeline), \
           "Lenghts of input arrays ({}, {}, {}, {}, {}, {}, {}, {}) should be equal!".format(len(temp_array), len(Tss), len(Tis), len(Tib), len(h_ib), len(h_is), len(h_ss), len(timeline))
    
    
    mesh_Z = np.array([levels]*len(temp_array))
    has_snow = (abs(h_ss - h_is) >= snow_thickness_threshold)
    
    # == Formation of mesh and temperature arrays for snow and ice ==
    filter_ice = (h_ib.reshape(-1, 1) < mesh_Z) & (mesh_Z < h_is.reshape(-1, 1))
    filter_snow = (h_is.reshape(-1, 1) < mesh_Z) & (mesh_Z < h_ss.reshape(-1, 1))
    
    Z_ice = [np.concatenate(([surf], line_ice[filt_ice], [base])) \
             for surf, line_ice, filt_ice, base \
             in zip(h_is, mesh_Z, filter_ice, h_ib)]
    Z_snow = [np.concatenate(([surf], line_snow[filt_snow], [base])) \
              for surf, line_snow, filt_snow, base \
              in zip(h_ss, mesh_Z, filter_snow, h_is)]
    
    temp_ice = [np.concatenate(([surf], line_ice[filt_ice], [base])) \
                for surf, line_ice, filt_ice, base \
                in zip(Tis, temp_array, filter_ice, Tib)]
    temp_snow = [np.concatenate(([surf], line_ice[filt_ice], [base])) \
                 for surf, line_ice, filt_ice, base \
                 in zip(Tss, temp_array, filter_snow, Tis)]
    
    # == Interpolating internal nodes into the new mesh ==
    sigma_ice_nodes = np.concatenate(([0.0], dsigma_ice.cumsum()))
    sigma_ice_centers = (sigma_ice_nodes[:-1] + sigma_ice_nodes[1:])/2
    Z_points_ice = [Z_ice_line[-1] + sigma_ice_centers*(Z_ice_line[0] - Z_ice_line[-1]) for Z_ice_line in Z_ice]
    T_points_ice = [np.interp(Z_points_ice_line, Z_ice_line[::-1], temp_ice_line[::-1]) \
                    for Z_points_ice_line, Z_ice_line, temp_ice_line \
                    in zip(Z_points_ice, Z_ice, temp_ice)]
    
    sigma_snow_nodes = np.concatenate(([0.0], dsigma_snow.cumsum()))
    sigma_snow_centers = (sigma_snow_nodes[:-1] + sigma_snow_nodes[1:])/2
    Z_points_snow = [Z_snow_line[-1] + sigma_snow_centers*(Z_snow_line[0] - Z_snow_line[-1]) \
                     for Z_snow_line, snow in zip(Z_snow, has_snow)]
    T_points_snow = [np.interp(Z_points_snow_line, Z_snow_line[::-1], temp_snow_line[::-1]) if snow \
                     else np.array([np.nan]*len(dsigma_snow)) \
                     for Z_points_snow_line, Z_snow_line, temp_snow_line, snow \
                     in zip(Z_points_snow, Z_snow, temp_snow, has_snow)]
    
    print(np.array(dsigma_ice*(h_is - h_ib).reshape(-1, 1)))
    print(np.array(dsigma_snow*(h_ss - h_is).reshape(-1, 1)))
    
    return Process(dzi_arr_init=dsigma_ice*(h_is - h_ib).reshape(-1, 1),
                   dzs_arr_init=dsigma_snow*(h_ss - h_is).reshape(-1, 1),
                   timeline_init=timeline,
                   temp_oi_arr_init=Tib,
                   temp_ice_arr_init=T_points_ice,
                   temp_is_arr_init=Tis,
                   temp_snow_arr_init=T_points_snow,
                   temp_sa_arr_init=[temp if snow else np.nan for temp, snow in zip(Tss, has_snow)],
                   rho_ice_arr_init=np.ones((len(timeline), len(dsigma_ice)))*rho_ice,
                   snow_filter_init=has_snow)

def secant_solver(f, x0, x1, tol=1e-9, max_it=1000):
    x_new = x1
    x_old = x0
    it = 0
    while abs(x_new - x_old) > tol and it < max_it:
        it += 1
        temp = x_new
        x_new = x_new - f(x_new)*(x_new-x_old)/(f(x_new) - f(x_old))
        x_old = temp
    
    if it == max_it:
        X = np.linspace(x0, x1, 100)
        plt.figure(figsize=(15, 10))
        plt.plot(X, [f(x) for x in X])
        plt.grid()
        plt.show()
        raise Exception("secant solver won't converge!")
        
    return x_new, it


def thomas_solver(T_under, T_diag, T_over, rhs):
    assert len(T_under) == len(T_diag) - 1
    assert len(T_over) == len(T_diag) - 1
    assert len(T_diag) == len(rhs)
    
    sol = np.zeros(len(T_diag));
    c_new = np.zeros(len(T_diag) - 1);
    d_new = np.zeros(len(T_diag));

    # forward iteration
    c_new[0] = T_over[0]/T_diag[0];
    for i in range(1, len(T_diag) - 1):
        c_new[i] = T_over[i]/(T_diag[i] - T_under[i-1]*c_new[i-1])

    # backward iteration
    d_new[0] = rhs[0]/T_diag[0]
    for i in range(1, len(T_diag)):
        d_new[i] = (rhs[i] - T_under[i-1]*d_new[i-1])/(T_diag[i] - T_under[i-1]*c_new[i-1])

    # final solution
    sol[len(T_diag) - 1] = d_new[len(T_diag)-1]
    for i in range(len(T_diag) - 2, -1, -1): 
        sol[i] = d_new[i] - c_new[i]*sol[i+1]

    return sol


def Update_dz(dz_cells_old, omega_ocn, omega_atm, time_step):
    dz_old = dz_cells_old.sum()
    return dz_cells_old - time_step*(omega_atm - omega_ocn)*dz_cells_old/dz_old


def compute_radiation_nodes(dzi, F_sw, dzs=None):
    dzi_extended = np.concatenate(([0], dzi[::-1]))
    if dzs is None:
        radiation_nodes_ice = F_sw*(1-albedo_i)*i0_i*np.exp(-kappa_i*np.cumsum(dzi_extended))
        return radiation_nodes_ice[::-1]
    else:
        dzs_extended = np.concatenate(([0], dzs[::-1]))
        radiation_nodes_snow = F_sw*(1-albedo_s)*i0_s*np.exp(-kappa_s*np.cumsum(dzs_extended))
        radiation_nodes_ice = F_sw*(1-albedo_s)*i0_s*np.exp(-kappa_i*np.cumsum(dzi_extended)\
                                                            - sum(kappa_s*dzs))
        return radiation_nodes_snow[::-1], radiation_nodes_ice[::-1]


def k_eff_node(k_low, k_high, dz_low, dz_high):
    return 2*k_low*k_high/(k_low*dz_high + k_high*dz_low)


def get_matrix_upwind(T_cells_prev, T_cells_old,
                      T_ice_ocn_new, T_ice_ocn_prev, T_ice_ocn_old,
                      T_ice_atm_new, T_ice_atm_prev, T_ice_atm_old,
                      omega_ice_ocn, omega_ice_atm,
                      dz_cells_new, dz_cells_old,
                      salinity_cells,
                      radiation_nodes,
                      E, c, k, rho,
                      time_step, is_snow=False):
    
    N = len(T_cells_old)
    
    # проверка идентичности размеров массивов температур
    assert len(T_cells_prev) == N
    
    # проверка размерности массивов толщин ячеек
    assert len(dz_cells_new) == N
    assert len(dz_cells_old) == N
    
    # проверка размерности массива радиации
    assert len(radiation_nodes) == (N + 1)
    
    # расчет коэффициентов интерполяции (в каждом узле будет 2 коэффициента,
    # в первом и последнем узле коэффициенты нулевые) 
    a = np.zeros(N+1)
    b = np.zeros(N+1)
    
    for i in range(1, N):
        a[i] = dz_cells_new[i]/(dz_cells_new[i-1] + dz_cells_new[i])
        b[i] = dz_cells_new[i-1]/(dz_cells_new[i-1] + dz_cells_new[i])
    
    
    # расчет приближений коэффициента теплопроводности, теплоемкости и температуры в узлы
    k_nodes = np.zeros(N+1)
    
    if not is_snow:
    
        for i in range(1, N):
            k_nodes[i] = a[i]*k(T_cells_prev[i-1], salinity_cells[i-1]) + \
                         b[i]*k(T_cells_prev[i], salinity_cells[i])

        k_nodes[0] = k(T_cells_prev[0], salinity_cells[0])
        k_nodes[-1] = k(T_cells_prev[-1], salinity_cells[-1])

        c_ice_ocn = c(T_ice_ocn_prev, T_ice_ocn_old, salinity_cells[0])
        c_ice_atm = c(T_ice_atm_prev, T_ice_atm_old, salinity_cells[-1])

        E_ice_ocn = E(T_ice_ocn_old, salinity_cells[0])
        E_ice_atm = E(T_ice_atm_old, salinity_cells[-1])
        
    else:
        
        for i in range(1, N):
            k_nodes[i] = a[i]*k(T_cells_prev[i-1]) + \
                         b[i]*k(T_cells_prev[i])

        k_nodes[0] = k(T_cells_prev[0])
        k_nodes[-1] = k(T_cells_prev[-1])

        c_ice_ocn = c(T_ice_ocn_prev, T_ice_ocn_old)
        c_ice_atm = c(T_ice_atm_prev, T_ice_atm_old)

        E_ice_ocn = E(T_ice_ocn_old)
        E_ice_atm = E(T_ice_atm_old)
    
    # расчет узловых значений omega, теплоемкости и энтальпии
    omega_nodes = np.zeros(N+1)
    
    for i in range(0, N+1): 
        omega_nodes[i] = omega_ice_ocn\
                       + (sum(dz_cells_new[:i])/sum(dz_cells_new))*(omega_ice_atm - omega_ice_ocn)
    
    # расчет значений энтальпии, теплоемкости в ячейках
    E_cells = np.zeros(N)
    c_cells = np.zeros(N)
    
    for i in range(0, N):
        
        if not is_snow: 
            E_cells[i] = E(T_cells_old[i], salinity_cells[i])
            c_cells[i] = c(T_cells_prev[i], T_cells_old[i], salinity_cells[i])
        else:
            E_cells[i] = E(T_cells_old[i])
            c_cells[i] = c(T_cells_prev[i], T_cells_old[i])
        
    ### сборка диагоналей матрицы и вектора правой части
    A = np.zeros(N)
    B = np.zeros(N)
    C = np.zeros(N)
    RHS = np.zeros(N)
    
    # первая строка 
    B[0] = rho*c_cells[0]*dz_cells_new[0]/time_step + \
    (2.0*k_nodes[1])/(dz_cells_new[0] + dz_cells_new[1]) + (2.0*k_nodes[0])/(dz_cells_new[0])
    
    C[0] = -(2.0*k_nodes[1])/(dz_cells_new[0] + dz_cells_new[1])
    
    RHS[0] = rho*c_cells[0]*dz_cells_new[0]*T_cells_old[0]/time_step - \
    rho*E_cells[0]*(dz_cells_new[0] - dz_cells_old[0])/time_step + \
    (radiation_nodes[1] - radiation_nodes[0]) + \
    (2.0*k_nodes[0]*T_ice_ocn_new)/(dz_cells_new[0])
    
    if (omega_nodes[0] >= 0):
        RHS[0] += rho*c_ice_ocn*T_ice_ocn_new*omega_nodes[0] - \
        rho*(c_ice_ocn*T_ice_ocn_old - E_ice_ocn)*omega_nodes[0]
    else:
        B[0] += -rho*c_cells[0]*omega_nodes[0]
        RHS[0] += -rho*(c_cells[0]*T_cells_old[0] - E_cells[0])*omega_nodes[0]
            
    if (omega_nodes[1] >= 0):
        B[0] += rho*c_cells[0]*omega_nodes[1]
        RHS[0] += rho*(c_cells[0]*T_cells_old[0] - E_cells[0])*omega_nodes[1]
    else:
        C[0] += rho*c_cells[1]*omega_nodes[1]
        RHS[0] += rho*(c_cells[1]*T_cells_old[1] - E_cells[1])*omega_nodes[1]
        
    # серединные строки
    for i in range(1, N-1):
        
        A[i] = -(2.0*k_nodes[i])/(dz_cells_new[i-1] + dz_cells_new[i])
        
        B[i] = rho*c_cells[i]*dz_cells_new[i]/time_step + \
        (2.0*k_nodes[i+1])/(dz_cells_new[i] + dz_cells_new[i+1]) + \
        (2.0*k_nodes[i])/(dz_cells_new[i-1] + dz_cells_new[i]) 
        
        C[i] = -(2.0*k_nodes[i+1])/(dz_cells_new[i] + dz_cells_new[i+1])
        
        RHS[i] = rho*c_cells[i]*dz_cells_new[i]*T_cells_old[i]/time_step - \
        rho*E_cells[i]*(dz_cells_new[i] - dz_cells_old[i])/time_step + \
        (radiation_nodes[i+1] - radiation_nodes[i])
        
        if (omega_nodes[i] >= 0):
            A[i] += -rho*c_cells[i-1]*omega_nodes[i]
            RHS[i] += -rho*(c_cells[i-1]*T_cells_old[i-1] - E_cells[i-1])*omega_nodes[i]
        else:
            B[i] += -rho*c_cells[i]*omega_nodes[i]
            RHS[i] += -rho*(c_cells[i]*T_cells_old[i] - E_cells[i])*omega_nodes[i]
            
        if (omega_nodes[i+1] >= 0):
            B[i] += rho*c_cells[i]*omega_nodes[i+1]
            RHS[i] += rho*(c_cells[i]*T_cells_old[i] - E_cells[i])*omega_nodes[i+1]
        else:
            C[i] += rho*c_cells[i+1]*omega_nodes[i+1]
            RHS[i] += rho*(c_cells[i+1]*T_cells_old[i+1] - E_cells[i+1])*omega_nodes[i+1]
        
    # последняя строка
    A[N-1] = -(2.0*k_nodes[N-1])/(dz_cells_new[N-2] + dz_cells_new[N-1])

    B[N-1] = rho*c_cells[N-1]*dz_cells_new[N-1]/time_step + \
    (2.0*k_nodes[N])/(dz_cells_new[N-1]) + (2.0*k_nodes[N-1])/(dz_cells_new[N-2] + dz_cells_new[N-1])
    
    RHS[N-1] = rho*c_cells[N-1]*T_cells_old[N-1]*dz_cells_new[N-1]/time_step - \
    rho*E_cells[N-1]*(dz_cells_new[N-1] - dz_cells_old[N-1])/time_step - \
    (radiation_nodes[N] - radiation_nodes[N-1]) + \
    (2.0*k_nodes[N]*T_ice_atm_new)/(dz_cells_new[N-1])
    
    if (omega_nodes[N-1] >= 0):
        A[N-1] += -rho*c_cells[N-2]*omega_nodes[N-1]
        RHS[N-1] += -rho*(c_cells[N-2]*T_cells_old[N-2] - E_cells[N-2])*omega_nodes[N-1]
    else:
        B[N-1] += -rho*c_cells[N-1]*omega_nodes[N-1]
        RHS[N-1] += -rho*(c_cells[N-1]*T_cells_old[N-1] - E_cells[N-1])*omega_nodes[N-1]
            
    if (omega_nodes[N] >= 0):
        B[N-1] += rho*c_cells[N-1]*omega_nodes[N]
        RHS[N-1] += rho*(c_cells[N-1]*T_cells_old[N-1] - E_cells[N-1])*omega_nodes[N]
    else:
        RHS[N-1] += -rho*c_ice_atm*T_ice_atm_new*omega_nodes[N] + \
        rho*(c_ice_atm*T_ice_atm_old - E_ice_atm)*omega_nodes[N]
        
    return A[1:], B, C[:-1], RHS


def concat_matrices(A_i, B_i, C_i, RHS_i,
                    A_s, B_s, C_s, RHS_s,
                    dzi, dzs,
                    k_i, k_s):
    
    # приближение второго порядка: T_i'z = alpha*T_is + beta*T_i-1/2 + gamma*T_i-3/2
    
    h1_i = dzi[-1]/2.0
    h2_i = dzi[-2]/2.0 + dzi[-1]
    
    h1_s = dzs[0]/2.0
    h2_s = dzi[1]/2.0 + dzi[0]
    
    alpha_i = (h1_i + h2_i)/(h1_i*h2_i)
    beta_i = -h2_i/(h1_i*(h2_i - h1_i))
    gamma_i = h1_i/(h2_i*(h2_i - h1_i))
    
    alpha_s = -(h1_s + h2_s)/(h1_s*h2_s)
    beta_s = h2_s/(h1_s*(h2_s - h1_s))
    gamma_s = -h1_s/(h2_s*(h2_s - h1_s))
    
    A_is = k_i*gamma_i
    B_is = k_i*beta_i
    C_is = k_i*alpha_i - k_s*alpha_s
    D_is = -k_s*beta_s
    E_is = -k_s*gamma_s
    RHS_is = 0
    
    # производим вычитание из промежуточной строки так, чтобы крайние элементы занулились
    
    # 1
    B_is -= A_is*B_i[-1]/A_i[-1]
    RHS_is -= A_is*RHS_i[-1]/A_i[-1]
    A_is = 0
    
    # 2
    D_is -= E_is*B_s[0]/C_s[0]
    RHS_is -= E_is*RHS_s[0]/C_s[0]
    E_is = 0
    
    # получаем трёхдиагональную матрицу, производим конкатенацию
    
    A_full = np.concatenate((A_i, [B_is, 0], A_s))
    B_full = np.concatenate((B_i, [C_is], B_s))
    C_full = np.concatenate((C_i, [0, D_is], C_s))
    RHS_full = np.concatenate((RHS_i, [RHS_is], RHS_s))
    
    return A_full, B_full, C_full, RHS_full


def W_from_BC(T_node, T_cells,
              dz_cells,
              salinity_cells,
              k, F, 
              L,
              is_surface=False, rho=rho_i, time_step=None):
    
    if is_surface:
    
        h1 = dz_cells[-1]/2.0
        h2 = dz_cells[-1] + dz_cells[-2]/2.0
        # возможно, тут стоит учесть солёность
        grad = lambda T: (T*(h2**2 - h1**2) - T_cells[-1]*h2**2 + T_cells[-2]*h1**2)/(h1*h2*(h2-h1)) 
        residue = (k*grad(T_node) - F)*time_step
        omega = 0
        melted_layers = 0

        for T, dz, S in zip(T_cells[::-1], dz_cells[::-1], salinity_cells[::-1]):
            # print('melted:', melted_layers, 'residue:', residue, 'layer enthalpy:', rho*L(T)*dz)
            if residue <= rho*L(T, S)*dz:
                residue -= rho*L(T, S)*dz
                omega += dz/time_step
                melted_layers += 1
            else:
                break

        omega += residue/(rho*L(T, S)*time_step)
        return omega
        
    else:
        
        #h1 = dz_cells[0]/2.0
        #h2 = dz_cells[0] + dz_cells[1]/2.0
        # grad = lambda T: (-T*(h2**2 - h1**2) + T_cells[0]*h2**2 - T_cells[1]*h1**2)/(h1*h2*(h2-h1))
        
        level = dz_cells[0]/2
        for i, (dz_prev, dz_curr) in enumerate(zip(dz_cells[:-1], dz_cells[1:]), start=1):
            level += (dz_prev + dz_curr)/2
            if level >= 0.4:
                break
            
        grad = (T_cells[i] - T_node)/level
        #omega = (k*grad + F)/(rho_i*L(T_cells[0], salinity_cells[0]))
        
        omega = 0
        if (k*grad + F) > 0:
            # basal melt
            omega = (k*grad + F)/(rho_i*L(T_cells[0], 4.0))
        else:
            # basal growth
            omega = (k*grad + F)/(rho_i*L(T_cells[0], 10.0))
            
        return omega, k*grad
    
    
def T_from_BC(T_cells, dz_cells, salinity_cells,
              k, omega,
              F, L,
              is_surface=False, rho=rho_i):
    
    # аппроксимация градиента на границе
    if is_surface:
        h1 = dz_cells[-1]/2.0
        h2 = dz_cells[-1] + dz_cells[-2]/2.0
        # возможно, тут стоит учесть солёность
        grad_su = lambda T: (T*(h2**2 - h1**2) - T_cells[-1]*h2**2 + T_cells[-2]*h1**2)/(h1*h2*(h2-h1))
        nonlin_func = lambda T: rho*L(T, salinity_cells[-1])*omega + F(T) - k*(grad_su(T))
        return secant_solver(nonlin_func, T_cells[-1]-10.0, Tf_i(salinity_cells[-1])+10.0)[0]
    
    else:
        h1 = dz_cells[0]/2.0
        h2 = dz_cells[0] + dz_cells[1]/2.0
        grad_su = lambda T: (-T*(h2**2 - h1**2) + T_cells[0]*h2**2 - T_cells[1]*h1**2)/(h1*h2*(h2-h1))
        nonlin_func = lambda T: rho_i*L(T, salinity_cells[0])*omega + F(T) - k*(grad_su(T))
        return secant_solver(nonlin_func, T_cells[0]-10.0, Tf_i(salinity_cells[0])+10.0)[0] 
    
    
def ice_freezing(Toi, Ti, Tia, F_atm, F_ocn, F_sw,
                 dzi, salinity,
                 N_pseudoiter,
                 time_step,
                 tol=1e-6,
                 err_func=lambda T_old, T_prev, T_new: \
                 np.linalg.norm(T_new - T_prev)/np.linalg.norm(T_old)
                ):
    
    Ti_old = Ti
    Ti_new = Ti
    Tia_old = Tia
    Tia_new = Tia
    
    dzi_old = dzi
    dzi_new = dzi
    
    salinity_cells = salinity
    
    Tia_values = [Tia_old]
    current_err = np.inf
    surface_err = np.inf
    prev_surface_err = surface_err
    
    for pseudoiter in range(N_pseudoiter):
        
        Ti_prev = Ti_new
        Tia_prev = Tia_new
        
        prev_surface_err = surface_err
        
        omega_io, cond = W_from_BC(Toi, Ti_prev,
                             dzi_new,
                             salinity_cells,
                             k_i(Toi, salinity_cells[0]),
                             F_ocn(Toi), 
                             L_i,
                             is_surface=False
                            )

        # оценка Tia
        Tia_new = T_from_BC(Ti_prev, dzi_new,
                            salinity_cells,
                            k_i(Tia_prev, salinity_cells[-1]), 0.0,
                            F_atm, L_i,
                            is_surface=True
                           )
        
        # пересчет толщин слоев
        dzi_new = Update_dz(dzi_old, omega_io, 0.0, time_step)
        
        # пересчёт радиации
        radiation_nodes_ice = compute_radiation_nodes(dzi_new, F_sw)

        # перессчет температур в ячейках
        A, B, C, RHS = get_matrix_upwind(Ti_prev, Ti_old,
                                         Toi, Toi, Toi,
                                         Tia_new, Tia_prev, Tia_old,
                                         omega_io, 0.0,
                                         dzi_new, dzi_old,
                                         salinity_cells,
                                         radiation_nodes_ice,
                                         E_i, c_i, k_i, rho_i,
                                         time_step)
        
        Ti_new = thomas_solver(A, B, C, RHS)

        # оценка ошибки
        prev_surface_err = surface_err
        surface_err = abs(Tia_new - Tia_prev)/abs(Tia_old)
        
        if (surface_err < prev_surface_err):
            Tia_values.append(Tia_new)
        else:
            Tia_new = sum(Tia_values)/len(Tia_values)
        
        T_full_old = np.append(Ti_old, Tia_old)
        T_full_prev = np.append(Ti_prev, Tia_prev)
        T_full_new = np.append(Ti_new, Tia_new)
        err = err_func(T_full_old, T_full_prev, T_full_new)
        curr_err = err
        
        if current_err < tol:
            break
    
    return Ti_new, Tia_new, dzi_new, cond


def ice_melting(Toi, Ti, Tia_old, F_atm, F_ocn, F_sw,
                dzi, salinity,
                N_pseudoiter,
                time_step,
                tol=1e-6,
                err_func=lambda T_old, T_prev, T_new: \
                np.linalg.norm(T_new - T_prev)/np.linalg.norm(T_old)):
    
    # инициализация текущих T
    Ti_old = Ti
    Ti_new = Ti
    
    dzi_old = dzi
    dzi_new = dzi
    
    salinity_cells = salinity
    
    current_err = np.inf
    
    for pseudoiter in range(N_pseudoiter):
    
        # инициализация текущих T
        Ti_prev = Ti_new

        # оценка omega на границе лед-океан
        omega_io, cond = W_from_BC(Toi, Ti_prev,
                             dzi_new,
                             salinity_cells,
                             k_i(Toi, salinity_cells[0]),
                             F_ocn(Toi), 
                             L_i,
                             is_surface=False)
        
        # оценка omega на границе лед-атмосфера
        omega_ia = W_from_BC(Tf_i(salinity_cells[-1]), Ti_prev,
                             dzi_new,
                             salinity_cells,
                             k_i(Tf_i(salinity_cells[-1]), salinity_cells[-1]),
                             F_atm(Tf_i(salinity_cells[-1])),
                             L_i,
                             is_surface=True,
                             time_step=time_step)
        
        # пересчет толщин слоев
        dzi_new = Update_dz(dzi_old, omega_io, omega_ia, time_step)
        
        # пересчёт радиации
        radiation_nodes_ice = compute_radiation_nodes(dzi_new, F_sw)

        # перессчет температур в ячейках
        A, B, C, RHS = get_matrix_upwind(Ti_prev, Ti_old,
                                         Toi, Toi, Toi,
                                         Tf_i(salinity_cells[-1]), Tf_i(salinity_cells[-1]), Tia_old,
                                         omega_io, omega_ia,
                                         dzi_new, dzi_old,
                                         salinity_cells,
                                         radiation_nodes_ice,
                                         E_i, c_i, k_i, rho_i,
                                         time_step)
        
        Ti_new = thomas_solver(A, B, C, RHS)

        # оценка ошибки
        err = err_func(Ti_old, Ti_prev, Ti_new)
        curr_err = err
        
        if current_err < tol:
            break
            
    return Ti_new, dzi_new, cond


def snow_ice_freezing(Toi, Ti, Ts, Tis, Tsa, Ta, F_atm, F_ocn, F_sw,
                      dzi, dzs, salinity,
                      N_pseudoiter,
                      time_step,
                      p,
                      tol=1e-6,
                      err_func=lambda T_old, T_prev, T_new: \
                      np.linalg.norm(T_new - T_prev)/np.linalg.norm(T_old)
                     ):
    
    Ti_old = Ti
    Ti_new = Ti
    Ts_old = Ts
    Ts_new = Ts
    Tis_old = Tis
    Tis_new = Tis
    Tsa_old = Tsa
    Tsa_new = Tsa
    
    dzi_old = dzi
    dzi_new = dzi
    dzs_old = dzs
    dzs_new = dzs
    
    salinity_cells = salinity
    
    Tsa_values = [Tsa_old]
    current_err = np.inf
    surface_err = np.inf
    prev_surface_err = surface_err
    
    omega_sa = -p*rho_w/rho_s if (Ta < 0.0) else 0.0
    
    for pseudoiter in range(N_pseudoiter):
        
        Ti_prev = Ti_new
        Ts_prev = Ts_new
        Tis_prev = Tis_new
        Tsa_prev = Tsa_new
        
        prev_surface_err = surface_err
        
        omega_io, cond = W_from_BC(Toi, Ti_prev,
                             dzi_new,
                             salinity_cells,
                             k_i(Toi, salinity_cells[0]),
                             F_ocn(Toi), 
                             L_i,
                             is_surface=False
                            )

        # оценка Tsa
        Tsa_new = T_from_BC(Ts_prev, dzs_new,
                            salinity_cells,
                            k_s(Tsa_prev, salinity_cells[-1]), omega_sa,
                            F_atm, L_s,
                            is_surface=True,
                            rho=rho_s
                           )
        
        # пересчет толщин слоев
        dzi_new = Update_dz(dzi_old, omega_io, 0.0, time_step)
        dzs_new = Update_dz(dzs_old, 0.0, omega_sa, time_step)
            
        # пересчёт радиации
        radiation_nodes_snow, radiation_nodes_ice = compute_radiation_nodes(dzi_new, F_sw, dzs_new)

        # перессчет температур в ячейках
        A_i, B_i, C_i, RHS_i = get_matrix_upwind(Ti_prev, Ti_old,
                                                 Toi, Toi, Toi,
                                                 Tis_new, Tis_prev, Tis_old,
                                                 omega_io, 0.0,
                                                 dzi_new, dzi_old,
                                                 salinity_cells,
                                                 radiation_nodes_ice,
                                                 E_i, c_i, k_i, rho_i,
                                                 time_step)
        
        A_s, B_s, C_s, RHS_s = get_matrix_upwind(Ts_prev, Ts_old,
                                                 Tis_new, Tis_prev, Tis_old,
                                                 Tsa_new, Tsa_prev, Tsa_old,
                                                 0.0, omega_sa,
                                                 dzs_new, dzs_old,
                                                 salinity_cells,
                                                 radiation_nodes_snow,
                                                 E_s, c_s, k_s, rho_s,
                                                 time_step, is_snow=True)
        
        A, B, C, RHS = concat_matrices(A_i, B_i, C_i, RHS_i,
                                       A_s, B_s, C_s, RHS_s,
                                       dzi_new, dzs_new,
                                       k_i(Tis_new, salinity_cells[-1]), k_s(Tis_new))
        
        T_vec_new = thomas_solver(A, B, C, RHS)
        
        Ti_new, Tis_new, Ts_new = T_vec_new[:len(Ti_old)],\
                                  T_vec_new[len(Ti_old)],\
                                  T_vec_new[len(Ti_old)+1:]

        # оценка ошибки
        prev_surface_err = surface_err
        surface_err = abs(Tsa_new - Tsa_prev)/(abs(Tsa_old) if Tsa_old != 0 else 1e-6)
        
        if (surface_err < prev_surface_err):
            Tsa_values.append(Tsa_new)
        else:
            Tsa_new = sum(Tsa_values)/len(Tsa_values)
        
        T_full_old = np.concatenate((Ti_old, [Tis_old], Ts_old, [Tsa_old]))
        T_full_prev = np.concatenate((Ti_prev, [Tis_prev], Ts_prev, [Tsa_prev]))
        T_full_new = np.concatenate((Ti_new, [Tis_new], Ts_new, [Tsa_new]))
        err = err_func(T_full_old, T_full_prev, T_full_new)
        curr_err = err
        
        if current_err < tol:
            break
            
    return Ti_new, Ts_new, Tis_new, Tsa_new, dzi_new, dzs_new, cond


def snow_melting(Toi, Ti, Ts, Tis, Tsa_old, Ta, F_atm, F_ocn, F_sw,
                 dzi, dzs, salinity,
                 N_pseudoiter,
                 time_step,
                 p,
                 tol=1e-6,
                 err_func=lambda T_old, T_prev, T_new: \
                 np.linalg.norm(T_new - T_prev)/np.linalg.norm(T_old)):
    
    # инициализация текущих T
    Ti_old = Ti
    Ti_new = Ti
    Ts_old = Ts
    Ts_new = Ts
    Tis_old = Tis
    Tis_new = Tis
    
    dzi_old = dzi
    dzi_new = dzi
    dzs_old = dzs
    dzs_new = dzs
    
    salinity_cells = salinity
    
    T_vec_old = np.concatenate((Ti_old, [Tis_old], Ts_old))
    T_vec_new = np.concatenate((Ti_new, [Tis_new], Ts_new))
    current_err = np.inf
    
    for pseudoiter in range(N_pseudoiter):
    
        # инициализация текущих T
        Ti_prev = Ti_new
        Ts_prev = Ts_new
        Tis_prev = Tis_new
        T_vec_prev = T_vec_new

        # оценка omega на границе лед-океан
        omega_io, cond = W_from_BC(Toi, Ti_prev,
                             dzi_new,
                             salinity_cells,
                             k_i(Toi, salinity_cells[0]),
                             F_ocn(Toi), 
                             L_i,
                             is_surface=False)
        
        # оценка omega на границе лед-атмосфера
        omega_sa = W_from_BC(0.0, Ts_prev,
                             dzs_new,
                             salinity_cells,
                             k_s(0.0, salinity_cells[-1]),
                             F_atm(0.0),
                             L_s,
                             is_surface=True,
                             rho=rho_s,
                             time_step=time_step)
        
        if Ta < 0:
            omega_sa -= p*rho_w/rho_s
        
        # пересчет толщин слоев
        dzi_new = Update_dz(dzi_old, omega_io, 0.0, time_step)
        dzs_new = Update_dz(dzs_old, 0.0, omega_sa, time_step)
        
        # пересчёт радиации
        radiation_nodes_snow, radiation_nodes_ice = compute_radiation_nodes(dzi_new, F_sw, dzs_new)

        # перессчет температур в ячейках
        A_i, B_i, C_i, RHS_i = get_matrix_upwind(Ti_prev, Ti_old,
                                                 Toi, Toi, Toi,
                                                 Tis_new, Tis_prev, Tis_old,
                                                 omega_io, 0.0,
                                                 dzi_new, dzi_old,
                                                 salinity_cells,
                                                 radiation_nodes_ice,
                                                 E_i, c_i, k_i, rho_i,
                                                 time_step)
        
        A_s, B_s, C_s, RHS_s = get_matrix_upwind(Ts_prev, Ts_old,
                                                 Tis_new, Tis_prev, Tis_old,
                                                 0.0, 0.0, Tsa_old,
                                                 0.0, omega_sa,
                                                 dzs_new, dzs_old,
                                                 salinity_cells,
                                                 radiation_nodes_snow,
                                                 E_s, c_s, k_s, rho_s,
                                                 time_step, is_snow=True)
        
        A, B, C, RHS = concat_matrices(A_i, B_i, C_i, RHS_i,
                                       A_s, B_s, C_s, RHS_s,
                                       dzi_new, dzs_new,
                                       k_i(Tis_new, salinity_cells[-1]), k_s(Tis_new))
        
        T_vec_new = thomas_solver(A, B, C, RHS)
        
        Ti_new, Tis_new, Ts_new = T_vec_new[:len(Ti_old)],\
                                  T_vec_new[len(Ti_old)],\
                                  T_vec_new[len(Ti_old)+1:]

        # оценка ошибки
        err = err_func(T_vec_old, T_vec_prev, T_vec_new)
        curr_err = err
        
        if current_err < tol:
            break
            
    return Ti_new, Ts_new, Tis_new, dzi_new, dzs_new, cond


def main_process(time_step, time_end,
                 N_pseudoiter,
                 Ti_init, Ts_init, Tis_init, Tsa_init,
                 dzi_init, dzs_init,
                 salinity,
                 snow_thickness_threshold,
                 Toi, Ta, p,
                 F_atm_ice, F_atm_snow,
                 F_sw,
                 F_ocn
                ):
    
    time = 0.0
    
    Ti_old = Ti_init
    Ti_new = Ti_init
    Ts_old = Ts_init
    Ts_new = Ts_init
    Tis_old = Tis_init
    Tis_new = Tis_init
    Tsa_old = Tsa_init
    Tsa_new = Tsa_init
    
    
    dzi_old = dzi_init
    dzi_new = dzi_init
    dzs_old = dzs_init
    dzs_new = dzs_init
    
    salinity_cells = salinity
    
    Fs = []
    
    Ns = len(dzs_init)
    
    process = Process([dzi_init], [dzs_init],
                      [time],
                      [Toi(time)], [Ti_init], [Tis_init], [Ts_init], [Tsa_init],
                      [[rho_i]*len(dzi_init)],
                      [sum(dzs_init) >= snow_thickness_threshold])
    
    while time < time_end:
        
        if time + time_step < time_end:
            time += time_step
        else:
            time_step = time_end - time
            time = time_end
            
        Ti_old, Ts_old, Tis_old, Tsa_old, dzi_old, dzs_old\
        = Ti_new, Ts_new, Tis_new, Tsa_new, dzi_new, dzs_new
        
        if sum(dzs_new) < snow_thickness_threshold:
            print('Time {:.1f} h.: Ice freezing...'.format(time/3600))
            Ti_new, Tis_new, dzi_new, cond = ice_freezing(Toi(time),
                                                    Ti_old,
                                                    Tis_old,
                                                    lambda T: F_atm_ice(T, time),
                                                    lambda T: F_ocn(T, time),
                                                    F_sw(time),
                                                    dzi_old,
                                                    salinity_cells,
                                                    N_pseudoiter,
                                                    time_step)
            
            if Tis_new >= Tf_i(salinity_cells[-1]):
                print('Time {:.1f} h.: Ice melting...'.format(time/3600))
                Tis_new = Tf_i(salinity_cells[-1])
                Ti_new, dzi_new, cond = ice_melting(Toi(time),
                                              Ti_old,
                                              Tis_old,
                                              lambda T: F_atm_ice(T, time),
                                              lambda T: F_ocn(T, time),
                                              F_sw(time),
                                              dzi_old,
                                              salinity_cells,
                                              N_pseudoiter,
                                              time_step)
                
            # присыпем чуть снега сверху, если он падает
            if Ta(time) < 0:
                if sum(dzs_old) != 0:
                    dzs_new = Update_dz(dzs_old, 0.0, -p(time)*rho_w/rho_s, time_step)
                else:
                    dzs_new = np.full(Ns, time_step*p(time)*rho_w/rho_s/Ns)
                    
                if sum(dzs_new) >= snow_thickness_threshold:
                    # инициализируем температуру появившегося снега
                    Tsa_new = Ta(time)
                    Ts_new = Tis_new + np.arange(0.5, Ns)/Ns * (Tsa_new - Tis_new)
                    
                else:
                    # задаём профиль и и интерфейс пустыми
                    Tsa_new = np.nan
                    Ts_new = np.array([np.nan]*Ns)
                    
            process.snow_presence_history = np.append(process.snow_presence_history, False)
            
        else:
            print('Time {:.1f} h.: Snow-ice freezing...'.format(time/3600))
            Ti_new, Ts_new, Tis_new, Tsa_new, dzi_new, dzs_new, cond =\
            snow_ice_freezing(Toi(time), Ti_old, Ts_old, Tis_old, Tsa_old, Ta(time),
                              lambda T: F_atm_snow(T, time),
                              lambda T: F_ocn(T, time),
                              F_sw(time),
                              dzi_old, dzs_old,
                              salinity_cells,
                              N_pseudoiter,
                              time_step,
                              p(time))
            
            if Tsa_new >= 0:
                print('Time {:.1f} h.: Snow-ice melting...'.format(time/3600))
                Tsa_new = 0
                Ti_new, Ts_new, Tis_new, dzi_new, dzs_new, cond =\
                snow_melting(Toi(time), Ti_old, Ts_old, Tis_old, Tsa_old, Ta(time),
                             lambda T: F_atm_snow(T, time),
                             lambda T: F_ocn(T, time),
                             F_sw(time),
                             dzi_old, dzs_old,
                             salinity_cells,
                             N_pseudoiter,
                             time_step,
                             p(time))
                
            process.snow_presence_history = np.append(process.snow_presence_history, True)
            
        

        process.ice_dz_history = np.append(process.ice_dz_history, [dzi_new.copy()], axis=0)
        process.snow_dz_history = np.append(process.snow_dz_history, [dzs_new.copy() if dzs_new.sum() > 0 else 0.0*dzs_new], axis=0)
        process.timeline = np.append(process.timeline, time)
        process.oi_temp_history = np.append(process.oi_temp_history, Toi(time))
        process.ice_temp_history = np.append(process.ice_temp_history, [Ti_new.copy()], axis=0)
        process.is_temp_history = np.append(process.is_temp_history, Tis_new)
        process.snow_temp_history = np.append(process.snow_temp_history, [Ts_new.copy()], axis=0)
        process.sa_temp_history = np.append(process.sa_temp_history, Tsa_new)
        process.ice_density_history = np.append(process.ice_density_history, [[rho_i]*len(dzi_init)],
                                                axis=0)
    
    return process