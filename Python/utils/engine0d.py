import utils.engine as ue
import numpy as np

S_i_const = 2.5
k_i_const = ue.k_i(-5.0, S_i_const)
k_s_const = ue.k0_s

def Update_thickness(thickness_old,
                     omega_down,
                     omega_up,
                     time_step):
    return thickness_old + time_step*(omega_up - omega_down)

def W_from_BC_0d(T_up,
                 T_down,
                 thickness,
                 k,
                 F,
                 L,
                 is_surface=False,
                 rho=ue.rho_i):
    
    if is_surface:
        
        T_melt = T_up + (T_down - T_up)/thickness*0.1
        grad = (T_up - T_down)/thickness
        omega = -(F - k*grad)/(rho*L(T_melt, 1.0))
        return omega
        
    else:
        
        grad = (T_up - T_down)/thickness
                    
        if (F + k*grad) > 0:
            # basal melt
            omega = (F + k*grad)/(rho*L(T_down, 4.0))
        else:
            # basal growth
            omega = (F + k*grad)/(rho*L(T_down, 10.0))
            
        return omega
    
def T_from_BC_0d(T_down,
                 thickness,
                 k,
                 omega_up,
                 F,
                 L,
                 is_surface=False,
                 rho=ue.rho_i):
    
    if is_surface:
        grad_su = lambda T: (T - T_down)/thickness
        nonlin_func = lambda T: rho*L(T, S_i_const)*omega_up + F(T) - k*(grad_su(T))
        return ue.secant_solver(nonlin_func, -50.0, 5.0)[0]
    
    else:
        # not needed
        return 0.0
    

def ice_freezing_0d(Toi,
                    Tia_init,
                    F_atm,
                    F_ocn,
                    thickness_i_init,
                    N_pseudoiter,
                    time_step,
                    tol=1e-12,
                    err_func=lambda T_old, T_prev, T_new: np.abs(T_new - T_prev)/(np.abs(T_old) + 5.0)
                    ):
    
    
    Tia_old = Tia_init
    Tia_new = Tia_init
    
    thickness_old = thickness_i_init
    thickness_new = thickness_i_init
    
    Tia_values = [Tia_old]
    surface_err = np.inf
    prev_surface_err = surface_err
    
    for pseudoiter in range(N_pseudoiter):
        
        Tia_prev = Tia_new
        thickness_prev = thickness_new
        
        # оценка omega на границе океан-лед
        omega_io = W_from_BC_0d(T_up=Tia_prev,
                                T_down=Toi,
                                thickness=thickness_prev,
                                k=k_i_const,
                                F=F_ocn(Toi), 
                                L=lambda T, S_i: -ue.L_i(T, S_i),
                                is_surface=False,
                                rho=ue.rho_i)
        

        # оценка Tia
        Tia_new = T_from_BC_0d(T_down=Toi,
                               thickness=thickness_prev,
                               k=k_i_const,
                               omega_up=0.0, 
                               F=F_atm,
                               L=lambda T, S_i: -ue.L_i(T, S_i),
                               is_surface=True)
        
        # пересчет толщины льда
        thickness_new = Update_thickness(thickness_old=thickness_old,
                                         omega_down=omega_io,
                                         omega_up=0.0,
                                         time_step=time_step)
        
        # оценка ошибки
        prev_surface_err = surface_err
        surface_err = err_func(T_old=Tia_old,
                               T_prev=Tia_prev,
                               T_new=Tia_new)
        
        if (surface_err < prev_surface_err):
            Tia_values.append(Tia_new)
        else:
            Tia_new = sum(Tia_values)/len(Tia_values)
        
        curr_err = err_func(T_old=Tia_old,
                            T_prev=Tia_prev,
                            T_new=Tia_new)
        
        if curr_err < tol:
            break
    
    return Tia_new, thickness_new

def ice_melting_0d(Toi,
                   F_atm,
                   F_ocn,
                   thickness_i_init,
                   N_pseudoiter,
                   time_step,
                   tol = 1e-12,
                   err_func=lambda h_old, h_prev, h_new: np.abs(h_new - h_prev)/(np.abs(h_old))):
    
    thickness_new, thickness_prev = thickness_i_init, thickness_i_init
    
    for pseudoiter in range(N_pseudoiter):
        
        thickness_prev = thickness_new
    
        # оценка omega на границе океан-лед
        omega_io = W_from_BC_0d(T_up=ue.Tf_i(1.0),
                                T_down=Toi,
                                thickness=thickness_prev,
                                k=k_i_const,
                                F=F_ocn(Toi), 
                                L=lambda T, S_i: -ue.L_i(T, S_i),
                                is_surface=False,
                                rho=ue.rho_i)

        # оценка omega на границе лед-атмосфера
        omega_ia = W_from_BC_0d(T_up=ue.Tf_i(1.0),
                                T_down=Toi,
                                thickness=thickness_prev,
                                k=k_i_const,
                                F=F_atm(ue.Tf_i(1.0)), 
                                L=lambda T, S_i: -ue.L_i(T, S_i),
                                is_surface=True,
                                rho=ue.rho_i)

         # пересчет толщины льда
        thickness_new = Update_thickness(thickness_old=thickness_i_init,
                                         omega_down=omega_io,
                                         omega_up=omega_ia,
                                         time_step=time_step)
        
                 
        # оценка ошибки
        err = err_func(h_old=thickness_i_init,
                        h_prev=thickness_prev,
                        h_new=thickness_new)
        
        if (err < tol):
            break

    return thickness_new

def snow_ice_freezing_0d(Toi,
                         Tis_init,
                         Tsa_init,
                         Ta,
                         F_atm,
                         F_ocn,
                         ice_thickness_init,
                         snow_thickness_init,
                         N_pseudoiter,
                         time_step,
                         p,
                         tol=1e-12,
                         err_func=lambda T_is_old, T_is_prev, T_is_new, T_sa_old, T_sa_prev, T_sa_new: \
                         np.abs(T_is_new - T_is_prev)/(np.abs(T_is_old)+ 5.0) + \
                         np.abs(T_sa_new - T_sa_prev)/(np.abs(T_sa_old)+ 5.0)):
    
    # инициализация массивов
    T_is_prev, T_is_new = Tis_init, Tis_init
    T_sa_prev, T_sa_new = Tsa_init, Tsa_init
    ice_thickness_new, ice_thickness_prev = ice_thickness_init, ice_thickness_init
    snow_thickness_new, snow_thickness_prev = snow_thickness_init, snow_thickness_init
    
    Tsa_values = [Tsa_init]
    surface_err, prev_surface_err = np.inf, np.inf
    
    # omega наверху определяется только осадками
    omega_sa = p*ue.rho_w/ue.rho_s if (Ta < 0.0) else 0.0
    
    for pseudoiter in range(N_pseudoiter):
        
        T_is_prev, T_sa_prev, ice_thickness_prev, snow_thickness_prev = \
        T_is_new, T_sa_new, ice_thickness_new, snow_thickness_new
        
        prev_surface_err = surface_err
        
        # оценка omega на границе океан-лед
        omega_io = W_from_BC_0d(T_up=T_is_prev,
                                T_down=Toi,
                                thickness=ice_thickness_prev,
                                k=k_i_const,
                                F=F_ocn(Toi), 
                                L=lambda T, S_i: -ue.L_i(T, S_i),
                                is_surface=False,
                                rho=ue.rho_i)
        
        # оценка Tsa
        T_sa_new = T_from_BC_0d(T_down=T_is_prev,
                                thickness=snow_thickness_prev,
                                k=k_s_const,
                                omega_up=omega_sa, 
                                F=F_atm,
                                L=lambda T, S_i: -ue.L_s(T, S_i),
                                is_surface=True,
                                rho=ue.rho_s)
        
        # пересчет толщины льда
        ice_thickness_new = Update_thickness(thickness_old=ice_thickness_init,
                                             omega_down=omega_io,
                                             omega_up=0.0,
                                             time_step=time_step)
        
        # пересчет толщины снега
        snow_thickness_new = Update_thickness(thickness_old=snow_thickness_init,
                                              omega_down=0.0,
                                              omega_up=omega_sa,
                                              time_step=time_step)
        
        # пересчет температуры интерфейса
        ratio = (k_i_const*snow_thickness_new)/(k_s_const*ice_thickness_new)
        T_is_new  = (1.0/(ratio + 1.0))*T_sa_new + ((ratio)/(ratio + 1.0))*Toi
        
        # оценка ошибки
        surface_err = err_func(T_is_old=Tis_init,
                               T_is_prev=T_is_prev,
                               T_is_new=T_is_new,
                               T_sa_old=Tsa_init,
                               T_sa_prev=T_sa_prev,
                               T_sa_new=T_sa_new)
                               
        if (surface_err < prev_surface_err):
            Tsa_values.append(T_is_new)
        else:
            T_sa_new = sum(Tsa_values)/len(Tsa_values)
        
        surface_err = err_func(T_is_old=Tis_init,
                               T_is_prev=T_is_prev,
                               T_is_new=T_is_new,
                               T_sa_old=Tsa_init,
                               T_sa_prev=T_sa_prev,
                               T_sa_new=T_sa_new)
        
        if surface_err < tol:
            break
            
    return T_is_new, T_sa_new, ice_thickness_new, snow_thickness_new

def snow_melting_0d(Toi,
                    Tis,
                    Ta,
                    F_atm,
                    F_ocn,
                    ice_thickness_init,
                    snow_thickness_init,
                    N_pseudoiter,
                    time_step,
                    p,
                    tol=1e-12,
                    err_func=lambda T_is_old, T_is_prev, T_is_new: \
                    np.abs(T_is_new - T_is_prev)/(np.abs(T_is_old) + 5.0)):
    
    # инициализация массивов
    T_is_old, T_is_prev, T_is_new = Tis, Tis, Tis
    ice_thickness_prev, ice_thickness_new = ice_thickness_init, ice_thickness_init
    snow_thickness_old, snow_thickness_new = snow_thickness_init, snow_thickness_init
    
    err = np.inf
    
    for pseudoiter in range(N_pseudoiter):
    
        T_is_prev = T_is_new
        ice_thickness_prev, snow_thickness_prev = ice_thickness_new, snow_thickness_new
                    
        # оценка omega на границе океан-лед
        omega_io = W_from_BC_0d(T_up=T_is_prev,
                                T_down=Toi,
                                thickness=ice_thickness_prev,
                                k=k_i_const,
                                F=F_ocn(Toi), 
                                L=lambda T, S_i: -ue.L_i(T, S_i),
                                is_surface=False,
                                rho=ue.rho_i)
        
        # оценка omega на границе снег-атмесфера
        omega_sa = W_from_BC_0d(T_up=0.0,
                                T_down=T_is_prev,
                                thickness=snow_thickness_prev,
                                k=k_s_const,
                                F=F_atm(0.0), 
                                L=lambda T, S_i: -ue.L_s(T, S_i),
                                is_surface=True,
                                rho=ue.rho_s)
        
        if Ta < 0:
            omega_sa += p*ue.rho_w/ue.rho_s
        
        # пересчет толщины льда
        ice_thickness_new = Update_thickness(thickness_old=ice_thickness_init,
                                             omega_down=omega_io,
                                             omega_up=0.0,
                                             time_step=time_step)
        # пересчет толщины снега
        snow_thickness_new = Update_thickness(thickness_old=snow_thickness_init,
                                              omega_down=0.0,
                                              omega_up=omega_sa,
                                              time_step=time_step)
        # пересчет температуры интерфейса
        ratio = (k_i_const*snow_thickness_new)/(k_s_const*ice_thickness_new)
        T_is_new  = (1.0/(ratio + 1.0))*0.0 + ((ratio)/(ratio + 1.0))*Toi
                    
        # оценка ошибки
        err = err_func(T_is_old=Tis,
                       T_is_prev=T_is_prev,
                       T_is_new=T_is_new)
        
        if err < tol:
            break
            
    return T_is_new, ice_thickness_new, snow_thickness_new

def make_linear_profile_ice(h_ice, T_ib, T_is, N_cells):
    dz_cells = np.ones(N_cells)*(h_ice/N_cells)
    T_cells = T_ib + np.arange(0.5, N_cells)/N_cells * (T_is - T_ib)
    return dz_cells, T_ib, T_cells, T_is

def main_process_0d(time_step,
                    time_end,
                    N_pseudoiter,
                    Tis_init,
                    Tsa_init,
                    ice_thickness_init,
                    snow_thickness_init,
                    snow_thickness_threshold,
                    Toi,
                    Ta,
                    p,
                    F_atm_ice,
                    F_atm_snow,
                    F_ocn,
                    N_cells_ice=5,
                    N_cells_snow=10):
    
    time = 0.0
    
    T_is_old, T_is_new, T_sa_old, T_sa_new = \
    Tis_init, Tis_init, Tsa_init, Tsa_init
    
    ice_thickness_old, ice_thickness_new, snow_thickness_old, snow_thickness_new = \
    ice_thickness_init, ice_thickness_init, snow_thickness_init, snow_thickness_init
    
    dz_cells_ice_init, T_ib_init, T_ice_cells_init, T_is_init = make_linear_profile_ice(h_ice=ice_thickness_init,
                                                                                        T_ib=Toi(0.0),
                                                                                        T_is=Tis_init,
                                                                                        N_cells=N_cells_ice)
    
    dz_cells_snow_init, T_is_init, T_snow_cells_init, T_sa_init = make_linear_profile_ice(h_ice=snow_thickness_init,
                                                                                          T_ib=Tis_init,
                                                                                          T_is=Tsa_init,
                                                                                          N_cells=N_cells_snow)
    
    process = ue.Process([dz_cells_ice_init], [dz_cells_snow_init],
                         [time],
                         [Toi(time)], [T_ice_cells_init], [T_is_init], [T_snow_cells_init], [T_sa_init],
                         [[ue.rho_i]*N_cells_ice],
                         [snow_thickness_init >= snow_thickness_threshold])
    
    while time < time_end:
        
        if time + time_step < time_end:
            time += time_step
        else:
            time_step = time_end - time
            time = time_end
            
        T_is_old, T_sa_old = T_is_new, T_sa_new
        ice_thickness_old, snow_thickness_old = ice_thickness_new, snow_thickness_new
        
        if snow_thickness_new < snow_thickness_threshold:
            print('Time {:.1f} h.: Ice freezing...'.format(time/3600))
            T_is_new, ice_thickness_new = ice_freezing_0d(Toi=Toi(time),
                                                          Tia_init=T_is_old,
                                                          F_atm=lambda T: F_atm_ice(T, time),
                                                          F_ocn=lambda T: F_ocn(T, time),
                                                          thickness_i_init=ice_thickness_old,
                                                          N_pseudoiter=N_pseudoiter,
                                                          time_step=time_step)
            
            if T_is_new >= ue.Tf_i(1.0):
                print('Time {:.1f} h.: Ice melting...'.format(time/3600))
                T_is_new = ue.Tf_i(1.0)
                ice_thickness_new = ice_melting_0d(Toi=Toi(time),
                                                   F_atm=lambda T: F_atm_ice(T, time),
                                                   F_ocn=lambda T: F_ocn(T, time),
                                                   thickness_i_init=ice_thickness_old,
                                                   N_pseudoiter=N_pseudoiter,
                                                   time_step=time_step)
                
            # присыпем чуть снега сверху, если он падает
            if Ta(time) < 0.0:
                
                snow_thickness_new = Update_thickness(thickness_old=snow_thickness_old,
                                                      omega_down=0.0,
                                                      omega_up=p(time)*ue.rho_w/ue.rho_s,
                                                      time_step=time_step)
                    
                if snow_thickness_new >= snow_thickness_threshold:
                    # инициализируем температуру появившегося снега
                    T_sa_new = Ta(time)
                    
                else:
                    # задаём профиль и и интерфейс пустыми
                    T_sa_new = np.nan
                    
            process.snow_presence_history = np.append(process.snow_presence_history, False)
            
        else:
            print('Time {:.1f} h.: Snow-ice freezing...'.format(time/3600))
            T_is_new, T_sa_new, ice_thickness_new, snow_thickness_new = \
            snow_ice_freezing_0d(Toi=Toi(time),
                                 Tis_init=T_is_old,
                                 Tsa_init=T_sa_old,
                                 Ta=Ta(time),
                                 F_atm=lambda T: F_atm_snow(T, time),
                                 F_ocn=lambda T: F_ocn(T, time),
                                 ice_thickness_init=ice_thickness_old,
                                 snow_thickness_init=snow_thickness_old,
                                 N_pseudoiter=N_pseudoiter,
                                 time_step=time_step,
                                 p=p(time))
            
            
            if T_sa_new >= 0.0:
                print('Time {:.1f} h.: Snow-ice melting...'.format(time/3600))
                T_sa_new = 0.0
                T_is_new, ice_thickness_new, snow_thickness_new = \
                snow_melting_0d(Toi=Toi(time),
                                Tis=T_is_old,
                                Ta=Ta(time),
                                F_atm=lambda T: F_atm_snow(T, time),
                                F_ocn=lambda T: F_ocn(T, time),
                                ice_thickness_init=ice_thickness_old,
                                snow_thickness_init=snow_thickness_old,
                                N_pseudoiter=N_pseudoiter,
                                time_step=time_step,
                                p=p(time))
                
            process.snow_presence_history = np.append(process.snow_presence_history, True)
            
        # вывод
        dz_cells_ice_out, T_ib_out, T_ice_cells_out, T_is_out = make_linear_profile_ice(h_ice=ice_thickness_new,
                                                                                        T_ib=Toi(time),
                                                                                        T_is=T_is_new,
                                                                                        N_cells=N_cells_ice)
        
        dz_cells_snow_out, T_is_out, T_snow_cells_out, T_sa_out = make_linear_profile_ice(h_ice=snow_thickness_new,
                                                                                          T_ib=T_is_new,
                                                                                          T_is=T_sa_new,
                                                                                          N_cells=N_cells_snow)
    
        
        process.ice_dz_history = np.append(process.ice_dz_history, [dz_cells_ice_out.copy()], axis=0)
        process.snow_dz_history = np.append(process.snow_dz_history, [dz_cells_snow_out.copy() if snow_thickness_new > 0.0 else 0.0*dz_cells_snow_out], axis=0)
        process.timeline = np.append(process.timeline, time)
        process.oi_temp_history = np.append(process.oi_temp_history, Toi(time))
        process.ice_temp_history = np.append(process.ice_temp_history, [T_ice_cells_out.copy()], axis=0)
        process.is_temp_history = np.append(process.is_temp_history, T_is_out)
        process.snow_temp_history = np.append(process.snow_temp_history, [T_snow_cells_out.copy()], axis=0)
        process.sa_temp_history = np.append(process.sa_temp_history, T_sa_out)
        process.ice_density_history = np.append(process.ice_density_history, [[ue.rho_i]*N_cells_ice], axis=0)
    
    return process
