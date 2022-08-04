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