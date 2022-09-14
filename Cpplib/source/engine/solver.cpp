#include "solver.hpp"

namespace icethermo
{
    template <typename NumType>
    ThermoSolver<NumType>::ThermoSolver(Mesh<NumType>* mesh_ice_,
                                        Mesh<NumType>* mesh_snow_,
                                        Dparam ice_rho_param_,
                                        Kparam ice_k_param_,
                                        Cparam ice_c_eff_param_,
                                        Eparam ice_E_param_,
                                        Lparam ice_L_param_,
                                        Dparam snow_rho_param_,
                                        Kparam snow_k_param_,
                                        Cparam snow_c_eff_param_,
                                        Eparam snow_E_param_,
                                        Lparam snow_L_param_)
    {
        this->mesh_ice = mesh_ice_;
        this->mesh_snow = mesh_snow_;
        this->ice_rho_param = ice_rho_param_;
        this->ice_k_param = ice_k_param_;
        this->ice_c_eff_param = ice_c_eff_param_;
        this->ice_E_param = ice_E_param_;
        this->ice_L_param = ice_L_param_;
        this->snow_rho_param = snow_rho_param_;
        this->snow_k_param = snow_k_param_;
        this->snow_c_eff_param = snow_c_eff_param_;
        this->snow_E_param = snow_E_param_;
        this->snow_L_param = snow_L_param_;
    }

    template<typename NumType>
    std::vector<NumType> ThermoSolver<NumType>::Update_dz(const std::vector<NumType>& dz_cells_old,
                                                          NumType omega_down,
                                                          NumType omega_up,
                                                          NumType time_step)
    {
        std::vector<NumType> dz_cells_new = dz_cells_old;
        NumType thickness = sum_vec(dz_cells_new);
        return dz_cells_new - (time_step*(omega_up - omega_down)/thickness)*dz_cells_new;
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::Update_dz(NumType dz_old,
                                             NumType omega_down,
                                             NumType omega_up,
                                             NumType time_step)
    {
        return dz_old + time_step*(omega_up - omega_down);
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::T_from_BC(const std::vector<NumType>& T_cells,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             NumType omega_value,
                                             FuncPtr<NumType> F,
                                             bool is_surface,
                                             ApproxOrder grad_approx_order,
                                             Dparam dparam,
                                             Kparam kparam,
                                             Lparam Lparam)
    {
        // surface gradient approximation
        if (is_surface)
        {
            FuncPtr<NumType> grad = [&](NumType T) {return 0.0;};

            if (grad_approx_order == ApproxOrder::first)
            {
                grad = [&](NumType T) {return (T - T_cells.back())/(0.5*dz_cells.back());};
            }
            else if (grad_approx_order == ApproxOrder::second)
            {
                int dz_size =  dz_cells.size();
                NumType h1 = 0.5*dz_cells.back();
                NumType h2 = dz_cells.back() + 0.5*dz_cells[dz_size-2];
                grad = [&](NumType T) {return (T*(h2*h2 - h1*h1) - T_cells[dz_size-1]*h2*h2 + T_cells[dz_size-2]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation orrder: first, second!");
            }

            // physical constants
            FuncPtr<NumType> rho = [&](NumType T){return Params<NumType>::Density(dparam, T, salinity_cells.back());};
            FuncPtr<NumType> k = [&](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells.back(), rho(T));};
            FuncPtr<NumType> L = [&](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells.back());};

            FuncPtr<NumType> nonlin_func = [&](NumType T){return F(T) - k(T)*grad(T) - rho(T)*L(T)*omega_value;};
            
            // solve nonlinear 1D equation
            return secant_solver<NumType>(nonlin_func, T_cells.back() - 10.0, GenConsts<NumType>::TempFusion(salinity_cells.back()) + 10.0).first;
        }
        else
        {
            // surface gradient approximation
            FuncPtr<NumType> grad = [&](NumType T) {return 0.0;};

            if (grad_approx_order == ApproxOrder::first)
            {
                grad = [&](NumType T) {return (T_cells[0] - T)/(0.5*dz_cells[0]);};
            }
            else if (grad_approx_order == ApproxOrder::second)
            {
                NumType h1 = 0.5*dz_cells[0];
                NumType h2 = dz_cells[0] + 0.5*dz_cells[1];
                grad = [&](NumType T) {return (-T*(h2*h2 - h1*h1) + T_cells[0]*h2*h2 - T_cells[1]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation orrder: first, second!");
            }

            // physical constants
            FuncPtr<NumType> rho = [&](NumType T){return Params<NumType>::Density(dparam, T, salinity_cells[0]);};
            FuncPtr<NumType> k = [&](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells[0], rho(T));};
            FuncPtr<NumType> L = [&](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells[0]);};

            // assemble nonlinear function for 1D solver
            FuncPtr<NumType> nonlin_func = [&](NumType T){return F(T) - k(T)*grad(T) - rho(T)*L(T)*omega_value;};
            
            // solve nonlinear 1D equation
            return secant_solver<NumType>(nonlin_func, T_cells.back() - 10.0, GenConsts<NumType>::TempFusion(salinity_cells.back()) + 10.0).first;
        }   
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>:: W_from_BC(NumType T_bnd,
                                              const std::vector<NumType>& T_cells,
                                              const std::vector<NumType>& dz_cells,
                                              const std::vector<NumType>& salinity_cells,
                                              FuncPtr<NumType> F,
                                              bool is_surface,
                                              ApproxOrder grad_approx_order,
                                              Dparam dparam,
                                              Kparam kparam,
                                              Lparam Lparam)
    {
        if (is_surface)
        {
            // surface gradient approximation
            FuncPtr<NumType> grad = [&](NumType T) {return 0.0;};

            if (grad_approx_order == ApproxOrder::first)
            {
                grad = [&](NumType T) {return (T - T_cells.back())/(0.5*dz_cells.back());};
            }
            else if (grad_approx_order == ApproxOrder::second)
            {
                int dz_size =  dz_cells.size();
                NumType h1 = 0.5*dz_cells.back();
                NumType h2 = dz_cells.back() + 0.5*dz_cells[dz_size-2];
                grad = [&](NumType T) {return (T*(h2*h2 - h1*h1) - T_cells[dz_size-1]*h2*h2 + T_cells[dz_size-2]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation orrder: first, second!");
            }

            // constants
            FuncPtr<NumType> rho = [&](NumType T){return Params<NumType>::Density(dparam, T, salinity_cells.back());};
            FuncPtr<NumType> k = [&](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells.back(), rho(T));};
            FuncPtr<NumType> L = [&](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells.back());};

            // calculate omega from boundary conditions
            return (F(T_bnd) - k(T_cells.back())*grad(T_bnd))/(rho(T_cells.back())*L(T_cells.back()));
        }
        else
        {
            // surface gradient approximation
            FuncPtr<NumType> grad = [&](NumType T) {return 0.0;};

            if (grad_approx_order == ApproxOrder::first)
            {
                grad = [&](NumType T) {return (T_cells[0] - T)/(0.5*dz_cells[0]);};
            }
            else if (grad_approx_order == ApproxOrder::second)
            {
                NumType h1 = 0.5*dz_cells[0];
                NumType h2 = dz_cells[0] + 0.5*dz_cells[1];
                grad = [&](NumType T) {return (-T*(h2*h2 - h1*h1) + T_cells[0]*h2*h2 - T_cells[1]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation orrder: first, second!");
            }

            // physical constants
            FuncPtr<NumType> rho = [&](NumType T){return Params<NumType>::Density(dparam, T, salinity_cells[0]);};
            FuncPtr<NumType> k = [&](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells[0], rho(T));};
            std::function<NumType(NumType, NumType)> L = [&](NumType T, NumType S){return Params<NumType>::FusionHeat(Lparam, T, S);};

            // calculate omegs from boundary conditions
            if (kparam == Kparam::BubblyBrine or kparam == Kparam::Untersteiner)
            {
                if ((F(T_bnd) - k(T_cells[0])*grad(T_bnd)) > 0)
                {
                    // basal growth with 10 psu salinity
                    return (F(T_bnd) - k(T_cells[0])*grad(T_bnd))/(rho(T_cells[0])*L(T_cells[0], 10.0));
                }
                else
                {
                    // basal melt with 4 psu salinity
                    return (F(T_bnd) - k(T_cells[0])*grad(T_bnd))/(rho(T_cells[0])*L(T_cells[0], 4.0));
                }
            }
            else
            {
                return (F(T_bnd) - k(T_cells[0])*grad(T_bnd))/(rho(T_cells[0])*L(T_cells[0], salinity_cells[0]));
            }
        }
    }

/*
    template<typename NumType>
    FourVecs<NumType> ThermoSolver<NumType>::Assemble_advdiff_martix_rhs(const std::vector<NumType>& T_cells_prev, const std::vector<NumType>& T_cells_old,
                                                                         NumType T_up_new, NumType T_up_prev, NumType T_up_old,
                                                                         NumType T_down_new, NumType T_down_prev, NumType T_down_old, 
                                                                         NumType omega_down, NumType omega_up, 
                                                                         const std::vector<NumType>& dz_cells_new, const std::vector<NumType>& dz_cells_old,
                                                                         const std::vector<NumType>& salinity_cells,
                                                                         const std::vector<NumType>& radiation_nodes,
                                                                         NumType time_step,
                                                                         Dparam dparam,
                                                                         Kparam kparam,
                                                                         Lparam Lparam,
                                                                         Eparam Eparam,
                                                                         Cparam cparam)
    {
        // get number of cells
        int N = T_cells_old.size();

        // check the size of input arrays
        if (T_cells_prev.size() != N)
        {
            THERMO_ERR("Wrong size of \'T_cells_prev\'!");
        }

        if (dz_cells_new.size() != N)
        {
            THERMO_ERR("Wrong size of \'dz_cells_new\'!")
        }

        if (dz_cells_old.size() != N)
        {
            THERMO_ERR("Wrong size of \'dz_cells_old\'!");
        }

        if (radiation_nodes.size() != N+1)
        {
            THERMO_ERR("Wrong size of \'radiation_nodes\'!");
        }

        // compute nodes effective thermal conductivity
        std::vector<NumType> nodes_eff_k(N+1);
        
        
        nodes_eff_k[0] = 2.0*Params<NumType>::Conductivity(kparam, T_cells_prev[0], salinity_cells[0], Params<NumType>::Density(dparam, T_cells_prev[0], salinity_cells[0]))/dz_cells_new[0];
        
        for (int i = 1; i < N; ++i)
        {
            NumType k_prev = Params<NumType>::Conductivity(kparam, T_cells_prev[i-1], salinity_cells[i-1], Params<NumType>::Density(dparam, T_cells_prev[i-1], salinity_cells[i-1]));
            NumType k_forw = Params<NumType>::Conductivity(kparam, T_cells_prev[i], salinity_cells[i], Params<NumType>::Density(dparam, T_cells_prev[i], salinity_cells[i]));
            nodes_eff_k[i] = 2.0*k_prev*k_forw/(k_prev*dz_cells_new[i] + k_forw*dz_cells_new[i-1]);
        }

        
        nodes_eff_k[N] = 2.0*Params<NumType>::Conductivity(kparam, T_cells_prev[N-1], salinity_cells[N-1], Params<NumType>::Density(dparam, T_cells_prev[N-1], salinity_cells[N-1]))/dz_cells_new[N-1];

        // compute effective heat capacity and enthalpy at the cells and interface nodes
        std::vector<NumType> cells_eff_c(N);
        std::vector<NumType> cells_E(N);


        NumType c_eff_down = Params<NumType>::EffCapacity(cparam, T_down_new, T_down_old, salinity_cells[0]);
        NumType E_down = Params<NumType>::Enthalpy(Eparam, T_down_old, salinity_cells[0]); 
        
        for (int i = 0; i < N; ++i)
        {
            cells_eff_c[i] = Params<NumType>::EffCapacity(cparam, T_cells_prev[i], T_cells_old[i], salinity_cells[i]);
            cells_E[i] = Params<NumType>::Enthalpy(Eparam, T_cells_old[i],  salinity_cells[i]);
        }

        NumType c_eff_up = Params<NumType>::EffCapacity(cparam, T_up_new, T_up_old, salinity_cells[N-1]);
        NumType E_up = Params<NumType>::Enthalpy(Eparam, T_up_old, salinity_cells[N-1]);

        // compute nodal values of omega
        std::vector<NumType> nodes_omega(N+1);

        nodes_omega[0] = omega_down;

        for (int i = 1; i < N; ++i)
        {
            nodes_omega[i] = omega_down + (sum_vec(dz_cells_new, 0, i)/sum_vec(dz_cells_new))*(omega_up - omega_down);
        }

        nodes_omega[N] = omega_up;
    

    
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

    }
*/
    // explicit instantiation
    template class ThermoSolver<float>;
    template class ThermoSolver<double>;
}