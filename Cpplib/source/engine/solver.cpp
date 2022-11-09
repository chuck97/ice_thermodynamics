#include "solver.hpp"

namespace icethermo
{
    template <typename NumType>
    ThermoSolver<NumType>::ThermoSolver(Mesh<NumType>* mesh_ice_,
                                        Mesh<NumType>* mesh_snow_,
                                        NumType time_step_,
                                        ApproxOrder grad_approx_order_,
                                        Kparam ice_k_param_,
                                        Cparam ice_c_eff_param_,
                                        Eparam ice_E_param_,
                                        Lparam ice_L_param_,
                                        Kparam snow_k_param_,
                                        Cparam snow_c_eff_param_,
                                        Eparam snow_E_param_,
                                        Lparam snow_L_param_)
    {
        this->mesh_ice = mesh_ice_;
        this->mesh_snow = mesh_snow_;
        this->time_step = time_step_;
        this->grad_approx_order = grad_approx_order_;
        this->ice_k_param = ice_k_param_;
        this->ice_c_eff_param = ice_c_eff_param_;
        this->ice_E_param = ice_E_param_;
        this->ice_L_param = ice_L_param_;
        this->snow_k_param = snow_k_param_;
        this->snow_c_eff_param = snow_c_eff_param_;
        this->snow_E_param = snow_E_param_;
        this->snow_L_param = snow_L_param_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateForcing(FuncPtr<NumType> F_up_,
                                              FuncPtr<NumType> F_down_,
                                              FuncPtr<NumType> F_sw_,
                                              FuncPtr<NumType> prec_rate_)
    {
        this->F_up = F_up_;
        this->F_down = F_down_;
        this->F_sw = F_sw_;
        this->prec_rate = prec_rate_;
    }

    template<typename NumType>
    std::vector<NumType> ThermoSolver<NumType>::Update_dz(const std::vector<NumType>& dz_cells_old,
                                                          NumType omega_down,
                                                          NumType omega_up)
    {
        NumType thickness = sum_vec(dz_cells_old);
        return dz_cells_old - (this->time_step*(omega_up - omega_down)/thickness)*dz_cells_old;
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::Update_dz(NumType dz_old,
                                             NumType omega_down,
                                             NumType omega_up)
    {
        return dz_old + this->time_step*(omega_up - omega_down);
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::T_from_BC(const std::vector<NumType>& T_cells,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             const std::vector<NumType>& rho_cells,
                                             NumType omega_value,
                                             bool is_ice,
                                             bool is_surface)
    {
        Kparam kparam = (is_ice) ? this->ice_k_param: this->snow_k_param;
        Lparam Lparam = (is_ice) ? this->ice_L_param: this->snow_L_param;

        // surface gradient approximation
        if (is_surface)
        {
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T - T_cells.back())/(0.5*dz_cells.back());};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                int dz_size =  dz_cells.size();
                NumType h1 = (NumType)0.5*dz_cells.back();
                NumType h2 = (NumType)dz_cells.back() + (NumType)0.5*dz_cells[dz_size-2];
                grad = [&T_cells, h1, h2, dz_size](NumType T) {return (T*(h2*h2 - h1*h1) - T_cells.back()*h2*h2 + T_cells[dz_size-2]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // physical constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells.back(), rho_cells.back());};
            FuncPtr<NumType> L = [&salinity_cells, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells.back());};

            FuncPtr<NumType> nonlin_func = [&rho_cells, omega_value, this, &k, &grad, &L](NumType T){return rho_cells.back()*L(T)*omega_value + this->F_up(T) - k(T)*grad(T);};
            
            // solve nonlinear 1D equation
            auto secant_res = secant_solver<NumType>(nonlin_func, T_cells.back() - (NumType)10.0, GenConsts<NumType>::TempFusion(salinity_cells.back() + (NumType)10.0));

            // handle possible errors in secant solver
            /*
            if (std::get<3>(secant_res))
            {
                return std::get<0>(secant_res);
            }
            else
            {
                if ((std::abs(std::get<1>(secant_res)/nonlin_func(T_cells.back() - 10.0)) < (NumType)ALLOWABLE_RELATIVE_1D_ERROR) or 
                    (std::get<0>(secant_res) > GenConsts<NumType>::TempFusion(salinity_cells.back())))
                {
                    return std::get<0>(secant_res);
                }
                else
                {
                    std::cout << std::get<0>(secant_res) << " " << std::get<1>(secant_res) << " " << std::get<2>(secant_res) << " " << std::get<3>(secant_res) << std::endl;
                    std::cout << std::get<1>(secant_res) << std::endl;
                    std::cout << std::abs(std::get<1>(secant_res)/nonlin_func(T_cells.back() - 10.0)) << std::endl;
                    std::cout << (NumType)ALLOWABLE_RELATIVE_1D_ERROR << std::endl;
                    THERMO_ERR("Impossible to solve boundary conditions on top interface!");
                }
            }
            */
           return std::get<0>(secant_res);
        }
        else
        {
            // surface gradient approximation
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T_cells[0] - T)/(0.5*dz_cells[0]);};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                NumType h1 = (NumType)0.5*dz_cells[0];
                NumType h2 = dz_cells[0] + (NumType)0.5*dz_cells[1];
                grad = [&T_cells, h1, h2](NumType T) {return (-T*(h2*h2 - h1*h1) + T_cells[0]*h2*h2 - T_cells[1]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // physical constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells[0], rho_cells[0]);};
            FuncPtr<NumType> L = [&salinity_cells, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells[0]);};

            // assemble nonlinear function for 1D solver
            FuncPtr<NumType> nonlin_func = [&rho_cells, omega_value, this, &k, &grad, &L](NumType T){return rho_cells[0]*L(T)*omega_value + this->F_down(T) - k(T)*grad(T);};
            
            // solve nonlinear 1D equation
            auto secant_res = secant_solver<NumType>(nonlin_func, T_cells[0] - (NumType)10.0, GenConsts<NumType>::TempFusion(salinity_cells[0]) + (NumType)10.0);
            
            // handle possible errors in secant solver
            /*
            if (std::get<3>(secant_res))
            {
                return std::get<0>(secant_res);
            }
            else
            {
                if (std::abs(std::get<1>(secant_res)/nonlin_func(T_cells[0] - 10.0)) < (NumType)ALLOWABLE_RELATIVE_1D_ERROR)
                {
                    return std::get<0>(secant_res);
                }
                else
                {
                    THERMO_ERR("Impossible to solve boundary conditions on bottom interface!");
                }
            }
            */
            return std::get<0>(secant_res);
        }   
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>:: W_from_BC(NumType T_bnd,
                                              const std::vector<NumType>& T_cells,
                                              const std::vector<NumType>& dz_cells,
                                              const std::vector<NumType>& salinity_cells,
                                              const std::vector<NumType>& rho_cells,
                                              bool is_ice,
                                              bool is_surface)
    {
        Kparam kparam = (is_ice) ? this->ice_k_param: this->snow_k_param;
        Lparam Lparam = (is_ice) ? this->ice_L_param: this->snow_L_param;

        if (is_surface)
        {
            // surface gradient approximation
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T - T_cells.back())/(0.5*dz_cells.back());};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                int dz_size =  dz_cells.size();
                NumType h1 = 0.5*dz_cells.back();
                NumType h2 = dz_cells.back() + 0.5*dz_cells[dz_size-2];
                grad = [&T_cells, h1, h2, dz_size](NumType T) {return (T*(h2*h2 - h1*h1) - T_cells.back()*h2*h2 + T_cells[dz_size-2]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells.back(), rho_cells.back());};
            FuncPtr<NumType> L = [&salinity_cells, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells.back());};

            // calculate omega from boundary conditions
            return (k(T_cells.back())*grad(T_bnd) - this->F_up(T_bnd))/(rho_cells.back()*L(T_cells.back()));
        }
        else
        {
            // surface gradient approximation
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T_cells[0] - T)/(0.5*dz_cells[0]);};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                NumType h1 = 0.5*dz_cells[0];
                NumType h2 = dz_cells[0] + 0.5*dz_cells[1];
                grad = [&T_cells, h1, h2](NumType T) {return (-T*(h2*h2 - h1*h1) + T_cells[0]*h2*h2 - T_cells[1]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // physical constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells[0], rho_cells[0]);};
            std::function<NumType(NumType, NumType)> L = [&salinity_cells, Lparam](NumType T, NumType S){return Params<NumType>::FusionHeat(Lparam, T, S);};

            // calculate omega from boundary conditions
            if (kparam == Kparam::BubblyBrine or kparam == Kparam::Untersteiner)
            {
                if ((k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd)) > 0)
                {
                    // basal melt with 4 psu salinity
                    return (k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd))/(rho_cells[0]*L(T_cells[0], 4.0));
                }
                else
                {
                    // basal growth with 10 psu salinity
                    return (k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd))/(rho_cells[0]*L(T_cells[0], 10.0));
                }
            }
            else
            {
                return (k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd))/(rho_cells[0]*L(T_cells[0], salinity_cells[0]));
            }
        }
    }


    template<typename NumType>
    FourVecs<NumType> ThermoSolver<NumType>::Assemble_advdiff_martix_rhs(const std::vector<NumType>& T_cells_prev,
                                                                         const std::vector<NumType>& T_cells_old,
                                                                         NumType T_up_new, 
                                                                         NumType T_up_prev, 
                                                                         NumType T_up_old,
                                                                         NumType T_down_new, 
                                                                         NumType T_down_prev, 
                                                                         NumType T_down_old, 
                                                                         NumType omega_down, 
                                                                         NumType omega_up, 
                                                                         const std::vector<NumType>& dz_cells_new, 
                                                                         const std::vector<NumType>& dz_cells_old,
                                                                         const std::vector<NumType>& salinity_cells,
                                                                         const std::vector<NumType>& rho_cells,
                                                                         const std::vector<NumType>& radiation_nodes,
                                                                         bool is_ice)
    {
        // get all parameterizations
        Kparam kparam = (is_ice) ? this->ice_k_param: this->snow_k_param;
        Cparam cparam = (is_ice) ? this->ice_c_eff_param: this->snow_c_eff_param;
        Eparam Eparam = (is_ice) ? this->ice_E_param: this->snow_E_param;

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

        // compute effective thermal conductivity at nodes
        std::vector<NumType> eff_k_nodes(N+1);
        
        
        eff_k_nodes[0] = Params<NumType>::Conductivity(kparam, T_cells_prev[0], salinity_cells[0], rho_cells[0])/dz_cells_new[0];
        
        for (int i = 1; i < N; ++i)
        {
            NumType k_prev = Params<NumType>::Conductivity(kparam, T_cells_prev[i-1], salinity_cells[i-1], rho_cells[i-1]);
            NumType k_forw = Params<NumType>::Conductivity(kparam, T_cells_prev[i], salinity_cells[i], rho_cells[i]);
            eff_k_nodes[i] = 2.0*k_prev*k_forw/(k_prev*dz_cells_new[i] + k_forw*dz_cells_new[i-1]);
        }

        eff_k_nodes[N] = Params<NumType>::Conductivity(kparam, T_cells_prev[N-1], salinity_cells[N-1], rho_cells[N-1])/dz_cells_new[N-1];

        // compute effective heat capacity and enthalpy at the cells and interface nodes
        std::vector<NumType> eff_c_cells(N);
        std::vector<NumType> E_cells(N);


        NumType eff_c_down = Params<NumType>::EffCapacity(cparam, T_down_prev, T_down_old, salinity_cells[0]);
        NumType E_down = Params<NumType>::Enthalpy(Eparam, T_down_old, salinity_cells[0]); 
        
        for (int i = 0; i < N; ++i)
        {
            eff_c_cells[i] = Params<NumType>::EffCapacity(cparam, T_cells_prev[i], T_cells_old[i], salinity_cells[i]);
            E_cells[i] = Params<NumType>::Enthalpy(Eparam, T_cells_old[i],  salinity_cells[i]);
        }

        NumType eff_c_up = Params<NumType>::EffCapacity(cparam, T_up_prev, T_up_old, salinity_cells[N-1]);
        NumType E_up = Params<NumType>::Enthalpy(Eparam, T_up_old, salinity_cells[N-1]);

        // compute nodal values of omega
        std::vector<NumType> omega_nodes(N+1);

        omega_nodes[0] = omega_down;

        for (int i = 1; i < N; ++i)
        {
            omega_nodes[i] = omega_down + (sum_vec(dz_cells_new, 0, i)/sum_vec(dz_cells_new))*(omega_up - omega_down);
        }

        omega_nodes[N] = omega_up;

        // construct diagonals of matrix and rhs vector
        std::vector<NumType> A(N);
        std::vector<NumType> B(N);
        std::vector<NumType> C(N);
        std::vector<NumType> RHS(N);

        // first row
        B[0] = rho_cells[0]*eff_c_cells[0]*dz_cells_new[0]/this->time_step + eff_k_nodes[1] + eff_k_nodes[0];

        C[0] = -eff_k_nodes[1];

        RHS[0] = rho_cells[0]*eff_c_cells[0]*dz_cells_new[0]*T_cells_old[0]/this->time_step - 
                  rho_cells[0]*E_cells[0]*(dz_cells_new[0] - dz_cells_old[0])/this->time_step +
                  (radiation_nodes[1] - radiation_nodes[0]) + 
                  eff_k_nodes[0]*T_down_new;
        
        if (omega_nodes[0] >= 0)
        {
            RHS[0] += rho_cells[0]*eff_c_down*T_down_new*omega_nodes[0] - 
                      rho_cells[0]*(eff_c_down*T_down_old - E_down)*omega_nodes[0];
        }
        else
        {
            B[0] += -rho_cells[0]*eff_c_cells[0]*omega_nodes[0];

            RHS[0] += -rho_cells[0]*(eff_c_cells[0]*T_cells_old[0] - E_cells[0])*omega_nodes[0];
        }

        if (omega_nodes[1] >= 0)
        {
            B[0] += rho_cells[0]*eff_c_cells[0]*omega_nodes[1];

            RHS[0] += rho_cells[0]*(eff_c_cells[0]*T_cells_old[0] - E_cells[0])*omega_nodes[1];
        }
        else
        {
            C[0] += rho_cells[0]*eff_c_cells[1]*omega_nodes[1];

            RHS[0] += rho_cells[0]*(eff_c_cells[1]*T_cells_old[1] - E_cells[1])*omega_nodes[1];
        }

        // middle rows
        for (int i = 1; i < N-1; ++i)
        {
            A[i] = -eff_k_nodes[i];

            B[i] = rho_cells[i]*eff_c_cells[i]*dz_cells_new[i]/this->time_step + eff_k_nodes[i+1] + eff_k_nodes[i];
            
            C[i] =  -eff_k_nodes[i+1];
            
            RHS[i] = rho_cells[i]*eff_c_cells[i]*dz_cells_new[i]*T_cells_old[i]/this->time_step -
                      rho_cells[i]*E_cells[i]*(dz_cells_new[i] - dz_cells_old[i])/this->time_step +
                      (radiation_nodes[i+1] - radiation_nodes[i]);

            if (omega_nodes[i] >= 0)
            {
                A[i] += -rho_cells[i]*eff_c_cells[i-1]*omega_nodes[i];

                RHS[i] += -rho_cells[i]*(eff_c_cells[i-1]*T_cells_old[i-1] - E_cells[i-1])*omega_nodes[i];
            }
            else
            {
                B[i] += -rho_cells[i]*eff_c_cells[i]*omega_nodes[i];

                RHS[i] += -rho_cells[i]*(eff_c_cells[i]*T_cells_old[i] - E_cells[i])*omega_nodes[i];
            }
            
            if (omega_nodes[i+1] >= 0)
            {
                B[i] += rho_cells[i]*eff_c_cells[i]*omega_nodes[i+1];

                RHS[i] += rho_cells[i]*(eff_c_cells[i]*T_cells_old[i] - E_cells[i])*omega_nodes[i+1];
            }
            else
            {
                C[i] += rho_cells[i]*eff_c_cells[i+1]*omega_nodes[i+1];

                RHS[i] += rho_cells[i]*(eff_c_cells[i+1]*T_cells_old[i+1] - E_cells[i+1])*omega_nodes[i+1];
            }
        }

        // last row
        A[N-1] = -eff_k_nodes[N-1];

        B[N-1] = rho_cells[N-1]*eff_c_cells[N-1]*dz_cells_new[N-1]/this->time_step + eff_k_nodes[N] + eff_k_nodes[N-1];

        RHS[N-1] = rho_cells[N-1]*eff_c_cells[N-1]*T_cells_old[N-1]*dz_cells_new[N-1]/this->time_step - 
                    rho_cells[N-1]*E_cells[N-1]*(dz_cells_new[N-1] - dz_cells_old[N-1])/this->time_step -
                    (radiation_nodes[N] - radiation_nodes[N-1]) +
                    eff_k_nodes[N]*T_up_new; 
        
        if (omega_nodes[N-1] >= 0)
        {
            A[N-1] += -rho_cells[N-1]*eff_c_cells[N-2]*omega_nodes[N-1];
            RHS[N-1] += -rho_cells[N-1]*(eff_c_cells[N-2]*T_cells_old[N-2] - E_cells[N-2])*omega_nodes[N-1];
        }
        else
        {
            B[N-1] += -rho_cells[N-1]*eff_c_cells[N-1]*omega_nodes[N-1];
            RHS[N-1] += -rho_cells[N-1]*(eff_c_cells[N-1]*T_cells_old[N-1] - E_cells[N-1])*omega_nodes[N-1];
        }

        if (omega_nodes[N] >= 0)
        {
            B[N-1] += rho_cells[N-1]*eff_c_cells[N-1]*omega_nodes[N];
            RHS[N-1] += rho_cells[N-1]*(eff_c_cells[N-1]*T_cells_old[N-1] - E_cells[N-1])*omega_nodes[N];
        }
        else
        {
            RHS[N-1] += -rho_cells[N-1]*eff_c_up*T_up_new*omega_nodes[N] + rho_cells[N-1]*(eff_c_up*T_up_old - E_up)*omega_nodes[N];
        }

        return {std::vector<NumType>{A.begin()+1, A.end()},
                B,
                std::vector<NumType>{C.begin(), C.end()-1},
                RHS};
    }

    template<typename NumType>
    std::pair<std::vector<NumType>, std::vector<NumType>> ThermoSolver<NumType>::Compute_radiation_nodes(NumType F_sw_value,
                                                                                                         NumType albedo_ice,
                                                                                                         NumType i0_ice,
                                                                                                         NumType kappa_ice,
                                                                                                         const std::vector<NumType>& dz_cells_ice,
                                                                                                         NumType albedo_snow,
                                                                                                         NumType i0_snow,
                                                                                                         NumType kappa_snow,
                                                                                                         const std::vector<NumType>& dz_cells_snow)
    {
        std::vector<NumType> dzi_extended = concatenate(std::vector<NumType>{0.0}, reverse_vec(dz_cells_ice));
        if (dz_cells_snow.size() == 0)
        {
            std::vector<NumType> radiation_nodes_ice = F_sw_value*((NumType)1.0 - albedo_ice)*i0_ice*exp_vec(-kappa_ice*cumsum(dzi_extended));
            return {reverse_vec(radiation_nodes_ice), std::vector<NumType>{0.0}};
        }
        else
        {
            std::vector<NumType> ones_array(dz_cells_ice.size() + 1);
            for (int i = 0; i < dz_cells_ice.size() + 1; ++i)
            {
                ones_array[i] = 1.0;
            }

            std::vector<NumType> dzs_extended = concatenate(std::vector<NumType>{0.0}, reverse_vec(dz_cells_snow));
            std::vector<NumType> radiation_nodes_snow = F_sw_value*((NumType)1.0 - albedo_snow)*i0_snow*exp_vec(-kappa_snow*cumsum(dzs_extended));
            std::vector<NumType> radiation_nodes_ice = F_sw_value*((NumType)1.0 - albedo_snow)*i0_snow*
                                                       exp_vec(-kappa_ice*cumsum(dzi_extended) - ones_array*sum_vec(kappa_snow*dz_cells_snow));
            return {reverse_vec(radiation_nodes_ice), reverse_vec(radiation_nodes_snow)};
        }
    }

    template<typename NumType>
    ThreeVecs<NumType> ThermoSolver<NumType>::sea_ice_freezing_1d(NumType T_ib,
                                                                  const std::vector<NumType>& T_cells,
                                                                  NumType T_is,
                                                                  const std::vector<NumType>& dz_cells,
                                                                  const std::vector<NumType>& salinity_cells,
                                                                  const std::vector<NumType>& rho_cells,
                                                                  bool is_radiation,
                                                                  int max_n_its,
                                                                  NumType tol)
    {
        std::vector<NumType> T_cells_new = T_cells;
        std::vector<NumType> T_cells_prev = T_cells;

        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType omega_ib = (NumType)0.0;

        std::vector<NumType> dz_cells_new = dz_cells;

        std::vector<NumType> radiation_nodes(T_cells.size() + 1);

        std::vector<NumType> T_is_history = {T_is};
        NumType current_err = std::numeric_limits<NumType>::max();
        NumType surface_err = std::numeric_limits<NumType>::max();
        NumType prev_surface_err = surface_err;

        std::vector<NumType> full_temp_vec(T_cells.size() + 1);
        std::vector<NumType> prev_temp_vec(T_cells.size() + 1);
        std::vector<NumType> old_temp_vec = concatenate<NumType>(T_cells, std::vector<NumType>{T_is});

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            ++npseudo;

            T_cells_prev = T_cells_new;
            T_is_prev = T_is_new;
            prev_surface_err = surface_err;

            // compute new value of omega at the base
            omega_ib = this->W_from_BC(T_ib,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       false);
            
            // compute new value of ice surface temperature
            T_is_new = this->T_from_BC(T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       0.0,
                                       true,
                                       true);

            // force the convergence of surface temperature
            surface_err = std::abs(T_is_new - T_is_prev)/std::abs(T_is);
            
            if (surface_err < prev_surface_err)
            {
                T_is_history.push_back(T_is_new);
            }
            else
            {
                T_is_new = sum_vec<NumType>(T_is_history)/T_is_history.size();
            }

            // recalculate ice thickness
            dz_cells_new = this->Update_dz(dz_cells,
                                           omega_ib,
                                           0.0);
            
            // recalculate radiation
            if (is_radiation)
            {
                radiation_nodes = Compute_radiation_nodes(F_sw(T_is_new),
                                                          IceConsts<NumType>::albedo_i,
                                                          IceConsts<NumType>::i0_i,
                                                          IceConsts<NumType>::kappa_i,
                                                          dz_cells_new).first;
            }

            // assemble matrix and rhs for ice 
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_cells_prev, T_cells,
                                                                T_is_new, T_is_prev, T_is,
                                                                T_ib, T_ib, T_ib, 
                                                                omega_ib, 0.0, 
                                                                dz_cells_new, dz_cells,
                                                                salinity_cells,
                                                                rho_cells,
                                                                radiation_nodes,
                                                                true);

            // solve linear system and update temperatures
            T_cells_new = thomas_solver<NumType>(std::get<0>(matrix_rhs),
                                                 std::get<1>(matrix_rhs),
                                                 std::get<2>(matrix_rhs),
                                                 std::get<3>(matrix_rhs));
            
            // evaluate the error
            prev_temp_vec = concatenate<NumType>(T_cells_prev, std::vector<NumType>{T_is_prev});
            
            full_temp_vec = concatenate<NumType>(T_cells_new, std::vector<NumType>{T_is_new}); 

            current_err = L2_norm(full_temp_vec - prev_temp_vec)/L2_norm(old_temp_vec);

            if (current_err < tol)
                break;
        }

        std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        return {T_cells_new, std::vector<NumType>{T_is_new}, dz_cells_new};
    }

    template<typename NumType>
    TwoVecs<NumType> ThermoSolver<NumType>::sea_ice_melting_1d(NumType T_ib,
                                                               NumType T_is,
                                                               const std::vector<NumType>& T_cells,
                                                               NumType T_is_old,
                                                               const std::vector<NumType>& dz_cells,
                                                               const std::vector<NumType>& salinity_cells,
                                                               const std::vector<NumType>& rho_cells,
                                                               bool is_radiation,
                                                               int max_n_its,
                                                               NumType tol)
    {
        std::vector<NumType> T_cells_new = T_cells;
        std::vector<NumType> T_cells_prev = T_cells;

        NumType omega_ib = 0.0;
        NumType omega_is = 0.0;

        std::vector<NumType> dz_cells_new = dz_cells;

        std::vector<NumType> radiation_nodes(T_cells.size() + 1);

        NumType current_err = std::numeric_limits<NumType>::max();

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            ++npseudo;

            T_cells_prev = T_cells_new;
            
            // compute new value of omega at the base
            omega_ib = this->W_from_BC(T_ib,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       false);

            // compute new value of omega at the base
            omega_is = this->W_from_BC(T_is,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       true);

            // recalculate ice thickness
            dz_cells_new = this->Update_dz(dz_cells,
                                           omega_ib,
                                           omega_is);
            
            // recalculate radiation
            if (is_radiation)
            {
                radiation_nodes = Compute_radiation_nodes(F_sw(T_is),
                                                          IceConsts<NumType>::albedo_i,
                                                          IceConsts<NumType>::i0_i,
                                                          IceConsts<NumType>::kappa_i,
                                                          dz_cells_new).first;
            }

            // assemble matrix and rhs for ice 
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_cells_prev, T_cells,
                                                                T_is, T_is, T_is_old,
                                                                T_ib, T_ib, T_ib, 
                                                                omega_ib, omega_is, 
                                                                dz_cells_new, dz_cells,
                                                                salinity_cells,
                                                                rho_cells,
                                                                radiation_nodes,
                                                                true);
            
            // solve linear system and update temperatures
            T_cells_new = thomas_solver<NumType>(std::get<0>(matrix_rhs),
                                                 std::get<1>(matrix_rhs),
                                                 std::get<2>(matrix_rhs),
                                                 std::get<3>(matrix_rhs));

            // evaluate the error
            current_err = L2_norm(T_cells_new - T_cells_prev)/L2_norm(T_cells);

            if (current_err < tol)
                break;
        }

        std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        return {T_cells_new, dz_cells_new};
    }

    // explicit instantiation
    template class ThermoSolver<float>;
    template class ThermoSolver<double>;
}