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
            FuncPtr<NumType> rho = [&](NumType T){return Params::Density(dparam, T, salinity_cells.back());};
            FuncPtr<NumType> k = [&](NumType T){return Params::Conductivity(kparam, T, salinity_cells.back(), rho(T));};
            FuncPtr<NumType> L = [&](NumType T){return Params::FusionHeat(Lparam, T, salinity_cells.back());};

            FuncPtr<NumType> nonlin_func = [&](NumType T){return F(T) - k(T)*grad(T) - rho(T)*L(T)*omega_value;};
            
            // solve nonlinear 1D equation
            return secant_solver<NumType>(nonlin_func, T_cells.back() - 10.0, GenConsts::TempFusion(salinity_cells.back()) + 10.0).first;
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
            FuncPtr<NumType> rho = [&](NumType T){return Params::Density(dparam, T, salinity_cells[0]);};
            FuncPtr<NumType> k = [&](NumType T){return Params::Conductivity(kparam, T, salinity_cells[0], rho(T));};
            FuncPtr<NumType> L = [&](NumType T){return Params::FusionHeat(Lparam, T, salinity_cells[0]);};

            // assemble nonlinear function for 1D solver
            FuncPtr<NumType> nonlin_func = [&](NumType T){return F(T) - k(T)*grad(T) - rho(T)*L(T)*omega_value;};
            
            // solve nonlinear 1D equation
            return secant_solver<NumType>(nonlin_func, T_cells.back() - 10.0, GenConsts::TempFusion(salinity_cells.back()) + 10.0).first;
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
            FuncPtr<NumType> rho = [&](NumType T){return Params::Density(dparam, T, salinity_cells.back());};
            FuncPtr<NumType> k = [&](NumType T){return Params::Conductivity(kparam, T, salinity_cells.back(), rho(T));};
            FuncPtr<NumType> L = [&](NumType T){return Params::FusionHeat(Lparam, T, salinity_cells.back());};

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
            FuncPtr<NumType> rho = [&](NumType T){return Params::Density(dparam, T, salinity_cells[0]);};
            FuncPtr<NumType> k = [&](NumType T){return Params::Conductivity(kparam, T, salinity_cells[0], rho(T));};
            std::function<NumType(NumType, NumType)> L = [&](NumType T, NumType S){return Params::FusionHeat(Lparam, T, S);};

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

    // explicit instantiation
    template class ThermoSolver<float>;
    template class ThermoSolver<double>;
}