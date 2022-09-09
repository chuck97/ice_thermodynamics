#include "solver.hpp"

namespace icethermo
{
    template <typename NumType>
    ThermoSolver<NumType>::ThermoSolver(Mesh<NumType>* mesh_ice_,
                                        Mesh<NumType>* mesh_snow_
                                        )
    {
        this->mesh_ice = mesh_ice_;
        this->mesh_snow = mesh_snow_;
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
                                             const std::vector<NumType> salinity_cells,
                                             NumType k_value,
                                             NumType omega_value,
                                             FuncPtr<NumType> F,
                                             bool is_surface,
                                             ApproxOrder grad_approx_order,
                                             NumType rho)
    {
        // аппроксимация градиента на границе
        if (is_surface)
        {
            /*
            h1 = dz_cells[-1]/2.0
            h2 = dz_cells[-1] + dz_cells[-2]/2.0
            # возможно, тут стоит учесть солёность
            grad_su = lambda T: (T*(h2**2 - h1**2) - T_cells[-1]*h2**2 + T_cells[-2]*h1**2)/(h1*h2*(h2-h1))
            nonlin_func = lambda T: rho*L(T, salinity_cells[-1])*omega + F(T) - k*(grad_su(T))
            return secant_solver(nonlin_func, T_cells[-1]-10.0, Tf_i(salinity_cells[-1])+10.0)[0]
            */
        }
        else
        {
            /*
            h1 = dz_cells[0]/2.0
            h2 = dz_cells[0] + dz_cells[1]/2.0
            grad_su = lambda T: (-T*(h2**2 - h1**2) + T_cells[0]*h2**2 - T_cells[1]*h1**2)/(h1*h2*(h2-h1))
            nonlin_func = lambda T: rho_i*L(T, salinity_cells[0])*omega + F(T) - k*(grad_su(T))
            return secant_solver(nonlin_func, T_cells[0]-10.0, Tf_i(salinity_cells[0])+10.0)[0] 
            */
        }   
        return (NumType)0;
    }

    // explicit instantiation
    template class ThermoSolver<float>;
    template class ThermoSolver<double>;
}