#include "solver.hpp"

namespace icethermo
{
    // virtual base class constructor
    template <typename NumType>
    ThermoSolver<NumType>::ThermoSolver(Mesh<NumType>* mesh_ice_,
                                        Mesh<NumType>* mesh_snow_
                                        )
    {
        this->mesh_ice = mesh_ice_;
        this->mesh_snow = mesh_snow_;
    }

    // 1D ice solver constructor
    template <typename NumType>
    Ice1D_Solver<NumType>::Ice1D_Solver(Mesh<NumType>* mesh_ice_,
                                        Dparam ice_rho_param_,
                                        Kparam ice_k_param_,
                                        Cparam ice_c_eff_param_,
                                        Eparam ice_E_param_,
                                        Lparam ice_L_param_):
        ThermoSolver<NumType>(mesh_ice_),
        ice_rho_param(ice_rho_param_),
        ice_k_param(ice_k_param_),
        ice_c_eff_param(ice_c_eff_param_),
        ice_E_param(ice_E_param_),
        ice_L_param(ice_L_param_)
    {
        // get main prognostic variables
        Ti_cells = this->mesh_ice->GetCellData("Ice cells temperature");
        dzi_cells = this->mesh_ice->GetCellsThickness();
        Si_cells = this->mesh_ice->GetCellData("Ice cells salinity");
        Ti_s = this->mesh_ice->GetSingleData("Ice surface temperature");
        Ti_b = this->mesh_ice->GetSingleData("Ice base temperature");

        // log
        std::cout << "1D ice solver class constructed!" << std::endl;
    } 

    // 1D ice solver evaluation
    template <typename NumType>
    void Ice1D_Solver<NumType>::Evaluate()
    {
        // to do
        std::cout << "1D ice solver evaluation!" << std::endl;
    }

    // explicit instantation 
    template class Ice1D_Solver<float>;
    template class Ice1D_Solver<double>;
}