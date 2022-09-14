#include "solver.hpp"

namespace icethermo
{
    // 1D ice solver constructor
    template <typename NumType>
    Ice1D_Solver<NumType>::Ice1D_Solver(Mesh<NumType>* mesh_ice_,
                                        Dparam ice_rho_param_,
                                        Kparam ice_k_param_,
                                        Cparam ice_c_eff_param_,
                                        Eparam ice_E_param_,
                                        Lparam ice_L_param_):
        ThermoSolver<NumType>(mesh_ice_,
                              NULL,
                              ice_rho_param_,
                              ice_k_param_,
                              ice_c_eff_param_,
                              ice_E_param_,
                              ice_L_param_)
        
    {
        // get main prognostic variables
        this->Ti_cells = this->mesh_ice->GetCellData("Ice cells temperature");
        this->dzi_cells = this->mesh_ice->GetCellsThickness();
        this->Si_cells = this->mesh_ice->GetCellData("Ice cells salinity");
        this->Ti_s = this->mesh_ice->GetSingleData("Ice surface temperature");
        this->Ti_b = this->mesh_ice->GetSingleData("Ice base temperature");

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