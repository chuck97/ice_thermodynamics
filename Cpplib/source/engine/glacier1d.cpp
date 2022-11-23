#include "solver.hpp"

namespace icethermo
{
    // 1D ice solver constructor
    template <typename NumType>
    Glacier1D_Solver<NumType>::Glacier1D_Solver(Mesh<NumType>* mesh_ice_,
                                                NumType time_step_,
                                                ApproxOrder grad_approx_order_,
                                                Kparam ice_k_param_,
                                                Cparam ice_c_eff_param_,
                                                Eparam ice_E_param_,
                                                Lparam ice_L_param_):
        ThermoSolver<NumType>(mesh_ice_,
                              NULL,
                              time_step_,
                              grad_approx_order_,
                              ice_k_param_,
                              ice_c_eff_param_,
                              ice_E_param_,
                              ice_L_param_)
    {
        // store input mesh
        this->UpdateMesh(mesh_ice_);

        // log
        std::cout << "1D glacier solver class constructed!" << std::endl;
    } 

    // 1D ice solver evaluation
    template <typename NumType>
    void Glacier1D_Solver<NumType>::Evaluate()
    {
        // recalculate temperatures in freezing mode
        auto freezing_values = this->glacier1d_freezing(*(this->Ti_b),
                                                        *(this->Ti_cells),
                                                        *(this->Ti_s),
                                                        *(this->dzi_cells),
                                                        std::vector<NumType>((*(this->Ti_cells)).size()),
                                                        *(this->rhoi_cells));
        
        // compute melting temperature for top layer
        NumType surface_fusion_temp = GenConsts<NumType>::TempFusion((NumType)0.0);
        
        // check if surface temperature exeeds melting point and recalculate in melting mode
        if ((std::get<2>(freezing_values))[0] >= surface_fusion_temp)
        {   
            // log mode
            std::cout << "1D-GLACIER MELTING MODE" << std::endl;

            // force surface temperature to melting point
            *(this->Ti_s) = surface_fusion_temp; 

            // recalculate mesh values
            auto melting_values = this->glacier1d_melting(*(this->Ti_b),
                                                         surface_fusion_temp,
                                                         *(this->Ti_cells),
                                                         *(this->Ti_s),
                                                         *(this->dzi_cells),
                                                         std::vector<NumType>((*(this->Ti_cells)).size()),
                                                         *(this->rhoi_cells));
            
            // update mesh values
            *(this->Ti_b) = (std::get<0>(melting_values))[0];
            *(this->Ti_cells) = std::get<1>(melting_values);
            *(this->dzi_cells) = std::get<2>(melting_values);
              
        }
        else
        {
            // log mode
            std::cout << "1D-GLACIER FREEZING MODE" << std::endl;

            // update mesh values
            *(this->Ti_b) = (std::get<0>(freezing_values))[0];
            *(this->Ti_cells) = std::get<1>(freezing_values);
            *(this->Ti_s) = (std::get<2>(freezing_values))[0];
            *(this->dzi_cells) = std::get<3>(freezing_values);
        }
    }

    template <typename NumType>
    void Glacier1D_Solver<NumType>::CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                                         Mesh<NumType>* snow_mesh)
    {
        if (!(ice_mesh->CheckCellsDataExistency("cells_temperature_array") &&
              ice_mesh->CheckCellsDataExistency("cells_density_array") &&
              ice_mesh->CheckSingleDataExistency("up_temperature") &&
              ice_mesh->CheckSingleDataExistency("down_temperature")))
        {
            THERMO_ERR((std::string)"Wrong ice mesh content, the following mesh variables should be created:\n"+
                       (std::string)"cells_temperature_array\n" +
                       (std::string)"cells_density_array\n" +
                       (std::string)"up_temperature\n" + 
                       (std::string)"down_temperature");
        }

        if (snow_mesh != NULL)
        {
            THERMO_ERR("Only ice mesh should be checked in /'Glacier1D_Solver/'");
        }
    }

    template <typename NumType>
    void Glacier1D_Solver<NumType>::UpdateMesh(Mesh<NumType>* mesh_ice_,
                                               Mesh<NumType>* mesh_snow_)
    {

        // check the costistency of ice mesh
        this->CheckMeshConsistency(mesh_ice_);

        if (mesh_ice_ != NULL)
        {
            this->mesh_ice = mesh_ice_;

            this->Ti_cells = this->mesh_ice->GetCellsData("cells_temperature_array");
            this->dzi_cells = this->mesh_ice->GetCellsThickness();
            this->rhoi_cells = this->mesh_ice->GetCellsData("cells_density_array");
            this->Ti_s = this->mesh_ice->GetSingleData("up_temperature");
            this->Ti_b = this->mesh_ice->GetSingleData("down_temperature"); 
        }

        // check consistency of snow mesh
        if (mesh_snow_ != NULL)
        {
            THERMO_ERR("Only ice mesh should be updated in /'Glacier1D_Solver/'");
        }
    }

    // explicit instantation 
    template class Glacier1D_Solver<float>;
    template class Glacier1D_Solver<double>;
}