#include "solver.hpp"

namespace icethermo
{
    // 1D ice solver with 0D snow constructor
    template <typename NumType>
    SeaIce1D_Snow0D_Solver<NumType>::SeaIce1D_Snow0D_Solver(Mesh<NumType>* mesh_ice_,
                                                            Mesh<NumType>* mesh_snow_,
                                                            NumType time_step_,
                                                            Kparam ice_k_param_,
                                                            Cparam ice_c_eff_param_,
                                                            Eparam ice_E_param_,
                                                            Lparam ice_L_param_,
                                                            Kparam snow_k_param_,
                                                            Lparam snow_L_param_):
        ThermoSolver<NumType>(mesh_ice_,
                              mesh_snow_,
                              time_step_,
                              ApproxOrder::first,
                              ice_k_param_,
                              ice_c_eff_param_,
                              ice_E_param_,
                              ice_L_param_,
                              snow_k_param_,
                              Cparam::FreshSnow,
                              Eparam::FreshSnow,
                              snow_L_param_)
    {
        // store input mesh
        this->UpdateMesh(mesh_ice_,
                         mesh_snow_);

        // log
        std::cout << "1D sea-ice with 0D snow solver class constructed!" << std::endl;
    } 


    // 1D ice with 0D snow solver evaluation
    template <typename NumType>
    void SeaIce1D_Snow0D_Solver<NumType>::Evaluate()
    {
        // force base temperature to freezing point
        *(this->Ti_b) = GenConsts<NumType>::TempFusion(*(this->So));

        /*
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        TODO
        */
       
        // recalculate temperatures in freezing mode
        auto seaice1d_snow0d_res = this->seaice1d_snow0d_freezing(*(this->Ti_b),
                                                                  *(this->Ti_cells),
                                                                  *(this->Ti_s),
                                                                  *(this->dzi_cells),
                                                                  *(this->Si_cells),
                                                                  *(this->rhoi_cells),
                                                                  *(this->Ts_s),
                                                                  (*(this->dzs_cells))[0],
                                                                  (*(this->rhos_cells))[0],
                                                                  this->prec_rate,
                                                                  this->atm_temp);

        auto freezing_ice = seaice1d_snow0d_res.first;
        auto freezing_snow = seaice1d_snow0d_res.second;

        // log mode
        std::cout << "ICE FREEZING WITH SNOW MODE" << std::endl;

        // update mesh values
        *(this->Ti_cells) = std::get<0>(freezing_ice);
        *(this->Ti_s) = std::get<1>(freezing_ice)[0];
        *(this->dzi_cells) = std::get<2>(freezing_ice);

        *(this->Ts_s) = std::get<0>(freezing_snow)[0];
        *(this->dzs_cells) = std::get<1>(freezing_snow);
        *(this->Ts_b) = *(this->Ti_s);
        (*(this->Ts_cells))[0] = (NumType)0.5*(*(this->Ts_s) + *(this->Ts_b));
        
    }

    template <typename NumType>
    void SeaIce1D_Snow0D_Solver<NumType>::CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                                               Mesh<NumType>* snow_mesh)
    {
        if (!(ice_mesh->CheckCellsDataExistency("cells_temperature_array") &&
              ice_mesh->CheckCellsDataExistency("cells_salinity_array") &&
              ice_mesh->CheckCellsDataExistency("cells_density_array") &&
              ice_mesh->CheckSingleDataExistency("up_temperature")))
        {
            THERMO_ERR((std::string)"Wrong ice mesh content, the following mesh variables should be created:\n"+
                       (std::string)"cells_temperature_array\n" +
                       (std::string)"cells_salinity_array\n" +
                       (std::string)"cells_density_array\n" +
                       (std::string)"up_temperature");
        }

        if (!(snow_mesh->CheckCellsDataExistency("cells_temperature_array") &&
              snow_mesh->CheckCellsDataExistency("cells_density_array") &&
              snow_mesh->CheckSingleDataExistency("up_temperature")))
        {
            THERMO_ERR((std::string)"Wrong snow mesh content, the following mesh variables should be created:\n"+
                       (std::string)"cells_temperature_array\n" +
                       (std::string)"cells_density_array\n" +
                       (std::string)"up_temperature");
        }
    }

    template <typename NumType>
    void SeaIce1D_Snow0D_Solver<NumType>::UpdateMesh(Mesh<NumType>* mesh_ice_,
                                                     Mesh<NumType>* mesh_snow_)
    {

        // check the costistency of ice mesh
        this->CheckMeshConsistency(mesh_ice_,
                                   mesh_snow_);

        // assign initial ice values
        if (mesh_ice_ != NULL)
        {
            this->mesh_ice = mesh_ice_;

            this->Ti_cells = this->mesh_ice->GetCellsData("cells_temperature_array");
            this->dzi_cells = this->mesh_ice->GetCellsThickness();
            this->Si_cells = this->mesh_ice->GetCellsData("cells_salinity_array");
            this->rhoi_cells = this->mesh_ice->GetCellsData("cells_density_array");
            this->Ti_s = this->mesh_ice->GetSingleData("up_temperature");

            if (!this->mesh_ice->CheckSingleDataExistency("ocean_salinity"))
            {
                this->So = this->mesh_ice->CreateSingleData("ocean_salinity", false);
                *(this->So) = (NumType)30.0;
            }

            if (!this->mesh_ice->CheckSingleDataExistency("down_temperature"))
            {
                this->Ti_b = this->mesh_ice->CreateSingleData("down_temperature");
                *(this->Ti_b) = GenConsts<NumType>::TempFusion(*(this->So));
            }
        }

        // assign initial snow values
        if (mesh_snow_ != NULL)
        {
            this->mesh_snow = mesh_snow_;

            this->Ts_cells = this->mesh_snow->GetCellsData("cells_temperature_array");
            this->dzs_cells = this->mesh_snow->GetCellsThickness();
            this->rhos_cells = this->mesh_snow->GetCellsData("cells_density_array");
            this->Ts_s = this->mesh_snow->GetSingleData("up_temperature");

            if (!this->mesh_snow->CheckSingleDataExistency("down_temperature"))
            {
                this->Ts_b = this->mesh_snow->CreateSingleData("down_temperature");
            }
        }
    }

    // explicit instantation 
    template class SeaIce1D_Snow0D_Solver<float>;
    template class SeaIce1D_Snow0D_Solver<double>;
}