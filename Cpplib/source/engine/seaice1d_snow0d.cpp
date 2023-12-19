#include "solver.hpp"

namespace icethermo
{
    // 1D ice solver with 0D snow constructor
    template <typename NumType>
    SeaIce1D_Snow0D_Solver<NumType>::SeaIce1D_Snow0D_Solver(Mesh<NumType>* mesh_ice_,
                                                            Mesh<NumType>* mesh_snow_,
                                                            NumType time_step_,
                                                            bool is_radiation_,
                                                            bool is_sublimation_,
                                                            bool is_verbose_,
                                                            Kparam ice_k_param_,
                                                            Cparam ice_c_eff_param_,
                                                            Eparam ice_E_param_,
                                                            Lparam ice_L_param_,
                                                            Aparam ice_albedo_param_,
                                                            Kparam snow_k_param_,
                                                            Lparam snow_L_param_,
                                                            Aparam snow_albedo_param_,
                                                            SnowIceTransition si_transition_mode_):
        SeaIce_Solver<NumType>(mesh_ice_,
                               mesh_snow_,
                               time_step_,
                               ApproxOrder::first,
                               is_radiation_,
                               is_sublimation_,
                               is_verbose_,
                               ice_k_param_,
                               ice_c_eff_param_,
                               ice_E_param_,
                               ice_L_param_,
                               ice_albedo_param_,
                               snow_k_param_,
                               Cparam::FreshSnow,
                               Eparam::FreshSnow,
                               snow_L_param_,
                               snow_albedo_param_)
    {
        // store snow->ice transition mode
        this->si_transition_mode = si_transition_mode_;
        
        // store input mesh
        this->UpdateMesh(mesh_ice_,
                         mesh_snow_);

        // log
        if (this->is_verbose)
        {
            std::cout << "1D sea-ice with 0D snow solver class constructed!" << std::endl;
        }
    } 

    // 1D ice with 0D snow solver evaluation
    template <typename NumType>
    void SeaIce1D_Snow0D_Solver<NumType>::Evaluate()
    {
        // force base temperature to freezing point
        *(this->Ti_b) = GenConsts<NumType>::TempFusion(*(this->So));

        // check if snow exists
        if (sum_vec(*(this->dzs_cells)) < (NumType)SNOW_THICKNESS_THRESHOLD)
        {
            // force snow temperatures to be the same as ice surface temp
            *(this->Ts_b) = *(this->Ti_s);
            *(this->Ts_s) = *(this->Ti_s);
            (*(this->Ts_cells))[0] = *(this->Ti_s);

            // computations without snow

            // recalculate temperatures in ice freezing mode
            auto freezing_values = this->seaice1d_freezing(*(this->Ti_b),
                                                           *(this->Ti_cells),
                                                           *(this->Ti_s),
                                                           *(this->dzi_cells),
                                                           *(this->Si_cells),
                                                           *(this->rhoi_cells));
            
            // compute melting temperature for top layer
            NumType surface_fusion_temp = GenConsts<NumType>::TempFusion((*(this->Si_cells)).back());

            // check if ice surface temperature exeeds melting point and recalculate in melting mode
            if ((std::get<1>(freezing_values))[0] >= surface_fusion_temp)
            {   
                // log mode
                if (this->is_verbose)
                {
                    std::cout << "1D-SEA-ICE MELTING MODE" << std::endl;
                }

                // force ice surface temperature to melting point
                *(this->Ti_s) = surface_fusion_temp; 

                // recompute temperatures in ice melting mode
                auto melting_values = this->seaice1d_melting(*(this->Ti_b),
                                                             surface_fusion_temp,
                                                             *(this->Ti_cells),
                                                             *(this->Ti_s),
                                                             *(this->dzi_cells),
                                                             *(this->Si_cells),
                                                             *(this->rhoi_cells));

                // update mesh values
                *(this->Ti_cells) = std::get<0>(melting_values);
                *(this->dzi_cells) = std::get<1>(melting_values);
            }
            else
            {
                // log mode
                if (this->is_verbose)
                {
                    std::cout << "1D-SEA-ICE FREEZING MODE" << std::endl;
                }

                // update mesh values
                *(this->Ti_cells) = std::get<0>(freezing_values);
                *(this->Ti_s) = (std::get<1>(freezing_values))[0];
                *(this->dzi_cells) = std::get<2>(freezing_values);
            }

            // if there are precipitations with less than zero atm temp add snow
            if (*(this->atm_temp) < (NumType)0.0)
            {
                // update snow thickness according to precipitation rate
                (*(this->dzs_cells))[0] = this->Update_dz_0D((*(this->dzs_cells))[0],
                                                              (NumType)0.0, 
                                                              -(*(this->prec_rate))*
                                                              WaterConsts<NumType>::rho_w/
                                                              SnowConsts<NumType>::rho_s);
                
                // if snow appeared initialize snow temperatures
                if (sum_vec(*(this->dzs_cells)) > (NumType)SNOW_THICKNESS_THRESHOLD)
                {
                    *(this->Ts_s) = *(this->atm_temp);
                    *(this->Ts_b) = *(this->Ti_s);
                    (*(this->Ts_cells))[0] = (NumType)0.5*(*(this->Ts_s) + *(this->Ts_b));
                }
            }
        }
        else
        {
            // computations with snow

            std::cout << (*(this->Ti_cells)).size() << std::endl;

            // recalculate temperatures in ice freezing with snow mode
            auto freezing_values = this->seaice1d_snow0d_freezing(*(this->Ti_b),
                                                                  *(this->Ti_cells),
                                                                  *(this->Ti_s),
                                                                  *(this->dzi_cells),
                                                                  *(this->Si_cells),
                                                                  *(this->rhoi_cells),
                                                                  *(this->Ts_s),
                                                                  (*(this->dzs_cells))[0],
                                                                  (*(this->rhos_cells))[0],
                                                                  *(this->prec_rate),
                                                                  *(this->atm_temp));

            auto freezing_ice = freezing_values.first;
            auto freezing_snow = freezing_values.second;

            // check if snow surface temperature exeeds melting point and recalculate in melting mode
            if ((std::get<0>(freezing_snow))[0] > (NumType)0.0)
            {
                // log mode
                if (this->is_verbose)
                {
                    std::cout << "1D-SEA-ICE WITH 0D-SNOW MELTING MODE" << std::endl;
                }

                // force snow surface temperature to melting point
                *(this->Ts_s) = (NumType)0.0; 

                auto melting_values = this->seaice1d_snow0d_melting(*(this->Ti_b),
                                                                    *(this->Ti_cells),
                                                                    *(this->Ti_s),
                                                                    *(this->dzi_cells),
                                                                    *(this->Si_cells),
                                                                    *(this->rhoi_cells),
                                                                    *(this->Ts_s),
                                                                    (*(this->dzs_cells))[0],
                                                                    (*(this->rhos_cells))[0],
                                                                    *(this->prec_rate),
                                                                    *(this->atm_temp));
                auto melting_ice = melting_values.first;
                auto melting_snow = melting_values.second;

                // update mesh values
                *(this->Ti_cells) = std::get<0>(melting_ice);
                *(this->Ti_s) = (std::get<1>(melting_ice))[0];
                *(this->dzi_cells) = std::get<2>(melting_ice);

                (*(this->dzs_cells))[0] = melting_snow;
                *(this->Ts_b) = *(this->Ti_s);
                (*(this->Ts_cells))[0] = (NumType)0.5*(*(this->Ts_b) + *(this->Ts_s));
            }
            else
            {
                // log mode
                if (this->is_verbose)
                {
                    std::cout << "1D-SEA-ICE WITH 0D-SNOW FREEZING MODE" << std::endl;
                }

                // update mesh values
                *(this->Ti_cells) = std::get<0>(freezing_ice);
                *(this->Ti_s) = std::get<1>(freezing_ice)[0];
                *(this->dzi_cells) = std::get<2>(freezing_ice);

                *(this->Ts_s) = std::get<0>(freezing_snow)[0];
                (*(this->dzs_cells))[0] = std::get<1>(freezing_snow)[0];
                *(this->Ts_b) = *(this->Ti_s);
                (*(this->Ts_cells))[0] = (NumType)0.5*(*(this->Ts_b) + *(this->Ts_s));
            }
        }
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
            else
            {
                this->So =  this->mesh_ice->GetSingleData("ocean_salinity");
            }

            if (!this->mesh_ice->CheckSingleDataExistency("down_temperature"))
            {
                this->Ti_b = this->mesh_ice->CreateSingleData("down_temperature");
                *(this->Ti_b) = GenConsts<NumType>::TempFusion(*(this->So));
            }
            else
            {
                this->Ti_b = this->mesh_ice->GetSingleData("down_temperature");
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
            else
            {
                this->Ts_b = this->mesh_snow->GetSingleData("down_temperature");
            }
        }
    }

    // explicit instantation 
    template class SeaIce1D_Snow0D_Solver<float>;
    template class SeaIce1D_Snow0D_Solver<double>;
}