#include "itinmcm0d.hpp"

using namespace icethermo;

ThermoModelsSet0D::ThermoModelsSet0D(double time_step_,           
                                     double min_ice_thick_,       
                                     int min_lon_ind_,            
                                     int max_lon_ind_,            
                                     int min_lat_ind_,            
                                     int max_lat_ind_,            
                                     double* init_ice_base_temp,  
                                     double* init_ice_surf_temp,  
                                     double* init_snow_surf_temp, 
                                     double* init_ice_thick,      
                                     double* init_snow_thick,     
                                     bool* water_marker,          
                                     bool is_verbose_
                                    ):
        time_step(time_step_),
        min_ice_thick(min_ice_thick_),
        min_lon_ind(min_lon_ind_),
        max_lon_ind(max_lon_ind_),
        min_lat_ind(min_lat_ind_),
        max_lat_ind(max_lat_ind_),
        is_verbose(is_verbose_)
{

    if (min_lon_ind_ > max_lon_ind_)
    {
        THERMO_ERR((std::string)"Thermodynamics constructor error: minimal longitude index is greater than maximal longitude index!");
    }

    if (min_lat_ind_ > max_lat_ind_)
    {
        THERMO_ERR((std::string)"Thermodynamics constructor error: minimal latitude index is greater than maximal latitude index!");
    }

    int lon_dim = max_lon_ind - min_lon_ind + 1;
    int lat_dim = max_lat_ind - min_lat_ind + 1;

    // store is_water array (and allocate do_compute array)
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;
        std::vector<bool> current_lat_is_water;
        std::vector<bool> current_lat_do_compute;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            bool current_is_water = water_marker[local_lat_ind*lon_dim + local_lon_ind];
            current_lat_is_water.push_back(current_is_water);
            current_lat_do_compute.push_back(current_is_water);
        }
        is_water.push_back(current_lat_is_water);
        do_compute.push_back(current_lat_do_compute);
    }

    // create 2d field of 0D meshes
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        std::vector<Mesh<double>*> current_lat_ice_meshes;
        std::vector<Mesh<double>*> current_lat_snow_meshes;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_thick = init_ice_thick[local_lat_ind*lon_dim + local_lon_ind];
            double current_snow_thick = init_snow_thick[local_lat_ind*lon_dim + local_lon_ind];

            current_lat_ice_meshes.push_back
            (
                new(Mesh<double>)(1, current_ice_thick)
            );

            current_lat_snow_meshes.push_back
            (
                new(Mesh<double>)(1, current_snow_thick)
            );
        }
        ice_meshes.push_back(current_lat_ice_meshes);
        snow_meshes.push_back(current_lat_snow_meshes);
    }

    if (is_verbose)
    {
        std::cout << "Meshes are created!\n";
    }

    // initialize 2d field of 0D meshes
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_base_temp = init_ice_base_temp[local_lat_ind*lon_dim + local_lon_ind];
            double current_ice_surf_temp = init_ice_surf_temp[local_lat_ind*lon_dim + local_lon_ind];
            double current_snow_surf_temp = init_snow_surf_temp[local_lat_ind*lon_dim + local_lon_ind];
            
            
            Mesh<double>* local_ice_mesh = ice_meshes[local_lat_ind][local_lon_ind];
            Mesh<double>* local_snow_mesh = snow_meshes[local_lat_ind][local_lon_ind];

            auto initial_ice_temp_cells = local_ice_mesh->CreateCellsData("cells_temperature_array");
            auto initial_ice_sal_cells = local_ice_mesh->CreateCellsData("cells_salinity_array");
            auto initial_ice_dens_cells = local_ice_mesh->CreateCellsData("cells_density_array");
            auto initial_ice_surf_temp = local_ice_mesh->CreateSingleData("up_temperature");
            auto initial_ice_base_temp = local_ice_mesh->CreateSingleData("down_temperature");

            auto initial_snow_temp_cells = local_snow_mesh->CreateCellsData("cells_temperature_array");
            auto initial_snow_dens_cells = local_snow_mesh->CreateCellsData("cells_density_array");
            auto initial_snow_surf_temp = local_snow_mesh->CreateSingleData("up_temperature");
            auto initial_snow_base_temp = local_snow_mesh->CreateSingleData("down_temperature");
            

            // initialize mandatory values
            (*initial_ice_base_temp) = current_ice_base_temp;
            (*initial_ice_temp_cells)[0] = 0.5*(current_ice_base_temp + current_ice_surf_temp);
            (*initial_ice_dens_cells)[0] = IceConsts<double>::rho_i;
            (*initial_ice_sal_cells)[0] = 2.5;
            (*initial_ice_surf_temp) = current_ice_surf_temp;

            (*initial_snow_temp_cells)[0] = 0.5*(current_ice_surf_temp + current_snow_surf_temp);
            (*initial_snow_dens_cells)[0] = SnowConsts<double>::rho_s;
            (*initial_snow_surf_temp) = current_snow_surf_temp;
            (*initial_snow_base_temp) = current_ice_surf_temp;
        }
    }

    if (is_verbose)
    {
        std::cout << "Meshes are initialized!\n";
    }

    // create 2d field of 0D solvers
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        std::vector<SeaIce0D_Snow0D_Solver<double>*> current_lat_solvers;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            Mesh<double>* local_ice_mesh = ice_meshes[local_lat_ind][local_lon_ind];
            Mesh<double>* local_snow_mesh = snow_meshes[local_lat_ind][local_lon_ind];

            current_lat_solvers.push_back
            (
                new(SeaIce0D_Snow0D_Solver<double>)(local_ice_mesh,
                                                    local_snow_mesh,
                                                    time_step,
                                                    false,
                                                    false)
            );
        }
        solvers.push_back(current_lat_solvers);
    }

    if (is_verbose)
    {
        std::cout << "Solvers are initialized!\n";
    }
}

void ThermoModelsSet0D::SetComputationMarker(bool* do_compute_,
                                             int min_lon_ind_,        
                                             int max_lon_ind_,        
                                             int min_lat_ind_,        
                                             int max_lat_ind_
                                             )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"SetComputationMarker error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"SetComputationMarker error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"SetComputationMarker error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"SetComputationMarker error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            do_compute[local_lat_ind][local_lon_ind] = (do_compute_[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] and
                                                       is_water[local_lat_ind][local_lon_ind]);
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d computational array!\n";
    }
}

void ThermoModelsSet0D::UpdateSwRadiation(double* sw_values,
                                          int min_lon_ind_,        
                                          int max_lon_ind_,        
                                          int min_lat_ind_,        
                                          int max_lat_ind_
                                          )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateSwRadiation error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateSwRadiation error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateSwRadiation error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateSwRadiation error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_sw_value = sw_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateShortWaveRadiation
            (
               [current_sw_value](double T){return current_sw_value;}
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d short-wave radiation array!\n";
    }
}

void ThermoModelsSet0D::UpdateLwRadiation(double* lw_values,
                                          int min_lon_ind_,        
                                          int max_lon_ind_,        
                                          int min_lat_ind_,        
                                          int max_lat_ind_
                                          )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateLwRadiation error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateLwRadiation error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateLwRadiation error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateLwRadiation error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_lw_value = lw_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateLongWaveRadiation
            (
               [current_lw_value](double T){return current_lw_value;}
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d long-wave radiation array!\n";
    }
}

void ThermoModelsSet0D::UpdateAirTemperature(double* ta_values,
                                                int min_lon_ind_,        
                                                int max_lon_ind_,        
                                                int min_lat_ind_,        
                                                int max_lat_ind_
                                                )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirTemperature error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirTemperature error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirTemperature error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirTemperature error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ta_value = ta_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateAtmosphereTemperature
            (
               current_ta_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d atmosphere temperature array!\n";
    }
}

void ThermoModelsSet0D::UpdatePrecipitationRate(double* pr_values,
                                                int min_lon_ind_,        
                                                int max_lon_ind_,        
                                                int min_lat_ind_,        
                                                int max_lat_ind_
                                                )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdatePrecipitationRate error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdatePrecipitationRate error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdatePrecipitationRate error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdatePrecipitationRate error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_pr_value = pr_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdatePrecipitationRate
            (
               current_pr_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d liquid precipitation rate array!\n";
    }
}

void ThermoModelsSet0D::UpdateAirPressure(double* ap_values,
                                          int min_lon_ind_,        
                                          int max_lon_ind_,        
                                          int min_lat_ind_,        
                                          int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirPressure error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirPressure error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirPressure error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirPressure error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ap_value = ap_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateAtmospherePressure
            (
               current_ap_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d atm pressure array!\n";
    }
}

void ThermoModelsSet0D::UpdateAirSpecificHumidity(double* sh_values,
                                                  int min_lon_ind_,        
                                                  int max_lon_ind_,        
                                                  int min_lat_ind_,        
                                                  int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirSpecificHumidity error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirSpecificHumidity error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirSpecificHumidity error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirSpecificHumidity error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_sh_value = sh_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateAirSpecificHumidity
            (
               current_sh_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d specific humidity array!\n";
    }
}

void ThermoModelsSet0D::UpdateAbsWindSpeed(double* ws_values,
                                           int min_lon_ind_,        
                                           int max_lon_ind_,        
                                           int min_lat_ind_,        
                                           int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAbsWindSpeed error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAbsWindSpeed error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAbsWindSpeed error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAbsWindSpeed error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ws_value = ws_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateAbsWindSpeed
            (
               current_ws_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d wind speed array!\n";
    }
}

void ThermoModelsSet0D::UpdateShCoeff(double* shc_values,
                                      int min_lon_ind_,        
                                      int max_lon_ind_,        
                                      int min_lat_ind_,        
                                      int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateShCoeff error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateShCoeff error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateShCoeff error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateShCoeff error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_shc_value = shc_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateShTransCoeff
            (
               current_shc_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d sensible heat coefficients array!\n";
    }
}

void ThermoModelsSet0D::UpdateLhCoeff(double* lhc_values,
                                      int min_lon_ind_,        
                                      int max_lon_ind_,        
                                      int min_lat_ind_,        
                                      int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateLhCoeff error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateLhCoeff error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateLhCoeff error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateLhCoeff error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_lhc_value = lhc_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateLhTransCoeff
            (
               current_lhc_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d latent heat coefficients array!\n";
    }
}

void ThermoModelsSet0D::AssembleTotalAtmFlux(int min_lon_ind_,        
                                             int max_lon_ind_,        
                                             int min_lat_ind_,        
                                             int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"AssembleTotalAtmFlux error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"AssembleTotalAtmFlux error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"AssembleTotalAtmFlux error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"AssembleTotalAtmFlux error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            solvers[local_lat_ind][local_lon_ind]->UpdateUpperFlux();
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d total atmosphere flux!\n";
    }
}

void ThermoModelsSet0D::UpdateOceanSalinity(double* os_values,
                                            int min_lon_ind_,        
                                            int max_lon_ind_,        
                                            int min_lat_ind_,        
                                            int max_lat_ind_
                                            )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanSalinity error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanSalinity error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanSalinity error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanSalinity error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_os_value = os_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateOceanSalinity
            (
               current_os_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d ocean salinity array!\n";
    }
}

void ThermoModelsSet0D::UpdateOceanFlux(double* of_values,
                                        int min_lon_ind_,        
                                        int max_lon_ind_,        
                                        int min_lat_ind_,        
                                        int max_lat_ind_
                                        )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanFlux error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanFlux error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanFlux error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanFlux error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_of_value = of_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[local_lat_ind][local_lon_ind]->UpdateLowerFlux
            (
               [current_of_value](double T){return current_of_value;}
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d ocean flux array!\n";
    }
}

void ThermoModelsSet0D::UpdateIceThickness(double* thick_values,
                                           int min_lon_ind_,        
                                           int max_lon_ind_,        
                                           int min_lat_ind_,        
                                           int max_lat_ind_
                                           )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateIceThickness error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateIceThickness error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateIceThickness error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateIceThickness error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_thickness = thick_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            (*ice_meshes[local_lat_ind][local_lon_ind]->GetCellsThickness())[0] = current_ice_thickness;
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d ice thickness array!\n";
    }
}

void ThermoModelsSet0D::UpdateSnowThickness(double* thick_values,
                                            int min_lon_ind_,        
                                            int max_lon_ind_,        
                                            int min_lat_ind_,        
                                            int max_lat_ind_
                                            )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateSnowThickness error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateSnowThickness error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateSnowThickness error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateSnowThickness error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_snow_thickness = thick_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            (*snow_meshes[local_lat_ind][local_lon_ind]->GetCellsThickness())[0] = current_snow_thickness;
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d snow thickness array!\n";
    }
}


void ThermoModelsSet0D::Evaluate(int min_lon_ind_,        
                                 int max_lon_ind_,        
                                 int min_lat_ind_,        
                                 int max_lat_ind_
                                 )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"Evaluate error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"Evaluate error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"Evaluate error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"Evaluate error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    // to do
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_thick = ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();

            if (is_water[local_lat_ind][local_lon_ind] and do_compute[local_lat_ind][local_lon_ind] and (current_ice_thick > min_ice_thick))
            {
                solvers[local_lat_ind][local_lon_ind]->Evaluate();
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Single step evaluation of thermodynamics solver!\n";
    }
}


void ThermoModelsSet0D::GetSurfaceTemperature(double* array,
                                              int min_lon_ind_,        
                                              int max_lon_ind_,        
                                              int min_lat_ind_,        
                                              int max_lat_ind_
                                              )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceTemperature error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceTemperature error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceTemperature error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceTemperature error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double temp_value;

            if (is_water[local_lat_ind][local_lon_ind] and 
                (ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness() > min_ice_thick)
               )
            {
                if (snow_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD)
                {
                    temp_value = *(snow_meshes[local_lat_ind][local_lon_ind]->GetSingleData("up_temperature"));
                }
                else
                {
                    temp_value = *(ice_meshes[local_lat_ind][local_lon_ind]->GetSingleData("up_temperature"));
                }
            }
            else
            {
                temp_value = NAN_TEMP_VALUE;
            }

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = temp_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d surface temperature array!\n";
    }
}

void ThermoModelsSet0D::GetIceThickness(double* array,
                                        int min_lon_ind_,        
                                        int max_lon_ind_,        
                                        int min_lat_ind_,        
                                        int max_lat_ind_
                                        )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetIceThickness error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetIceThickness error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetIceThickness error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetIceThickness error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double thick_value;

            if (is_water[local_lat_ind][local_lon_ind] and 
                (ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness() > min_ice_thick))
            {
                thick_value = ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();
            }
            else
            {
                thick_value = NAN_THICK_VALUE;
            }

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = thick_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d ice thickness array!\n";
    }
}

void ThermoModelsSet0D::GetSnowThickness(double* array,
                                         int min_lon_ind_,        
                                         int max_lon_ind_,        
                                         int min_lat_ind_,        
                                         int max_lat_ind_
                                         )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetSnowThickness error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetSnowThickness error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetSnowThickness error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetSnowThickness error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double thick_value;

            if (is_water[local_lat_ind][local_lon_ind] and 
                (snow_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD))
            {
                thick_value = snow_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();
            }
            else
            {
                thick_value = NAN_THICK_VALUE;
            }

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = thick_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d snow thickness array!\n";
    }
}

void ThermoModelsSet0D::GetIsIce(bool* array,
                                 int min_lon_ind_,        
                                 int max_lon_ind_,        
                                 int min_lat_ind_,        
                                 int max_lat_ind_
                                 )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetIsIce error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetIsIce error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetIsIce error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetIsIce error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            bool is_ice_value;

            if (is_water[local_lat_ind][local_lon_ind] and 
                (ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness() > min_ice_thick))
            {
                is_ice_value = true;
            }
            else
            {
                is_ice_value = false;
            }

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = is_ice_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d is_ice boolean array!\n";
    }
}

void ThermoModelsSet0D::GetIsSnow(bool* array,
                                  int min_lon_ind_,        
                                  int max_lon_ind_,        
                                  int min_lat_ind_,        
                                  int max_lat_ind_
                                 )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetIsSnow error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetIsSnow error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetIsSnow error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetIsSnow error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            bool is_snow_value;

            if (is_water[local_lat_ind][local_lon_ind] and 
                (snow_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD))
            {
                is_snow_value = true;
            }
            else
            {
                is_snow_value = false;
            }

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = is_snow_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d is_snow boolean array!\n";
    }
}


////////////////////////////////////////////////////
//                                                //
//             Relization for interfaces          //
//                                                //
////////////////////////////////////////////////////


void* InitThermodynamics(double time_step,
                         double min_ice_thick,
                         int min_lon_ind,
                         int max_lon_ind,
                         int min_lat_ind,
                         int max_lat_ind,
                         double* init_ice_base_temp,
                         double* init_ice_surf_temp,
                         double* init_snow_surf_temp,
                         double* init_ice_thick,
                         double* init_snow_thick,
                         bool* water_marker,
                         bool is_verbose
                         )
{
    ThermoModelsSet0D* ptr = new ThermoModelsSet0D(time_step,
                                                   min_ice_thick,
                                                   min_lon_ind,
                                                   max_lon_ind,
                                                   min_lat_ind,
                                                   max_lat_ind,
                                                   init_ice_base_temp,
                                                   init_ice_surf_temp,
                                                   init_snow_surf_temp,
                                                   init_ice_thick,
                                                   init_snow_thick,
                                                   water_marker,
                                                   is_verbose);
    
    return (void*)ptr;
}

void FinalizeThermodynamics(void* obj)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    delete ptr;
}

void SetComputationMarker(void* obj,
                          bool* do_compute,
                          int min_lon_ind,        
                          int max_lon_ind,        
                          int min_lat_ind,        
                          int max_lat_ind
                          )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->SetComputationMarker(do_compute,
                              min_lon_ind,        
                              max_lon_ind,        
                              min_lat_ind,        
                              max_lat_ind);
}

void UpdateSwRadiation(void* obj,
                       double* sw_values,
                       int min_lon_ind,        
                       int max_lon_ind,        
                       int min_lat_ind,        
                       int max_lat_ind
                       )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateSwRadiation(sw_values,
                           min_lon_ind,
                           max_lon_ind,
                           min_lat_ind,
                           max_lat_ind);
}

void UpdateLwRadiation(void* obj,
                       double* lw_values,
                       int min_lon_ind,        
                       int max_lon_ind,        
                       int min_lat_ind,        
                       int max_lat_ind
                       )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateLwRadiation(lw_values,
                           min_lon_ind,
                           max_lon_ind,
                           min_lat_ind,
                           max_lat_ind);
}

void UpdateAirTemperature(void* obj,
                          double* ta_values,
                          int min_lon_ind,        
                          int max_lon_ind,        
                          int min_lat_ind,        
                          int max_lat_ind
                          )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateAirTemperature(ta_values,
                              min_lon_ind,
                              max_lon_ind,
                              min_lat_ind,
                              max_lat_ind);
}

void UpdatePrecipitationRate(void* obj,
                             double* pr_values,
                             int min_lon_ind,        
                             int max_lon_ind,        
                             int min_lat_ind,        
                             int max_lat_ind
                             )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdatePrecipitationRate(pr_values,
                                 min_lon_ind,
                                 max_lon_ind,
                                 min_lat_ind,
                                 max_lat_ind);
}

void UpdateAirPressure(void* obj,
                       double* ap_values,
                       int min_lon_ind,        
                       int max_lon_ind,        
                       int min_lat_ind,        
                       int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateAirPressure(ap_values,
                           min_lon_ind,
                           max_lon_ind,
                           min_lat_ind,
                           max_lat_ind);
}

void UpdateAirSpecificHumidity(void* obj,
                               double* sh_values,
                               int min_lon_ind,        
                               int max_lon_ind,        
                               int min_lat_ind,        
                               int max_lat_ind
                               )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateAirSpecificHumidity(sh_values,
                                   min_lon_ind,
                                   max_lon_ind,
                                   min_lat_ind,
                                   max_lat_ind);
}

void UpdateAbsWindSpeed(void* obj,
                        double* ws_values,
                        int min_lon_ind,        
                        int max_lon_ind,        
                        int min_lat_ind,        
                        int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateAbsWindSpeed(ws_values,
                            min_lon_ind,
                            max_lon_ind,
                            min_lat_ind,
                            max_lat_ind);
}

void UpdateShCoeff(void* obj,
                   double* shc_values,
                   int min_lon_ind,        
                   int max_lon_ind,        
                   int min_lat_ind,        
                   int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateShCoeff(shc_values,
                       min_lon_ind,
                       max_lon_ind,
                       min_lat_ind,
                       max_lat_ind);
}

void UpdateLhCoeff(void* obj,
                   double* lhc_values,
                   int min_lon_ind,        
                   int max_lon_ind,        
                   int min_lat_ind,        
                   int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateLhCoeff(lhc_values,
                       min_lon_ind,
                       max_lon_ind,
                       min_lat_ind,
                       max_lat_ind);
}

void AssembleTotalAtmFlux(void* obj,
                          int min_lon_ind,        
                          int max_lon_ind,        
                          int min_lat_ind,        
                          int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->AssembleTotalAtmFlux(min_lon_ind,
                              max_lon_ind,
                              min_lat_ind,
                              max_lat_ind);
}

void UpdateOceanSalinity(void* obj,
                         double* os_values,
                         int min_lon_ind,        
                         int max_lon_ind,        
                         int min_lat_ind,        
                         int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateOceanSalinity(os_values,
                             min_lon_ind,
                             max_lon_ind,
                             min_lat_ind,
                             max_lat_ind);
}

void UpdateOceanFlux(void* obj,
                     double* of_values,
                     int min_lon_ind,        
                     int max_lon_ind,        
                     int min_lat_ind,        
                     int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateOceanFlux(of_values,
                         min_lon_ind,
                         max_lon_ind,
                         min_lat_ind,
                         max_lat_ind);
}

void UpdateIceThickness(void* obj,
                        double* thick_values,
                        int min_lon_ind,        
                        int max_lon_ind,        
                        int min_lat_ind,        
                        int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateIceThickness(thick_values,
                            min_lon_ind,
                            max_lon_ind,
                            min_lat_ind,
                            max_lat_ind);
}

void UpdateSnowThickness(void* obj,
                         double* thick_values,
                         int min_lon_ind,        
                         int max_lon_ind,        
                         int min_lat_ind,        
                         int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->UpdateSnowThickness(thick_values,
                             min_lon_ind,
                             max_lon_ind,
                             min_lat_ind,
                             max_lat_ind);
}

void Evaluate(void* obj,
              int min_lon_ind,        
              int max_lon_ind,        
              int min_lat_ind,        
              int max_lat_ind)
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->Evaluate(min_lon_ind,
                  max_lon_ind,
                  min_lat_ind,
                  max_lat_ind);
}

void GetSurfaceTemperature(void* obj,
                           double* array,
                           int min_lon_ind,        
                           int max_lon_ind,        
                           int min_lat_ind,        
                           int max_lat_ind
                           )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->GetSurfaceTemperature(array,
                               min_lon_ind,        
                               max_lon_ind,        
                               min_lat_ind,        
                               max_lat_ind);
}

void GetIceThickness(void* obj,
                     double* array,
                     int min_lon_ind,        
                     int max_lon_ind,        
                     int min_lat_ind,        
                     int max_lat_ind
                     )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->GetIceThickness(array,
                         min_lon_ind,        
                         max_lon_ind,        
                         min_lat_ind,        
                         max_lat_ind);
}

void GetSnowThickness(void* obj,
                      double* array,
                      int min_lon_ind,        
                      int max_lon_ind,        
                      int min_lat_ind,        
                      int max_lat_ind
                      )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->GetSnowThickness(array,
                          min_lon_ind,        
                          max_lon_ind,        
                          min_lat_ind,        
                          max_lat_ind);
}

void GetIsSnow(void* obj,
               bool* array,
               int min_lon_ind,        
               int max_lon_ind,        
               int min_lat_ind,        
               int max_lat_ind
             )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->GetIsSnow(array,
                   min_lon_ind,        
                   max_lon_ind,        
                   min_lat_ind,        
                   max_lat_ind);
}

void GetIsIce(void* obj,
              bool* array,
              int min_lon_ind,        
              int max_lon_ind,        
              int min_lat_ind,        
              int max_lat_ind
              )
{
    ThermoModelsSet0D* ptr = (ThermoModelsSet0D*)obj;
    ptr->GetIsIce(array,
                  min_lon_ind,        
                  max_lon_ind,        
                  min_lat_ind,        
                  max_lat_ind);
}