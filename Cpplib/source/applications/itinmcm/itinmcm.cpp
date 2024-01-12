#include "itinmcm.hpp"

using namespace icethermo;

ThermoModelsSet::ThermoModelsSet(double time_step_,  
                                 int num_ice_layers_,         
                                 double min_ice_thick_, 
                                 double min_snow_thick_,        
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
        num_ice_layers(num_ice_layers_),
        min_ice_thick(min_ice_thick_),
        min_snow_thick(min_snow_thick_),
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
}

ThermoModelsSetIce0dSnow0d::ThermoModelsSetIce0dSnow0d(double time_step_,         
                                                       double min_ice_thick_,
                                                       double min_snow_thick_,       
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
        ThermoModelsSet::ThermoModelsSet(time_step_,  
                                         1,         
                                         min_ice_thick_, 
                                         min_snow_thick_,      
                                         min_lon_ind_,            
                                         max_lon_ind_,            
                                         min_lat_ind_,            
                                         max_lat_ind_,            
                                         init_ice_base_temp,  
                                         init_ice_surf_temp,  
                                         init_snow_surf_temp, 
                                         init_ice_thick,      
                                         init_snow_thick,     
                                         water_marker,          
                                         is_verbose_)
{
    int lon_dim = max_lon_ind - min_lon_ind + 1;
    int lat_dim = max_lat_ind - min_lat_ind + 1;

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
            (*initial_ice_sal_cells)[0] = 0.0;
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

        std::vector<icethermo::SeaIce_Solver<double>*> current_lat_solvers;

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
                                                    min_ice_thick,
                                                    min_snow_thick,
                                                    false,
                                                    is_verbose,
                                                    Kparam::FreshIce,
                                                    Lparam::FreshIce,
                                                    Aparam::MeltingFreezingIce,
                                                    Kparam::FreshSnow,
                                                    Lparam::FreshSnow,
                                                    Aparam::MeltingFreezingSnow
                                                    )
            );
        }
        solvers.push_back(current_lat_solvers);
    }

    if (is_verbose)
    {
        std::cout << "Solvers are initialized!\n";
    }
}

ThermoModelsSetIce1dSnow0d::ThermoModelsSetIce1dSnow0d(double time_step_,
                                                       int num_ice_layers_,
                                                       double min_ice_thick_,
                                                       double min_snow_thick_,
                                                       int min_lon_ind_,
                                                       int max_lon_ind_,
                                                       int min_lat_ind_,
                                                       int max_lat_ind_,
                                                       double* init_ice_base_temp,
                                                       double* init_ice_surf_temp,
                                                       double* init_snow_surf_temp,
                                                       double* init_ice_thick,
                                                       double* init_snow_thick,
                                                       double* init_ice_base_salinity,
                                                       double* init_ice_surface_salinity,
                                                       bool* water_marker,
                                                       bool is_verbose_
                                                       ):
        ThermoModelsSet::ThermoModelsSet(time_step_,  
                                         num_ice_layers_,         
                                         min_ice_thick_,
                                         min_snow_thick_,       
                                         min_lon_ind_,            
                                         max_lon_ind_,            
                                         min_lat_ind_,            
                                         max_lat_ind_,            
                                         init_ice_base_temp,  
                                         init_ice_surf_temp,  
                                         init_snow_surf_temp, 
                                         init_ice_thick,      
                                         init_snow_thick,     
                                         water_marker,          
                                         is_verbose_)
{
    int lon_dim = max_lon_ind - min_lon_ind + 1;
    int lat_dim = max_lat_ind - min_lat_ind + 1;

    // create 2d field of 0D meshes for snow and 1d meshes for ice
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
                new(Mesh<double>)(num_ice_layers, current_ice_thick)
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

    // initialize 2d field of 0D meshes for snow and 1d meshes for ice
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_base_temp = init_ice_base_temp[local_lat_ind*lon_dim + local_lon_ind];
            double current_ice_surf_temp = init_ice_surf_temp[local_lat_ind*lon_dim + local_lon_ind];
            double current_snow_surf_temp = init_snow_surf_temp[local_lat_ind*lon_dim + local_lon_ind];
            double current_ice_base_salinity = init_ice_base_salinity[local_lat_ind*lon_dim + local_lon_ind];
            double current_ice_surface_salinity = init_ice_surface_salinity[local_lat_ind*lon_dim + local_lon_ind];
            
            
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
            int n_ice_cells = num_ice_layers;
            (*initial_ice_base_temp) = current_ice_base_temp;

            // initialize mandatory values for ice
            for (int i = 0; i < n_ice_cells; ++i)
            {
                (*initial_ice_temp_cells)[i] = current_ice_base_temp + 1.0*(i + 0.5)/(n_ice_cells)*(current_ice_surf_temp - current_ice_base_temp);
                (*initial_ice_sal_cells)[i] = current_ice_base_salinity + 1.0*(i + 0.5)/(n_ice_cells)*(current_ice_surface_salinity - current_ice_base_salinity);
                (*initial_ice_dens_cells)[i] = IceConsts<double>::rho_i;
            }

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

    // create 2d field of solvers
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        std::vector<icethermo::SeaIce_Solver<double>*> current_lat_solvers;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            Mesh<double>* local_ice_mesh = ice_meshes[local_lat_ind][local_lon_ind];
            Mesh<double>* local_snow_mesh = snow_meshes[local_lat_ind][local_lon_ind];

            current_lat_solvers.push_back
            (
                new(SeaIce1D_Snow0D_Solver<double>)(local_ice_mesh,
                                                    local_snow_mesh,
                                                    time_step,
                                                    min_ice_thick,
                                                    min_snow_thick,
                                                    true,
                                                    false,
                                                    is_verbose,
                                                    ice_cond_param,
                                                    ice_eff_capacity_param,
                                                    ice_enthalpy_param,
                                                    ice_heat_of_fusion_param,
                                                    ice_albedo_paaram,
                                                    snow_cond_param,
                                                    snow_heat_of_fusion_param,
                                                    snow_albedo_param
                                                    )
            );
        }
        solvers.push_back(current_lat_solvers);
    }
    if (is_verbose)
    {
        std::cout << "Solvers are initialized!\n";
    }
}

void ThermoModelsSet::SetComputationMarker(bool* do_compute_,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            do_compute[global_lat_ind][global_lon_ind] = (do_compute_[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] and
                                                       is_water[local_lat_ind][local_lon_ind]);
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d computational array!\n";
    }
}

void ThermoModelsSet::UpdateSwRadiation(double* sw_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_sw_value = sw_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateShortWaveRadiation
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

void ThermoModelsSet::UpdateLwRadiation(double* lw_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_lw_value = lw_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateLongWaveRadiation
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

void ThermoModelsSet::UpdateAirTemperature(double* ta_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_ta_value = ta_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateAtmosphereTemperature
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

void ThermoModelsSet::UpdatePrecipitationRate(double* pr_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_pr_value = pr_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdatePrecipitationRate
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

void ThermoModelsSet::UpdateAirPressure(double* ap_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_ap_value = ap_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateAtmospherePressure
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

void ThermoModelsSet::UpdateAirDensity(double* rho_values,
                                       int min_lon_ind_,        
                                       int max_lon_ind_,        
                                       int min_lat_ind_,        
                                       int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirDensity error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAirDensity error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirDensity error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAirDensity error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_rho_value = rho_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateAirDensity
            (
               current_rho_value
            );
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d atm density array!\n";
    }
}

void ThermoModelsSet::UpdateAirSpecificHumidity(double* sh_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_sh_value = sh_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateAirSpecificHumidity
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

void ThermoModelsSet::UpdateAbsWindSpeed(double* ws_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_ws_value = ws_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateAbsWindSpeed
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

void ThermoModelsSet::UpdateShCoeff(double* shc_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_shc_value = shc_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateShTransCoeff
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

void ThermoModelsSet::UpdateLhCoeff(double* lhc_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_lhc_value = lhc_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateLhTransCoeff
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

void ThermoModelsSet::AssembleTotalAtmFlux(int min_lon_ind_,        
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            solvers[global_lat_ind][global_lon_ind]->UpdateUpperFlux();
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d total atmosphere flux!\n";
    }
}

void ThermoModelsSet::UpdateOceanSalinity(double* os_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_os_value = os_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateOceanSalinity
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

void ThermoModelsSet::UpdateOceanFlux(double* of_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_of_value = of_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            solvers[global_lat_ind][global_lon_ind]->UpdateLowerFlux
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

void ThermoModelsSet::UpdateIceThickness(double* thick_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_ice_thickness = thick_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            for (int i = 0; i < num_ice_layers; ++i)
            {
                (*ice_meshes[global_lat_ind][global_lon_ind]->GetCellsThickness())[i] = current_ice_thickness/num_ice_layers;
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d ice thickness array!\n";
    }
}

void ThermoModelsSet::UpdateSnowThickness(double* thick_values,
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_snow_thickness = thick_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];

            (*snow_meshes[global_lat_ind][global_lon_ind]->GetCellsThickness())[0] = current_snow_thickness;
        }
    }

    if (is_verbose)
    {
        std::cout << "Updated 2d snow thickness array!\n";
    }
}

void ThermoModelsSetIce0dSnow0d::AddPrecipitation(int min_lon_ind_,        
                                                  int max_lon_ind_,        
                                                  int min_lat_ind_,        
                                                  int max_lat_ind_
                                                  )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"AddPrecipitation error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"AddPrecipitation error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"AddPrecipitation error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"AddPrecipitation error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    // to do
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_ice_thick = ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();

            if (is_water[global_lat_ind][global_lon_ind] and do_compute[global_lat_ind][global_lon_ind] and (current_ice_thick > min_ice_thick))
            {
                solvers[global_lat_ind][global_lon_ind]->AddPrecipitation();
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Added precipitation to snow thickness!\n";
    }
}


void ThermoModelsSet::Evaluate(int min_lon_ind_,        
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
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double current_ice_thick = ice_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();

            if (is_water[global_lat_ind][global_lon_ind] and do_compute[global_lat_ind][global_lon_ind] and (current_ice_thick > min_ice_thick))
            {
                solvers[global_lat_ind][global_lon_ind]->Evaluate();
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Single step evaluation of thermodynamics solver!\n";
    }
}


void ThermoModelsSet::GetSurfaceTemperature(double* array,
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
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double temp_value;

            if (is_water[global_lat_ind][global_lon_ind] and 
                (ice_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > min_ice_thick)
               )
            {
                if (snow_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD)
                {
                    temp_value = *(snow_meshes[global_lat_ind][global_lon_ind]->GetSingleData("up_temperature"));
                }
                else
                {
                    temp_value = *(ice_meshes[global_lat_ind][global_lon_ind]->GetSingleData("up_temperature"));
                }
            }
            else
            {
                temp_value = NAN_TEMP_VALUE;
            }

            array[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] = temp_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d surface temperature array!\n";
    }
}

void ThermoModelsSet::GetIceThickness(double* array,
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
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double thick_value;

            if (is_water[global_lat_ind][global_lon_ind] and 
                (ice_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > min_ice_thick))
            {
                thick_value = ice_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness();
            }
            else
            {
                thick_value = NAN_THICK_VALUE;
            }

            array[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] = thick_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d ice thickness array!\n";
    }
}

void ThermoModelsSet::GetSnowThickness(double* array,
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
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            double thick_value;

            if (is_water[global_lat_ind][global_lon_ind] and 
                (snow_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD))
            {
                thick_value = snow_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness();
            }
            else
            {
                thick_value = NAN_THICK_VALUE;
            }

            array[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] = thick_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d snow thickness array!\n";
    }
}

void ThermoModelsSet::GetIsIce(bool* array,
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
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            bool is_ice_value;

            if (is_water[global_lat_ind][global_lon_ind] and 
                (ice_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > min_ice_thick))
            {
                is_ice_value = true;
            }
            else
            {
                is_ice_value = false;
            }

            array[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] = is_ice_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d is_ice boolean array!\n";
    }
}

void ThermoModelsSet::GetIsSnow(bool* array,
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
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            bool is_snow_value;

            if (is_water[global_lat_ind][global_lon_ind] and 
                (snow_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD))
            {
                is_snow_value = true;
            }
            else
            {
                is_snow_value = false;
            }

            array[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind] = is_snow_value;
        }
    }

    if (is_verbose)
    {
        std::cout << "Output of 2d is_snow boolean array!\n";
    }
}

void ThermoModelsSetIce1dSnow0d::UpdateOceanIceMassFlux(double* flux_values,
                                                        int min_lon_ind_,        
                                                        int max_lon_ind_,        
                                                        int min_lat_ind_,        
                                                        int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanIceMassFlux error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanIceMassFlux error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanIceMassFlux error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateOceanIceMassFlux error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            if (is_water[global_lat_ind][global_lon_ind])
            {
                double current_mass_flux = flux_values[local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind];
                solvers[global_lat_ind][global_lon_ind]->UpdateOceanIceMassFlux(current_mass_flux);
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Update of 2d ocean->ice mass flux array!\n";
    }

}

void ThermoModelsSetIce1dSnow0d::UpdateTemperatureProfile(double* temp_values,
                                                          int min_lon_ind_,        
                                                          int max_lon_ind_,        
                                                          int min_lat_ind_,        
                                                          int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateTemperatureProfile error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateTemperatureProfile error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateTemperatureProfile error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateTemperatureProfile error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            std::vector<double> ice_temp_profile;
            double ice_surface_temp;
            double snow_surface_temp;

            if (is_water[global_lat_ind][global_lon_ind])
            {
                if (snow_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD)
                {
                    // snow is available

                    // get ice temperature profile + snow surface temp
                    for (int i = 0; i < num_ice_layers; ++i)
                    {
                        ice_temp_profile.push_back
                        (
                            temp_values[i*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                        + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                        + local_lon_ind]
                        );
                    }

                    // get snow surface temp
                    snow_surface_temp = temp_values[num_ice_layers*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                                    + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                                    + local_lon_ind];
                
                    // calculate ice-snow interface temp
                    double ice_temp = ice_temp_profile[num_ice_layers - 1];
                    double ice_thick = 0.5*(*(ice_meshes[global_lat_ind][global_lon_ind]->GetCellsThickness()))[num_ice_layers - 1];
                    double snow_thick = (*(snow_meshes[global_lat_ind][global_lon_ind]->GetCellsThickness()))[0];
                
                    double k_i = Params<double>::Conductivity
                    (
                        ice_cond_param,
                        ice_temp,
                        (*(ice_meshes[global_lat_ind][global_lon_ind]->GetCellsData("cells_salinity_array")))[num_ice_layers - 1],
                        (*(ice_meshes[global_lat_ind][global_lon_ind]->GetCellsData("cells_density_array")))[num_ice_layers - 1]
                    );

                    double k_s = Params<double>::Conductivity
                    (
                        snow_cond_param,
                        snow_surface_temp,
                        0.0,
                        (*(snow_meshes[global_lat_ind][global_lon_ind]->GetCellsData("cells_density_array")))[0]
                    );

                    double ratio = (k_i*snow_thick)/(k_s*ice_thick);
                    ice_surface_temp  = (1.0/(ratio + 1.0))*snow_surface_temp + ((ratio)/(ratio + 1.0))*ice_temp;
                }
                else
                {
                    // no snow

                    // get ice temperature profile + ice surface temp
                    for (int i = 0; i < num_ice_layers; ++i)
                    {
                        ice_temp_profile.push_back
                        (
                            temp_values[i*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                        + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                        + local_lon_ind]
                        );
                    }

                    // get ice surface temp
                    ice_surface_temp = temp_values[num_ice_layers*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                                    + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                                    + local_lon_ind];

                    snow_surface_temp = ice_surface_temp;
                }

                (*ice_meshes[global_lat_ind][global_lon_ind]->GetCellsData("cells_temperature_array")) = ice_temp_profile;
                (*ice_meshes[global_lat_ind][global_lon_ind]->GetSingleData("up_temperature")) = ice_surface_temp;
                (*snow_meshes[global_lat_ind][global_lon_ind]->GetSingleData("down_temperature")) = ice_surface_temp;
                (*snow_meshes[global_lat_ind][global_lon_ind]->GetCellsData("cells_temperature_array"))[0] = 0.5*(ice_surface_temp + snow_surface_temp);
                (*snow_meshes[global_lat_ind][global_lon_ind]->GetSingleData("up_temperature")) = snow_surface_temp;
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Update of 3d temperature profile array!\n";
    }
}

void ThermoModelsSetIce1dSnow0d::GetTemperatureProfile(double* temp_values,
                                                       int min_lon_ind_,        
                                                       int max_lon_ind_,        
                                                       int min_lat_ind_,        
                                                       int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetTemperatureProfile error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetTemperatureProfile error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetTemperatureProfile error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetTemperatureProfile error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;
        int global_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            int global_lon_ind = lon_ind - min_lon_ind;

            if (is_water[global_lat_ind][global_lon_ind])
            {
                for (int i = 0; i < num_ice_layers; ++i)
                {
                    if (ice_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > min_ice_thick)
                    {
                        temp_values[i*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                    + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                    + local_lon_ind] = 
                        (*ice_meshes[global_lat_ind][global_lon_ind]->GetCellsData("cells_temperature_array"))[i];
                    }
                    else
                    {
                        temp_values[i*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                    + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                    + local_lon_ind] = 
                        NAN_TEMP_VALUE;
                    }
                }

                if (snow_meshes[global_lat_ind][global_lon_ind]->GetTotalThickness() > SNOW_THICKNESS_THRESHOLD)
                {
                    // snow is available
                    temp_values[num_ice_layers*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                + local_lon_ind] = 
                    *snow_meshes[global_lat_ind][global_lon_ind]->GetSingleData("up_temperature");
                }
                else
                {
                    temp_values[num_ice_layers*(max_lon_ind_ - min_lon_ind_ + 1)*(max_lat_ind_ - min_lat_ind_ + 1) 
                                + local_lat_ind*(max_lon_ind_ - min_lon_ind_ + 1) 
                                + local_lon_ind] = 
                    NAN_TEMP_VALUE;
                }
            }
        }
    }

    if (is_verbose)
    {
        std::cout << "Recieve of 3d temperature profile array!\n";
    }
}


////////////////////////////////////////////////////
//                                                //
//      Relization for interface functions        //
//                                                //
////////////////////////////////////////////////////


void* InitThermodynamics0d(double time_step,
                           double min_ice_thick,
                           double min_snow_thick,
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
    ThermoModelsSetIce0dSnow0d* ptr = new ThermoModelsSetIce0dSnow0d(time_step,
                                                                     min_ice_thick,
                                                                     min_snow_thick,
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

void* InitThermodynamics1d(double time_step,
                           int num_ice_layers,
                           double min_ice_thick,
                           double min_snow_thick,
                           int min_lon_ind,
                           int max_lon_ind,
                           int min_lat_ind,
                           int max_lat_ind,
                           double* init_ice_base_temp,
                           double* init_ice_surf_temp,
                           double* init_snow_surf_temp,
                           double* init_ice_thick,
                           double* init_snow_thick,
                           double* init_ice_base_sal,
                           double* init_ice_surface_sal,
                           bool* water_marker,
                           bool is_verbose
                           )
{
    ThermoModelsSetIce1dSnow0d* ptr = new ThermoModelsSetIce1dSnow0d(time_step,
                                                                     num_ice_layers,
                                                                     min_ice_thick,
                                                                     min_snow_thick,
                                                                     min_lon_ind,
                                                                     max_lon_ind,
                                                                     min_lat_ind,
                                                                     max_lat_ind,
                                                                     init_ice_base_temp,
                                                                     init_ice_surf_temp,
                                                                     init_snow_surf_temp,
                                                                     init_ice_thick,
                                                                     init_snow_thick,
                                                                     init_ice_base_sal,
                                                                     init_ice_surface_sal,
                                                                     water_marker,
                                                                     is_verbose);
    
    return (void*)ptr;
}

void FinalizeThermodynamics0d(void* obj)
{
    ThermoModelsSetIce0dSnow0d* ptr = (ThermoModelsSetIce0dSnow0d*)obj;
    delete ptr;
}

void FinalizeThermodynamics1d(void* obj)
{
    ThermoModelsSetIce1dSnow0d* ptr = (ThermoModelsSetIce1dSnow0d*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateAirSpecificHumidity(sh_values,
                                   min_lon_ind,
                                   max_lon_ind,
                                   min_lat_ind,
                                   max_lat_ind);
}

void UpdateAirDensity(void* obj,
                      double* rho_values,
                      int min_lon_ind,        
                      int max_lon_ind,        
                      int min_lat_ind,        
                      int max_lat_ind)
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateAirDensity(rho_values,
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateOceanFlux(of_values,
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateSnowThickness(thick_values,
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
    
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateIceThickness(thick_values,
                            min_lon_ind,
                            max_lon_ind,
                            min_lat_ind,
                            max_lat_ind);
}

void AddPrecipitation(void* obj,
                      int min_lon_ind,        
                      int max_lon_ind,        
                      int min_lat_ind,        
                      int max_lat_ind)
{
    ThermoModelsSetIce0dSnow0d* ptr = (ThermoModelsSetIce0dSnow0d*)obj;
    ptr->AddPrecipitation(min_lon_ind,
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
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
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->GetIsIce(array,
                  min_lon_ind,        
                  max_lon_ind,        
                  min_lat_ind,        
                  max_lat_ind);
}

void UpdateOceanIceMassFlux(void* obj,
                            double* flux_values,
                            int min_lon_ind_,        
                            int max_lon_ind_,        
                            int min_lat_ind_,        
                            int max_lat_ind_)
{
    ThermoModelsSetIce1dSnow0d* ptr = (ThermoModelsSetIce1dSnow0d*)obj;
    ptr->UpdateOceanIceMassFlux(flux_values,
                                min_lon_ind_,
                                max_lon_ind_,
                                min_lat_ind_,
                                max_lat_ind_);
}

void UpdateTemperatureProfile(void* obj,
                              double* array,
                              int min_lon_ind_,        
                              int max_lon_ind_,        
                              int min_lat_ind_,        
                              int max_lat_ind_)
{
    ThermoModelsSetIce1dSnow0d* ptr = (ThermoModelsSetIce1dSnow0d*)obj;
    ptr->UpdateTemperatureProfile(array,
                                  min_lon_ind_,
                                  max_lon_ind_,
                                  min_lat_ind_,
                                  max_lat_ind_);
}

void GetTemperatureProfile(void* obj,
                           double* array,
                           int min_lon_ind_,        
                           int max_lon_ind_,        
                           int min_lat_ind_,        
                           int max_lat_ind_)
{
    ThermoModelsSetIce1dSnow0d* ptr = (ThermoModelsSetIce1dSnow0d*)obj;
    ptr->GetTemperatureProfile(array,
                               min_lon_ind_,
                               max_lon_ind_,
                               min_lat_ind_,
                               max_lat_ind_);
}