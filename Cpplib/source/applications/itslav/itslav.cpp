#include "itslav.hpp"

using namespace icethermo;

// constructor 
ThermoModelsSet::ThermoModelsSet(double time_step_,
                                 int num_ice_cells_,
                                 double min_ice_thick_,
                                 int min_lon_ind_,
                                 int max_lon_ind_,
                                 int min_lat_ind_,
                                 int max_lat_ind_,
                                 double* init_base_temp,
                                 double* init_surf_temp,
                                 double* init_ice_thick,
                                 bool is_verbose_
                                 ):
        time_step(time_step_),
        num_ice_cells(num_ice_cells_),
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

    // create 2d field of 1D meshes
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        std::vector<Mesh<double>*> current_lat_meshes;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_thick = init_ice_thick[local_lat_ind*lon_dim + local_lon_ind];

            current_lat_meshes.push_back
            (
                new(Mesh<double>)(num_ice_cells, current_ice_thick)
            );
        }
        current_meshes.push_back(current_lat_meshes);
    }

    if (is_verbose)
    {
        std::cout << "Meshes are created!\n";
    }

    // initialize 2d field of 1D meshes (linear temperature profiles from base to surface temperature, constant densities)
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            double current_ice_thick = init_ice_thick[local_lat_ind*lon_dim + local_lon_ind];

            double current_base_temp = init_base_temp[local_lat_ind*lon_dim + local_lon_ind];

            double current_surf_temp = init_surf_temp[local_lat_ind*lon_dim + local_lon_ind];
            
            Mesh<double>* local_mesh = current_meshes[local_lat_ind][local_lon_ind];

            auto initial_temp_cells = local_mesh->CreateCellsData("cells_temperature_array");
            auto initial_dens_cells = local_mesh->CreateCellsData("cells_density_array");
            auto initial_ice_surf_temp = local_mesh->CreateSingleData("up_temperature");
            auto initial_ice_base_temp = local_mesh->CreateSingleData("down_temperature");

            // initialize mandatory values
            (*initial_ice_base_temp) = current_base_temp;

            for (int i = 0; i < local_mesh->GetCellsNum(); ++i)
            {
                (*initial_temp_cells)[i] = current_base_temp + (i + 0.5)/(local_mesh->GetCellsNum())*(current_surf_temp - current_base_temp);
                (*initial_dens_cells)[i] = IceConsts<double>::rho_i;
            }

            (*initial_ice_surf_temp) = current_surf_temp;
        }
    }

    if (is_verbose)
    {
        std::cout << "Meshes are initialized!\n";
    }

    // create 2d field of 1D solvers
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        std::vector<Glacier1D_Solver<double>*> current_lat_solvers;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;
            Mesh<double>* local_mesh = current_meshes[local_lat_ind][local_lon_ind];

            current_lat_solvers.push_back
            (
                new(Glacier1D_Solver<double>)(local_mesh,
                                              time_step,
                                              true,
                                              true,
                                              false,
                                              ApproxOrder::second,
                                              Kparam::FreshIce,
                                              Cparam::FreshIce,
                                              Eparam::FreshIce,
                                              Lparam::FreshIce)
            );
        }
        solvers.push_back(current_lat_solvers);
    }

    if (is_verbose)
    {
        std::cout << "Solvers are initialized!\n";
    }

    // make copy of current meeshes to store
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind_;

        std::vector<Mesh<double>*> current_lat_meshes;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind_;

            current_lat_meshes.push_back
            (
                new(Mesh<double>)(*(current_meshes[local_lat_ind][local_lon_ind]))
            );
        }
        stored_meshes.push_back(current_lat_meshes);
    }

    if (is_verbose)
    {
        std::cout << "Meshes are stored!\n";
    }
}

// realization of update total atmospheric flux
void ThermoModelsSet::UpdateAtmFlux(double* atm_flux_values,
                                    int min_lon_ind_,        
                                    int max_lon_ind_,        
                                    int min_lat_ind_,        
                                    int max_lat_ind_
                                    )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAtmFlux error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateAtmFlux error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAtmFlux error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateAtmFlux error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double current_flux_value = atm_flux_values[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_];

            solvers[local_lat_ind][local_lon_ind]->UpdateUpperFlux
            (
               [current_flux_value](double T){return current_flux_value;}
            );
        }
    }
}

// realization of update short-wave radiation
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
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double current_sw_value = sw_values[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_];

            solvers[local_lat_ind][local_lon_ind]->UpdateShortWaveRadiation
            (
                [current_sw_value](double temp) 
                {
                    return current_sw_value;
                }
            );     
        }
    }
}

// realization of update latent heat flux
void ThermoModelsSet::UpdateLatentHeatFlux(double* lh_values,
                                           int min_lon_ind_,        
                                           int max_lon_ind_,        
                                           int min_lat_ind_,        
                                           int max_lat_ind_
                                           )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateLatentHeatFlux error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"UpdateLatentHeatFlux error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateLatentHeatFlux error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"UpdateLatentHeatFlux error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double current_lh_value = lh_values[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_];

            solvers[local_lat_ind][local_lon_ind]->UpdateLatentHeatFlux
            (
                [current_lh_value](double temp) 
                {
                    return current_lh_value;
                }
            );     
        }
    }
}

// realization of one-step evaluation
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

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;

            double current_ice_thick = current_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();

            if (current_ice_thick > min_ice_thick)
            {
                solvers[local_lat_ind][local_lon_ind]->Evaluate();
            }
        }
    }
}

// realization of get 2d-array of surface temperatures
void ThermoModelsSet::GetIceSurfaceTemperature(double* array,
                                               int min_lon_ind_,        
                                               int max_lon_ind_,        
                                               int min_lat_ind_,        
                                               int max_lat_ind_
                                               )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetIceSurfaceTemperature error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetIceSurfaceTemperature error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetIceSurfaceTemperature error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetIceSurfaceTemperature error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double current_surf_temp = *(current_meshes[local_lat_ind][local_lon_ind]->GetSingleData("up_temperature"));

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = current_surf_temp;
        }
    }
}

// realization of get 2d-array of ice thicknesses
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
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double current_thick = current_meshes[local_lat_ind][local_lon_ind]->GetTotalThickness();

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = current_thick;
        }
    }
}

// realization of get 2d-array of conductive flux at the ice surface
void ThermoModelsSet::GetSurfaceConductiveFlux(double* array,
                                               int min_lon_ind_,        
                                               int max_lon_ind_,        
                                               int min_lat_ind_,        
                                               int max_lat_ind_
                                               )
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceConductiveFlux error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceConductiveFlux error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceConductiveFlux error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"GetSurfaceConductiveFlux error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;
        int local_lat_ind_ = lat_ind - min_lat_ind_;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;
            int local_lon_ind_ = lon_ind - min_lon_ind_;

            double surf_temp = *(current_meshes[local_lat_ind][local_lon_ind]->GetSingleData("up_temperature"));
            double top_layer_temp = (*(current_meshes[local_lat_ind][local_lon_ind]->GetCellsData("cells_temperature_array"))).back();
            double top_layer_size =  (*(current_meshes[local_lat_ind][local_lon_ind]->GetCellsThickness())).back();
            double k_value = Params<double>::Conductivity(k_param, top_layer_temp, 0.0, IceConsts<double>::rho_i);

            double cond_flux = (top_layer_size < REAL_MIN_VAL(double)) ? 0.0 : -k_value*(surf_temp - top_layer_temp)/(0.5*top_layer_size);

            array[local_lat_ind_*(max_lon_ind_ - min_lon_ind_ + 1) + local_lon_ind_] = cond_flux;
        }
    }
}

// realization of save current ModelsSet current_state to stored_state 
void ThermoModelsSet::StoreState(int min_lon_ind_,        
                                 int max_lon_ind_,        
                                 int min_lat_ind_,        
                                 int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"StoreState error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"StoreState error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"StoreState error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"StoreState error: maximal latitude index is greater than maximal mesh latitude index!");
    }
    
    // store 2d field of 1D meshes
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;

            icethermo::Mesh<double>* stored_mesh_ptr = stored_meshes[local_lat_ind][local_lon_ind];
            delete stored_mesh_ptr;
            stored_mesh_ptr = new(icethermo::Mesh<double>)(*current_meshes[local_lat_ind][local_lon_ind]);
        }
    }
}

void ThermoModelsSet::RestoreState(int min_lon_ind_,        
                                   int max_lon_ind_,        
                                   int min_lat_ind_,        
                                   int max_lat_ind_)
{
    if (min_lon_ind_ < min_lon_ind)
    {
        THERMO_ERR((std::string)"RestoreState error: minimal longitude index is less than minimal mesh longitude index!");
    }

    if (max_lon_ind_ > max_lon_ind)
    {
        THERMO_ERR((std::string)"RestoreState error: maximal longitude index is greater than maximal mesh longitude index!");
    }

    if (min_lat_ind_ < min_lat_ind)
    {
        THERMO_ERR((std::string)"RestoreState error: minimal latitude index is less than minimal mesh latitude index!");
    }

    if (max_lat_ind_ > max_lat_ind)
    {
        THERMO_ERR((std::string)"RestoreState error: maximal latitude index is greater than maximal mesh latitude index!");
    }

    // restore 2d field of 1D meshes and update solvers
    for (int lat_ind = min_lat_ind_; lat_ind < max_lat_ind_ + 1; ++lat_ind)
    {
        int local_lat_ind = lat_ind - min_lat_ind;

        for (int lon_ind = min_lon_ind_; lon_ind < max_lon_ind_ + 1; ++lon_ind)
        {
            int local_lon_ind = lon_ind - min_lon_ind;

            icethermo::Mesh<double>* current_mesh_ptr = current_meshes[local_lat_ind][local_lon_ind];
            delete current_mesh_ptr;
            current_mesh_ptr = new(icethermo::Mesh<double>)(*stored_meshes[local_lat_ind][local_lon_ind]);
            solvers[local_lat_ind][local_lon_ind]->UpdateMesh(current_mesh_ptr);
        }
    }
}

// realization of initialization of thermodynamics
void* InitThermodynamics(double time_step,
                         int num_ice_cells,
                         double min_ice_thick,
                         int min_lon_ind,
                         int max_lon_ind,
                         int min_lat_ind,
                         int max_lat_ind,
                         double* init_base_temp,
                         double* init_surf_temp,
                         double* init_ice_thick
                         )
{
    ThermoModelsSet* ptr = new ThermoModelsSet(time_step,
                                               num_ice_cells,
                                               min_ice_thick,
                                               min_lon_ind,
                                               max_lon_ind,
                                               min_lat_ind,
                                               max_lat_ind,
                                               init_base_temp,
                                               init_surf_temp,
                                               init_ice_thick);
    
    return (void*)ptr;
}

// realization of finalization of thermodynamics
void FinalizeThermodynamics(void* obj)
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    delete ptr;
}


// realization of update atmosphere flux
void UpdateAtmFlux(void* obj,
                   double* atm_flux_values,
                   int min_lon_ind,        
                   int max_lon_ind,        
                   int min_lat_ind,        
                   int max_lat_ind
                   )
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateAtmFlux(atm_flux_values,
                       min_lon_ind,        
                       max_lon_ind,        
                       min_lat_ind,        
                       max_lat_ind);
}

// realization of update short-wave radiation
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

// realization of update latent heat flux
void UpdateLatentHeatFlux(void* obj,
                          double* lh_values,
                          int min_lon_ind,        
                          int max_lon_ind,        
                          int min_lat_ind,        
                          int max_lat_ind)
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->UpdateLatentHeatFlux(lh_values,
                              min_lon_ind,        
                              max_lon_ind,        
                              min_lat_ind,        
                              max_lat_ind);
}

// realization of one-step evaluation
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

// realization of return ice surface temperature
void GetIceSurfaceTemperature(void* obj,
                              double* array,
                              int min_lon_ind,        
                              int max_lon_ind,        
                              int min_lat_ind,        
                              int max_lat_ind)
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->GetIceSurfaceTemperature(array,
                                  min_lon_ind,        
                                  max_lon_ind,        
                                  min_lat_ind,        
                                  max_lat_ind);
}

// realization of return ice thickness
void GetIceThickness(void* obj,
                     double* array,
                     int min_lon_ind,        
                     int max_lon_ind,        
                     int min_lat_ind,        
                     int max_lat_ind)
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->GetIceThickness(array,
                         min_lon_ind,        
                         max_lon_ind,        
                         min_lat_ind,        
                         max_lat_ind);
}

// realization of return ice surface conductive heat flux
void GetSurfaceConductiveFlux(void* obj,
                              double* array,
                              int min_lon_ind,        
                              int max_lon_ind,        
                              int min_lat_ind,        
                              int max_lat_ind)
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->GetSurfaceConductiveFlux(array,
                                  min_lon_ind,        
                                  max_lon_ind,        
                                  min_lat_ind,        
                                  max_lat_ind);
}

// realization of storing
void StoreState(void* obj,
                int min_lon_ind,        
                int max_lon_ind,        
                int min_lat_ind,        
                int max_lat_ind
                )
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->StoreState(min_lon_ind,        
                    max_lon_ind,        
                    min_lat_ind,        
                    max_lat_ind);
}

// realization of restoring
void RestoreState(void* obj,
                  int min_lon_ind,        
                  int max_lon_ind,        
                  int min_lat_ind,        
                  int max_lat_ind
                  )
{
    ThermoModelsSet* ptr = (ThermoModelsSet*)obj;
    ptr->RestoreState(min_lon_ind,        
                      max_lon_ind,        
                      min_lat_ind,        
                      max_lat_ind);
}