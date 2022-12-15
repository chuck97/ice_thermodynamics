#pragma once
#include "icethermo.hpp"

/*
2D mesh for 1D ice thermodinamic models
------------------------------------------------------------------------------
|  (lon_min, lat_min)     (lon_min+1, lat_min)     ...   (lon_max, lat_min)  |
| (lon_min, lat_min+1)   (lon_min+1, lat_min+1)    ...  (lon_max, lat_min+1) |
|         ...                     ...              ...          ...          |
|  (lon_min, lat_max)      (lon_min+1, lat_max)    ...   (lon_max, lat_max)  |
------------------------------------------------------------------------------
*/


class ThermoModelsSet
{
public:
    // constructor
    ThermoModelsSet(double time_step_,        // time step in seconds
                    int num_ice_cells,        // number of vertical ice cells for every 1D model
                    double min_ice_thick_,    // minimal ice thickness (meters) - if thickness is less than this value, model would not be evaluated
                    int min_lon_ind,          // minimal longitude index
                    int max_lon_ind,          // maximal longitude index (inclusive)
                    int min_lat_ind,          // minimal latitude index 
                    int max_lat_ind,          // maximal latitude index (inclusive)
                    double* init_base_temp,   // initial ice base temperatures array
                    double* init_surf_temp,   // initial ice surface temperatures array
                    double* init_ice_thick,   // initial ice thickness array
                    bool is_verbose_ = false  // verbose output? (default - true)
                    );


    // update total atmospheric flux
    void UpdateAtmFlux(double* atm_flux_values,
                       int min_lon_ind_,        
                       int max_lon_ind_,        
                       int min_lat_ind_,        
                       int max_lat_ind_
                       );

    // update short-wave radiation
    void UpdateSwRadiation(double* sw_values,
                           int min_lon_ind_,        
                           int max_lon_ind_,        
                           int min_lat_ind_,        
                           int max_lat_ind_
                           );

    // update latent heat flux
    void UpdateLatentHeatFlux(double* lh_values,
                              int min_lon_ind_,        
                              int max_lon_ind_,        
                              int min_lat_ind_,        
                              int max_lat_ind_
                              );

    // one-step evaluation
    void Evaluate(int min_lon_ind_,        
                  int max_lon_ind_,        
                  int min_lat_ind_,        
                  int max_lat_ind_
                  );

    // get 2d-array of surface temperatures
    void GetIceSurfaceTemperature(double* array,
                                  int min_lon_ind_,        
                                  int max_lon_ind_,        
                                  int min_lat_ind_,        
                                  int max_lat_ind_
                                  );

    // get 2d-array of ice thicknesses
    void GetIceThickness(double* array,
                         int min_lon_ind_,        
                         int max_lon_ind_,        
                         int min_lat_ind_,        
                         int max_lat_ind_
                         );

    // get 2d-array of conductive flux at the ice surface
    void GetSurfaceConductiveFlux(double* array,
                                  int min_lon_ind_,        
                                  int max_lon_ind_,        
                                  int min_lat_ind_,        
                                  int max_lat_ind_
                                  );

    // destructor
    ~ThermoModelsSet();

private:
    int time_step;
    double min_ice_thick;
    int lon_dim, lat_dim;
    int min_lon_ind, max_lon_ind;
    int min_lat_ind, max_lat_ind;
    bool is_verbose;
    std::vector<std::vector<icethermo::Mesh<double>*>> meshes;
    std::vector<std::vector<icethermo::Glacier1D_Solver<double>*>> solvers;

    icethermo::Kparam k_param = icethermo::Kparam::FreshIce;
};

// Interface functions to call from FORTRAN
extern "C" 
{
   // initialization of thermodynamics
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
                            );

    // finalization of thermodynamics
    void FinalizeThermodynamics(void* obj);

    // update atmosphere flux
    void UpdateAtmFlux(void* obj,
                       double* atm_flux_values,
                       int min_lon_ind,        
                       int max_lon_ind,        
                       int min_lat_ind,        
                       int max_lat_ind
                       );
    
    // update short-wave radiation
    void UpdateSwRadiation(void* obj,
                           double* sw_values,
                           int min_lon_ind,        
                           int max_lon_ind,        
                           int min_lat_ind,        
                           int max_lat_ind
                           );
    
    // update latent heat flux
    void UpdateLatentHeatFlux(void* obj,
                              double* lh_values,
                              int min_lon_ind,        
                              int max_lon_ind,        
                              int min_lat_ind,        
                              int max_lat_ind
                              );
    
    // one-step evaluation
    void Evaluate(void* obj,
                  int min_lon_ind,        
                  int max_lon_ind,        
                  int min_lat_ind,        
                  int max_lat_ind
                  );

    // return ice surface temperature
    void GetIceSurfaceTemperature(void* obj,
                                  double* array,
                                  int min_lon_ind,        
                                  int max_lon_ind,        
                                  int min_lat_ind,        
                                  int max_lat_ind
                                  );
    
    // return ice thickness
    void GetIceThickness(void* obj,
                         double* array,
                         int min_lon_ind,        
                         int max_lon_ind,        
                         int min_lat_ind,        
                         int max_lat_ind
                         );
    
    // return ice surface conductive heat flux
    void GetSurfaceConductiveFlux(void* obj,
                                  double* array,
                                  int min_lon_ind,        
                                  int max_lon_ind,        
                                  int min_lat_ind,        
                                  int max_lat_ind
                                  );
}