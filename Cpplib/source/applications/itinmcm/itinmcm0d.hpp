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
class ThermoModelsSet0D
{
public:
    // constructor
    ThermoModelsSet0D(double time_step_,           // time step in seconds
                      double min_ice_thick_,       // minimal ice thickness (meters) - if thickness is less than this value, model would not be evaluated
                      int min_lon_ind_,            // minimal longitude index
                      int max_lon_ind_,            // maximal longitude index (inclusive)
                      int min_lat_ind_,            // minimal latitude index 
                      int max_lat_ind_,            // maximal latitude index (inclusive)
                      double* init_ice_base_temp,  // initial ice base temperatures array
                      double* init_ice_surf_temp,  // initial ice surface temperatures array (2d)
                      double* init_snow_surf_temp, // initial snow surface temperatures array (2d)
                      double* init_ice_thick,      // initial ice thickness array (2d)
                      double* init_snow_thick,     // initial ice thickness array (2d)
                      bool* water_marker,          // marker of water array (2d)
                      bool is_verbose_ = true      // verbose output? (default = true)
                      );

    // set 2d-computation marker
    void SetComputationMarker(bool* do_compute_,
                              int min_lon_ind_,        
                              int max_lon_ind_,        
                              int min_lat_ind_,        
                              int max_lat_ind_
                              );

    // update descending short-wave radiation
    void UpdateSwRadiation(double* sw_values,
                           int min_lon_ind_,        
                           int max_lon_ind_,        
                           int min_lat_ind_,        
                           int max_lat_ind_
                           );
    
    // update descending long-wave radiation
    void UpdateLwRadiation(double* lw_values,
                           int min_lon_ind_,        
                           int max_lon_ind_,        
                           int min_lat_ind_,        
                           int max_lat_ind_
                           );
    
    // update atmosphere temperature
    void UpdateAirTemperature(double* ta_values,
                              int min_lon_ind_,        
                              int max_lon_ind_,        
                              int min_lat_ind_,        
                              int max_lat_ind_
                             );

    // update liquid precipitation rate
    void UpdatePrecipitationRate(double* pr_values,
                                 int min_lon_ind_,        
                                 int max_lon_ind_,        
                                 int min_lat_ind_,        
                                 int max_lat_ind_
                                 );
    
    // update atmospheric pressure
    void UpdateAirPressure(double* ap_values,
                           int min_lon_ind_,        
                           int max_lon_ind_,        
                           int min_lat_ind_,        
                           int max_lat_ind_);
    
    // update air specific humidity
    void UpdateAirSpecificHumidity(double* sh_values,
                                   int min_lon_ind_,        
                                   int max_lon_ind_,        
                                   int min_lat_ind_,        
                                   int max_lat_ind_);
    
    // update air specific humidity
    void UpdateAbsWindSpeed(double* ws_values,
                            int min_lon_ind_,        
                            int max_lon_ind_,        
                            int min_lat_ind_,        
                            int max_lat_ind_);
    
    // update surface sensible heat transfer coefficient
    void UpdateShCoeff(double* shc_values,
                       int min_lon_ind_,        
                       int max_lon_ind_,        
                       int min_lat_ind_,        
                       int max_lat_ind_);

    // update surface latent heat transport coefficient
    void UpdateLhCoeff(double* lhc_values,
                       int min_lon_ind_,        
                       int max_lon_ind_,        
                       int min_lat_ind_,        
                       int max_lat_ind_);

    // update ocean salinity
    void UpdateOceanSalinity(double* os_values,
                             int min_lon_ind_,        
                             int max_lon_ind_,        
                             int min_lat_ind_,        
                             int max_lat_ind_
                             );

    // update ocean heat flux
    void UpdateOceanFlux(double* of_values,
                         int min_lon_ind_,        
                         int max_lon_ind_,        
                         int min_lat_ind_,        
                         int max_lat_ind_
                        );
    
    // update ice thickness
    void UpdateIceThickness(double* thick_values,
                            int min_lon_ind_,        
                            int max_lon_ind_,        
                            int min_lat_ind_,        
                            int max_lat_ind_
                            );

    // update snow thickness
    void UpdateSnowThickness(double* thick_values,
                             int min_lon_ind_,        
                             int max_lon_ind_,        
                             int min_lat_ind_,        
                             int max_lat_ind_
                             );

    // single-step evaluation (computation of temps and thicks)
    void Evaluate(int min_lon_ind_,        
                  int max_lon_ind_,        
                  int min_lat_ind_,        
                  int max_lat_ind_
                  );

    // recieve 2d-array of surface temperatures (ice or snow)
    void GetSurfaceTemperature(double* array,
                               int min_lon_ind_,        
                               int max_lon_ind_,        
                               int min_lat_ind_,        
                               int max_lat_ind_
                               );
    
    // recieve 2d-array (boolean) of snow presence
    void GetIsSnow(bool* array,
                   int min_lon_ind_,        
                   int max_lon_ind_,        
                   int min_lat_ind_,        
                   int max_lat_ind_
                   );
    
    // get 2d-array (boolean) of ice presence
    void GetIsIce(bool* array,
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
    
    // get 2d-array of snow thicknesses
    void GetSnowThickness(double* array,
                          int min_lon_ind_,        
                          int max_lon_ind_,        
                          int min_lat_ind_,        
                          int max_lat_ind_
                          );

private:

    double time_step;

    int lon_dim;
    int lat_dim;

    int min_lon_ind;
    int max_lon_ind;
    int min_lat_ind;
    int max_lat_ind;

    double min_ice_thick;
    bool is_verbose;

    std::vector<std::vector<icethermo::Mesh<double>*>> ice_meshes;
    std::vector<std::vector<icethermo::Mesh<double>*>> snow_meshes;

    std::vector<std::vector<bool>> is_water;
    std::vector<std::vector<bool>> do_compute;

    std::vector<std::vector<icethermo::SeaIce0D_Snow0D_Solver<double>*>> solvers;
};

// Interface functions to call from FORTRAN
extern "C" 
{
   // initialization of thermodynamics
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
                            bool is_verbose = true
                            );

    // finalization of thermodynamics
    void FinalizeThermodynamics(void* obj);

    // setup 2d computation marker
    void SetComputationMarker(void* obj,
                              bool* do_compute,
                              int min_lon_ind,        
                              int max_lon_ind,        
                              int min_lat_ind,        
                              int max_lat_ind
                              );

    // update descending 2d short-wave radiation
    void UpdateSwRadiation(void* obj,
                           double* sw_values,
                           int min_lon_ind,        
                           int max_lon_ind,        
                           int min_lat_ind,        
                           int max_lat_ind
                           );
    
    // update descending 2d long-wave radiation
    void UpdateLwRadiation(void* obj,
                           double* lw_values,
                           int min_lon_ind,        
                           int max_lon_ind,        
                           int min_lat_ind,        
                           int max_lat_ind
                           );

    // update 2d air temperature
    void UpdateAirTemperature(void* obj,
                              double* ta_values,
                              int min_lon_ind,        
                              int max_lon_ind,        
                              int min_lat_ind,        
                              int max_lat_ind
                              );

    // update 2d liquid precipitation rate
    void UpdatePrecipitationRate(void* obj,
                                 double* pr_values,
                                 int min_lon_ind,        
                                 int max_lon_ind,        
                                 int min_lat_ind,        
                                 int max_lat_ind
                                 );
    
    // update 2d air pressure
    void UpdateAirPressure(void* obj,
                           double* ap_values,
                           int min_lon_ind,        
                           int max_lon_ind,        
                           int min_lat_ind,        
                           int max_lat_ind
                           );

    // update 2d air specific humidity
    void UpdateAirSpecificHumidity(void* obj,
                                   double* sh_values,
                                   int min_lon_ind,        
                                   int max_lon_ind,        
                                   int min_lat_ind,        
                                   int max_lat_ind
                                   );

    // update 2d air wind speed
    void UpdateAbsWindSpeed(void* obj,
                            double* ws_values,
                            int min_lon_ind,        
                            int max_lon_ind,        
                            int min_lat_ind,        
                            int max_lat_ind);
    
    // update 2d sensible heat coefficient
    void UpdateShCoeff(void* obj,
                       double* shc_values,
                       int min_lon_ind,        
                       int max_lon_ind,        
                       int min_lat_ind,        
                       int max_lat_ind);

    // update 2d latent heat coefficient
    void UpdateLhCoeff(void* obj,
                       double* lhc_values,
                       int min_lon_ind,        
                       int max_lon_ind,        
                       int min_lat_ind,        
                       int max_lat_ind);

    // update 2d ocean salinity
    void UpdateOceanSalinity(void* obj,
                             double* os_values,
                             int min_lon_ind,        
                             int max_lon_ind,        
                             int min_lat_ind,        
                             int max_lat_ind
                             );
    
    // update 2d ocean heat flux
    void UpdateOceanFlux(void* obj,
                         double* of_values,
                         int min_lon_ind_,        
                         int max_lon_ind_,        
                         int min_lat_ind_,        
                         int max_lat_ind_
                        );
    
    // update 2d ice thickness
    void UpdateIceThickness(void* obj,
                            double* thick_values,
                            int min_lon_ind_,        
                            int max_lon_ind_,        
                            int min_lat_ind_,        
                            int max_lat_ind_
                            );
    
    // update 2d snow thickness
    void UpdateSnowThickness(void* obj,
                             double* thick_values,
                             int min_lon_ind_,        
                             int max_lon_ind_,        
                             int min_lat_ind_,        
                             int max_lat_ind_
                             );
    
    // single-step evaluation (computation of temps and thicks)
    void Evaluate(void* obj,
                  int min_lon_ind,        
                  int max_lon_ind,        
                  int min_lat_ind,        
                  int max_lat_ind
                  );

    // recieve 2d-array of surface temperatures
    void GetSurfaceTemperature(void* obj,
                               double* array,
                               int min_lon_ind,        
                               int max_lon_ind,        
                               int min_lat_ind,        
                               int max_lat_ind
                               );
    
    // recieve 2d ice thickness
    void GetIceThickness(void* obj,
                         double* array,
                         int min_lon_ind,        
                         int max_lon_ind,        
                         int min_lat_ind,        
                         int max_lat_ind
                         );
    
    // recieve 2d snow thickness
    void GetSnowThickness(void* obj,
                         double* array,
                         int min_lon_ind,        
                         int max_lon_ind,        
                         int min_lat_ind,        
                         int max_lat_ind
                         );
    
    // recieve 2d-array (boolean) of snow presence
    void GetIsSnow(bool* array,
                   int min_lon_ind_,        
                   int max_lon_ind_,        
                   int min_lat_ind_,        
                   int max_lat_ind_
                   );
    
    // get 2d-array (boolean) of ice presence
    void GetIsIce(bool* array,
                  int min_lon_ind_,        
                  int max_lon_ind_,        
                  int min_lat_ind_,        
                  int max_lat_ind_
                  );
}