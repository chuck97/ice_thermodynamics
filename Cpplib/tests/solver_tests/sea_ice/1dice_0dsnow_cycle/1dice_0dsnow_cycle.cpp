#include "icethermo.hpp"
#include <string>

using namespace icethermo;

template <typename NumType>
NumType mod(NumType lhs, NumType rhs)
{
    return lhs - std::floor(lhs/rhs)*rhs;
}

// "annual" atm temperature and precipitation rate
template<typename NumType>
std::pair<NumType, NumType> annual_temp_prec(NumType winter_duration, // duration of winter (seconds)
                                             NumType spring_duration, // duration of spring (seconds)
                                             NumType summer_duration, // duration of summer (seconds)
                                             NumType autumn_duration, // duration of autumn (seconds)
                                             NumType lowest_temp,     // lowest atm temperature (deg Cel)
                                             NumType highest_temp,    // highest atm temperature (deg Cel)
                                             NumType lowest_p,        // lowest precipitation rate (m s-1) - not needed in this test
                                             NumType highest_p,       // highest precipitation rate (m s-1) - not needed in this test
                                             NumType time)            // current time (seconds) 
{
    NumType t = mod(time, winter_duration + spring_duration + summer_duration + autumn_duration);

    if (t < winter_duration)
    {
        return {lowest_temp, highest_p};
    }
    else if (t < (winter_duration+ spring_duration))
    {
        return {lowest_temp + (t - winter_duration)*(highest_temp - lowest_temp)/spring_duration,
                highest_p + (t - winter_duration)*(lowest_p - highest_p)/spring_duration};
    }
    else if (t < (winter_duration + spring_duration + summer_duration))
    {
        return {highest_temp, lowest_p};
    }
    else
    {
        return {highest_temp + (t - (winter_duration + spring_duration + summer_duration))*(lowest_temp - highest_temp)/autumn_duration,
                lowest_p + (t - (winter_duration + spring_duration + summer_duration))*(highest_p - lowest_p)/autumn_duration};
    }
}

// non-stationary atmospheric flux function
template<typename NumType>
NumType atm_flux(NumType temp, NumType time, NumType atm_temp)
{
    // simple sensible heat flux parameterization
    return AirConsts<NumType>::rho_a *
           AirConsts<NumType>::cp_a *
           GenConsts<NumType>::C_sh *
           (NumType)15.0 * 
           (atm_temp - temp);
}

// model launcher
template <typename NumType>
void run_model(NumType time_step,
               int num_steps,
               int output_frequency,
               int num_cells,
               NumType initial_ice_thickness,
               NumType initial_snow_thickness,
               NumType initial_interface_temperature,
               NumType initial_snow_surface_temperature,
               NumType winter_duration, 
               NumType spring_duration, 
               NumType summer_duration, 
               NumType autumn_duration, 
               NumType lowest_temp,     
               NumType highest_temp,
               NumType lowest_prec_rate,     
               NumType highest_prec_rate,    
               const std::string& output_ice_prefix,
               const std::string& output_snow_prefix,
               Kparam conductivity_ice_parameterization,
               Cparam eff_capacity_ice_parameterization,
               Eparam enthalpy_ice_parameterization,
               Lparam fusion_heat_ice_parameterization,
               Aparam albedo_ice_perameterization,
               Kparam conductivity_snow_parameterization,
               Lparam fusion_heat_snow_parameterization,
               Aparam albedo_snow_parameterization)
{

    // configure constants from .json file
    Configure<NumType>("/home/chuck/Desktop/ice_thermodynamics/Cpplib/configs/test_config.json");

    // create uniform sigma-mesh for 1d ice
    Mesh<NumType>* ice_mesh = new(Mesh<NumType>)(num_cells, initial_ice_thickness);

    // create mandatory initial values for ice
    auto initial_ice_temp_cells = ice_mesh->CreateCellsData("cells_temperature_array");
    auto initial_ice_sal_cells = ice_mesh->CreateCellsData("cells_salinity_array");
    auto initial_ice_dens_cells = ice_mesh->CreateCellsData("cells_density_array");
    auto initial_ice_surf_temp = ice_mesh->CreateSingleData("up_temperature");
    int n_ice_cells = ice_mesh->GetCellsNum();

    // create single-cell sigma-mesh for 0d snow
    Mesh<NumType>* snow_mesh = new(Mesh<NumType>)(1, initial_snow_thickness);

    // create mandatory initial values for snow
    auto initial_snow_temp_cells = snow_mesh->CreateCellsData("cells_temperature_array");
    auto initial_snow_dens_cells = snow_mesh->CreateCellsData("cells_density_array");
    auto initial_snow_surf_temp = snow_mesh->CreateSingleData("up_temperature");

    // freezing point of salty water with 30 psu
    NumType fusion_temp = Params<NumType>::TempFusion(30.0);

    // initialize mandatory values for ice
    for (int i = 0; i < n_ice_cells; ++i)
    {
       (*initial_ice_temp_cells)[i] = fusion_temp + (NumType)1.0*(i + (NumType)0.5)/(n_ice_cells)*(initial_interface_temperature - fusion_temp);
       (*initial_ice_sal_cells)[i] = (NumType)4.0 + (NumType)1.0*(i + (NumType)0.5)/(n_ice_cells)*((NumType)1.0 - (NumType)4.0);
       (*initial_ice_dens_cells)[i] = IceConsts<NumType>::rho_i;
    }
    (*initial_ice_surf_temp) = initial_interface_temperature;

    // initialize mandatory values for snow
    (*initial_snow_temp_cells)[0] = (NumType)0.5*(initial_snow_surface_temperature + initial_interface_temperature);
    (*initial_snow_dens_cells)[0] = SnowConsts<NumType>::rho_s;
    (*initial_snow_surf_temp) = initial_snow_surface_temperature;

    // create 1d ice with 0d snow solver class
    SeaIce1D_Snow0D_Solver<NumType> thermo_solver(ice_mesh,
                                                  snow_mesh,
                                                  time_step,
                                                  0.001,
                                                  0.001,
                                                  true,
                                                  true,
                                                  true,
                                                  conductivity_ice_parameterization,
                                                  eff_capacity_ice_parameterization,
                                                  enthalpy_ice_parameterization,
                                                  fusion_heat_ice_parameterization,
                                                  albedo_ice_perameterization,
                                                  conductivity_snow_parameterization,
                                                  fusion_heat_snow_parameterization,
                                                  albedo_snow_parameterization);

    // save initial state to file
#ifdef USE_JSON_OUTPUT
    ice_mesh->SaveJSON((std::string)"./" + output_ice_prefix, 0);
    snow_mesh->SaveJSON((std::string)"./" + output_snow_prefix, 0);
#endif

    // time stepping
    for (int step_num = 1; step_num < num_steps + 1; ++step_num)
    {
        // update atmospheric flux
        auto curr_forc = annual_temp_prec(winter_duration,
                                          spring_duration,
                                          summer_duration,
                                          autumn_duration,
                                          lowest_temp,
                                          highest_temp,
                                          lowest_prec_rate,
                                          highest_prec_rate,
                                          step_num*time_step);
        thermo_solver.UpdateUpperFlux
        (
            [&step_num,
             &time_step,
             &curr_forc](NumType temp)
            {
                return atm_flux<NumType>(temp, step_num*time_step, curr_forc.first);
            }
        );

        // update ocean salinity to 30 psu
        thermo_solver.UpdateOceanSalinity((NumType)30.0);

        // update atmosphere temperature
        thermo_solver.UpdateAtmosphereTemperature(curr_forc.first);

        // update precipitation rate
        thermo_solver.UpdatePrecipitationRate(curr_forc.second);

        // evaluation of thermo solver
        thermo_solver.Evaluate();
        
        // write mesh to file
        if (step_num % output_frequency == 0)
        {
#ifdef USE_JSON_OUTPUT
            ice_mesh->SaveJSON((std::string)"./" + output_ice_prefix, step_num);
            snow_mesh->SaveJSON((std::string)"./" + output_snow_prefix, step_num);
#endif
        }
    }

    // delete mesh
    delete ice_mesh; 
    delete snow_mesh;
}

// main function
int main()
{
    // model launcher (you can choose float or double)
    run_model<double>(3600.0,                        // time step (seconds) 
                      1,                          // number of time steps 
                      1,                             // output frequency N (every N-th step would be written to file) 
                      5,                             // number of uniform sigma-cells in ice 
                      2.0,                           // initial ice thickness (meters) 
                      0.1,                           // initial snow thickness (meters)
                      -15.0,                         // initial ice-snow interface temperature (deg Cel) 
                      -20.0,                         // initial snow surface temperature (deg Cel) 
                      350.0*3600.0,                  // winter duration (seconds)
                      200.0*3600.0,                  // spring duration (seconds) 
                      150.0*3600.0,                  // summer duration (seconds) 
                      200.0*3600.0,                  // autumn duration (seconds) 
                      -30.0,                         // lowest atm temperature (deg C)
                      10.0,                          // highest atm temperature (deg C)    
                      0.0,                           // lowest prec rate (m s-1)
                      5e-8,                          // highest prec rate (m s-1)
                      "ice_mesh",                    // output ice prefix
                      "snow_mesh",                   // output snow prefix
                      Kparam::BubblyBrine,           // ice conductivity parameterization
                      Cparam::SeaIce,                // ice effective capacity parameterization
                      Eparam::SeaIce,                // ice enthalpy parameterization
                      Lparam::SeaIce,                // ice heat of fusion parameterization
                      Aparam::MeltingFreezingIce,    // ice albedo parameterization
                      Kparam::FreshSnow,             // snow conductivity parameterization
                      Lparam::FreshSnow,             // snow heat of fusion parameterization
                      Aparam::MeltingFreezingSnow);  // snow albedo parameterization
}