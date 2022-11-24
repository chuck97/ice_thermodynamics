#include "icethermo.hpp"
#include <string>

using namespace icethermo;

// linearly warming atmosphere
template<typename NumType>
NumType T_atm(NumType time, NumType melt_time_end)
{
    return (time < melt_time_end) ? (NumType)-20.0 + (NumType)30.0*time/melt_time_end : (NumType)10.0;
}

// non-stationary atmospheric flux function
template<typename NumType>
NumType atm_flux(NumType temp, NumType time, NumType melt_time_end)
{
    // simple sensible heat flux parameterization
    return AirConsts<NumType>::rho_a *
           AirConsts<NumType>::cp_a *
           GenConsts<NumType>::C_sh *
           (NumType)15.0 * 
           ((NumType)T_atm(time, melt_time_end) - temp);
}

// model launcher
template <typename NumType>
void run_model(NumType time_step,
               int num_steps,
               int output_frequency,
               int num_cells,
               NumType precipitation_rate,
               NumType initial_ice_thickness,
               NumType initial_snow_thickness,
               NumType initial_interface_temperature,
               NumType initial_snow_surface_temperature,
               const std::string& output_ice_prefix,
               const std::string& output_snow_prefix,
               Kparam conductivity_ice_parameterization,
               Cparam eff_capacity_ice_parameterization,
               Eparam enthalpy_ice_parameterization,
               Lparam fusion_heat_ice_parameterization,
               Kparam conductivity_snow_parameterization,
               Lparam fusion_heat_snow_parameterization)
{
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
    NumType fusion_temp = GenConsts<NumType>::TempFusion(30.0);

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

    // assign melt time end in hours
    NumType melt_time_end = (NumType)200.0*time_step;

    // create 1d ice with 0d snow solver class
    SeaIce1D_Snow0D_Solver<NumType> thermo_solver(ice_mesh,
                                                  snow_mesh,
                                                  time_step,
                                                  true,
                                                  true,
                                                  conductivity_ice_parameterization,
                                                  eff_capacity_ice_parameterization,
                                                  enthalpy_ice_parameterization,
                                                  fusion_heat_ice_parameterization,
                                                  conductivity_snow_parameterization,
                                                  fusion_heat_snow_parameterization);

    // save initial state to file
    ice_mesh->SaveJSON((std::string)"./" + output_ice_prefix, 0);
    snow_mesh->SaveJSON((std::string)"./" + output_snow_prefix, 0);

    // time stepping
    for (int step_num = 1; step_num < num_steps + 1; ++step_num)
    {
        // update atmospheric flux (it is not necessary in stationary case)
        thermo_solver.UpdateUpperFlux
        (
            [&melt_time_end, &step_num, &time_step](NumType temp)
            {
                return atm_flux<NumType>(temp, step_num*time_step, melt_time_end);
            }
        );

        // update ocean salinity to 30 psu
        thermo_solver.UpdateOceanSalinity((NumType)30.0);

        // update atmosphere temperature
        thermo_solver.UpdateAtmosphereTemperature(T_atm(step_num*time_step, melt_time_end));

        // update precipitation rate
        thermo_solver.UpdatePrecipitationRate(precipitation_rate);

        // evaluation of thermo solver
        thermo_solver.Evaluate();
        
        // write mesh to file
        if (step_num % output_frequency == 0)
        {
            ice_mesh->SaveJSON((std::string)"./" + output_ice_prefix, step_num);
            snow_mesh->SaveJSON((std::string)"./" + output_snow_prefix, step_num);
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
    run_model<double>(3600.0,              // time step (seconds)
                      500,                 // number of time steps
                      1,                   // output frequency N (every N-th step would be written to file)
                      20,                  // number of uniform sigma-cells in ice
                      5e-8,                // precipitation rate (m s-1)
                      2.0,                 // initial ice thickness (meters)
                      0.1,                 // initial snow thickness (meters)
                      -15.0,               // initial ice-snow interface temperature (deg Cel)
                      -20.0,               // initial snow surface temperature (deg Cel)
                      "ice_mesh",          // output ice prefix
                      "snow_mesh",         // output snow prefix
                      Kparam::BubblyBrine, // ice conductivity parameterization
                      Cparam::SeaIce,      // ice effective capacity parameterization
                      Eparam::SeaIce,      // ice enthalpy parameterization
                      Lparam::SeaIce,      // ice heat of fusion parameterization
                      Kparam::FreshSnow,   // snow conductivity parameterization
                      Lparam::FreshSnow);  // snow heat of fusion parameterization
}