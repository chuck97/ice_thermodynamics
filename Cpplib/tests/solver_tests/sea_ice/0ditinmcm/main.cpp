#include "icethermo.hpp"
#include <string>

using namespace icethermo;

// model launcher
template <typename NumType>
void run_model(NumType time_step,
               int num_steps,
               NumType initial_ice_thickness,
               NumType initial_snow_thickness,
               NumType initial_interface_temperature,
               NumType initial_snow_surface_temperature,   
               NumType precipitation_test,
               NumType salinity_ocean_test,
               NumType sw_radiation_test,
               NumType lw_radiation_test,
               NumType atm_temp_test,
               NumType spec_humid_test,
               NumType sh_coeff_test,
               NumType lh_coeff_test,
               NumType abs_wind_speed_test,
               NumType atm_press_test,
               NumType ocean_flux_test)
{
    // create 1-layer sigma-mesh for 0d ice
    Mesh<NumType>* ice_mesh = new(Mesh<NumType>)(1, initial_ice_thickness);

    // create mandatory initial values for ice
    auto initial_ice_temp_cells = ice_mesh->CreateCellsData("cells_temperature_array");
    auto initial_ice_sal_cells = ice_mesh->CreateCellsData("cells_salinity_array");
    auto initial_ice_dens_cells = ice_mesh->CreateCellsData("cells_density_array");
    auto initial_ice_surf_temp = ice_mesh->CreateSingleData("up_temperature");

    // create 1-layer sigma-mesh for 0d snow
    Mesh<NumType>* snow_mesh = new(Mesh<NumType>)(1, initial_snow_thickness);

    // create mandatory initial values for snow
    auto initial_snow_temp_cells = snow_mesh->CreateCellsData("cells_temperature_array");
    auto initial_snow_dens_cells = snow_mesh->CreateCellsData("cells_density_array");
    auto initial_snow_surf_temp = snow_mesh->CreateSingleData("up_temperature");

    // freezing point of salty water with given psu
    NumType fusion_temp = GenConsts<NumType>::TempFusion(salinity_ocean_test);

    // initialize mandatory values for 0d-ice
    (*initial_ice_temp_cells)[0] = (NumType)0.5*(fusion_temp + initial_interface_temperature);
    (*initial_ice_sal_cells)[0] = 2.5;
    (*initial_ice_dens_cells)[0] = IceConsts<NumType>::rho_i;
    (*initial_ice_surf_temp) = initial_interface_temperature;

    // initialize mandatory values for 0d-snow
    (*initial_snow_temp_cells)[0] = (NumType)0.5*(initial_snow_surface_temperature + initial_interface_temperature);
    (*initial_snow_dens_cells)[0] = SnowConsts<NumType>::rho_s;
    (*initial_snow_surf_temp) = initial_snow_surface_temperature;


    // create 0d ice with 0d snow solver class
    SeaIce0D_Snow0D_Solver<NumType> thermo_solver(ice_mesh,
                                                  snow_mesh,
                                                  time_step,
                                                  true,
                                                  true,
                                                  Kparam::BubblyBrine, 
                                                  Cparam::SeaIce,      
                                                  Eparam::SeaIce,      
                                                  Lparam::SeaIce,      
                                                  Kparam::FreshSnow,   
                                                  Lparam::FreshSnow);
    // save initial state to file
#ifdef USE_JSON_OUTPUT
    ice_mesh->SaveJSON((std::string)"./" + "ice", 0);
    snow_mesh->SaveJSON((std::string)"./" + "snow", 0);
#endif    

    // print initial state
    std::cout << "INITIAL VALUES " << std::endl;
    std::cout << "ICE THICKNESS: " << ice_mesh->GetTotalThickness() << std::endl;
    std::cout << "SNOW THICKNESS: " << snow_mesh->GetTotalThickness() << std::endl;
    std::cout << "SNOW SURFACE TEMPERATURE: " << *(snow_mesh->GetSingleData("up_temperature")) << std::endl;
    std::cout << std::endl;

    // time stepping
    for (int step_num = 1; step_num < num_steps + 1; ++step_num)
    {
        // update atmosphere temperature
        thermo_solver.UpdateAtmosphereTemperature(atm_temp_test);

        // update atmosphere pressure
        thermo_solver.UpdateAtmospherePressure(atm_press_test);

        // update precipitation rate
        thermo_solver.UpdatePrecipitationRate(precipitation_test);

        // update atm humidity
        thermo_solver.UpdateAirSpecificHumidity(spec_humid_test);

        // update wind speed
        thermo_solver.UpdateAbsWindSpeed(abs_wind_speed_test);

        // update downwards long-wawe radiation
        thermo_solver.UpdateLongWaveRadiation
        (
            [lw_radiation_test](NumType T){return lw_radiation_test;}
        );

        // update downwards short-wawe radiation
        thermo_solver.UpdateShortWaveRadiation
        (
            [sw_radiation_test](NumType T){return sw_radiation_test;}
        );

        // update sh coeff
        thermo_solver.UpdateShTransCoeff(sh_coeff_test);

        // update lh coeff
        thermo_solver.UpdateLhTransCoeff(lh_coeff_test);

        // update default atm flux parameterization
        thermo_solver.UpdateUpperFlux();

        // updatee ocean flux
        thermo_solver.UpdateLowerFlux
        (
            [ocean_flux_test](NumType T){return ocean_flux_test;}
        );

        // evaluation of thermo solver
        thermo_solver.Evaluate();

        // write mesh to file
#ifdef USE_JSON_OUTPUT
            ice_mesh->SaveJSON((std::string)"./" + "ice", step_num);
            snow_mesh->SaveJSON((std::string)"./" + "snow", step_num);
#endif
        
        // print current state
        std::cout << "STEP " << step_num << std::endl;
        std::cout << "ICE THICKNESS: " << ice_mesh->GetTotalThickness() << std::endl;
        std::cout << "SNOW THICKNESS: " << snow_mesh->GetTotalThickness() << std::endl;
        std::cout << "SNOW SURFACE TEMPERATURE: " << *(snow_mesh->GetSingleData("up_temperature")) << std::endl;
        std::cout << std::endl;

        // update default atm flux parameterization
        thermo_solver.UpdateUpperFlux();
    }

    // delete mesh
    delete ice_mesh; 
    delete snow_mesh;
}

// main function
int main()
{
    // model launcher (you can choose float or double)
    run_model<double>(1800.0,                             // time step (seconds) 
                      1,                                  // number of steps 
                      2.76812768470955,                   // initial ice thickness (meters) 
                      0.205363345324322,                  // initial snow thickness (meters)
                      -5.0,                               // initial ice-snow interface temperature (deg Cel) 
                      -10,                                // initial snow surface temperature (deg Cel) 
                      1.852127370040868E-009,             // precipitation rate in water equivalent (m s-1)
                      28.8696956634521,                   // ocean salinity (psu)
                      0.0,                                // short-wave radiation from atm (W m-2)
                      160.093663854041,                   // long-wave radiation from atm (W m-2)
                      -38.1966667175293,                  // atmosphere temperature (deg Cel)
                      9.191176650347188E-002,             // atmosphere specific humidity (g kg-1)
                      1.249765278771520E-003,             // sensible heat tranport coefficient (-)
                      1.249765278771520E-003,             // latent heat tranport coefficient (-)
                      2.37398432193731,                   // abs wind speed (m s-1)
                      100510.515625000,                   // atmosphere pressure (Pa)
                      0.0                                 // total flux from ocean (W m-2)
                      );
}