#include "icethermo.hpp"
#include <string>
#include <cmath>

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

template<typename NumType>
void run_model(NumType time_step,
               int num_steps,
               int output_frequency,
               int num_cells,
               NumType initial_thickness,
               NumType winter_duration, 
               NumType spring_duration, 
               NumType summer_duration, 
               NumType autumn_duration, 
               NumType lowest_temp,     
               NumType highest_temp,    
               const std::string& output_prefix,
               ApproxOrder grad_approx_order,
               Kparam conductivity_parameterization,
               Cparam eff_capacity_parameterization,
               Eparam enthalpy_parameterization,
               Lparam fusion_heat_parameterization)
{
    // create uniform sigma-mesh
    Mesh<NumType>* ice_mesh = new(Mesh<NumType>)(num_cells, initial_thickness);

    // create mandatory initial values on mesh
    auto initial_temp_cells = ice_mesh->CreateCellsData("cells_temperature_array");
    auto initial_sal_cells = ice_mesh->CreateCellsData("cells_salinity_array");
    auto initial_dens_cells = ice_mesh->CreateCellsData("cells_density_array");
    auto initial_ice_surf_temp = ice_mesh->CreateSingleData("up_temperature");
    int n_cells = ice_mesh->GetCellsNum();

    // salty water with 30 psu freezing temperature
    NumType fusion_temp = GenConsts<NumType>::TempFusion(30.0);

    // initialize mandatory values
    for (int i = 0; i < n_cells; ++i)
    {
       (*initial_temp_cells)[i] = fusion_temp + (NumType)1.0*(i + 0.5)/(n_cells)*((NumType)-10.0 - fusion_temp);
       (*initial_sal_cells)[i] = 1.0 + 1.0*(i + 0.5)/(n_cells)*((NumType)4.0 - (NumType)1.0);
       (*initial_dens_cells)[i] = IceConsts<NumType>::rho_i;
    }

    (*initial_ice_surf_temp) = (NumType)-10.0;

    // create 1dice solver class
    SeaIce1D_Solver<NumType> thermo_solver(ice_mesh,
                                           time_step, 
                                           grad_approx_order,
                                           conductivity_parameterization,
                                           eff_capacity_parameterization,
                                           enthalpy_parameterization,
                                           fusion_heat_parameterization);
    
    // update ocean salinity to 30 psu
    thermo_solver.UpdateOceanSalinity((NumType)30.0);

    // time stepping
    for (int step_num = 0; step_num < num_steps + 1; ++step_num)
    {
        // update atmospheric flux
        thermo_solver.UpdateForcing
        (
            [&step_num,
             &time_step,
             &winter_duration, 
             &spring_duration, 
             &summer_duration, 
             &autumn_duration,
             &lowest_temp,
             &highest_temp](NumType temp)
            {
                return atm_flux<NumType>(temp, step_num*time_step, annual_temp_prec(winter_duration,
                                                                                    spring_duration,
                                                                                    summer_duration,
                                                                                    autumn_duration,
                                                                                    lowest_temp,
                                                                                    highest_temp,
                                                                                    (NumType)0.0,
                                                                                    (NumType)0.0,
                                                                                    step_num*time_step).first);
            }
        );

        // evaluation of thermo solver
        thermo_solver.Evaluate();
        
        // write mesh to file
        if (step_num % output_frequency == 0)
        {
            ice_mesh->SaveJSON(std::string("./") + output_prefix, step_num);
        }
    }

    // delete mesh
    delete ice_mesh; 
}

int main()
{
    // model launcher (you can choose float or double)
    run_model<double>(3600.0,                       // time step (seconds)
                      2000,                         // number of time steps
                      1,                            // output frequency N (every N-th step would be written to file)
                      10,                           // number of uniform sigma-cells
                      2.0,                          // initial ice thickness (meters)
                      350.0*3600.0,                 // winter duration (seconds)
                      200.0*3600.0,                 // spring duration (seconds) 
                      150.0*3600.0,                 // summer duration (seconds) 
                      200.0*3600.0,                 // autumn duration (seconds) 
                      -30.0,                        // lowest atm temperature (deg C)     
                      10.0,                         // highest atm temperature (deg C)    
                      "ice_cycle",                  // output ice prefix
                      ApproxOrder::second,          // gradient approximation order
                      Kparam::BubblyBrine,          // conductivity parameterization
                      Cparam::SeaIce,               // effective capacity parameterization
                      Eparam::SeaIce,               // enthalpy parameterization
                      Lparam::SeaIce);              // heat of fusion parameterization
}