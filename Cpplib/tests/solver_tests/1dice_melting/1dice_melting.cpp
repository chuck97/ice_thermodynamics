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

template<typename NumType>
void run_model(NumType time_step,
               int num_steps,
               int output_frequency,
               int num_cells,
               NumType initial_thickness,
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

    // freezing point of salty water with 30 psu
    NumType fusion_temp = GenConsts<NumType>::TempFusion(30.0);

    // initialize mandatory values
    for (int i = 0; i < n_cells; ++i)
    {
       (*initial_temp_cells)[i] = fusion_temp + (NumType)1.0*(i + 0.5)/(n_cells)*((NumType)-10.0 - fusion_temp);
       (*initial_sal_cells)[i] = 1.0 + 1.0*(i + 0.5)/(n_cells)*((NumType)4.0 - (NumType)1.0);
       (*initial_dens_cells)[i] = IceConsts<NumType>::rho_i;
    }

    (*initial_ice_surf_temp) = (NumType)-10.0;

    // assign melt time end in hours
    NumType melt_time_end = (NumType)200.0*time_step;

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
            [&melt_time_end, &step_num, &time_step](NumType temp)
            {
                return atm_flux<NumType>(temp, step_num*time_step, melt_time_end);
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
    run_model<double>(3600.0,                    // time step (seconds)
                      500,                       // number of time steps
                      1,                         // output frequency N (every N-th step would be written to file)
                      10,                        // number of uniform sigma-cells
                      2.0,                       // initial ice thickness (meters)
                      "ice_melting",             // output ice prefix
                      ApproxOrder::second,       // gradient approximation order
                      Kparam::BubblyBrine,       // conductivity parameterization
                      Cparam::SeaIce,            // effective capacity parameterization
                      Eparam::SeaIce,            // enthalpy parameterization
                      Lparam::SeaIce);           // heat of fusion parameterization
}