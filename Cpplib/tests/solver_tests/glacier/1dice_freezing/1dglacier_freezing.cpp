#include "icethermo.hpp"
#include <string>

using namespace icethermo;

// stationary atmospheric heat flux
template<typename NumType>
NumType atm_flux(NumType temp)
{
    // simple sensible heat flux parameterization
    return AirConsts<NumType>::rho_a*
           AirConsts<NumType>::cp_a*
           GenConsts<NumType>::C_sh*
           (NumType)15.0*
           ((NumType)-1.0);
}

// model launcher
template <typename NumType>
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
    auto initial_dens_cells = ice_mesh->CreateCellsData("cells_density_array");
    auto initial_ice_surf_temp = ice_mesh->CreateSingleData("up_temperature");
    auto initial_ice_base_temp = ice_mesh->CreateSingleData("down_temperature");
    int n_cells = ice_mesh->GetCellsNum();

    // freezing point of salty water with 30 psu
    NumType fusion_temp = GenConsts<NumType>::TempFusion(30.0);

    // initialize mandatory values
    for (int i = 0; i < n_cells; ++i)
    {
       (*initial_temp_cells)[i] = 0.0 + (NumType)1.0*(i + (NumType)0.5)/(n_cells)*((NumType)-5.0 - 0.0);
       (*initial_dens_cells)[i] = IceConsts<NumType>::rho_i;
    }

    (*initial_ice_surf_temp) = (NumType)-5.0;
    (*initial_ice_base_temp) = (NumType)0.0;

    // create test 1d ice solver class
    Glacier1D_Solver<NumType> thermo_solver(ice_mesh,
                                            time_step,
                                            true,
                                            true,
                                            true,
                                            grad_approx_order,
                                            conductivity_parameterization,
                                            eff_capacity_parameterization,
                                            enthalpy_parameterization,
                                            fusion_heat_parameterization);

    // time stepping
    for (int step_num = 0; step_num < num_steps + 1; ++step_num)
    {
        // update atmospheric flux (it is not necessary in stationary case)
        thermo_solver.UpdateUpperFlux(atm_flux<NumType>);

        // evaluation of thermo solver
        thermo_solver.Evaluate();
        
        // write mesh to file
        if (step_num % output_frequency == 0)
        {
#ifdef USE_JSON_OUTPUT
            ice_mesh->SaveJSON((std::string)"./" + output_prefix, step_num);
#endif
        }
    }

    // delete mesh
    delete ice_mesh; 
}

// main function
int main()
{
    // model launcher (you can choose float or double)
    run_model<double>(3600.0,                // time step (seconds)
                      1000,                  // number of time steps
                      1,                     // output frequency N (every N-th step would be written to file)
                      10,                    // number of uniform sigma-cells
                      1.0,                   // initial ice thickness (meters)
                      "glacier_freezing",    // output ice prefix
                      ApproxOrder::second,   // gradient approximation order
                      Kparam::FreshIce,      // conductivity parameterization
                      Cparam::FreshIce,      // effective capacity parameterization
                      Eparam::FreshIce,      // enthalpy parameterization
                      Lparam::FreshIce);     // heat of fusion parameterization
}