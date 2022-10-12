#include "icethermo.hpp"

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
           ((NumType)-30.0 - temp);
}

int main()
{
    // create test ice mesh with 10 uniform cells with total size of 2 meters
    Mesh<float>* ice_mesh = new(Mesh<float>)(10, 2.0);

    // create mandatory initial values on mesh
    auto initial_temp_cells = ice_mesh->CreateCellsData("cells_temperature_array");
    auto initial_sal_cells = ice_mesh->CreateCellsData("cells_salinity_array");
    auto initial_dens_cells = ice_mesh->CreateCellsData("cells_density_array");
    auto initial_ice_surf_temp = ice_mesh->CreateSingleData("up_temperature");
    auto initial_ice_base_temp = ice_mesh->CreateSingleData("down_temperature");
    auto initial_ocn_salinity = ice_mesh->CreateSingleData("ocean_salinity", false);
    int n_cells = ice_mesh->GetCellsNum();

    float fusion_temp = GenConsts<float>::TempFusion(30.0);

    // initialize mandatory values
    for (int i = 0; i < n_cells; ++i)
    {
       (*initial_temp_cells)[i] = fusion_temp + 1.0f*(i + 0.5f)/(n_cells)*(-10.0f - fusion_temp);
       (*initial_sal_cells)[i] = 1.0f + 1.0f*(i + 0.5f)/(n_cells)*(4.0f - 1.0f);
       (*initial_dens_cells)[i] = IceConsts<float>::rho_i;
    }

    (*initial_ice_surf_temp) = -10.0f;
    (*initial_ice_base_temp) = fusion_temp;
    (*initial_ocn_salinity) = 30.0;

    // create test 1d ice solver class
    SeaIce1D_Solver<float> thermo_solver(ice_mesh, 3600.0, ApproxOrder::second);

    // time stepping
    for (int step_num = 0; step_num < 251; ++step_num)
    {
        // update atmospheric flux (it is not necessary in stationary case)
        thermo_solver.UpdateForcing(atm_flux<float>);

        // evaluation of thermo solver
        thermo_solver.Evaluate();
        
        // write mesh to file
        if (step_num % 10 == 0)
        {
            ice_mesh->SaveJSON("./ice_freezing", step_num);
            std::cout << "ice thickness: " << ice_mesh->GetTotalThickness() << std::endl;
        }
    }

    // delete mesh
    delete ice_mesh; 
}