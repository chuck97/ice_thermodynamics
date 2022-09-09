#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    // create test ice mesh with 10 uniform cells with total size of 2m
    Mesh<float>* ice_mesh = new(Mesh<float>)(10, 2.0);

    // make mandatory initial values
    auto initial_temp_cells = ice_mesh->CreateCellData("Ice cells temperature");
    auto initial_sal_cells = ice_mesh->CreateCellData("Ice cells salinity");
    auto initial_ice_surf_temp = ice_mesh->CreateSingleData("Ice surface temperature");
    auto initial_ice_base_temp = ice_mesh->CreateSingleData("Ice base temperature");
    int n_cells = ice_mesh->GetCellsNum();

    // fill initial values
    for (int i = 0; i < ice_mesh->GetCellsNum(); ++i)
    {
       (*initial_temp_cells)[i] = -10.0f;
       (*initial_sal_cells)[i] = 2.5f;
    }

    (*initial_ice_surf_temp) = -10.0f;
    (*initial_ice_base_temp) = -10.0f;

    // create test 1d ice solver class
    Ice1D_Solver<float> thermo_solver(ice_mesh);

    thermo_solver.Evaluate();

    // delete mesh
    delete ice_mesh; 
}