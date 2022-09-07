#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    // ### examples of Mesh class ###

    // constructor with given total thickness (default is 10 cells)
    Mesh<float> mesh1(1.0f);
    
    // how to create cells data (new variables should be references!!)
    auto& cells_temp =  mesh1.CreateCellData("cells_temperature");
    auto& cells_capacity = mesh1.CreateCellData("cells_capacity", true);
    auto& cells_enthalpy = mesh1.CreateCellData("cells_enthalpy", false);

    // one can modify cell data
    cells_temp[0] = 1.0f; cells_temp[1] = 2.0f; cells_temp[2] = 3.0f; 

    // how to create nodes data (new variables should be references!!)
    auto& nodes_k = mesh1.CreateNodeData("nodes_k");
    auto& nodes_enthalpy = mesh1.CreateNodeData("nodes_enthalpy", false);

    // one can modify nodes data
    nodes_k[0] = -5.0f; nodes_k.back() = -3.0f;

    // one can delete cell or node data 
    mesh1.DeleteCellData("cells_enthalpy");
    mesh1.DeleteCellData("nodes_enthalpy");

    // its is better to avoid this, but one can get another reference to created data
    auto& another_cells_temp = mesh1.GetCellData("cells_temperature");
    auto& another_nodes_k = mesh1.GetNodeData("nodes_k");
    another_cells_temp[0] = -5.0f;
    another_nodes_k[0] *= 2.0f;

    // one can get array with cell thicknesses
    std::cout << "current cell thickness array: " << mesh1.GetCellsThickness() << std::endl;

    // one can get total thickness
    std::cout << "current total cell thickness: " << mesh1.GetTotalThickness() << std::endl;

    // one could manually mute and unmute variables
    mesh1.MuteCellData("cells_temperature");
    mesh1.UnmuteCellData("cells_temperature");

    mesh1.MuteNodeData("nodes_enthalpy");
    mesh1.UnmuteNodeData("nodes_enthalpy");
    
    // one can save mesh to .txt file
    mesh1.SaveTXT("./mesh");

    // one can save mesh to .txt file with postfix number (relevant for time series)
    mesh1.SaveTXT("./mesh", 1488);

    // ### examples of another Mesh class constructor ###

    // construct uniform mesh with given cells num and total thickness
    Mesh<double> mesh2(15, 1.0);
    mesh2.SaveTXT("./mesh2");

    // construct arbitrary mesh with given unit segment partition and total thickness
    Mesh<double> mesh3({0.5, 0.5}, 5.0);
    mesh3.SaveTXT("./mesh3");
    
    // wrong constructor (it should be unit segment partition 0.5 + 0.4 != 1.0)
    Mesh<float> mesh4({0.5, 0.4}, 5.0);
    mesh4.SaveTXT("./mesh4");
}