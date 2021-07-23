#include "thermo_view.h"
#include "thomas_solver.h"
#include "output_operators.h"
#include "ice_thermo.h"

using namespace std;
using namespace INMOST;

int main()
{
	/*
	ThermoMesh t_m(5);
	t_m.SetWidth(0.2);
	t_m.AssignCellThickness({0.1, 0.2, 0.3, 0.4, 0.5});
	t_m.AssignCellData("temperature", {1.0, 2.0, 3.0, 4.0, 5.0});
	t_m.AssignNodeData("pressure", {-7.0, -8.0, -9.0, -10.0, -11.0, -12.0});
	t_m.PlotVertical("./vertical.vtu");
	t_m.PlotSigma("./sigma.vtu");
	*/

	IceThermo it(10, 0.0, 5.0, 1.0);
	it.SetTimeStep(1.0);
	it.WriteVTU("./", 0);
}
