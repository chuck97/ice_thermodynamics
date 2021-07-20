#include "thermo_view.h"
#include "thomas_solver.h"
#include "output_operators.h"

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


	/*
	1 0 0  x0     10
	0 1 0  x1  =  20
	0 0 1  x2     30

	sol: [10 20 30]^T
	*/
	vector<double> ld = {0.0, 0.0};
	vector<double> c =  {1.0, 1.0, 1.0};
	vector<double> ru = {0.0, 0.0};

	vector<double> rhs = {10.0, 20.0, 30.0};
	
	cout << thomas_solver({ld, c, ru}, rhs) << endl;

	/*
	1 3  0  x0     10
	2 1 -5  x1  =  20
	0 2  1  x2     30

	sol: [-80 30 -30]^T
	*/
	ld = {2.0, 2.0};
	c  = {1.0, 1.0, 1.0};
	ru = {3.0, -5.0};

	rhs = {10.0, 20.0, 30.0};
	
	cout << thomas_solver({ld, c, ru}, rhs) << endl;

	/*
	 1 2 x0   =  10
	-3 4 x1      20

	sol: [0 5]^T
	*/
	ld = {-3};
	c  = {1.0, 4.0};
	ru = {2.0};

	rhs = {10.0, 20.0};
	
	cout << thomas_solver({ld, c, ru}, rhs) << endl;

	/*
	1 3  0  0  0  x0     10
	2 1 -5  0  0  x1     20
	0 2  1  2  0  x2  =  30
	0 0  3  6  4  x3     40
	0 0  0 -7 -2  x5    -50

	sol: [-170 60 -60 -15 77.5]^T
	*/
	ld = {2.0, 2.0, 3.0, -7.0};
	c  = {1.0, 1.0, 1.0, 6.0, -2.0};
	ru = {3.0, -5.0, 2.0, 4.0};

	rhs = {10.0, 20.0, 30.0, 40.0, -50.0};
	
	cout << thomas_solver({ld, c, ru}, rhs) << endl;
}
