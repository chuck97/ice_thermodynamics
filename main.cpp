#include "ice_thermo.h"
#include <cmath>

using namespace std;
using namespace INMOST;

double AirTemp(double t)
{
	return -15.0 + sin(M_PI*t/(2.0*24.0*3600.0));
}

double BottomFlux(double t)
{
	return 0.0;
}

double SurfaceFlux(double t)
{
	double rho_air = 1.28;
	double c_air = 1.01e3;
	double C_tr_air = 1.0e-3;
	double abs_ua = 5.0;
	double T_surf = 0.0;

	return (rho_air*c_air*C_tr_air*abs_ua*(AirTemp(t) - T_surf));
}

int main()
{
	double t = 0.0;
	double time_step = 3600.0;
	IceThermo it(10, -0.5, AirTemp(t), 1.0);
	it.SetTimeStep(time_step);
	it.WriteVTU("./", 0);

	for (unsigned int n_step = 1; n_step < 240; ++n_step)
	{
		t+= time_step;
		it.UpdateFluxes(BottomFlux(t), SurfaceFlux(t));
		it.Evaluate();
		it.WriteVTU("./", n_step);
	}
}
