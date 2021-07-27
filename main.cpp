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
    double abs_ua = 1.0;
    double T_surf = C_to_K(0.0);

    return (rho_air*c_air*C_tr_air*abs_ua*(C_to_K(AirTemp(t)) - T_surf));
}

int main()
{
	double t = 0.0;
	double time_step = 3600.0;
    IceThermo proc(10, C_to_K(0.0), C_to_K(AirTemp(t)), 1.0);
    proc.SetTimeStep(time_step);
    proc.WriteVTU("./", 0);

    for (unsigned int n_step = 1; n_step < 10; ++n_step)
	{
        t += time_step;
        proc.UpdateFluxes(BottomFlux(t), SurfaceFlux(t));
        proc.Evaluate();
        proc.WriteVTU("./", n_step);
	}
}
