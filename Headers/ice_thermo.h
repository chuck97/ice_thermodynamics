#pragma once

#include <vector>
#include "thomas_solver.h"
#include "thermo_view.h"

#define IS_COEFFS_SIGMA true
#define ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}

class IceThermo
{
public:
    IceThermo(unsigned int n_cells,
              double T_b,
              double T_su,
              double h);

    IceThermo(unsigned int n_cells,
              const std::vector<double>& T_init_cells,
              const std::vector<double>& dz_init_cells);

    void SetTimeStep(double time_step);
    void UpdateFluxes(double flux_boundary, double flux_surface);
    void Evaluate();
    void WriteVTU();

private:
    double dt_;
    double F_b;
    double F_su;
    const unsigned int n_cells_;
    const unsigned int n_nodes_;
    std::vector<double> T_cells_old_;
    std::vector<double> T_cells_last_;
    std::vector<double> T_cells_new_;
    std::vector<double> T_nodes_old_;
    std::vector<double> T_nodes_last_;
    std::vector<double> T_nodes_new_;
    std::vector<double> dz_cells_old_;
    std::vector<double> dz_cells_new_;
    std::vector<double> w_nodes_;
    std::vector<double> c_cells_;
    std::vector<double> c_nodes_;
    std::vector<double> k_nodes_;
    std::vector<double> E_cells_;
    std::vector<double> E_nodes_;
    std::vector<double> dR_cells_;
    std::vector<double> a_nodes_;
    std::vector<double> b_nodes_;
    std::vector<double> S_nodes_;

    void UpdateDeltaZ();
    void UpdateCoefs();
    void UpdateT();
    void UpdateW();
};

struct Consts
{
    const double rho_ice;
    const double c0_ice;
    const double cw_water;
    const double L0_ice;
    const double mu_ice;
    const double k0_ice;
    const double k1_ice;
};
