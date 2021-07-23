#pragma once

#include <vector>
#include <iostream>
#include <sstream>
#include "thomas_solver.h"
#include "thermo_view.h"

#define IS_COEFFS_SIGMA true
#define ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}

struct Consts
{
    const double rho_ice = 917.0;
    const double c0_ice = 2108.0;
    const double cw_water = 4184.0;
    const double L0_ice = 333500.0;
    const double mu_ice = 0.054;
    const double k0_ice = 2.03;
    const double k1_ice = 0.117;
    const double sal_ice = 0.005;
};

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
    void WriteVTU(const std::string& path,
                  unsigned int step_num) const;

private:
    double dt_ = 0.0;
    double F_b_ = 0.0;
    double F_su_ = 0.0;
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
    std::vector<double> S_cells_;

    Consts consts;

    void UpdateDeltaZ();
    void UpdateCoeffs();
    void UpdateT();
    void UpdateW();

    double CalculateE(double T, double S);
    double CalculateK(double T, double S);
};
