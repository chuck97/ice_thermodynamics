#pragma once

#include <vector>
#include <iostream>
#include <sstream>
#include "vec_operators.h"
#include "vec_norms.h"
#include "thomas_solver.h"
#include "thermo_view.h"
#include "inmost.h"

#define IS_COEFFS_SIGMA false
#define N_PSEUDOITS 10

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
    const double sal_ice = 5.0;
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
    double dt_ = 3600.0;
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

    void RecalculateCoeffs(const std::vector<double>& dz_array);

    double CalculateE(double T, double S);
    double CalculateC(double T_last, double T_old, double S);
    double CalculateK(double T, double S);

    void InterpolateScalarNodes(const std::vector<double>& scalar_cells_array,
                                      std::vector<double>& scalar_nodes_array);

    void UpdateDeltaZ();

    void UpdateT();

    //TODO: combine functions for updating parameter in one

    void UpdateC(const std::vector<double>& T_old_array,
                 const std::vector<double>& T_new_array,
                 const std::vector<double>& S_array,
                       std::vector<double>& C_array);

    void UpdateE(const std::vector<double>& T_array,
                 const std::vector<double>& S_array,
                       std::vector<double>& E_array);

    void UpdateK(const std::vector<double>& T_array,
                 const std::vector<double>& S_array,
                       std::vector<double>& K_array);

    void UpdateW(const std::vector<double>& T_cells_array,
                 const std::vector<double>& S_cells_array,
                       std::vector<double>& w_nodes_array);

};

double C_to_K(double T);
double K_to_C(double T);
std::vector<double> C_to_K(const std::vector<double>& T_C_array);
std::vector<double> K_to_C(const std::vector<double>& T_K_array);
