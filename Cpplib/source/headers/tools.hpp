#pragma once
#include <vector>
#include <utility>
#include <string>
#include "defines.hpp"


typedef double(*FuncPtr)(double);

// 1-D secant solver
std::pair<double, int> secant_solver(FuncPtr func,
                                     double x_left, double x_right,
                                     double tol = 1e-9, int max_it = 1000);

// SLAE Thomas solver
std::vector<double> thomas_solver(const std::vector<double>& T_under,
                                  const std::vector<double>& T_diag,
                                  const std::vector<double>& T_over,
                                  const std::vector<double>& rhs);