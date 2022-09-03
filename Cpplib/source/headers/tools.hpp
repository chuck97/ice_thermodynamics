#pragma once
#include <vector>
#include <utility>
#include <string>
#include "defines.hpp"

// https://stackoverflow.com/questions/14848924/how-to-define-typedef-of-function-pointer-which-has-template-arguments
template <typename NumType>
using FuncPtr = NumType(*)(NumType);

// 1-D secant solver
template <typename NumType>
std::pair<NumType, int> secant_solver(FuncPtr<NumType> func,
                                      NumType x_left, NumType x_right,
                                      NumType tol = 1e-9, int max_it = 1000);

// SLAE Thomas solver
template <typename NumType>
std::vector<NumType> thomas_solver(const std::vector<NumType>& T_under,
                                   const std::vector<NumType>& T_diag,
                                   const std::vector<NumType>& T_over,
                                   const std::vector<NumType>& rhs);
