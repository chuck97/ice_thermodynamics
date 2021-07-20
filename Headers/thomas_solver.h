#pragma once
#include <vector>
#include <iostream>

#define ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}

// lhs is vector of 3 diagonals: ld, c, ru
std::vector<double> thomas_solver(const std::vector<std::vector<double>>& lhs,
                                  const std::vector<double>& rhs);