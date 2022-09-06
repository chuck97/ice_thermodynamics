#include <iostream>
#include <cassert>
#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    std::vector<float> under_float1 = {1.0f, 1.0f};
    std::vector<float> diag_float1 = {2.0f, 2.0f, 2.0f};
    std::vector<float> over_float1 = {1.0f, 1.0f};
    std::vector<float> rhs_float1 = {0.0f, 2.0f, 0.0f};
    std::vector<float> exact_sol_float1 = {-1.0, 2.0, -1.0};
    std::vector<float> sol_float1 = thomas_solver(under_float1, diag_float1, over_float1, rhs_float1);
    std::cout << "exact1: " << exact_sol_float1 << ", computed1: " << sol_float1 << std::endl;

    std::vector<double> under_double1 = {-1.0, -1.0, -1.0, -1.0};
    std::vector<double> diag_double1 = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> over_double1 = {-1.0, -1.0, -1.0, -1.0};
    std::vector<double> rhs_double1 = {-3.0, 2.0, 0.0, 2.0, -3.0};
    std::vector<double> exact_sol_double1 = {-1.0, 1.0, 1.0, 1.0, -1.0};
    std::vector<double> sol_double1 = thomas_solver(under_double1, diag_double1, over_double1, rhs_double1);
    std::cout << "exact2: " << exact_sol_double1 << ", computed2: " << sol_double1 << std::endl;

    std::vector<float> under_float2 = {3.0f};
    std::vector<float> diag_float2 = {1.0f, 4.0f};
    std::vector<float> over_float2 = {-2.0f};
    std::vector<float> rhs_float2 = {0.0f, -5.0f};
    std::vector<float> exact_sol_float2 = {-1.0f, -0.5f};
    std::vector<float> sol_float2 = thomas_solver(under_float2, diag_float2, over_float2, rhs_float2);
    std::cout << "exact3: " << exact_sol_float2 << ", computed3: " << sol_float2 << std::endl;

    std::vector<double> under_double2 = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> diag_double2 = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> over_double2 = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> rhs_double2 = {1.0, -1.0, 1.0, -1.0, 1.0};
    std::vector<double> exact_sol_double2 = {1.0/1.0, -1.0/2.0, 1.0/3.0, -1.0/4.0, 1.0/5.0};
    std::vector<double> sol_double2 = thomas_solver(under_double2, diag_double2, over_double2, rhs_double2);
    std::cout << "exact4: " << exact_sol_double2 << ", computed4: " << sol_double2 << std::endl;
    
    std::cout << "Thomas solver tests finished!\n";
}