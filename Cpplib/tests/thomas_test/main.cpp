#include <iostream>
#include "icethermo.hpp"

int main()
{
    std::vector<float> under_float = {1.0f, 1.0f};
    std::vector<float> diag_float = {2.0f, 2.0f, 2.0f};
    std::vector<float> over_float = {1.0f, 1.0f};
    std::vector<float> rhs_float = {0.0f, 2.0f, 0.0f};

    std::vector<double> under_double = {-1.0, -1.0, -1.0, -1.0};
    std::vector<double> diag_double = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> over_double = {-1.0, -1.0, -1.0, -1.0};
    std::vector<double> rhs_double = {-3.0, 2.0, 0.0, 2.0, -3.0};

    std::vector<float> sol_float = thomas_solver<float>(under_float, diag_float, over_float, rhs_float);
    std::vector<double> sol_double = thomas_solver<double>(under_double, diag_double, over_double, rhs_double);

    // sol_float = [-1.0, 2.0, -1.0]
    std::cout << "sol_float = [";
    for (auto it = sol_float.begin(); it != sol_float.end(); it++)
    {
        std::cout << *it;
        if (std::next(it) != sol_float.end())
        {
            std::cout << ", ";
        }
    }
    std::cout << ']' << std::endl;

    // sol_double = [-1.0, 1.0, 1.0, 1.0, -1.0]
    std::cout << "sol_double = [";
    for (auto it = sol_double.begin(); it != sol_double.end(); it++)
    {
        std::cout << *it;
        if (std::next(it) != sol_double.end())
        {
            std::cout << ", ";
        }
    }
    std::cout << ']' << std::endl;
}
