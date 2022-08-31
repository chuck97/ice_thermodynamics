#include "tools.hpp"

std::pair<double, int> secant_solver(FuncPtr func,
                                     double x_left, double x_right,
                                     double tol, int max_it)
{
    double x_new = x_right;
    double x_old = x_left;
    int it = 0;

    while(std::abs(x_new - x_old) > tol && it < max_it)
    {
        it++;
        double temp = x_new;
        x_new -= func(x_new)*(x_new - x_old) / (func(x_new) - func(x_old));
        x_old = temp;
    }

    if (it == max_it)
    {
        THERMO_ERR("Secant solver doesn't converge! Exceed " 
                   + std::to_string(max_it)
                   + " iterations.");
    }

    return {x_new, it};
}