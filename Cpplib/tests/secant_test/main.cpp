#include <iostream>
#include <cmath>
#include "icethermo.hpp"

double square_func(double x)
{
    return x*x - 5.0;
}

double sin_func(double x)
{
    return std::sin(x);
}

int main()
{
    auto square_res = secant_solver(square_func, 0.0, 4.0);
    auto sin_res = secant_solver(sin_func, 2.0, 4.0);

    std::cout << "square sol: " << square_res.first << " , iters: " << square_res.second << std::endl;
    std::cout << "sin sol: " << sin_res.first << " , iters: " << sin_res.second << std::endl;
    return 0;
}