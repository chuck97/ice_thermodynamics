#include <iostream>
#include <cmath>
#include "icethermo.hpp"

using namespace icethermo;

double square_func(double x)
{
    return x*x - 5.0;
}

float sin_func(float x)
{
    return std::sin(x);
}

int main()
{
    auto square_res = secant_solver<double>(square_func, 0.0, 4.0, 1e-1);
    auto sin_res = secant_solver<float>(sin_func, 2.0f, 4.0f);

    std::cout << "square sol: " << square_res.first << " , iters: " << square_res.second << std::endl;
    std::cout << "sin sol: " << sin_res.first << " , iters: " << sin_res.second << std::endl;
    return 0;
}