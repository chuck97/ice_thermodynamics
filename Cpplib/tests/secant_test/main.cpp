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

    std::cout << "square sol: " << std::get<0>(square_res) << " , iters: " << std::get<2>(square_res) << std::endl;
    std::cout << "sin sol: " << std::get<0>(sin_res) << " , iters: " << std::get<2>(sin_res) << std::endl;
    return 0;
}