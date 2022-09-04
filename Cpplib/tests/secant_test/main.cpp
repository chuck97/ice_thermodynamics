#include <iostream>
#include <cmath>
#include "icethermo.hpp"

float square_func_float(float x)
{
    return x*x - 5.0;
}

float sin_func_float(float x)
{
    return std::sin(x);
}

double square_func_double(double x)
{
    return x*x - 5.0;
}

double sin_func_double(double x)
{
    return std::sin(x);
}

int main()
{
    auto square_res_float = secant_solver(square_func_float, 0.0f, 4.0f);
    auto sin_res_float = secant_solver(sin_func_float, 2.0f, 4.0f);
    auto square_res_double = secant_solver(square_func_double, 0.0, 4.0);
    auto sin_res_double = secant_solver(sin_func_double, 2.0, 4.0);

    std::cout << "square sol float: " << square_res_float.first << " , iters: " << square_res_float.second << std::endl;
    std::cout << "sin sol float: " << sin_res_float.first << " , iters: " << sin_res_float.second << std::endl;
    std::cout << "square sol double: " << square_res_double.first << " , iters: " << square_res_double.second << std::endl;
    std::cout << "sin sol double: " << sin_res_double.first << " , iters: " << sin_res_double.second << std::endl;
    return 0;
}
