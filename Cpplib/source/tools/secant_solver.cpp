#include "tools.hpp"

namespace icethermo
{
    template <typename NumType>
    std::pair<NumType, int> secant_solver(FuncPtr<NumType> func,
                                          NumType x_left, NumType x_right,
                                          NumType tol, int max_it)
    {
        
        NumType x_new = x_right;
        NumType x_old = x_left;
        int it = 0;
    
        while(std::abs(x_new - x_old) > tol && it < max_it)
        {
            it++;
            NumType temp = x_new;
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
    
    // explicit instantation for storing template functions declaration into .hpp file
    template std::pair<float, int> secant_solver(FuncPtr<float>, float, float, float, int);
    
    template std::pair<double, int> secant_solver(FuncPtr<double>, double, double, double, int);
}