#include "tools.hpp"

namespace icethermo
{
    template <typename NumType>
    std::tuple<NumType, NumType, int, bool> secant_solver(FuncPtr<NumType> func,
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

            // handle zero denominator
            if (std::abs(func(x_new) - func(x_old)) < REAL_MIN_VAL(NumType))
            {
                return {x_new, func(x_new), it, false};
            }
            else
            {
                x_new -= func(x_new)*(x_new - x_old) / (func(x_new) - func(x_old));
            }

            x_old = temp;
        }
        
        return {x_new, func(x_new), it, (it != max_it) ? true : false};
    }
    
    // explicit instantation for storing template functions declaration into .hpp file
    template std::tuple<float, float, int, bool> secant_solver(FuncPtr<float>, float, float, float, int);
    
    template std::tuple<double, double, int, bool> secant_solver(FuncPtr<double>, double, double, double, int);
}