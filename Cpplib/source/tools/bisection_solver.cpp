#include "tools.hpp"

namespace icethermo
{
    template <typename NumType>
    std::tuple<NumType, NumType, int, bool> bisection_solver(FuncPtr<NumType> func,
                                                             NumType x_left, NumType x_right,
                                                             NumType tol, int max_it)
    {
        NumType x_RIGHT = x_right;
        NumType x_LEFT = x_left;
        NumType x_MIDDLE = 0.0;
        NumType f_MIDDLE = (NumType)0.0;
        int iteration = 0;
    
        for (int it = 0; it < max_it; ++it)
        {
            iteration++;
            x_MIDDLE = (NumType)(0.5*(x_RIGHT + x_LEFT));
            f_MIDDLE = (NumType)func(x_MIDDLE);

            if (f_MIDDLE <= 0)
            {
                x_RIGHT = x_MIDDLE;
            }
            else
            {
                x_LEFT = x_MIDDLE;
            }
            std::cout << x_MIDDLE << " " << f_MIDDLE << std::endl;
            if (std::abs(x_RIGHT - x_LEFT) < tol)
            {
                break;
            }
        }
        std::cout << "result:" << x_MIDDLE << std::endl;
        return {x_MIDDLE, f_MIDDLE, iteration, (iteration != max_it) ? true: false};
    }
    
    // explicit instantation for storing template functions declaration into .hpp file
    template std::tuple<float, float, int, bool> bisection_solver(FuncPtr<float>, float, float, float, int);
    
    template std::tuple<double, double, int, bool> bisection_solver(FuncPtr<double>, double, double, double, int);
}