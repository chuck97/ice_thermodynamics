#include "tools.hpp"

namespace icethermo
{
    template <typename NumType>
    std::vector<NumType> thomas_solver(const std::vector<NumType>& under,
                                       const std::vector<NumType>& diag,
                                       const std::vector<NumType>& over,
                                       const std::vector<NumType>& rhs)
    {
        if (diag.size() != under.size() + 1
         || diag.size() != over.size() + 1
         || diag.size() != rhs.size())
        {
            THERMO_ERR("Sizes of under-diagonal (" + std::to_string(under.size())
                      + "), diagonal (" + std::to_string(diag.size())
                      + "), over-diagonal (" + std::to_string(over.size())
                      + ") and right-hand side ("  + std::to_string(rhs.size())
                      + ") arrays are not fit the matrix!")
        }

        std::vector<NumType> c_new(diag.size());
        std::vector<NumType> d_new(diag.size());
        std::vector<NumType> sol(diag.size());

        // forward iteration
        c_new[0] = over[0]/diag[0];
        for (int i = 1; i < c_new.size(); i++)
        {
            c_new[i] = over[i]/(diag[i] - under[i-1]*c_new[i-1]);
        }

        // backward iteration
        d_new[0] = rhs[0]/diag[0];
        for (int i = 1; i < d_new.size(); i++)
        {
            d_new[i] = (rhs[i] - under[i-1]*d_new[i-1])/(diag[i] - under[i-1]*c_new[i-1]);
        }

        // final solution
        sol.back() = d_new.back();
        for (int i = d_new.size() - 2; i >= 0; i--)
        {
            sol[i] = d_new[i] - c_new[i]*sol[i+1];
        }

        return sol;
    }

    // explicit instantation for storing template functions declaration into .hpp file
    template std::vector<float> thomas_solver(const std::vector<float>& under,
                                              const std::vector<float>& diag,
                                              const std::vector<float>& over,
                                              const std::vector<float>& rhs);

    template std::vector<double> thomas_solver(const std::vector<double>& under,
                                               const std::vector<double>& diag,
                                               const std::vector<double>& over,
                                               const std::vector<double>& rhs);
}