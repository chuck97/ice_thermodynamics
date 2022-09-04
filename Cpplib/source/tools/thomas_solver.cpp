#include "tools.hpp"

template <typename NumType>
std::vector<NumType> thomas_solver(const std::vector<NumType>& T_under,
                                   const std::vector<NumType>& T_diag,
                                   const std::vector<NumType>& T_over,
                                   const std::vector<NumType>& rhs)
{
    if (T_diag.size() != T_under.size() + 1
     || T_diag.size() != T_over.size() + 1
     || T_diag.size() != rhs.size())
    {
        THERMO_ERR("Sizes of under-diagonal (" + std::to_string(T_under.size())
                  + "), diagonal (" + std::to_string(T_diag.size())
                  + "), over-diagonal (" + std::to_string(T_over.size())
                  + ") and right-hand side ("  + std::to_string(rhs.size())
                  + ") arrays are not fit the matrix!")
    }

    std::vector<NumType> c_new(T_diag.size());
    std::vector<NumType> d_new(T_diag.size());
    std::vector<NumType> sol(T_diag.size());

    // forward iteration
    c_new[0] = T_over[0]/T_diag[0];
    for (int i = 0; i < c_new.size(); i++)
    {
        c_new[i] = T_over[i]/(T_diag[i] - T_under[i-1]*c_new[i-1]);
    }

    // backward iteration
    d_new[0] = rhs[0]/T_diag[0];
    for (int i = 0; i < d_new.size(); i++)
    {
        d_new[i] = (rhs[i] - T_under[i-1]*d_new[i-1])/(T_diag[i] - T_under[i-1]*c_new[i-1]);
    }

    // final solution
    sol.back() = d_new.back();
    for (int i = d_new.size() - 1; i >= 0; i--)
    {
        sol[i] = d_new[i] - c_new[i]*sol[i+1];
    }

    return sol;
}

// explicit instantation for storing template functions declaration into .hpp file

template
std::vector<float> thomas_solver(const std::vector<float>& T_under,
                                 const std::vector<float>& T_diag,
                                 const std::vector<float>& T_over,
                                 const std::vector<float>& rhs);

template
std::vector<double> thomas_solver(const std::vector<double>& T_under,
                                  const std::vector<double>& T_diag,
                                  const std::vector<double>& T_over,
                                  const std::vector<double>& rhs);
