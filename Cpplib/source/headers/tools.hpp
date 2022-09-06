#pragma once
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <tuple>
#include "defines.hpp"

namespace icethermo
{
    template <typename NumType>
    using FuncPtr = NumType(*)(NumType);

    template <typename NumType>
    using FourVecs = std::tuple<std::vector<NumType>, std::vector<NumType>, std::vector<NumType>, std::vector<NumType>>;

    // 1-D secant solver
    template <typename NumType>
    std::pair<NumType, int> secant_solver(FuncPtr<NumType> func,
                                          NumType x_left, NumType x_right,
                                          NumType tol = 1e-9, int max_it = 1000);

    // SLAE Thomas solver
    template <typename NumType>
    std::vector<NumType> thomas_solver(const std::vector<NumType>& under,
                                       const std::vector<NumType>& diag,
                                       const std::vector<NumType>& over,
                                       const std::vector<NumType>& rhs);

    template <typename NumType>
    FourVecs<NumType> concat_matrices(const std::vector<NumType>& under_first,
                                      const std::vector<NumType>& diag_first,
                                      const std::vector<NumType>& over_first,
                                      const std::vector<NumType>& rhs_first,
                                      const std::vector<NumType>& under_second,
                                      const std::vector<NumType>& diag_second,
                                      const std::vector<NumType>& over_second,
                                      const std::vector<NumType>& rhs_second,
                                      const std::vector<NumType>& linker_);
}
