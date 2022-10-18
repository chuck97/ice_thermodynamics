#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include "defines.hpp"

namespace icethermo
{
    template <typename NumType> std::vector<NumType> operator+(const std::vector<NumType>& vec1,
                                                               const std::vector<NumType>& vec2);

    template <typename NumType> std::vector<NumType> operator-(const std::vector<NumType>& vec1,
                                                               const std::vector<NumType>& vec2);
    
    template <typename NumType> std::vector<NumType> operator*(const std::vector<NumType>& vec1,
                                                               const std::vector<NumType>& vec2);

    template <typename NumType> NumType L2_norm(const std::vector<NumType>& vec);

    template <typename NumType> std::vector<NumType> operator*(const std::vector<NumType>& vec, NumType scal);
    template <typename NumType> std::vector<NumType> operator*(NumType scal, const std::vector<NumType>& vec);

    template <typename NumType> std::vector<NumType> concatenate(const std::vector<NumType>& vec1, const std::vector<NumType>& vec2);
    template <typename NumType> std::vector<NumType> concatenate(const std::vector<std::vector<NumType>>& vecs);

    template <typename NumType> std::ostream& operator<<(std::ostream& os, const std::vector<NumType>& vec);

    template <typename NumType> NumType sum_vec(const std::vector<NumType>& vec);
    template <typename NumType> NumType sum_vec(const std::vector<NumType>& vec, int start, int end);
}