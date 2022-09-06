#include "tools.hpp"

template <typename NumType>
FourVecs<NumType> concat_matrices(const std::vector<NumType>& under_first,
                                  const std::vector<NumType>& diag_first,
                                  const std::vector<NumType>& over_first,
                                  const std::vector<NumType>& rhs_first,
                                  const std::vector<NumType>& under_second,
                                  const std::vector<NumType>& diag_second,
                                  const std::vector<NumType>& over_second,
                                  const std::vector<NumType>& rhs_second,
                                  std::vector<NumType>& linker)
{
    if (linker.size() != 3 || linker.size() != 5)
    {
        THERMO_ERR("Linker size (" + std::to_string(linker.size()) + ") have the wrong size!")
    }

    NumType linker_rhs = 0;

    if (linker.size() == 5)
    {
        linker[1] -= linker[0] * diag_first.back()/under_first.back();
        linker_rhs -= linker[0] * rhs_first.back()/under_first.back();

        linker[3] -= linker[4] * diag_second[0]/over_second[0];
        linker_rhs -= linker[4] * rhs_second[0]/over_second[0];
    }

    std::vector<NumType> under_concat(under_first);
    under_concat.insert(under_concat.end(), {linker[1], 0});
    under_concat.insert(under_concat.end(), under_second);

    std::vector<NumType> diag_concat(diag_first);
    diag_concat.push_back(linker[2]);
    diag_concat.insert(diag_concat.end(), diag_second);

    std::vector<NumType> over_concat(over_first);
    over_concat.insert(over_concat.end(), {0, linker[3]});
    over_concat.insert(over_concat.end(), over_second);

    std::vector<NumType> rhs_concat(rhs_first);
    rhs_concat.push_back(linker_rhs);
    rhs_concat.insert(rhs_second.begin(), rhs_second);

    return {under_concat, diag_concat, over_concat, rhs_concat};
}

template <typename NumType>
std::vector<NumType> operator+(const std::vector<NumType>& vec1,
                               const std::vector<NumType>& vec2)
{
    std::vector<NumType> vec_sum(vec1);
    for (int i = 0; i < vec2.size(); i++)
    {
        vec_sum[i] += vec2[i];
    }

    return vec_sum;
}

template <typename NumType>
std::vector<NumType> operator-(const std::vector<NumType>& vec1,
                               const std::vector<NumType>& vec2)
{
    std::vector<NumType> vec_sum(vec1);
    for (int i = 0; i < vec2.size(); i++)
    {
        vec_sum[i] -= vec2[i];
    }

    return vec_sum;
}

template <typename NumType>
NumType L2(const std::vector<NumType>& vec)
{
    NumType vec_norm(0.0);
    for (NumType el: vec)
    {
        vec_norm += el*el;
    }

    return std::sqrt(vec_norm);
}

template
FourVecs<float> concat_matrices(const std::vector<float>& under_first,
                                const std::vector<float>& diag_first,
                                const std::vector<float>& over_first,
                                const std::vector<float>& under_second,
                                const std::vector<float>& diag_second,
                                const std::vector<float>& over_second,
                                std::vector<float>& linker);

template
FourVecs<double> concat_matrices(const std::vector<double>& under_first,
                                 const std::vector<double>& diag_first,
                                 const std::vector<double>& over_first,
                                 const std::vector<double>& under_second,
                                 const std::vector<double>& diag_second,
                                 const std::vector<double>& over_second,
                                 std::vector<double>& linker);

template
std::vector<float> operator+(const std::vector<float>& vec1,
                               const std::vector<float>& vec2);

template
std::vector<double> operator+(const std::vector<double>& vec1,
                               const std::vector<double>& vec2);

template
std::vector<float> operator-(const std::vector<float>& vec1,
                               const std::vector<float>& vec2);

template
std::vector<double> operator-(const std::vector<double>& vec1,
                               const std::vector<double>& vec2);

template
float L2(const std::vector<float>& vec);

template
double L2(const std::vector<double>& vec);
