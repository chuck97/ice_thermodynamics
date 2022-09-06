#include "matvec.hpp"
#include "tools.hpp"


namespace icethermo
{
    template <typename NumType> std::vector<NumType> concatenate(const std::vector<std::vector<NumType>>& vecs)
    {
        std::vector<NumType> res;

        for (int i = 0; i < vecs.size(); i++)
        {
            res.insert
            (
                res.end(),
                vecs[i].begin(),
                vecs[i].end()
            );
        }
        return res;
    }

    template <typename NumType>
    FourVecs<NumType> concat_matrices(const std::vector<NumType>& under_first,
                                      const std::vector<NumType>& diag_first,
                                      const std::vector<NumType>& over_first,
                                      const std::vector<NumType>& rhs_first,
                                      const std::vector<NumType>& under_second,
                                      const std::vector<NumType>& diag_second,
                                      const std::vector<NumType>& over_second,
                                      const std::vector<NumType>& rhs_second,
                                      const std::vector<NumType>& linker_)
    {
        auto linker = linker_;

        if (linker.size() != 3 && linker.size() != 5)
        {
            THERMO_ERR("Linker size (" + std::to_string(linker.size()) + ") have the wrong size!")
        }

        NumType linker_rhs = 0.0;

        if (linker.size() == 5)
        {
            linker[1] -= linker[0] * diag_first.back()/under_first.back();
            linker_rhs -= linker[0] * rhs_first.back()/under_first.back();

            linker[3] -= linker[4] * diag_second[0]/over_second[0];
            linker_rhs -= linker[4] * rhs_second[0]/over_second[0];

            linker = std::vector<NumType>(linker.begin() + 1, linker.end() - 1);
        }

        std::vector<NumType> under_concat = concatenate<NumType>({under_first, {linker[0], 0.0}, under_second});
        std::vector<NumType> diag_concat = concatenate<NumType>({diag_first, {linker[1]}, diag_second});
        std::vector<NumType> over_concat = concatenate<NumType>({over_first, {0.0, linker[2]}, over_second});

        std::vector<NumType> rhs_concat = concatenate<NumType>({rhs_first, {linker_rhs}, rhs_second});

        return {under_concat, diag_concat, over_concat, rhs_concat};
    }

    template <typename NumType>
    std::vector<NumType> operator+(const std::vector<NumType>& vec1,
                                   const std::vector<NumType>& vec2)
    {
        if (vec1.size() != vec2.size())
        {
            THERMO_ERR("can't compute the sum of two vectors with differrent size!");
        }

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
        if (vec1.size() != vec2.size())
        {
            THERMO_ERR("can't compute the diff of two vectors with differrent size!");
        }

        std::vector<NumType> vec_diff(vec1);
        for (int i = 0; i < vec2.size(); i++)
        {
            vec_diff[i] -= vec2[i];
        }

        return vec_diff;
    }

    template <typename NumType>
    NumType operator*(const std::vector<NumType>& vec1,
                      const std::vector<NumType>& vec2)
    {
        if (vec1.size() != vec2.size())
        {
            THERMO_ERR("can't compute scalar product of two vectors with differrent size!");
        }

        NumType res = 0.0;
        for (int i = 0; i < vec2.size(); i++)
        {
            res += vec1[i]*vec2[i];
        }
        return res;
    }

    template <typename NumType>
    NumType L2_norm(const std::vector<NumType>& vec)
    {                              
        NumType res;

        for (int i = 0; i < vec.size(); i++)
        {
            res += vec[i]*vec[i];
        }
        return std::sqrt(res);
    }

    template <typename NumType>
    std::vector<NumType> operator*(const std::vector<NumType>& vec, NumType scal)
    {
        auto res_vec = vec;
        for (int i = 0 ; i < vec.size(); i++)
        {
            res_vec[i] *= scal;
        }
        return res_vec;
    }

    template <typename NumType> 
    std::vector<NumType> operator*(NumType scal, const std::vector<NumType>& vec)
    {
        auto res_vec = vec;
        for (int i = 0 ; i < vec.size(); i++)
        {
            res_vec[i] *= scal;
        }
        return res_vec;
    }

    template <typename NumType> 
    std::ostream& operator<<(std::ostream& os, const std::vector<NumType>& vec)
    {
        os << "[";
        for (auto it = vec.begin(); it != vec.end(); it++)
        {
            os << *it;
            if (std::next(it) != vec.end())
                os << ", ";
        }
        os << "]";
        return os;
    }


    // explicit instantaion
    template std::vector<int> concatenate(const std::vector<std::vector<int>>& vecs);
    template std::vector<float> concatenate(const std::vector<std::vector<float>>& vecs);
    template std::vector<double> concatenate(const std::vector<std::vector<double>>& vecs);


    template FourVecs<float> concat_matrices(const std::vector<float>& under_first,
                                             const std::vector<float>& diag_first,
                                             const std::vector<float>& over_first,
                                             const std::vector<float>& rhs_first,
                                             const std::vector<float>& under_second,
                                             const std::vector<float>& diag_second,
                                             const std::vector<float>& over_second,
                                             const std::vector<float>& rhs_second,
                                             const std::vector<float>& linker_);

    template FourVecs<double> concat_matrices(const std::vector<double>& under_first,
                                              const std::vector<double>& diag_first,
                                              const std::vector<double>& over_first,
                                              const std::vector<double>& rhs_first,
                                              const std::vector<double>& under_second,
                                              const std::vector<double>& diag_second,
                                              const std::vector<double>& over_second,
                                              const std::vector<double>& rhs_second,
                                              const std::vector<double>& linker_);

    template std::vector<float> operator+(const std::vector<float>& vec1,
                                          const std::vector<float>& vec2);

    template std::vector<double> operator+(const std::vector<double>& vec1,
                                           const std::vector<double>& vec2);

    template std::vector<float> operator-(const std::vector<float>& vec1,
                                          const std::vector<float>& vec2);

    template std::vector<double> operator-(const std::vector<double>& vec1,
                                           const std::vector<double>& vec2);

    template float operator*(const std::vector<float>& vec1,
                             const std::vector<float>& vec2);

    template double operator*(const std::vector<double>& vec1,
                              const std::vector<double>& vec2);

    template float L2_norm(const std::vector<float>& vec);
    template double L2_norm(const std::vector<double>& vec);

    template std::vector<float> operator*(const std::vector<float>& vec, float scal);
    template std::vector<double> operator*(const std::vector<double>& vec, double scal);

    template std::vector<float> operator*(float scal, const std::vector<float>& vec);
    template std::vector<double> operator*(double scal, const std::vector<double>& vec);

    template std::ostream& operator<<(std::ostream& os, const std::vector<float>& vec);
    template std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec);
}