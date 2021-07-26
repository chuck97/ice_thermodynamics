#pragma once
#include "vec_norms.h"

using namespace std;

double vec_C_norm(const vector<double>& vec)
{
    double current_max = abs(vec[0]);

    for (unsigned int i = 1; i < vec.size(); ++i)
    {
        current_max = max(current_max, abs(vec[i]));
    }
    return current_max;
}

double vec_L2_norm(const vector<double>& vec)
{
    double sum_sqr = 0.0;

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        sum_sqr += vec[i]*vec[i];
    }
    return sqrt(sum_sqr);
}