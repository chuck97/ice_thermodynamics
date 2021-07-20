#include "thomas_solver.h"

using namespace std;

// lhs is vector of 3 diagonals: ld, c, ru
vector<double> thomas_solver(const vector<vector<double>>& lhs,
                             const vector<double>& rhs)
{
    // check if lhs 3-diagonal
    if (lhs.size() != 3)
    {
        ERR("lhs must be tridiagonal matrix: {ld, c, ru}");
    }

    vector<double> ld = lhs[0];
    vector<double> c =  lhs[1];
    vector<double> ru = lhs[2];


    //  check that c.size = ld.size+1 = ru.size+1
    if (c.size() != (ld.size() + 1) or
        c.size() != (ru.size() + 1))
    {
        ERR("c.size must be equal to (ld.size+1) and (ru.size+1)");
    } 

    // check that c.size = rhs.size
    if (c.size() != rhs.size())
    {
        ERR("c.size must be equal to rhs.size");
    }

    vector<double> sol(c.size());
    vector<double> c_new(c.size()-1);
    vector<double> d_new(c.size());

    // forward iteration
    c_new[0] = ru[0]/c[0];
    for (int i = 1; i<(c.size()-1); ++i)
    {
        c_new[i] = ru[i]/(c[i] - ld[i-1]*c_new[i-1]);
    }

    // backward iteration
    d_new[0] = rhs[0]/c[0];
    for (int i = 1; i<(c.size()); ++i)
    {
        d_new[i] = (rhs[i] - ld[i-1]*d_new[i-1])/(c[i] - ld[i-1]*c_new[i-1]);
    }

    // final solution
    sol[c.size()-1] = d_new[c.size()-1];
    for (int i = c.size()-2; i >= 0; --i)
    {
        sol[i] = d_new[i] - c_new[i]*sol[i+1];
    }

    return sol;
}