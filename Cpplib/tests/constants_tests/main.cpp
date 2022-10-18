#include <iostream>
#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    std::cout << Params<double>::Conductivity(Kparam::Untersteiner, -5.0, 20.0, IceConsts<double>::rho_i) << std::endl;
    std::cout << GenConsts<float>::mu << std::endl;
    return 0;
}