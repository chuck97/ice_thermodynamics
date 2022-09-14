#include <iostream>
#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    std::cout << Params::Conductivity(Kparam::Untersteiner, -5.0, 20.0, IceConsts::rho_i) << std::endl;
    std::cout << GenConsts::mu << std::endl;
    return 0;
}