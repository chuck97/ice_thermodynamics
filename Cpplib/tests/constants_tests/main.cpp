#include <iostream>
#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    std::cout << IceInfo::Conductivity(Kparam::Untersteiner, -5.0, 20.0) << std::endl;
    std::cout << IceInfo::rho_i << std::endl;
    return 0;
}