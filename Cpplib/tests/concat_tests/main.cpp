#include <iostream>
#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    auto result_double = concat_matrices(std::vector<double>{1.0, 2.0, 3.0},
                                         std::vector<double>{1.0, 2.0, 3.0, 4.0},
                                         std::vector<double>{1.0, 2.0, 3.0},
                                         std::vector<double>{5.0, 5.0, 5.0, 5.0},
                                         std::vector<double>{-1.0, -2.0, -3.0},
                                         std::vector<double>{-1.0, -2.0, -3.0, -4.0},
                                         std::vector<double>{-1.0, -2.0, -3.0},
                                         std::vector<double>{6.0, 6.0, 6.0, 6.0},
                                         std::vector<double>{7.0, 8.0, 9.0});

    std::cout << std::get<0>(result_double) << "\n" 
              << std::get<1>(result_double) << "\n" 
              << std::get<2>(result_double) << "\n"
              << std::get<3>(result_double) << "\n";
}