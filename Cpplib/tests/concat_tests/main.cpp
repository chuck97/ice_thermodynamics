#include <iostream>
#include "icethermo.hpp"

using namespace icethermo;

int main()
{
    auto ordinary_double = concat_matrices(std::vector<double>{1.0, 2.0, 3.0},
                                           std::vector<double>{1.0, 2.0, 3.0, 4.0},
                                           std::vector<double>{1.0, 2.0, 3.0},
                                           std::vector<double>{5.0, 5.0, 5.0, 5.0},
                                           std::vector<double>{-1.0, -2.0, -3.0},
                                           std::vector<double>{-1.0, -2.0, -3.0, -4.0},
                                           std::vector<double>{-1.0, -2.0, -3.0},
                                           std::vector<double>{6.0, 6.0, 6.0, 6.0},
                                           std::vector<double>{7.0, 8.0, 9.0});

    std::cout << "Ordinary double vecs: \n"
              << std::get<0>(ordinary_double) << "\n"
              << std::get<1>(ordinary_double) << "\n"
              << std::get<2>(ordinary_double) << "\n"
              << std::get<3>(ordinary_double) << "\n";

    auto ordinary_float = concat_matrices(std::vector<float>{1.0, 2.0, 3.0},
                                          std::vector<float>{1.0, 2.0, 3.0, 4.0},
                                          std::vector<float>{1.0, 2.0, 3.0},
                                          std::vector<float>{5.0, 5.0, 5.0, 5.0},
                                          std::vector<float>{-1.0, -2.0, -3.0},
                                          std::vector<float>{-1.0, -2.0, -3.0, -4.0},
                                          std::vector<float>{-1.0, -2.0, -3.0},
                                          std::vector<float>{6.0, 6.0, 6.0, 6.0},
                                          std::vector<float>{7.0, 8.0, 9.0});

    std::cout << "Ordinary float vecs: \n"
              << std::get<0>(ordinary_float) << "\n"
              << std::get<1>(ordinary_float) << "\n"
              << std::get<2>(ordinary_float) << "\n"
              << std::get<3>(ordinary_float) << "\n";

    auto secord_double = concat_matrices(std::vector<double>{1.0, 2.0},
                                         std::vector<double>{1.0, 2.0, 3.0},
                                         std::vector<double>{2.0, 3.0},
                                         std::vector<double>{6.0, 6.0, 6.0},
                                         std::vector<double>{3.0, 4.0},
                                         std::vector<double>{3.0, 4.0, 5.0},
                                         std::vector<double>{4.0, 5.0},
                                         std::vector<double>{6.0, 6.0, 6.0},
                                         std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0});

    std::cout << "Second-order double vecs: \n"
              << std::get<0>(secord_double) << "\n"
              << std::get<1>(secord_double) << "\n"
              << std::get<2>(secord_double) << "\n"
              << std::get<3>(secord_double) << "\n";

    auto secord_float = concat_matrices(std::vector<float>{1.0, 2.0},
                                        std::vector<float>{1.0, 2.0, 3.0},
                                        std::vector<float>{2.0, 3.0},
                                        std::vector<float>{6.0, 6.0, 6.0},
                                        std::vector<float>{3.0, 4.0},
                                        std::vector<float>{3.0, 4.0, 5.0},
                                        std::vector<float>{4.0, 5.0},
                                        std::vector<float>{6.0, 6.0, 6.0},
                                        std::vector<float>{1.0, 2.0, 3.0, 4.0, 5.0});

    std::cout << "Second-order float vecs: \n"
              << std::get<0>(secord_float) << "\n"
              << std::get<1>(secord_float) << "\n"
              << std::get<2>(secord_float) << "\n"
              << std::get<3>(secord_float) << "\n";
}
