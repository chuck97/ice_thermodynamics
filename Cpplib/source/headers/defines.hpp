#pragma once
#include <iostream>

#define THERMO_ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}
#define REAL_MIN_VAL 1e-30