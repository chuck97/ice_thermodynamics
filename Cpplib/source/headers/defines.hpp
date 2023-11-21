#pragma once
#include <iostream>
#include <limits>

namespace icethermo
{
    #define THERMO_ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}
    #define REAL_MIN_VAL(NumType) std::numeric_limits<NumType>::min()*1e10  
    #define ALLOWABLE_RELATIVE_1D_ERROR 1e-1
    #define SNOW_THICKNESS_THRESHOLD 1e-3
    #define NAN_TEMP_VALUE -1000000.0
    #define NAN_THICK_VALUE -1000000.0
}
