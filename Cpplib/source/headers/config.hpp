#pragma once
#include "constants.hpp"
#include "defines.hpp"
#include <fstream>
#include <string>
#include <map>

#ifdef USE_JSON_OUTPUT
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#endif

namespace icethermo
{

    // ! nemelist for general constants in config (.json) file
    static const std::map<std::string, std::string> GenConstNames = 
    {
        {"mu", "Fusion temperature-salinity constant (deg C psu-1)"},
        {"sigma", "Stephan-Boltzman constant (W m-2 K-4)"},
        {"T0", "Zero Celsium (K)"},
        {"emissivity", "Emissivity of ice and snow (-)"},
        {"C_sh", "Default sensible heat transport coefficient (-)"},
        {"C_lh", "Default latent heat transport coefficient (-)"},
        {"L_s", "Latent heat of sublimation (J kg-1)"}
    };


    // ! pointer to class with configured constants
    extern void* config_consts;

    // ! class with configured constants
    template<typename NumType>
    struct Config
    {
        struct gen_consts
        {
            NumType mu;
            NumType sigma;
            NumType T0;
            NumType emissivity;
            NumType C_sh;
            NumType C_lh;
            NumType L_s;
        };

        gen_consts GenConsts;
    };

    // ! configure constants using .json file
    template<typename NumType>
    void Configure(const char* filepath);

    // ! parse .json configuration file and fill out constants class
    template<typename NumType>
    void ParseJson(const char* filepath, Config<NumType>* config_ptr);

    // ! check out if configuration has already been made
    bool Configured();

    // ! get pointer to class with configured constants
    template<typename NumType>
    Config<NumType>* GetConfigConsts();
}