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

    static const std::map<std::string, std::string> AirConstNames = 
    {
        {"rho_a", "Air density (kg m-3)"},
        {"cp_a", "Heat capacity of air (J kg-1 K-1)"}
    };

    static const std::map<std::string, std::string> WaterConstNames = 
    {
        {"rho_w", "Water density (kg m-3)"},
        {"cp_w", "Heat capacity of fresh water (J kg-1 K-1)"}
    };

    static const std::map<std::string, std::string> IceConstNames = 
    {
        {"c0_i", "Heat capacity of freshwater ice (J kg-1 K-1)"},
        {"rho_i", "Density of freshwater ice (kg m-3)"},
        {"L0_i", "Latent heat of fusion of freshwater ice (J kg-1)"},
        {"kappa_i", "Bulk radiation extinction coefficient of ice (m-1)"},
        {"albedo_i", "Default ice albedo (-)"},
        {"i0_i", "Surface radiation transmission coefficient of ice (-)"},
        {"k0_i", "Thermal conductivity of freshwater ice (W m-1 K-1)"},
        {"albedo_dry_i", "Albedo of dry ice (-)"},
        {"albedo_wet_i", "Albedo of wet ice (-)"}
    };

    static const std::map<std::string, std::string> SnowConstNames = 
    {
        {"c0_s", "Heat capacity of snow (J kg-1 K-1)"},
        {"rho_s", "Density of snow (kg m-3)"},
        {"L_f0", "Latent heat of fusion of snow (J kg-1)"},
        {"kappa_s", "Bulk radiation extinction coefficient of snow (m-1)"},
        {"albedo_s", "Default snow albedo (-)"},
        {"i0_s", "Surface radiation transmission coefficient of snow (-)"},
        {"k0_s", "Thermal conductivity of snow (W m-1 K-1)"},
        {"albedo_dry_s", "Albedo of dry snow (-)"},
        {"albedo_wet_s", "Albedo of wet snow (-)"}
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

        struct air_consts
        {
            NumType rho_a;
            NumType cp_a;
        };

        struct water_consts
        {
            NumType rho_w;
            NumType cp_w;
        };

        struct ice_consts
        {
            NumType c0_i;
            NumType rho_i;
            NumType L0_i;
            NumType kappa_i;
            NumType albedo_i;
            NumType i0_i;
            NumType k0_i;
            NumType albedo_dry_i;
            NumType albedo_wet_i;
        };

        struct snow_consts
        {
            NumType c0_s;
            NumType rho_s;
            NumType L_f0;
            NumType kappa_s;
            NumType albedo_s;
            NumType i0_s;
            NumType k0_s;
            NumType albedo_dry_s;
            NumType albedo_wet_s;
        };

        gen_consts GenConsts;
        air_consts AirConsts;
        water_consts WaterConsts;
        ice_consts IceConsts;
        snow_consts SnowConsts;
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