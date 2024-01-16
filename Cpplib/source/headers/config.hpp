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

    // ! nemelist for constants in config (.json) file
    static const std::map<std::string, std::string> GenConstNames = 
    {
        {"mu", "Fusion temperature-salinity constant (deg C psu-1)"},
        {"sigma", "Stephan-Boltzman constant (W m-2 K-4)"},
        {"T0", "Zero Celsium (K)"},
        {"emissivity", "Emissivity of ice and snow (-)"},
        {"C_sh", "Default sensible heat transport coefficient (-)"},
        {"C_lh", "Default latent heat transport coefficient (-)"},
        {"L_sublim", "Latent heat of sublimation (J kg-1)"}
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
        {"L0_s", "Latent heat of fusion of snow (J kg-1)"},
        {"kappa_s", "Bulk radiation extinction coefficient of snow (m-1)"},
        {"albedo_s", "Default snow albedo (-)"},
        {"i0_s", "Surface radiation transmission coefficient of snow (-)"},
        {"k0_s", "Thermal conductivity of snow (W m-1 K-1)"},
        {"albedo_dry_s", "Albedo of dry snow (-)"},
        {"albedo_wet_s", "Albedo of wet snow (-)"}
    };

    // ! types of root-finding solvers
    enum class Solver1dType
    {
        bisection,
        secant, 
        ridders
    };

    // ! namelist for 1d root-finding solvers
    static const std::map<std::string, Solver1dType> Solver1dNameToType = 
    {
        {"Bisection", Solver1dType::bisection},
        {"Secant", Solver1dType::secant},
        {"Ridders", Solver1dType::ridders}
    };


    // ! pointer to class with configured constants
    extern void* config_consts;

    // ! class with configured constants
    template<typename NumType>
    struct Config
    {
        struct gen_consts
        {
            NumType mu = icethermo::GenConsts<NumType>::mu;
            NumType sigma = icethermo::GenConsts<NumType>::sigma;
            NumType T0 = icethermo::GenConsts<NumType>::T0;
            NumType emissivity = icethermo::GenConsts<NumType>::emissivity;
            NumType C_sh = icethermo::GenConsts<NumType>::C_sh;
            NumType C_lh = icethermo::GenConsts<NumType>::C_lh;
            NumType L_sublim = icethermo::GenConsts<NumType>::L_sublim;
        };

        struct air_consts
        {
            NumType rho_a = icethermo::AirConsts<NumType>::rho_a;
            NumType cp_a = icethermo::AirConsts<NumType>::cp_a;
        };

        struct water_consts
        {
            NumType rho_w = icethermo::WaterConsts<NumType>::rho_w;
            NumType cp_w = icethermo::WaterConsts<NumType>::cp_w;
        };

        struct ice_consts
        {
            NumType c0_i = icethermo::IceConsts<NumType>::c0_i;
            NumType rho_i = icethermo::IceConsts<NumType>::rho_i;
            NumType L0_i = icethermo::IceConsts<NumType>::L0_i;
            NumType kappa_i = icethermo::IceConsts<NumType>::kappa_i;
            NumType albedo_i = icethermo::IceConsts<NumType>::albedo_i;
            NumType i0_i = icethermo::IceConsts<NumType>::i0_i;
            NumType k0_i = icethermo::IceConsts<NumType>::k0_i;
            NumType albedo_dry_i = icethermo::IceConsts<NumType>::albedo_dry_i;
            NumType albedo_wet_i = icethermo::IceConsts<NumType>::albedo_wet_i;
        };

        struct snow_consts
        {
            NumType c0_s = icethermo::SnowConsts<NumType>::c0_s;
            NumType rho_s = icethermo::SnowConsts<NumType>::rho_s;
            NumType L0_s = icethermo::SnowConsts<NumType>::L0_s;
            NumType kappa_s = icethermo::SnowConsts<NumType>::kappa_s;
            NumType albedo_s = icethermo::SnowConsts<NumType>::albedo_s;
            NumType i0_s = icethermo::SnowConsts<NumType>::i0_s;
            NumType k0_s = icethermo::SnowConsts<NumType>::k0_s;
            NumType albedo_dry_s = icethermo::SnowConsts<NumType>::albedo_dry_s;
            NumType albedo_wet_s = icethermo::SnowConsts<NumType>::albedo_wet_s;
        };

        struct solver_1d
        {
            Solver1dType type = Solver1dType::bisection;
            int max_nits = MAX_1D_SOLVER_ITS;
            NumType residual = NONLIN_SOLVER_ACCUR;
        };

        struct solver_relaxation
        {
            bool force_surf_conv = true;
            int max_nits = MAX_RELAXATION_ITS;
            NumType inc_error = RELAXATION_SOLVER_ACCUR;
        };

        struct nans
        {
            NumType temp_nan = NAN_TEMP_VALUE;
            NumType thick_nan = NAN_THICK_VALUE;
        };

        gen_consts GenConsts;
        air_consts AirConsts;
        water_consts WaterConsts;
        ice_consts IceConsts;
        snow_consts SnowConsts;
        solver_1d Solver1d;
        solver_relaxation SolverRelaxation;
        nans NaNs;
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