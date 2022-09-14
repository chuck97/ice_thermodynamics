#pragma once
#include "defines.hpp"
#include <vector>

namespace icethermo
{
    // thermal conductivity parametrization enum
    enum class Kparam
    {
        Untersteiner,
        BubblyBrine,
        FreshIce,
        FreshSnow
    };

    // density parametrization enum
    enum class Dparam
    {
        SeaIce,
        FreshIce,
        FreshSnow
    };

    // Effective heat capacity parametrization enum
    enum class Cparam
    {
        SeaIce,
        FreshIce,
        FreshSnow
    };

    // Enthalpy parametrization enum
    enum class Eparam
    {
        SeaIce,
        FreshIce,
        FreshSnow
    };

    // Specific heat of fusion parametrization enum
    enum class Lparam
    {
        SeaIce,
        FreshIce,
        FreshSnow
    };

    struct GenConsts
    {
        inline static double TempFusion(double S) {return -mu*S;}; 

        static constexpr double mu = 0.054;
        static constexpr double sigma = 5.67e-8;
        static constexpr double T0 = 273.16;
        static constexpr double emissivity = 0.99;
        static constexpr double C_sh = 1e-3;
        static constexpr double C_lh = 1e-3;
    };

    struct AirConsts
    {
        static constexpr double rho_a = 1.28;
        static constexpr double cp_a = 1.01e3;
        static constexpr double P_surf = 1013.25;
    };

    struct WaterConsts
    {
        static constexpr double c_w =  3.99e3;
        static constexpr double c_pw =  4.17e3;
        static constexpr double rho_w =  1023.0;
        static constexpr double c1_w =  17.27;
        static constexpr double c2_w =  35.86;
    };

    struct IceConsts
    {
        static constexpr double c0_i = 2.06e3;
        static constexpr double rho_i = 917.0;
        static constexpr double L0_i = 3.34e5;
        static constexpr double kappa_i = 1.5;
        static constexpr double albedo_i = 0.6;
        static constexpr double i0_i = 0.15;
        static constexpr double c1_i = 21.87;
        static constexpr double c2_i = 7.66;
        static constexpr double k0_i = 2.03;
    };

    struct SnowConsts
    {
        static constexpr double c0_s = 2.06e3;
        static constexpr double rho_s = 330.0;
        static constexpr double k0_s = 0.31;
        static constexpr double L_f0 = 3.33e5;
        static constexpr double L_s0 = 2.83e6;
        static constexpr double kappa_s = 10;
        static constexpr double albedo_s = 0.8;
        static constexpr double i0_s = 0.08;
    };

    struct Params
    {
        static double Density(Dparam param, double T, double S);
        static double EffCapacity(Cparam param, double T, double T_old, double S);
        static double Enthalpy(Eparam param, double T, double S);
        static double FusionHeat(Lparam param, double T, double S);
        static double Conductivity(Kparam param, double T, double S, double rho);
    };
}