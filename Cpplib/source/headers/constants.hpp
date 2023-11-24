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

    // albedo parameterizations
    enum class Aparam
    {
        ConstantIce,
        ConstantSnow,
        MeltingFreezingIce,
        MeltingFreezingSnow
    };

    template<typename NumType>
    struct GenConsts
    {
        inline static NumType TempFusion(double S) {return -mu*S;}; 

        static constexpr NumType mu = 0.0575;
        static constexpr NumType sigma = 5.67e-8;
        static constexpr NumType T0 = 273.15;
        static constexpr NumType emissivity = 0.99;
        static constexpr NumType C_sh = 1e-3;
        static constexpr NumType C_lh = 1e-3;
        static constexpr NumType L_s = 2.834e6;
        static constexpr NumType si_trans_time_scale = 1e7;
    };

    template<typename NumType>
    struct AirConsts
    {
        static constexpr NumType rho_a = 1.28;
        static constexpr NumType cp_a = 1.004e3;
        static constexpr NumType P_surf = 1013.25;
    };

    template<typename NumType>
    struct WaterConsts
    {
        static constexpr NumType c_w =  3.99e3;
        static constexpr NumType c_pw =  4.17e3;
        static constexpr NumType rho_w =  1026.0;
        static constexpr NumType c1_w =  21.87;
        static constexpr NumType c2_w =  7.66;
    };

    template<typename NumType>
    struct IceConsts
    {
        static constexpr NumType c0_i = 2.06e3;
        static constexpr NumType rho_i = 917.0;
        static constexpr NumType L0_i = 3.34e5;
        static constexpr NumType kappa_i = 1.5;
        static constexpr NumType albedo_i = 0.6;
        static constexpr NumType i0_i = 0.17;
        static constexpr NumType c1_i = 21.87;
        static constexpr NumType c2_i = 7.66;
        static constexpr NumType k0_i = 2.04;
        static constexpr NumType albedo_dry_i = 0.8;
        static constexpr NumType albedo_wet_i = 0.65;
    };

    template<typename NumType>
    struct SnowConsts
    {
        static constexpr NumType c0_s = 2.06e3;
        static constexpr NumType rho_s = 320.0;
        static constexpr NumType k0_s = 0.31;
        static constexpr NumType L_f0 = 3.33e5;
        static constexpr NumType kappa_s = 10;
        static constexpr NumType albedo_s = 0.8;
        static constexpr NumType i0_s = 0.0;
        static constexpr NumType albedo_dry_s = 0.8;
        static constexpr NumType albedo_wet_s = 0.65;
    };

    template<typename NumType>
    struct Params
    {
        static NumType Density(Dparam param, NumType T, NumType S);
        static NumType EffCapacity(Cparam param, NumType T, NumType T_old, NumType S);
        static NumType Enthalpy(Eparam param, NumType T, NumType S);
        static NumType FusionHeat(Lparam param, NumType T, NumType S);
        static NumType Conductivity(Kparam param, NumType T, NumType S, NumType rho);
        static NumType Albedo(Aparam param, NumType T, NumType h, NumType Tfsurf);
    };
}