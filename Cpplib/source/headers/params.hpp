#pragma once
#include "constants.hpp"
#include "config.hpp"

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
    struct Params
    {
        static NumType TempFusion(double S); 
        static NumType Density(Dparam param, NumType T, NumType S);
        static NumType EffCapacity(Cparam param, NumType T, NumType T_old, NumType S);
        static NumType Enthalpy(Eparam param, NumType T, NumType S);
        static NumType FusionHeat(Lparam param, NumType T, NumType S);
        static NumType Conductivity(Kparam param, NumType T, NumType S, NumType rho);
        static NumType Albedo(Aparam param, NumType T, NumType h, NumType Tfsurf);
    };
}