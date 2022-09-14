#include "constants.hpp"

using namespace icethermo;

double Params::Density(Dparam param, double T, double S)
{
    if (param == Dparam::SeaIce)
    {
        return IceConsts::rho_i;
    }
    else if (param == Dparam::FreshIce)
    {
        return IceConsts::rho_i;
    }
    else if (param == Dparam::FreshSnow)
    {
        return SnowConsts::rho_s;
    }
    else
    {
        THERMO_ERR("Available Density parametrizations: SeaIce, FreshIce, FreshSnow!");
    }
}

double Params::EffCapacity(Cparam param, double T, double T_old, double S)
{
    if (param == Cparam::SeaIce)
    {
        double Tf_i = GenConsts::TempFusion(S);
        return IceConsts::c0_i - IceConsts::L0_i*Tf_i/(T*T_old);
    }
    else if (param == Cparam::FreshIce)
    {
        return IceConsts::c0_i;
    }
    else if (param == Cparam::FreshSnow)
    {
        return SnowConsts::c0_s;
    }
    else
    {
        THERMO_ERR("Available Effective Heat Capacity parametrizations: SeaIce, FreshIce, FreshSnow!");
    }
}

double Params::Enthalpy(Eparam param, double T, double S)
{
    if (param == Eparam::SeaIce)
    {
        double Tf_i = GenConsts::TempFusion(S);
        double c_w = WaterConsts::c_w;
        return IceConsts::c0_i*(T - Tf_i) - IceConsts::L0_i*(1.0 - Tf_i/T) + c_w*Tf_i;
    }
    else if (param == Eparam::FreshIce)
    {
        return IceConsts::c0_i*T - IceConsts::L0_i;
    }
    else if (param == Eparam::FreshSnow)
    {
        return SnowConsts::c0_s*T - SnowConsts::L_f0;
    }
    else
    {
        THERMO_ERR("Available Enthalpy parametrizations: SeaIce, FreshIce, FreshSnow!");
    }
    
}

double Params::FusionHeat(Lparam param, double T, double S)
{
    if (param == Lparam::SeaIce)
    {
        double Tf_i = GenConsts::TempFusion(S);
        double c_w = WaterConsts::c_w;
        return IceConsts::c0_i*(Tf_i - T) + IceConsts::L0_i*(1.0 + Tf_i/T);
    }
    else if (param == Lparam::FreshIce)
    {
        return -IceConsts::c0_i*T + IceConsts::L0_i;
    }
    else if (param == Lparam::FreshSnow)
    {
        return -SnowConsts::c0_s*T + SnowConsts::L_f0;
    }
    else
    {
        THERMO_ERR("Available Fusion Heat parametrizations: SeaIce, FreshIce, FreshSnow!");
    }
}

double Params::Conductivity(Kparam param, double T, double S, double rho)
{
    if (param == Kparam::Untersteiner)
    {
        if (std::abs(T) < REAL_MIN_VAL)
        {
            return 0.0;
        }
        else
        {
            return 2.03 + 0.1172 * (S/T);
        }
        
    }
    else if (param == Kparam::BubblyBrine)
    {
        if (std::abs(T) < REAL_MIN_VAL)
        {
            return 0.0;
        }
        else
        {
            return rho/IceConsts::rho_i * (2.11 - 0.011*T + (0.09*S/T));
        }
    }
    else if (param == Kparam::FreshIce)
    {
        return IceConsts::k0_i;
    }
    else if (param == Kparam::FreshSnow)
    {
        return SnowConsts::k0_s;
    }
    else
    {
        THERMO_ERR("Available Thermal Conductivity parametrizations: Untersteiner, BubblyBrine, FreshIce, FreshSnow!");
    }
}