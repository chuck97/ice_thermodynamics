#include "constants.hpp"

using namespace icethermo;

double IceInfo::Density(Dparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Dparam::SeaIce)
    {
        return rho_i;
    }
    else if (param == Dparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Ice Density parametrizations: SeaIce, FreshSnow, Custom!");
    }
    return rho_i;
}

double IceInfo::EffCapacity(Cparam param, double T, double T_old, double S, std::vector<double> params, CustomEffFuncPtr fptr)
{
    if (param == Cparam::SeaIce)
    {
        double Tf_i = GenInfo::TempFusion(S);
        return c0_i - L0_i*Tf_i/(T*T_old);
    }
    else if (param == Cparam::Custom)
    {
        return fptr(T, T_old, S, params);
    }
    else
    {
        THERMO_ERR("Available Ice Effective Heat Capacity parametrizations: SeaIce, FreshSnow, Custom!");
    }
}

double IceInfo::Enthalpy(Eparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Eparam::SeaIce)
    {
        double Tf_i = GenInfo::TempFusion(S);
        double c_w = WaterInfo::c_w;
        return c0_i*(T - Tf_i) - L0_i*(1.0 - Tf_i/T) + c_w*Tf_i;
    }
    else if (param == Eparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Ice Enthalpy parametrizations: SeaIce, FreshSnow, Custom!");
    }
    
}

double IceInfo::FusionHeat(Lparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Lparam::SeaIce)
    {
        double Tf_i = GenInfo::TempFusion(S);
        double c_w = WaterInfo::c_w;
        return c0_i*(Tf_i - T) + L0_i*(1.0 + Tf_i/T);
    }
    else if (param == Lparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Ice Fusion Heat parametrizations: SeaIce, FreshSnow, Custom!");
    }
}

double IceInfo::Conductivity(Kparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
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
            double rho = Density(Dparam::SeaIce, T, S);
            return rho/rho_i * (2.11 - 0.011*T + (0.09*S/T));
        }
    }
    else if (param == Kparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Ice Thermal Conductivity parametrizations: Untersteiner, BubblyBrine, Custom!");
    }
}

double SnowInfo::Density(Dparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Dparam::FreshSnow)
    {
        return rho_s;
    }
    else if (param == Dparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Snow Density parametrizations: FreshSnow, Custom!");
    }
}

double SnowInfo::EffCapacity(Cparam param, double T, double T_old, double S, std::vector<double> params, CustomEffFuncPtr fptr)
{
    if (param == Cparam::FreshSnow)
    {
        return c0_s;
    }
    else if (param == Cparam::Custom)
    {
        return fptr(T, T_old, S, params);
    }
    else
    {
        THERMO_ERR("Available Snow Effective Capacity parametrizations: FreshSnow, Custom!");
    }
}

double SnowInfo::Enthalpy(Eparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Eparam::FreshSnow)
    {
        return  c0_s*T - L_f0;
    }
    else if (param == Eparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Snow Enthalpy parametrizations: FreshSnow, Custom!");
    }
}

double SnowInfo::FusionHeat(Lparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Lparam::FreshSnow)
    {
        return -c0_s*T + L_f0;
    }
    else if (param == Lparam::Custom) 
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Snow Fusion Heat parametrizations: FreshSnow, Custom!");
    }
}

double SnowInfo::Conductivity(Kparam param, double T, double S, std::vector<double> params, CustomFuncPtr fptr)
{
    if (param == Kparam::FreshSnow)
    {
        return k0_s;
    }
    else if (param == Kparam::Custom)
    {
        return fptr(T, S, params);
    }
    else
    {
        THERMO_ERR("Available Snow Thermal Conductivity parametrizations: FreshSnow, Custom!");
    }
}