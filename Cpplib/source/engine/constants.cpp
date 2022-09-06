#include "constants.hpp"

using namespace icethermo;

double IceInfo::IceDensity(double T, double S)
{
    return rho_i;
}

double IceInfo::EffCapacity(double T, double T_old, double S)
{
    double Tf_i = GenInfo::TempFusion(S);
    return c0_i - L0_i*Tf_i/(T*T_old);
}

double IceInfo::Enthalpy(double T, double S)
{
    double Tf_i = GenInfo::TempFusion(S);
    double c_w = WaterInfo::c_w;
    return c0_i*(T - Tf_i) - L0_i*(1.0 - Tf_i/T) + c_w*Tf_i;
}

double IceInfo::FusionHeat(double T, double S)
{
    double Tf_i = GenInfo::TempFusion(S);
    double c_w = WaterInfo::c_w;
    return c0_i*(Tf_i - T) + L0_i*(1.0 + Tf_i/T);
}

double IceInfo::Conductivity(Kparam param, double T, double S)
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
            double rho = IceDensity(T, S);
            return rho/rho_i * (2.11 - 0.011*T + (0.09*S/T));
        }
    }
    else
    {
        THERMO_ERR("Only Untersteiner and bubbly-brine thermal conductivity parametrizations available!");
    }
}

double SnowInfo::EffCapacity(double T, double T_old)
{
    return c0_s;
}

double SnowInfo::Enthalpy(double T)
{
    return c0_s*T - L_f0;
}

double SnowInfo::FusionHeat(double T)
{
    return -c0_s*T + L_f0;
}

double SnowInfo::Conductivity(double T)
{
    return k0_s;
}