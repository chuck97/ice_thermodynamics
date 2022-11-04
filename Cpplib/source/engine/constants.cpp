#include "constants.hpp"

namespace icethermo
{
    template<typename NumType>
    NumType Params<NumType>::Density(Dparam param, NumType T, NumType S)
    {
        if (param == Dparam::SeaIce)
        {
            return IceConsts<NumType>::rho_i;
        }
        else if (param == Dparam::FreshIce)
        {
            return IceConsts<NumType>::rho_i;
        }
        else if (param == Dparam::FreshSnow)
        {
            return SnowConsts<NumType>::rho_s;
        }
        else
        {
            THERMO_ERR("Available Density parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::EffCapacity(Cparam param, NumType T, NumType T_old, NumType S)
    {
        if (param == Cparam::SeaIce)
        {
            NumType Tf_i = GenConsts<NumType>::TempFusion(S);
            if ((std::abs(T) < REAL_MIN_VAL(NumType)) or (std::abs(T_old) < REAL_MIN_VAL(NumType)))
            {
                return IceConsts<NumType>::c0_i;
            }
            else
            {
                return IceConsts<NumType>::c0_i - IceConsts<NumType>::L0_i*Tf_i/(T*T_old);
            }
        }
        else if (param == Cparam::FreshIce)
        {
            return IceConsts<NumType>::c0_i;
        }
        else if (param == Cparam::FreshSnow)
        {
            return SnowConsts<NumType>::c0_s;
        }
        else
        {
            THERMO_ERR("Available Effective Heat Capacity parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::Enthalpy(Eparam param, NumType T, NumType S)
    {
        if (param == Eparam::SeaIce)
        {
            NumType Tf_i = GenConsts<NumType>::TempFusion(S);
            NumType c_w = WaterConsts<NumType>::c_w;
            if (std::abs(T) < REAL_MIN_VAL(NumType))
            {
                return IceConsts<NumType>::c0_i*(T - Tf_i) + c_w*Tf_i;;
            }
            else
            {
                return IceConsts<NumType>::c0_i*(T - Tf_i) - IceConsts<NumType>::L0_i*(1.0 - Tf_i/T) + c_w*Tf_i;
            }
        }
        else if (param == Eparam::FreshIce)
        {
            return IceConsts<NumType>::c0_i*T - IceConsts<NumType>::L0_i;
        }
        else if (param == Eparam::FreshSnow)
        {
            return SnowConsts<NumType>::c0_s*T - SnowConsts<NumType>::L_f0;
        }
        else
        {
            THERMO_ERR("Available Enthalpy parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::FusionHeat(Lparam param, NumType T, NumType S)
    {
        if (param == Lparam::SeaIce)
        {
            NumType Tf_i = GenConsts<NumType>::TempFusion(S);
            NumType c_w = WaterConsts<NumType>::c_w;
            if (std::abs(T) < REAL_MIN_VAL(NumType))
            {
                return IceConsts<NumType>::c0_i*T - IceConsts<NumType>::L0_i;
            }
            else if (std::abs(T - Tf_i) < REAL_MIN_VAL(NumType))
            {
                return -IceConsts<NumType>::L0_i;
            }
            else
            {
                return IceConsts<NumType>::c0_i*(T - Tf_i) - IceConsts<NumType>::L0_i*(1.0 - Tf_i/T);
            }
        }
        else if (param == Lparam::FreshIce)
        {
            return IceConsts<NumType>::c0_i*T - IceConsts<NumType>::L0_i;
        }
        else if (param == Lparam::FreshSnow)
        {
            return SnowConsts<NumType>::c0_s*T - SnowConsts<NumType>::L_f0;
        }
        else
        {
            THERMO_ERR("Available Fusion Heat parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::Conductivity(Kparam param, NumType T, NumType S, NumType rho)
    {
        if (param == Kparam::Untersteiner)
        {
            if (std::abs(T) < REAL_MIN_VAL(NumType))
            {
                return 2.03;
            }
            else
            {
                return 2.03 + 0.1172 * (S/T);
            }

        }
        else if (param == Kparam::BubblyBrine)
        {
            if (std::abs(T) < REAL_MIN_VAL(NumType))
            {
                return rho/IceConsts<NumType>::rho_i*2.11;
            }
            else
            {
                return rho/IceConsts<NumType>::rho_i*(2.11 - 0.011*T + (0.09*S/T));
            }
        }
        else if (param == Kparam::FreshIce)
        {
            return IceConsts<NumType>::k0_i;
        }
        else if (param == Kparam::FreshSnow)
        {
            return SnowConsts<NumType>::k0_s;
        }
        else
        {
            THERMO_ERR("Available Thermal Conductivity parametrizations: Untersteiner, BubblyBrine, FreshIce, FreshSnow!");
        }
    }

    // explicit instantiation
    template struct GenConsts<float>;
    template struct GenConsts<double>;

    template struct AirConsts<float>;
    template struct AirConsts<double>;

    template struct WaterConsts<float>;
    template struct WaterConsts<double>;

    template struct IceConsts<float>;
    template struct IceConsts<double>;

    template struct SnowConsts<float>;
    template struct SnowConsts<double>;

    template struct Params<float>;
    template struct Params<double>;
}