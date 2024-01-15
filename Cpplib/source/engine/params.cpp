#include "params.hpp"

namespace icethermo
{
    template<typename NumType>
    NumType Params<NumType>::TempFusion(double S)
    {
        NumType mu = (Configured()) ? GetConfigConsts<NumType>()->GenConsts.mu : GenConsts<NumType>::mu;
        return -mu*S;
    }

    template<typename NumType>
    NumType Params<NumType>::Density(Dparam param, NumType T, NumType S)
    {
        NumType r_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.rho_i : IceConsts<NumType>::rho_i;
        NumType r_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.rho_s : SnowConsts<NumType>::rho_s;

        if (param == Dparam::SeaIce)
        {
            return r_i;
        }
        else if (param == Dparam::FreshIce)
        {
            return r_i;
        }
        else if (param == Dparam::FreshSnow)
        {
            return r_s;
        }
        else
        {
            THERMO_ERR("Available Density parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::EffCapacity(Cparam param, NumType T, NumType T_old, NumType S)
    {
        NumType c0i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.c0_i : IceConsts<NumType>::c0_i;
        NumType c0s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.c0_s : SnowConsts<NumType>::c0_s;
        NumType L0i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.L0_i : IceConsts<NumType>::L0_i;

        if (param == Cparam::SeaIce)
        {
            NumType Tf_i = Params<NumType>::TempFusion(S);
            if ((std::abs(T) < REAL_MIN_VAL(NumType)) or (std::abs(T_old) < REAL_MIN_VAL(NumType)))
            {
                return c0i;
            }
            else
            {
                return c0i - L0i*Tf_i/(T*T_old);
            }
        }
        else if (param == Cparam::FreshIce)
        {
            return c0i;
        }
        else if (param == Cparam::FreshSnow)
        {
            return c0s;
        }
        else
        {
            THERMO_ERR("Available Effective Heat Capacity parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::Enthalpy(Eparam param, NumType T, NumType S)
    {
        NumType c0i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.c0_i : IceConsts<NumType>::c0_i;
        NumType c0s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.c0_s : SnowConsts<NumType>::c0_s;
        NumType L0i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.L0_i : IceConsts<NumType>::L0_i;
        NumType Lf0 = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.L_f0 :SnowConsts<NumType>::L_f0;

        if (param == Eparam::SeaIce)
        {
            NumType Tf_i = Params<NumType>::TempFusion(S);
            NumType c_w = WaterConsts<NumType>::c_w;
            if (std::abs(T) < REAL_MIN_VAL(NumType))
            {
                return c0i*(T - Tf_i) + c_w*Tf_i;;
            }
            else
            {
                return c0i*(T - Tf_i) - L0i*(1.0 - Tf_i/T) + c_w*Tf_i;
            }
        }
        else if (param == Eparam::FreshIce)
        {
            return c0i*T - L0i;
        }
        else if (param == Eparam::FreshSnow)
        {
            return c0s*T - Lf0;
        }
        else
        {
            THERMO_ERR("Available Enthalpy parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::FusionHeat(Lparam param, NumType T, NumType S)
    {
        NumType c0i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.c0_i : IceConsts<NumType>::c0_i;
        NumType c0s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.c0_s : SnowConsts<NumType>::c0_s;
        NumType L0i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.L0_i : IceConsts<NumType>::L0_i;
        NumType Lf0 = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.L_f0 :SnowConsts<NumType>::L_f0;

        if (param == Lparam::SeaIce)
        {
            NumType Tf_i = Params<NumType>::TempFusion(S);
            NumType c_w = WaterConsts<NumType>::c_w;
            if (std::abs(T) < REAL_MIN_VAL(NumType))
            {
                return c0i*T - L0i;
            }
            else if (std::abs(T - Tf_i) < REAL_MIN_VAL(NumType))
            {
                return -L0i;
            }
            else
            {
                return c0i*(T - Tf_i) - L0i*(1.0 - Tf_i/T);
            }
        }
        else if (param == Lparam::FreshIce)
        {
            return c0i*T - L0i;
        }
        else if (param == Lparam::FreshSnow)
        {
            return c0s*T - Lf0;
        }
        else
        {
            THERMO_ERR("Available Fusion Heat parametrizations: SeaIce, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::Conductivity(Kparam param, NumType T, NumType S, NumType rho)
    {
        NumType r_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.rho_i : IceConsts<NumType>::rho_i;
        NumType k0_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.k0_i : IceConsts<NumType>::k0_i;
        NumType k0_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.k0_s : SnowConsts<NumType>::k0_s;

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
                return rho/r_i*2.11;
            }
            else
            {
                return rho/r_i*(2.11 - 0.011*T + (0.09*S/T));
            }
        }
        else if (param == Kparam::FreshIce)
        {
            return k0_i;
        }
        else if (param == Kparam::FreshSnow)
        {
            return k0_s;
        }
        else
        {
            THERMO_ERR("Available Thermal Conductivity parametrizations: Untersteiner, BubblyBrine, FreshIce, FreshSnow!");
        }
    }

    template<typename NumType>
    NumType Params<NumType>::Albedo(Aparam param, NumType T, NumType h, NumType Tfsurf)
    {
        NumType albedo_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.albedo_i : IceConsts<NumType>::albedo_i;
        NumType albedo_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.albedo_s : SnowConsts<NumType>::albedo_s;
        NumType albedo_dry_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.albedo_dry_i : IceConsts<NumType>::albedo_dry_i;
        NumType albedo_wet_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.albedo_wet_i : IceConsts<NumType>::albedo_wet_i;
        NumType albedo_dry_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.albedo_dry_s : SnowConsts<NumType>::albedo_dry_s;
        NumType albedo_wet_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.albedo_wet_s : SnowConsts<NumType>::albedo_wet_s;
        

        if (param == Aparam::ConstantIce)
        {
            return albedo_i;
        }
        else if (param == Aparam::ConstantSnow)
        {
            return albedo_s;
        }
        else if (param == Aparam::MeltingFreezingIce)
        {
            if (std::abs(T - Tfsurf) <= REAL_MIN_VAL(NumType))
            {
                return albedo_wet_i;
            }
            else
            {
                return albedo_dry_i;
            }
        }
        else if (param == Aparam::MeltingFreezingSnow)
        {
            if (std::abs(T - Tfsurf) <= REAL_MIN_VAL(NumType))
            {
                return albedo_wet_s;
            }
            else
            {
                return albedo_dry_s;
            }
        }
        else
        {
            THERMO_ERR("Available Albedo parametrizations: ConstantIce, ConstantSnow, MeltingFreezingIce, MeltingFreezingSnow!");
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