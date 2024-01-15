#include "config.hpp"

namespace icethermo
{

// ! define foo pointer to configured consts class
void* config_consts = NULL;

// ! configure constants using .json file
template<typename NumType>
void Configure(const char* filepath)
{
    if (config_consts != NULL)
    {
        delete (Config<NumType>*)config_consts;
    }

    Config<NumType>* config_ptr = new Config<NumType>;
    #ifdef USE_JSON_OUTPUT
        // ? read configuration file and fill out Config object
        ParseJson<NumType>(filepath, config_ptr);
    #else
        THERMO_ERR("Can't configure constants without nlohmann::json library!");
    #endif
        
    config_consts = config_ptr;
}

#ifdef USE_JSON_OUTPUT
// ! parse .json configuration file and fill out constants class
template<typename NumType>
void ParseJson(const char* filepath, Config<NumType>* config_ptr)
{
    std::ifstream ifs(filepath);
    if (!ifs.is_open())
    {
        std::string str_filename(filepath);
        THERMO_ERR("failed to open file \"" + str_filename + "\"!");
    }

    nlohmann::json j_file = nlohmann::json::parse(ifs);

    // ? parse general constants
    if (!j_file["General constants"].empty())
    {
        nlohmann::json j_gen = j_file["General constants"];

        if (!j_gen[GenConstNames.at("mu")].empty())
        {
            config_ptr->GenConsts.mu = (NumType)j_gen[GenConstNames.at("mu")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << GenConstNames.at("mu") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_gen[GenConstNames.at("sigma")].empty())
        {
            config_ptr->GenConsts.sigma = (NumType)j_gen[GenConstNames.at("sigma")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << GenConstNames.at("sigma") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_gen[GenConstNames.at("T0")].empty())
        {
            config_ptr->GenConsts.T0 = (NumType)j_gen[GenConstNames.at("T0")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << GenConstNames.at("T0") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_gen[GenConstNames.at("emissivity")].empty())
        {
            config_ptr->GenConsts.emissivity = (NumType)j_gen[GenConstNames.at("emissivity")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << GenConstNames.at("emissivity") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_gen[GenConstNames.at("C_sh")].empty())
        {
            config_ptr->GenConsts.C_sh = (NumType)j_gen[GenConstNames.at("C_sh")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << GenConstNames.at("C_sh") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_gen[GenConstNames.at("C_lh")].empty())
        {
            config_ptr->GenConsts.C_lh = (NumType)j_gen[GenConstNames.at("C_lh")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << GenConstNames.at("C_lh") << "\' not found, used default value instead!" << std::endl;
        }
    }

    // ? parse Air constants
    if (!j_file["Air constants"].empty())
    {
        nlohmann::json j_air = j_file["Air constants"];
        if (!j_air[AirConstNames.at("rho_a")].empty())
        {
            config_ptr->AirConsts.rho_a = (NumType)j_air[AirConstNames.at("rho_a")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << AirConstNames.at("rho_a") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_air[AirConstNames.at("cp_a")].empty())
        {
            config_ptr->AirConsts.cp_a = (NumType)j_air[AirConstNames.at("cp_a")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << AirConstNames.at("cp_a") << "\' not found, used default value instead!" << std::endl;
        }
    }

    // ? parse Water constants 
    if (!j_file["Water constants"].empty())
    {
        nlohmann::json j_water = j_file["Water constants"];

        if (!j_water[WaterConstNames.at("rho_w")].empty())
        {
            config_ptr->WaterConsts.rho_w = (NumType)j_water[WaterConstNames.at("rho_w")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << WaterConstNames.at("rho_w") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_water[WaterConstNames.at("cp_w")].empty())
        {
            config_ptr->WaterConsts.cp_w = (NumType)j_water[WaterConstNames.at("cp_w")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << WaterConstNames.at("cp_w") << "\' not found, used default value instead!" << std::endl;
        }
    }


    // ? parse Ice constants 
    if (!j_file["Ice constants"].empty())
    {
        nlohmann::json j_ice = j_file["Ice constants"];

        if (!j_ice[IceConstNames.at("c0_i")].empty())
        {
            config_ptr->IceConsts.c0_i = (NumType)j_ice[IceConstNames.at("c0_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("c0_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("rho_i")].empty())
        {
            config_ptr->IceConsts.rho_i = (NumType)j_ice[IceConstNames.at("rho_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("rho_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("L0_i")].empty())
        {
            config_ptr->IceConsts.L0_i = (NumType)j_ice[IceConstNames.at("L0_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("L0_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("kappa_i")].empty())
        {
            config_ptr->IceConsts.kappa_i = (NumType)j_ice[IceConstNames.at("kappa_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("kappa_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("albedo_i")].empty())
        {
            config_ptr->IceConsts.albedo_i = (NumType)j_ice[IceConstNames.at("albedo_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("albedo_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("i0_i")].empty())
        {
            config_ptr->IceConsts.i0_i = (NumType)j_ice[IceConstNames.at("i0_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("i0_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("k0_i")].empty())
        {
            config_ptr->IceConsts.k0_i = (NumType)j_ice[IceConstNames.at("k0_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("k0_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("albedo_dry_i")].empty())
        {
            config_ptr->IceConsts.albedo_dry_i = (NumType)j_ice[IceConstNames.at("albedo_dry_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("albedo_dry_i") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_ice[IceConstNames.at("albedo_wet_i")].empty())
        {
            config_ptr->IceConsts.albedo_wet_i = (NumType)j_ice[IceConstNames.at("albedo_wet_i")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << IceConstNames.at("albedo_wet_i") << "\' not found, used default value instead!" << std::endl;
        }
    }

    // ? parse Snow constants 
    if (!j_file["Snow constants"].empty())
    {
        nlohmann::json j_snow = j_file["Snow constants"];

        if (!j_snow[SnowConstNames.at("c0_s")].empty())
        {
            config_ptr->SnowConsts.c0_s = (NumType)j_snow[SnowConstNames.at("c0_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("c0_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("rho_s")].empty())
        {
            config_ptr->SnowConsts.rho_s = (NumType)j_snow[SnowConstNames.at("rho_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("rho_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("L0_s")].empty())
        {
            config_ptr->SnowConsts.L0_s = (NumType)j_snow[SnowConstNames.at("L0_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("L0_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("kappa_s")].empty())
        {
            config_ptr->SnowConsts.kappa_s = (NumType)j_snow[SnowConstNames.at("kappa_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("kappa_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("albedo_s")].empty())
        {
            config_ptr->SnowConsts.albedo_s = (NumType)j_snow[SnowConstNames.at("albedo_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("albedo_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("i0_s")].empty())
        {
            config_ptr->SnowConsts.i0_s = (NumType)j_snow[SnowConstNames.at("i0_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("i0_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("k0_s")].empty())
        {
            config_ptr->SnowConsts.k0_s = (NumType)j_snow[SnowConstNames.at("k0_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("k0_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("albedo_dry_s")].empty())
        {
            config_ptr->SnowConsts.albedo_dry_s = (NumType)j_snow[SnowConstNames.at("albedo_dry_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("albedo_dry_s") << "\' not found, used default value instead!" << std::endl;
        }

        if (!j_snow[SnowConstNames.at("albedo_wet_s")].empty())
        {
            config_ptr->SnowConsts.albedo_wet_s = (NumType)j_snow[SnowConstNames.at("albedo_wet_s")];
        }
        else
        {
            std::cout << "IceThermo Warning! \'" << SnowConstNames.at("albedo_wet_s") << "\' not found, used default value instead!" << std::endl;
        }
    }

    // ? parse 1D root solver parameters
    if (!j_file["1D root solver"].empty())
    {
        nlohmann::json j_solver1d = j_file["1D root solver"];

        if (!j_solver1d["Type"].empty())
        {
            for (auto key : j_solver1d["Type"])
            {
                if (j_solver1d["Type"][key])
                {
                    config_ptr->Solver1d.type = Solver1dNameToType.at(key);
                    break;
                }
            }
        }
        else
        {
            std::cout << "IceThermo Warning! 1D root solver::\"Type\" not found, used default value instead!" << std::endl;
        }

        if (!j_solver1d["Maximum number of iterations"].empty())
        {
            config_ptr->Solver1d.max_nits = j_solver1d["Maximum number of iterations"];
        }
        else
        {
            std::cout << "IceThermo Warning! 1D root solver::\"Maximum number of iterations\" not found, used default value instead!" << std::endl;
        }

        if (!j_solver1d["Residual accuracy"].empty())
        {
            config_ptr->Solver1d.residual = (NumType)j_solver1d["Residual accuracy"];
        }
        else
        {
            std::cout << "IceThermo Warning! 1D root solver::\"Residual accuracy\" not found, used default value instead!" << std::endl;
        }
    }

    
    // TODO: parse Relaxation solver parameters

    // TODO: parse NaN values

    // TODO: parse Ice parameterizations

    // TODO: parse Snow parameterizations 
    
}
#endif

// ! check out if configuration has already been made
bool Configured()
{
    return (config_consts != NULL);
}

template<typename NumType>
Config<NumType>* GetConfigConsts()
{
    if (config_consts != NULL)
    {
        return (Config<NumType>*)config_consts;
    }
    else
    {
        THERMO_ERR("Configuration hasn't been done yet, can't get configured constants!");
    }
}

// ! explicit instantiation
template struct Config<float>;
template struct Config<double>;

template void Configure<float>(const char* filepath);
template void Configure<double>(const char* filepath);

#ifdef USE_JSON_OUTPUT
template void ParseJson(const char* filepath, Config<float>* config_ptr);
template void ParseJson(const char* filepath, Config<double>* config_ptr);
#endif

template Config<float>* GetConfigConsts();
template Config<double>* GetConfigConsts();

}