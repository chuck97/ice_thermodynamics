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
        // TODO: read configuration file and fill out Config object
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
    nlohmann::json j_gen = j_file["General constants"];

    if (!j_gen[GenConstNames.at("mu")].empty())
    {
        config_ptr->GenConsts.mu = (NumType)j_gen[GenConstNames.at("mu")];
    }
    else
    {
        std::cout << "Warning! \'" << GenConstNames.at("mu") << "\' not found, used default value instead!" << std::endl;
        config_ptr->GenConsts.mu = GenConsts<NumType>::sigma;
    }

    if (!j_gen[GenConstNames.at("sigma")].empty())
    {
        config_ptr->GenConsts.sigma = (NumType)j_gen[GenConstNames.at("sigma")];
    }
    else
    {
        std::cout << "Warning! \'" << GenConstNames.at("sigma") << "\' not found, used default value instead!" << std::endl;
        config_ptr->GenConsts.sigma = GenConsts<NumType>::sigma;
    }

    if (!j_gen[GenConstNames.at("T0")].empty())
    {
        config_ptr->GenConsts.T0 = (NumType)j_gen[GenConstNames.at("T0")];
    }
    else
    {
        std::cout << "Warning! \'" << GenConstNames.at("T0") << "\' not found, used default value instead!" << std::endl;
        config_ptr->GenConsts.T0 = GenConsts<NumType>::T0;
    }

    if (!j_gen[GenConstNames.at("emissivity")].empty())
    {
        config_ptr->GenConsts.emissivity = (NumType)j_gen[GenConstNames.at("emissivity")];
    }
    else
    {
        std::cout << "Warning! \'" << GenConstNames.at("emissivity") << "\' not found, used default value instead!" << std::endl;
        config_ptr->GenConsts.emissivity = GenConsts<NumType>::emissivity;
    }

    if (!j_gen[GenConstNames.at("C_sh")].empty())
    {
        config_ptr->GenConsts.C_sh = (NumType)j_gen[GenConstNames.at("C_sh")];
    }
    else
    {
        std::cout << "Warning! \'" << GenConstNames.at("C_sh") << "\' not found, used default value instead!" << std::endl;
        config_ptr->GenConsts.C_sh = GenConsts<NumType>::C_sh;
    }

    if (!j_gen[GenConstNames.at("C_lh")].empty())
    {
        config_ptr->GenConsts.C_lh = (NumType)j_gen[GenConstNames.at("C_lh")];
    }
    else
    {
        std::cout << "Warning! \'" << GenConstNames.at("C_lh") << "\' not found, used default value instead!" << std::endl;
        config_ptr->GenConsts.C_lh = GenConsts<NumType>::C_lh;
    }
    
    // TODO: parse Air constants

    // TODO: parse Water constants 

    // TODO: parse Ice constants 

    // TODO: parse Snow constants 

    // TODO: parse Ice parameterizations

    // TODO: parse Snow parameterizations 

    // TODO: parse 1D root solver parameters

    // TODO: parse Relaxation solver parameters

    // TODO: parse NaN values
    
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

template void ParseJson(const char* filepath, Config<float>* config_ptr);
template void ParseJson(const char* filepath, Config<double>* config_ptr);

template Config<float>* GetConfigConsts();
template Config<double>* GetConfigConsts();

}