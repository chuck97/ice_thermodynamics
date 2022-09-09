#pragma once

#include "mesh.hpp"
#include "constants.hpp"

namespace icethermo
{
    template <typename NumType>
    class ThermoSolver
    {
    public:
        ThermoSolver(Mesh<NumType>* mesh_ice_,
                     Mesh<NumType>* mesh_snow_ = NULL
                     );

        virtual void Evaluate() = 0;

    
    protected:
    // ice and snow mesh
        Mesh<NumType>* mesh_ice;
        Mesh<NumType>* mesh_snow;
        
    };

    template<typename NumType>
    class Ice1D_Solver : public ThermoSolver<NumType>
    {
    public:
        Ice1D_Solver(Mesh<NumType>* mesh_ice_,
                     Dparam ice_rho_param_ = Dparam::SeaIce,
                     Kparam ice_k_param_ = Kparam::Untersteiner,
                     Cparam ice_c_eff_param_ = Cparam::SeaIce,
                     Eparam ice_E_param_ = Eparam::SeaIce,
                     Lparam ice_L_param_ = Lparam::SeaIce
                     );

        void Evaluate() override;

    private:
        // configuration variables
        Dparam ice_rho_param;
        Kparam ice_k_param;
        Cparam ice_c_eff_param;
        Eparam ice_E_param;
        Lparam ice_L_param;

        // mandatory prognostic mesh variables
        std::shared_ptr<std::vector<NumType>> Ti_cells;
        std::shared_ptr<std::vector<NumType>> dzi_cells;
        std::shared_ptr<std::vector<NumType>> Si_cells;
        std::shared_ptr<NumType> Ti_s, Ti_b;

        // auxilary mesh variables

    };

/*    
    template <typename NumType>
    class Ice1D_Snow1D_Solver : public Solver<NumType>
    {
    public:
        Ice1D_Snow1D_Solver(Mesh<NumType>* mesh_ice_,
                            Mesh<NumType>* mesh_snow_,
                            Dparam ice_rho_param = Dparam::SeaIce,
                            Kparam ice_k_param = Kparam::Untersteiner,
                            Cparam ice_effc_param = Cparam::SeaIce,
                            Eparam ice_E_param = Eparam::SeaIce,
                            Lparam ice_L_param = Lparam::SeaIce,
                            Kparam snow_rho_param = Dparam::FreshSnow,
                            Kparam snow_k_param = Kparam::FreshSnow,
                            Cparam snow_effc_param = Cparam::FreshSnow,
                            Eparam snow_E_param = Eparam::FreshSnow,
                            Lparam snow_L_param = Lparam::FreshSnow
                            );

    private:

        // mandatory prognostic mesh variables
        bool is_snow;
        std::vector<NumType>& Ti_cells, Ts_cells;
        std::vector<NumType>& dzi_cells, dzs_cells;
        std::vector<NumType>& Si_cells, Ss_cells;
        NumType& Tis_interface, T_surface, Ti_bottom;

        // auxilary mesh variables


    };

    template <typename NumType>
    class Ice1D_Snow0D_Solver : public Solver<NumType>
    {

    };

    template <typename NumType>
    class Ice0D_Snow0D_Solver : public Solver<NumType>
    {

    };
*/
}