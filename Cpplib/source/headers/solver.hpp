#pragma once

#include <vector>
#include <functional>

#include "defines.hpp"
#include "mesh.hpp"
#include "constants.hpp"
#include "matvec.hpp"
#include "tools.hpp"

namespace icethermo
{
    // order of approximation enum
    enum class ApproxOrder
    {
        first, 
        second
    };

    template <typename NumType>
    class ThermoSolver
    {
    public:
        // class constructor using initial mesh
        ThermoSolver(Mesh<NumType>* mesh_ice_,
                     Mesh<NumType>* mesh_snow_,
                     Dparam ice_rho_param_ = Dparam::FreshIce,
                     Kparam ice_k_param_ = Kparam::FreshIce,
                     Cparam ice_c_eff_param_ = Cparam::FreshIce,
                     Eparam ice_E_param_ = Eparam::FreshIce,
                     Lparam ice_L_param_ = Lparam::FreshIce,
                     Dparam snow_rho_param_ = Dparam::FreshSnow,
                     Kparam snow_k_param_ = Kparam::FreshSnow,
                     Cparam snow_c_eff_param_ = Cparam::FreshSnow,
                     Eparam snow_E_param_ = Eparam::FreshSnow,
                     Lparam snow_L_param_ = Lparam::FreshSnow
                     );

        // virtual Evaluation function
        virtual void Evaluate() = 0;

    
    protected:
        // ice and snow mesh
        Mesh<NumType>* mesh_ice;
        Mesh<NumType>* mesh_snow;

        // mandatory prognostic variables
        std::shared_ptr<std::vector<NumType>> Ti_cells;
        std::shared_ptr<std::vector<NumType>> dzi_cells;
        std::shared_ptr<std::vector<NumType>> Si_cells;
        std::shared_ptr<NumType> Ti_s, Ti_b;

        // parametrizations
        Dparam ice_rho_param;
        Kparam ice_k_param;
        Cparam ice_c_eff_param;
        Eparam ice_E_param;
        Lparam ice_L_param;

        Dparam snow_rho_param;
        Kparam snow_k_param;
        Cparam snow_c_eff_param;
        Eparam snow_E_param;
        Lparam snow_L_param;

    protected:

        // ##### auxilary functions #####

        // update 1D thickness
        std::vector<NumType> Update_dz(const std::vector<NumType>& dz_cells_old,
                                       NumType omega_down,
                                       NumType omega_up,
                                       NumType time_step);

        // update 0D thickness
        NumType Update_dz(NumType dz_old,
                          NumType omega_down,
                          NumType omega_up,
                          NumType time_step);

        // find temperature consistent with boundary conditions
        NumType T_from_BC(const std::vector<NumType>& T_cells,
                          const std::vector<NumType>& dz_cells,
                          const std::vector<NumType>& salinity_cells,
                          NumType omega_value,
                          FuncPtr<NumType> F,
                          bool is_surface,
                          ApproxOrder grad_approx_order = ApproxOrder::first,
                          Dparam dparam = Dparam::FreshIce,
                          Kparam kparam = Kparam::FreshIce,
                          Lparam Lparam = Lparam::FreshIce);
        
        // find omega value consistent with boundary conditions
        NumType W_from_BC(NumType T_bnd,
                          const std::vector<NumType>& T_cells,
                          const std::vector<NumType>& dz_cells,
                          const std::vector<NumType>& salinity_cells,
                          FuncPtr<NumType> F,
                          bool is_surface,
                          ApproxOrder grad_approx_order,
                          Dparam dparam = Dparam::FreshIce,
                          Kparam kparam = Kparam::FreshIce,
                          Lparam Lparam = Lparam::FreshIce);
        
        // construct tridiagonal discrete advection-diffusion matrix and rhs
        FourVecs<NumType> Assemble_advdiff_martix_rhs(const std::vector<NumType>& T_cells_prev, const std::vector<NumType>& T_cells_old,
                                                      NumType T_up_new, NumType T_up_prev, NumType T_up_old,
                                                      NumType T_down_new, NumType T_down_prev, NumType T_down_old, 
                                                      NumType omega_down, NumType omega_up, 
                                                      const std::vector<NumType>& dz_cells_new, const std::vector<NumType>& dz_cells_old,
                                                      const std::vector<NumType>& salinity_cells,
                                                      const std::vector<NumType>& radiation_nodes,
                                                      NumType time_step,
                                                      Dparam dparam = Dparam::FreshIce,
                                                      Kparam kparam = Kparam::FreshIce,
                                                      Eparam Eparam = Eparam::FreshIce,
                                                      Cparam cparam = Cparam::FreshIce);

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