#pragma once

#include <vector>
#include <functional>
#include <limits>

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
    
    // base thermodynamics solver class
    template <typename NumType>
    class ThermoSolver
    {
    public:
        // class constructor using initial mesh
        ThermoSolver(Mesh<NumType>* mesh_ice_,
                     Mesh<NumType>* mesh_snow_,
                     NumType time_step_,
                     ApproxOrder grad_approx_order_ = ApproxOrder::first,
                     Kparam ice_k_param_ = Kparam::FreshIce,
                     Cparam ice_c_eff_param_ = Cparam::FreshIce,
                     Eparam ice_E_param_ = Eparam::FreshIce,
                     Lparam ice_L_param_ = Lparam::FreshIce,
                     Kparam snow_k_param_ = Kparam::FreshSnow,
                     Cparam snow_c_eff_param_ = Cparam::FreshSnow,
                     Eparam snow_E_param_ = Eparam::FreshSnow,
                     Lparam snow_L_param_ = Lparam::FreshSnow
                     );

        // virtual Evaluation function
        virtual void Evaluate() = 0; 
    
    protected:
        // time step
        NumType time_step;

        // ice and snow mesh
        Mesh<NumType>* mesh_ice;
        Mesh<NumType>* mesh_snow;

        // upwards, downwards and short-wave radiation forcing
        FuncPtr<NumType> F_up;
        FuncPtr<NumType> F_down;
        FuncPtr<NumType> F_sw;

        // mandatory ice prognostic variables
        std::shared_ptr<std::vector<NumType>> Ti_cells;
        std::shared_ptr<std::vector<NumType>> dzi_cells;
        std::shared_ptr<std::vector<NumType>> Si_cells;
        std::shared_ptr<std::vector<NumType>> rhoi_cells;
        std::shared_ptr<NumType> Ti_s, Ti_b, So;

        // mandatory snow prognostic variables
        std::shared_ptr<std::vector<NumType>> Ts_cells;
        std::shared_ptr<std::vector<NumType>> dzs_cells;
        std::shared_ptr<std::vector<NumType>> rhos_cells;
        std::shared_ptr<NumType> Ts_s, Ts_b;

        // parametrizations
        Kparam ice_k_param;
        Cparam ice_c_eff_param;
        Eparam ice_E_param;
        Lparam ice_L_param;

        Kparam snow_k_param;
        Cparam snow_c_eff_param;
        Eparam snow_E_param;
        Lparam snow_L_param;

        // the order of gradient approximation
        ApproxOrder grad_approx_order;

    protected:
        // functions for checking mesh consistency
        virtual void CheckMeshConsistency(Mesh<NumType>* mesh_ice,
                                          Mesh<NumType>* mesh_snow) = 0;

    public:
        // update forcing
        void UpdateForcing(FuncPtr<NumType> F_up_ = [](NumType a){return (NumType)0.0;},
                           FuncPtr<NumType> F_down_ = [](NumType a){return (NumType)0.0;},
                           FuncPtr<NumType> F_sw_ = [](NumType a){return (NumType)0.0;});

        // update mesh
        virtual void UpdateMesh(Mesh<NumType>* mesh_ice_,
                                Mesh<NumType>* mesh_snow_) = 0;
    
    protected:

        // ##### auxilary functions #####

        // update 1D thickness
        std::vector<NumType> Update_dz(const std::vector<NumType>& dz_cells_old,
                                       NumType omega_down,
                                       NumType omega_up);

        // update 0D thickness
        NumType Update_dz(NumType dz_old,
                          NumType omega_down,
                          NumType omega_up);

        // find temperature consistent with boundary conditions
        NumType T_from_BC(const std::vector<NumType>& T_cells,
                          const std::vector<NumType>& dz_cells,
                          const std::vector<NumType>& salinity_cells,
                          const std::vector<NumType>& rho_cells,
                          NumType omega_value,
                          bool is_ice,
                          bool is_surface);
        
        // find omega value consistent with boundary conditions
        NumType W_from_BC(NumType T_bnd,
                          const std::vector<NumType>& T_cells,
                          const std::vector<NumType>& dz_cells,
                          const std::vector<NumType>& salinity_cells,
                          const std::vector<NumType>& rho_cells,
                          bool is_ice,
                          bool is_surface);
        
        // construct tridiagonal discrete advection-diffusion matrix and rhs
        FourVecs<NumType> Assemble_advdiff_martix_rhs(const std::vector<NumType>& T_cells_prev,
                                                      const std::vector<NumType>& T_cells_old,
                                                      NumType T_up_new,
                                                      NumType T_up_prev,
                                                      NumType T_up_old,
                                                      NumType T_down_new,
                                                      NumType T_down_prev,
                                                      NumType T_down_old, 
                                                      NumType omega_down,
                                                      NumType omega_up, 
                                                      const std::vector<NumType>& dz_cells_new,
                                                      const std::vector<NumType>& dz_cells_old,
                                                      const std::vector<NumType>& salinity_cells,
                                                      const std::vector<NumType>& rho_cells,
                                                      const std::vector<NumType>& radiation_nodes,
                                                      bool is_ice);
        
        // ice freezing mode (for 1d profile)
        ThreeVecs<NumType> sea_ice_freezing_1d(NumType T_ib,
                                               const std::vector<NumType>& T_cells,
                                               NumType T_is,
                                               const std::vector<NumType>& dz_cells,
                                               const std::vector<NumType>& salinity_cells,
                                               const std::vector<NumType>& rho_cells,
                                               bool is_radiation = true,
                                               int max_n_its = 50,
                                               NumType tol = 1e-6);
    };

    template<typename NumType>
    class SeaIce1D_Solver : public ThermoSolver<NumType>
    {
    public:
        SeaIce1D_Solver(Mesh<NumType>* mesh_ice_,
                        NumType time_step,
                        ApproxOrder grad_approx_order_ = ApproxOrder::first,
                        Kparam ice_k_param_ = Kparam::Untersteiner,
                        Cparam ice_c_eff_param_ = Cparam::FreshIce,
                        Eparam ice_E_param_ = Eparam::FreshIce,
                        Lparam ice_L_param_ = Lparam::FreshIce
                        );

        void Evaluate() override;

    private:
        void CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                  Mesh<NumType>* snow_mesh = NULL) override;
    
    public:
        void UpdateMesh(Mesh<NumType>* mesh_ice_,
                        Mesh<NumType>* mesh_snow_ = 0) override;
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