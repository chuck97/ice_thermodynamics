#pragma once

#include <vector>
#include <functional>
#include <limits>

#include "defines.hpp"
#include "mesh.hpp"
#include "constants.hpp"
#include "params.hpp"
#include "matvec.hpp"
#include "tools.hpp"
#include "config.hpp"

namespace icethermo
{
    // order of approximation enum
    enum class ApproxOrder
    {
        first, 
        second
    };

    // snow-ice transition strategy
    enum class SnowIceTransition
    {
        None,
        SnowAging
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
                     NumType min_ice_thick_ = (NumType)0.0,
                     NumType min_snow_thick_ = (NumType)0.0,
                     ApproxOrder grad_approx_order_ = ApproxOrder::first,
                     bool is_radiation = false,
                     bool is_sublimation_ = false,
                     bool is_verbose_ = false,
                     Kparam ice_k_param_ = Kparam::FreshIce,
                     Cparam ice_c_eff_param_ = Cparam::FreshIce,
                     Eparam ice_E_param_ = Eparam::FreshIce,
                     Lparam ice_L_param_ = Lparam::FreshIce,
                     Aparam ice_albedo_param_ = Aparam::ConstantIce,
                     Kparam snow_k_param_ = Kparam::FreshSnow,
                     Cparam snow_c_eff_param_ = Cparam::FreshSnow,
                     Eparam snow_E_param_ = Eparam::FreshSnow,
                     Lparam snow_L_param_ = Lparam::FreshSnow,
                     Aparam snow_albedo_param_ = Aparam::ConstantSnow,
                     SnowIceTransition si_transition_mode_ = SnowIceTransition::None);

        // destructor
        ~ThermoSolver();

        // virtual Evaluation function
        virtual void Evaluate() = 0; 
    
    protected:
        // time step
        NumType time_step;

        // treshold values for thicknesses
        NumType min_ice_thick;
        NumType min_snow_thick;

        // ice and snow mesh
        Mesh<NumType>* mesh_ice = NULL;
        Mesh<NumType>* mesh_snow = NULL;


        // upwards, downwards, short-wave radiation, latent heat flux, precipitation rate, atmosphere temperature
        
        // total atmosphere flux
        FuncPtr<NumType> F_up = [](NumType T){return 0.0;};

        // total ground/ocean flux
        FuncPtr<NumType> F_down = [](NumType T){return 0.0;};

        // descending short-wave radiation
        FuncPtr<NumType> F_sw = [](NumType T){return 0.0;};

        // descending long-wave radiation
        FuncPtr<NumType> F_lw = [](NumType T){return 0.0;};

        // atmosphere latent heat flux
        FuncPtr<NumType> F_lh = [](NumType T){return 0.0;};
        
        // atmosphere sensible heat flux
        FuncPtr<NumType> F_sh = [](NumType T){return 0.0;};

        // ascending long-wave ice/snow emitting
        FuncPtr<NumType> F_lwi = [](NumType T){return 0.0;};

        // precipitation heat flux
        FuncPtr<NumType> F_P = [](NumType T){return 0.0;};


        std::shared_ptr<NumType> prec_rate = NULL;      // in snow-water equavalent (mm/s)
        std::shared_ptr<NumType> atm_temp = NULL;       // in deg Cel
        std::shared_ptr<NumType> atm_press = NULL;      // in Pascals
        std::shared_ptr<NumType> atm_humid = NULL;      // in g/kg
        std::shared_ptr<NumType> abs_wind_speed = NULL; // m/s
        std::shared_ptr<NumType> sh_trans_coeff = NULL; // -
        std::shared_ptr<NumType> lh_trans_coeff = NULL; // -
        std::shared_ptr<NumType> atm_dens = NULL;       // kg/m^3

        // mandatory ice prognostic variables
        std::shared_ptr<std::vector<NumType>> Ti_cells = NULL;
        std::shared_ptr<std::vector<NumType>> dzi_cells = NULL;
        std::shared_ptr<std::vector<NumType>> Si_cells = NULL;
        std::shared_ptr<std::vector<NumType>> rhoi_cells = NULL;
        std::shared_ptr<NumType> Ti_s = NULL;
        std::shared_ptr<NumType> Ti_b = NULL;

        // mandatory snow prognostic variables
        std::shared_ptr<std::vector<NumType>> Ts_cells = NULL;
        std::shared_ptr<std::vector<NumType>> dzs_cells = NULL;
        std::shared_ptr<std::vector<NumType>> rhos_cells = NULL;
        std::shared_ptr<NumType> Ts_s = NULL;
        std::shared_ptr<NumType> Ts_b= NULL;

        // parametrizations
        Kparam ice_k_param;
        Cparam ice_c_eff_param;
        Eparam ice_E_param;
        Lparam ice_L_param;
        Aparam ice_albedo_param;

        Kparam snow_k_param;
        Cparam snow_c_eff_param;
        Eparam snow_E_param;
        Lparam snow_L_param;
        Aparam snow_albedo_param;

        // switchers for radiation and sublimation
        bool is_radiation;
        bool is_sublimation;

        // the order of gradient approximation
        ApproxOrder grad_approx_order;

        // snow->ice transition mode
        SnowIceTransition si_transition_mode;

        // verbose output?
        bool is_verbose;

    protected:
        // functions for checking mesh consistency
        virtual void CheckMeshConsistency(Mesh<NumType>* mesh_ice,
                                          Mesh<NumType>* mesh_snow) = 0;

    public:
        // update upper (atmosphere) flux
        void UpdateUpperFlux(FuncPtr<NumType> F_up_);
        void UpdateUpperFlux(); // assign default total flux

        // update lower (ocean/soil flux)
        void UpdateLowerFlux(FuncPtr<NumType> F_down_);

        // update short-wave radiation (used to compute penetrating radiation)
        void UpdateShortWaveRadiation(FuncPtr<NumType> F_sw_);

        // update long-wave radiation from atmosphere
        void UpdateLongWaveRadiation(FuncPtr<NumType> F_lw_);

        // update latent heat flux 
        void UpdateLatentHeatFlux(FuncPtr<NumType> F_lh_);
        void UpdateLatentHeatFlux(); // assign default latent heat flux parameterization

        // update sensible heat flux
        void UpdateSensibleHeatFlux(FuncPtr<NumType> F_sh_);
        void UpdateSensibleHeatFlux(); // assign default sensible heat flux parameterization

        // update sensible heat flux
        void UpdatePrecipitationHeatFlux(FuncPtr<NumType> F_p_);
        void UpdatePrecipitationHeatFlux(); // assign default precipitation heat flux parameterization

        // update sensible heat flux
        void UpdateEmittingHeatFlux(FuncPtr<NumType> F_lwi_);
        void UpdateEmittingHeatFlux(); // assign default surface emitting heat flux parameterization

        // update precipitation rate (m s-1)
        void UpdatePrecipitationRate(NumType prec_rate_mm_sm1_);

        // update atmosphere temperature (deg C)
        void UpdateAtmosphereTemperature(NumType atm_temp_);

        // update atmosphere pressure (Pascals = N/m^2)
        void UpdateAtmospherePressure(NumType atm_press_);

        // update atmosphere specific humidity (g/kg)
        void UpdateAirSpecificHumidity(NumType atm_humid_);

        // update absolute wind speed (m/s)
        void UpdateAbsWindSpeed(NumType abs_wind_speed_);

        // update sensible heat transfer coefficient (-)
        void UpdateShTransCoeff(NumType sh_coeff_);

        // update latent heat transfer coefficient (-)
        void UpdateLhTransCoeff(NumType lh_coeff_);

        // update atmosphere density
        void UpdateAirDensity(NumType atm_dens_);

        // add precipitation
        virtual void AddPrecipitation();
        
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
        NumType Update_dz_0D(NumType dz_old,
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

        // find temperature consistent with boundary conditions (0D)
        NumType T_from_BC_0D(FuncPtr<NumType> T_op_interface,
                             NumType thickness,
                             NumType k_value,
                             NumType sal_value,
                             NumType omega_value,
                             NumType rho_value,
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

        // find omega value consistent with boundary conditions (0D)
        NumType W_from_BC_0D(NumType T_bnd,
                             NumType T_op_interface,
                             NumType thickness,
                             NumType k_value,
                             NumType sal_value,
                             NumType rho_value,
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
        
        // compute nodal radiation values (output - pair {radiation nodes ice, radiation nodes snow (optional)})
        std::pair<std::vector<NumType>, std::vector<NumType>> Compute_radiation_nodes(NumType F_sw_value,
                                                                                      NumType albedo_ice,
                                                                                      NumType i0_ice,
                                                                                      NumType kappa_ice,
                                                                                      const std::vector<NumType>& dz_cells_ice,
                                                                                      NumType albedo_snow = 0.0,
                                                                                      NumType i0_snow = 0.0,
                                                                                      NumType kappa_snow = 0.0,
                                                                                      const std::vector<NumType>& dz_cells_snow = std::vector<NumType>{0.0});
        
        
        // !! glacier freezing mode (for 1d profile) !!
        /* 
            output - 1d-ice values
            
            1d-ice values:
                new ice base temperature
                vector of new ice cells temperatures
                new ice surface temperature
                vector of new ice cells thicknesses
        */
        FourVecs<NumType> glacier1d_freezing(NumType T_ib,
                                             const std::vector<NumType>& T_cells,
                                             NumType T_is,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             const std::vector<NumType>& rho_cells,
                                             int max_n_its = MAX_RELAXATION_ITS,
                                             NumType tol = 1e-6);

        // !! galcier melting mode (for 1d profile) !!
        /* 
            output - 1d-ice values
            
            1d-ice values:
                the ice base temperature
                vector of new ice cells temperatures
                vector of new ice cells thicknesses
        */
        ThreeVecs<NumType> glacier1d_melting(NumType T_ib,
                                             NumType T_is,
                                             const std::vector<NumType>& T_cells,
                                             NumType T_is_old,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             const std::vector<NumType>& rho_cells,
                                             int max_n_its = MAX_RELAXATION_ITS,
                                             NumType tol = RELAXATION_SOLVER_ACCUR);

    };

    //   ########################################################
    //   ############# SEA ICE THERMO SOLVER  CLASS #############
    //   ########################################################
    template<typename NumType>
    class SeaIce_Solver : public ThermoSolver<NumType>
    {
    public:
        // constructor
        SeaIce_Solver(Mesh<NumType>* mesh_ice_,
                      Mesh<NumType>* mesh_snow_,
                      NumType time_step_,
                      NumType min_ice_thick_,
                      NumType min_snow_thick_,
                      ApproxOrder grad_approx_order_ = ApproxOrder::first,
                      bool is_radiation = true,
                      bool is_sublimation_ = true,
                      bool is_verbose_ = true,
                      Kparam ice_k_param_ = Kparam::FreshIce,
                      Cparam ice_c_eff_param_ = Cparam::FreshIce,
                      Eparam ice_E_param_ = Eparam::FreshIce,
                      Lparam ice_L_param_ = Lparam::FreshIce,
                      Aparam ice_albedo_param_ = Aparam::ConstantIce,
                      Kparam snow_k_param_ = Kparam::FreshSnow,
                      Cparam snow_c_eff_param_ = Cparam::FreshSnow,
                      Eparam snow_E_param_ = Eparam::FreshSnow,
                      Lparam snow_L_param_ = Lparam::FreshSnow,
                      Aparam snow_albedo_param_ = Aparam::ConstantSnow,
                      SnowIceTransition si_transition_mode_ = SnowIceTransition::None);
    
        // virtual Evaluation function
        virtual void Evaluate() = 0; 
    
    protected:
        std::shared_ptr<NumType> So = NULL; // psu
        std::shared_ptr<NumType> omega_ib_expl_mass = NULL; // m s-1

    public:
        // update ocean salinity
        void UpdateOceanSalinity(NumType ocn_sal_);

        // update explicit mass flux from ocean to ice
        void UpdateOceanIceMassFlux(NumType mass_flux_);

        // !! sea-ice freezing mode (for 1d profile) !!
        /* 
            output - 1d-ice values
            
            1d-ice values:
                vector of new ice cells temperatures
                new ice surface temperature
                vector of new ice cells thicknesses
        */
        ThreeVecs<NumType> seaice1d_freezing(NumType T_ib,
                                             const std::vector<NumType>& T_cells,
                                             NumType T_is,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             const std::vector<NumType>& rho_cells,
                                             int max_n_its = MAX_RELAXATION_ITS,
                                             NumType tol = RELAXATION_SOLVER_ACCUR);
        
        
        // !! sea-ice melting mode (for 1d profile) !!
        /* 
            output - 1d-ice values
            
            1d-ice values:
                vector of new ice cells temperatures
                vector of new ice cells thicknesses
        */
        TwoVecs<NumType> seaice1d_melting(NumType T_ib,
                                          NumType T_is,
                                          const std::vector<NumType>& T_cells,
                                          NumType T_is_old,
                                          const std::vector<NumType>& dz_cells,
                                          const std::vector<NumType>& salinity_cells,
                                          const std::vector<NumType>& rho_cells,
                                          int max_n_its = MAX_RELAXATION_ITS,
                                          NumType tol = RELAXATION_SOLVER_ACCUR);

        // !! sea-ice freezing mode with snow (for 1d-ice and 0d-snow profile) !!
        /* 
            output - pair of 1d-ice and 0d-snow values
            
            1d-ice values:
                vector of new ice cells temperatures
                new temperature of ice-snow interface
                vector of new ice cell thicknesses
            
            0d-snow values:
                new snow surface temperature
                new snow thickness
        */
        std::pair<ThreeVecs<NumType>, TwoVecs<NumType>> seaice1d_snow0d_freezing(NumType T_ib,
                                                                                 const std::vector<NumType>& T_i_cells,
                                                                                 NumType T_is,
                                                                                 const std::vector<NumType>& dz_i_cells,
                                                                                 const std::vector<NumType>& salinity_i_cells,
                                                                                 const std::vector<NumType>& rho_i_cells,
                                                                                 NumType T_ss,
                                                                                 NumType thickness_s,
                                                                                 NumType rho_s,
                                                                                 NumType precipitation_rate,
                                                                                 NumType atm_temperature,
                                                                                 int max_n_its = MAX_RELAXATION_ITS,
                                                                                 NumType tol = RELAXATION_SOLVER_ACCUR);
        
        // !! snow melting mode with sea-ice (for 1d-ice and 0d-snow profile) !!
        /* 
            output - pair of 1d-ice and 0d-snow values
            
            1d-ice values:
                vector of new ice cells temperatures
                new temperature of ice-snow interface
                vector of new ice cell thicknesses
            
            0d-snow values:
                new snow thickness
        */
        std::pair<ThreeVecs<NumType>, NumType> seaice1d_snow0d_melting(NumType T_ib,
                                                                       const std::vector<NumType>& T_i_cells,
                                                                       NumType T_is,
                                                                       const std::vector<NumType>& dz_i_cells,
                                                                       const std::vector<NumType>& salinity_i_cells,
                                                                       const std::vector<NumType>& rho_i_cells,
                                                                       NumType T_ss_old,
                                                                       NumType thickness_s,
                                                                       NumType rho_s,
                                                                       NumType precipitation_rate,
                                                                       NumType atm_temperature,
                                                                       int max_n_its = MAX_RELAXATION_ITS,
                                                                       NumType tol = RELAXATION_SOLVER_ACCUR);

        // !! sea-ice freezing mode with snow (for 0d-ice and 0d-snow profile) !!
        /* 
            output - pair of 0d-ice and 0d-snow values
            
            0d-ice values:
                new temperature of ice-snow interface
                new ice thicknesses
            
            0d-snow values:
                new snow surface temperature
                new snow thickness
        */
        std::pair<std::pair<NumType, NumType>, std::pair<NumType, NumType>> seaice0d_snow0d_freezing(NumType T_ib,
                                                                                                     NumType T_is,
                                                                                                     NumType thickness_i,
                                                                                                     NumType salinity_i,
                                                                                                     NumType rho_i,
                                                                                                     NumType T_ss,
                                                                                                     NumType thickness_s,
                                                                                                     NumType rho_s,
                                                                                                     NumType precipitation_rate,
                                                                                                     NumType atm_temperature,
                                                                                                     int max_n_its = MAX_RELAXATION_ITS,
                                                                                                     NumType tol = 1e-6);
        
        
        // !! snow melting mode with sea-ice (for 0d-ice and 0d-snow profile) !!
        /* 
            output - pair of 0d-ice and 0d-snow values
            
            0d-ice values:
                new temperature of ice-snow interface
                new ice thickness
            
            0d-snow values:
                new snow thickness
        */
        std::pair<std::pair<NumType, NumType>, NumType> seaice0d_snow0d_melting(NumType T_ib,
                                                                                NumType T_is,
                                                                                NumType thickness_i,
                                                                                NumType salinity_i,
                                                                                NumType rho_i,
                                                                                NumType T_ss,
                                                                                NumType thickness_s,
                                                                                NumType rho_s,
                                                                                NumType precipitation_rate,
                                                                                NumType atm_temperature,
                                                                                int max_n_its = MAX_RELAXATION_ITS,
                                                                                NumType tol = 1e-6);

        // !! seaice freezing mode (for 0d profile) !!
        /* 
            output - 0d-ice values
            
            0d-ice values:
                new ice surface temperature
                new ice thickness
        */

        std::pair<NumType, NumType> seaice0d_freezing(NumType T_ib,
                                                      NumType T_is,
                                                      NumType thickness_i,
                                                      NumType salinity_i,
                                                      NumType rho_i,
                                                      int max_n_its = MAX_RELAXATION_ITS,
                                                      NumType tol = 1e-6);

        // !! seaice melting mode (for 0d profile) !!
        /* 
            output - 0d-ice values
            
            0d-ice values:
                new ice thickness
        */
        NumType seaice0d_melting(NumType T_ib,
                                 NumType T_is,
                                 NumType thickness_i,
                                 NumType salinity_i,
                                 NumType rho_i,
                                 int max_n_its = MAX_RELAXATION_ITS,
                                 NumType tol = 1e-6);
    };



    //   #############################################
    //   ############# 1D SEA ICE  CLASS #############
    //   #############################################

    template<typename NumType>
    class SeaIce1D_Solver : public SeaIce_Solver<NumType>
    {
    public:
        // constructor
        SeaIce1D_Solver(Mesh<NumType>* mesh_ice_,
                        NumType time_step_,
                        NumType min_ice_thick_,
                        bool is_radiation = false,
                        bool is_sublimation = false,
                        bool is_verbose_ = false,
                        ApproxOrder grad_approx_order_ = ApproxOrder::first,
                        Kparam ice_k_param_ = Kparam::FreshIce,
                        Cparam ice_c_eff_param_ = Cparam::FreshIce,
                        Eparam ice_E_param_ = Eparam::FreshIce,
                        Lparam ice_L_param_ = Lparam::FreshIce,
                        Aparam ice_albedo_param_ = Aparam::ConstantIce);

        // one-step evaluation
        void Evaluate() override;

    private:
        void CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                  Mesh<NumType>* snow_mesh = NULL) override;
    
    public:
        void UpdateMesh(Mesh<NumType>* mesh_ice_,
                        Mesh<NumType>* mesh_snow_ = NULL) override;

    };


    //   ######################################################
    //   ############# 1D SEA ICE + 0D SNOW CLASS #############
    //   ######################################################

    template<typename NumType>
    class SeaIce1D_Snow0D_Solver : public SeaIce_Solver<NumType>
    {
    public:
        SeaIce1D_Snow0D_Solver(Mesh<NumType>* mesh_ice_,
                               Mesh<NumType>* mesh_snow_,
                               NumType time_step_,
                               NumType min_ice_thick_,
                               NumType min_snow_thick_,
                               bool is_radiation_ = false,
                               bool is_sublimation_ = false,
                               bool is_verbose_ = false,
                               Kparam ice_k_param_ = Kparam::FreshIce,
                               Cparam ice_c_eff_param_ = Cparam::FreshIce,
                               Eparam ice_E_param_ = Eparam::FreshIce,
                               Lparam ice_L_param_ = Lparam::FreshIce,
                               Aparam ice_albedo_param_ = Aparam::ConstantIce,
                               Kparam snow_k_param_ = Kparam::FreshSnow,
                               Lparam snow_L_param_ = Lparam::FreshSnow,
                               Aparam snow_albedo_param_ = Aparam::ConstantSnow,
                               SnowIceTransition si_transition_mode_ = SnowIceTransition::None);
        
        // one-step evaluation
        void Evaluate() override;

    private:
        void CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                  Mesh<NumType>* snow_mesh) override;
    
    public:
        void UpdateMesh(Mesh<NumType>* mesh_ice_,
                        Mesh<NumType>* mesh_snow_) override;

    };

    //   ######################################################
    //   ############# 0D SEA ICE + 0D SNOW CLASS #############
    //   ######################################################

    template<typename NumType>
    class SeaIce0D_Snow0D_Solver : public SeaIce_Solver<NumType>
    {
    public:
        SeaIce0D_Snow0D_Solver(Mesh<NumType>* mesh_ice_,
                               Mesh<NumType>* mesh_snow_,
                               NumType time_step_,
                               NumType min_ice_thick_,
                               NumType min_snow_thick_,
                               bool is_sublimation_ = false,
                               bool is_verbose_ = false,
                               Kparam ice_k_param_ = Kparam::FreshIce,
                               Lparam ice_L_param_ = Lparam::FreshIce,
                               Aparam ice_albedo_param_ = Aparam::ConstantIce,
                               Kparam snow_k_param_ = Kparam::FreshSnow,
                               Lparam snow_L_param_ = Lparam::FreshSnow,
                               Aparam snow_albedo_param_ = Aparam::ConstantSnow,
                               SnowIceTransition si_transition_mode_ = SnowIceTransition::None);
        
        // one-step evaluation
        void Evaluate() override;

    private:
        void CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                  Mesh<NumType>* snow_mesh) override;
    
    public:
        void UpdateMesh(Mesh<NumType>* mesh_ice_,
                        Mesh<NumType>* mesh_snow_) override;

    };

    //   ######################################################
    //   ################## 1D GLACIER CLASS ##################
    //   ######################################################
    template<typename NumType>
    class Glacier1D_Solver : public ThermoSolver<NumType>
    {
    public:
        // constructors
        Glacier1D_Solver();

        Glacier1D_Solver(Mesh<NumType>* mesh_ice_,
                         NumType time_step_,
                         NumType min_ice_thick_,
                         bool is_radiation_ = false,
                         bool is_sublimation_ = false,
                         bool is_verbose_ = false,
                         ApproxOrder grad_approx_order_ = ApproxOrder::first,
                         Kparam ice_k_param_ = Kparam::FreshIce,
                         Cparam ice_c_eff_param_ = Cparam::FreshIce,
                         Eparam ice_E_param_ = Eparam::FreshIce,
                         Lparam ice_L_param_ = Lparam::FreshIce,
                         Aparam ice_albedo_param_ = Aparam::ConstantIce);

        // one evaluation step
        void Evaluate() override;

    private:
        void CheckMeshConsistency(Mesh<NumType>* ice_mesh,
                                  Mesh<NumType>* snow_mesh = NULL) override;
    
    public:
        void UpdateMesh(Mesh<NumType>* mesh_ice_,
                        Mesh<NumType>* mesh_snow_ = NULL) override;
                        
    };
 
}