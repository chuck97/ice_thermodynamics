#include "solver.hpp"

namespace icethermo
{
    template <typename NumType>
    SeaIce_Solver<NumType>::SeaIce_Solver(Mesh<NumType>* mesh_ice_,
                                          Mesh<NumType>* mesh_snow_,
                                          NumType time_step_,
                                          NumType min_ice_thick_,
                                          NumType min_snow_thick_,
                                          ApproxOrder grad_approx_order_,
                                          bool is_radiation_,
                                          bool is_sublimation_,
                                          bool is_verbose_,
                                          Kparam ice_k_param_,
                                          Cparam ice_c_eff_param_,
                                          Eparam ice_E_param_,
                                          Lparam ice_L_param_,
                                          Aparam ice_albedo_param_,
                                          Kparam snow_k_param_,
                                          Cparam snow_c_eff_param_,
                                          Eparam snow_E_param_,
                                          Lparam snow_L_param_,
                                          Aparam snow_albedo_param_,
                                          SnowIceTransition si_transition_mode_):
        ThermoSolver<NumType>(mesh_ice_,
                              mesh_snow_,
                              time_step_,
                              min_ice_thick_,
                              min_snow_thick_,
                              grad_approx_order_,
                              is_radiation_,
                              is_sublimation_,
                              is_verbose_,
                              ice_k_param_,
                              ice_c_eff_param_,
                              ice_E_param_,
                              ice_L_param_,
                              ice_albedo_param_,
                              snow_k_param_,
                              snow_c_eff_param_,
                              snow_E_param_,
                              snow_L_param_,
                              snow_albedo_param_,
                              si_transition_mode_)
    {
        this->So = std::make_shared<NumType>(30.0); // psu
        this->omega_ib_expl_mass = std::make_shared<NumType>(0.0); // m s-1
    }

    // Update ocean salinity
    template<typename NumType>
    void SeaIce_Solver<NumType>::UpdateOceanSalinity(NumType ocn_sal_)
    {
        *(this->So) = ocn_sal_;
    }

    // Update ocean salinity
    template<typename NumType>
    void SeaIce_Solver<NumType>::UpdateOceanIceMassFlux(NumType mass_flux_)
    {
        *(this->omega_ib_expl_mass) = mass_flux_;
    }

    // 1d sea-ice freezing mode
    template<typename NumType>
    ThreeVecs<NumType> SeaIce_Solver<NumType>::seaice1d_freezing(NumType T_ib,
                                                                 const std::vector<NumType>& T_cells,
                                                                 NumType T_is,
                                                                 const std::vector<NumType>& dz_cells,
                                                                 const std::vector<NumType>& salinity_cells,
                                                                 const std::vector<NumType>& rho_cells,
                                                                 int max_n_its,
                                                                 NumType tol)
    {
        std::vector<NumType> T_cells_new = T_cells;
        std::vector<NumType> T_cells_prev = T_cells;

        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType omega_ib = (NumType)0.0;

        std::vector<NumType> dz_cells_new = dz_cells;

        std::vector<NumType> radiation_nodes(T_cells.size() + 1);

        std::vector<NumType> T_is_history = {T_is};
        NumType current_err = std::numeric_limits<NumType>::max();
        NumType surface_err = std::numeric_limits<NumType>::max();
        NumType prev_surface_err = surface_err;

        std::vector<NumType> full_temp_vec(T_cells.size() + 1);
        std::vector<NumType> prev_temp_vec(T_cells.size() + 1);
        std::vector<NumType> old_temp_vec = concatenate<NumType>(T_cells, std::vector<NumType>{T_is});

        NumType omega_is = (NumType)0.0; 

        NumType kappa_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.kappa_i : IceConsts<NumType>::kappa_i;
        NumType i0_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.i0_i : IceConsts<NumType>::i0_i;

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;
            T_cells_prev = T_cells_new;
            T_is_prev = T_is_new;
            prev_surface_err = surface_err;

            // compute new value of omega at the base
            omega_ib = this->W_from_BC(T_ib,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       false);

            // add sublimation
            NumType Ls = (Configured()) ? GetConfigConsts<NumType>()->GenConsts.L_sublim : GenConsts<NumType>::L_sublim;
            omega_is = (this->is_sublimation) ?
                -this->F_lh(T_is_new)/(rho_cells.back()*Ls) :
                (NumType)0.0;
            
            // compute new value of ice surface temperature
            T_is_new = this->T_from_BC(T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       0.0,
                                       true,
                                       true);

            // force the convergence of surface temperature
            if ((Configured()) ? GetConfigConsts<NumType>()->SolverRelaxation.force_surf_conv : true)
            {
                surface_err = std::abs(T_is_new - T_is_prev)/(std::abs(T_is) + (NumType)0.1);

                if (surface_err < prev_surface_err)
                {
                    T_is_history.push_back(T_is_new);
                }
                else
                {
                    T_is_new = sum_vec<NumType>(T_is_history)/T_is_history.size();
                }
            }

            // recalculate ice thickness
            dz_cells_new = this->Update_dz(dz_cells,
                                           omega_ib + (*(this->omega_ib_expl_mass)),
                                           omega_is);
            
            // recalculate radiation
            if (this->is_radiation)
            {
                NumType h_i = sum_vec<NumType>(dz_cells_new);
                NumType Tfi = Params<NumType>::TempFusion((*(this->Si_cells)).back());
                auto albedo_param = this->ice_albedo_param;
                FuncPtr<NumType> alb_i = [albedo_param, h_i, Tfi] (NumType T)
                {
                    return Params<NumType>::Albedo(albedo_param, T, h_i, Tfi);
                };

                radiation_nodes = this->Compute_radiation_nodes(this->F_sw(T_is_new),
                                                                alb_i(T_is_new),
                                                                i0_i,
                                                                kappa_i,
                                                                dz_cells_new).first;
            }

            // assemble matrix and rhs for ice 
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_cells_prev, T_cells,
                                                                T_is_new, T_is_prev, T_is,
                                                                T_ib, T_ib, T_ib, 
                                                                omega_ib, omega_is, 
                                                                dz_cells_new, dz_cells,
                                                                salinity_cells,
                                                                rho_cells,
                                                                radiation_nodes,
                                                                true);

            // solve linear system and update temperatures
            T_cells_new = thomas_solver<NumType>(std::get<0>(matrix_rhs),
                                                 std::get<1>(matrix_rhs),
                                                 std::get<2>(matrix_rhs),
                                                 std::get<3>(matrix_rhs));
            
            // evaluate the error
            prev_temp_vec = concatenate<NumType>(T_cells_prev, std::vector<NumType>{T_is_prev});
            
            full_temp_vec = concatenate<NumType>(T_cells_new, std::vector<NumType>{T_is_new}); 

            current_err = L2_norm(full_temp_vec - prev_temp_vec)/(L2_norm(old_temp_vec) + (NumType)0.1);

            if (current_err < tol)
                break;
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }
        
        return {T_cells_new, std::vector<NumType>{T_is_new}, dz_cells_new};
    }

    // 1d sea-ice freezing mode
    template<typename NumType>
    std::pair<NumType, NumType> SeaIce_Solver<NumType>::seaice0d_freezing(NumType T_ib,
                                                                          NumType T_is,
                                                                          NumType thickness_i,
                                                                          NumType salinity_i,
                                                                          NumType rho_i,
                                                                          int max_n_its,
                                                                          NumType tol)
    {
        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType thickness_old = thickness_i;
        NumType thickness_new = thickness_i;
        NumType thickness_prev = thickness_i;

        NumType omega_ib = (NumType)0.0;
        NumType omega_is = (NumType)0.0;

        std::vector<NumType> T_is_history = {T_is};

        NumType current_err = std::numeric_limits<NumType>::max();
        NumType surface_err = std::numeric_limits<NumType>::max();

        NumType prev_surface_err = surface_err; 
        

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_is_prev = T_is_new;
            thickness_prev = thickness_new;
            prev_surface_err = surface_err;

            // recalculate conductivity of ice
            NumType k_i = Params<NumType>::Conductivity(this->ice_k_param, salinity_i, 0.5*(T_ib + T_is_new), rho_i);

            // compute new value of omega at the base
            omega_ib = this->W_from_BC_0D(T_ib,
                                          T_is_prev,
                                          thickness_prev,
                                          k_i,
                                          salinity_i,
                                          rho_i,
                                          true,
                                          false);

            // add sublimation
            omega_is = (NumType)0.0;
            
            // compute new value of ice surface temperature
            T_is_new = this->T_from_BC_0D([T_ib](NumType T){return T_ib;},
                                          thickness_prev,
                                          k_i,
                                          salinity_i,
                                          omega_is,
                                          rho_i,
                                          true,
                                          true);

            // force the convergence of surface temperature
            surface_err = std::abs(T_is_new - T_is_prev)/(std::abs(T_is) + (NumType)0.1);
            
            //if (surface_err <= prev_surface_err)
            //{
            //    T_is_history.push_back(T_is_new);
            //}
            //else
            //{
            //    T_is_new = sum_vec<NumType>(T_is_history)/T_is_history.size();
            //}

            // recalculate ice thickness
            thickness_new = this->Update_dz_0D(thickness_old,
                                               omega_ib,
                                               omega_is);

            //if (surface_err < tol)
            //{
            //    break;
            //}
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << surface_err << std::endl;
        }
        
        return {T_is_new, thickness_new};
    }

    // 1d sea-ice melting mode
    template<typename NumType>
    TwoVecs<NumType> SeaIce_Solver<NumType>::seaice1d_melting(NumType T_ib,
                                                              NumType T_is,
                                                              const std::vector<NumType>& T_cells,
                                                              NumType T_is_old,
                                                              const std::vector<NumType>& dz_cells,
                                                              const std::vector<NumType>& salinity_cells,
                                                              const std::vector<NumType>& rho_cells,
                                                              int max_n_its,
                                                              NumType tol)
    {
        std::vector<NumType> T_cells_new = T_cells;
        std::vector<NumType> T_cells_prev = T_cells;

        NumType omega_ib = 0.0;
        NumType omega_is = 0.0;

        std::vector<NumType> dz_cells_new = dz_cells;

        std::vector<NumType> radiation_nodes(T_cells.size() + 1);

        NumType current_err = std::numeric_limits<NumType>::max();

        NumType omega_is_sublim;

        NumType kappa_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.kappa_i : IceConsts<NumType>::kappa_i;
        NumType i0_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.i0_i : IceConsts<NumType>::i0_i;

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_cells_prev = T_cells_new;
            
            // compute new value of omega at the base
            omega_ib = this->W_from_BC(T_ib,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       false);

            // compute new value of omega at the base
            omega_is = this->W_from_BC(T_is,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       true);
            
            // add sublimation
            NumType Ls = (Configured()) ? GetConfigConsts<NumType>()->GenConsts.L_sublim : GenConsts<NumType>::L_sublim;
            omega_is_sublim = (this->is_sublimation) ?
                -this->F_lh(T_is)/(rho_cells.back()*Ls) :
                (NumType)0.0;

            // recalculate ice thickness
            dz_cells_new = this->Update_dz(dz_cells,
                                           omega_ib + (*(this->omega_ib_expl_mass)),
                                           omega_is + omega_is_sublim);
            
            // recalculate radiation
            if (this->is_radiation)
            {
                NumType h_i = sum_vec<NumType>(dz_cells_new);
                NumType Tfi = Params<NumType>::TempFusion((*(this->Si_cells)).back());
                auto albedo_param = this->ice_albedo_param;
                FuncPtr<NumType> alb_i = [albedo_param, h_i, Tfi] (NumType T)
                {
                    return Params<NumType>::Albedo(albedo_param, T, h_i, Tfi);
                };

                radiation_nodes = this->Compute_radiation_nodes(this->F_sw(T_is),
                                                                alb_i(T_is),
                                                                i0_i,
                                                                kappa_i,
                                                                dz_cells_new).first;
            }

            // assemble matrix and rhs for ice 
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_cells_prev, T_cells,
                                                                T_is, T_is, T_is_old,
                                                                T_ib, T_ib, T_ib, 
                                                                omega_ib, omega_is + omega_is_sublim, 
                                                                dz_cells_new, dz_cells,
                                                                salinity_cells,
                                                                rho_cells,
                                                                radiation_nodes,
                                                                true);
            
            // solve linear system and update temperatures
            T_cells_new = thomas_solver<NumType>(std::get<0>(matrix_rhs),
                                                 std::get<1>(matrix_rhs),
                                                 std::get<2>(matrix_rhs),
                                                 std::get<3>(matrix_rhs));

            // evaluate the error
            current_err = L2_norm(T_cells_new - T_cells_prev)/(L2_norm(T_cells) + (NumType)0.1);

            //if (current_err < tol)
            //    break;
        }
        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }
        return {T_cells_new, dz_cells_new};
    }

    // 0d sea-ice melting mode
    template<typename NumType>
    NumType SeaIce_Solver<NumType>::seaice0d_melting(NumType T_ib,
                                                     NumType T_is,
                                                     NumType thickness_i,
                                                     NumType salinity_i,
                                                     NumType rho_i,
                                                     int max_n_its,
                                                     NumType tol)
    {

        NumType omega_ib = 0.0;
        NumType omega_is = 0.0;

        NumType thickness_new = thickness_i;
        NumType thickness_prev = thickness_i;
        NumType thickness_old = thickness_i;

        NumType current_err = std::numeric_limits<NumType>::max();

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {   
            npseudo++;

            thickness_prev = thickness_new;

            // recalculate conductivity of ice
            NumType k_i = Params<NumType>::Conductivity(this->ice_k_param, salinity_i, 0.5*(T_ib + T_is), rho_i);

            // compute new value of omega at the base
            omega_ib = this->W_from_BC_0D(T_ib,
                                          T_is,
                                          thickness_prev,
                                          k_i,
                                          salinity_i,
                                          rho_i,
                                          true,
                                          false);

            // compute new value of omega at the base
            omega_is = this->W_from_BC_0D(T_is,
                                          T_ib,
                                          thickness_prev,
                                          k_i,
                                          salinity_i,
                                          rho_i,
                                          true,
                                          true);

            // recalculate ice thickness
            thickness_new = this->Update_dz_0D(thickness_old,
                                               omega_ib,
                                               omega_is);
            

            // evaluate the error
            current_err = std::abs(thickness_new - thickness_prev)/std::abs(thickness_old + 0.1);

            //if (current_err < tol)
            //    break;
        }
        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }
        return thickness_new;
    }

    template<typename NumType>
    std::pair<std::pair<NumType, NumType>, std::pair<NumType, NumType>> SeaIce_Solver<NumType>::seaice0d_snow0d_freezing
            (NumType T_ib,
             NumType T_is,
             NumType thickness_i,
             NumType salinity_i,
             NumType rho_i,
             NumType T_ss,
             NumType thickness_s,
             NumType rho_s,
             NumType precipitation_rate,
             NumType atm_temperature,
             int max_n_its,
             NumType tol)
    {
        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType T_ss_new = T_ss;
        NumType T_ss_prev = T_ss;

        NumType thickness_s_old = thickness_s;
        NumType thickness_s_new = thickness_s;
        NumType thickness_s_prev = thickness_s;

        NumType thickness_i_old = thickness_i;
        NumType thickness_i_new = thickness_i;
        NumType thickness_i_prev = thickness_i;

        NumType omega_ib = (NumType)0.0;

        std::vector<NumType> T_ss_history = {T_ss};

        NumType surface_err = std::numeric_limits<NumType>::max();
        NumType current_err = std::numeric_limits<NumType>::max();
        NumType prev_surface_err = surface_err; 

        // setup interface temperature as a function of snow surface temperature


        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_is_prev = T_is_new;
            T_ss_prev = T_ss_new;
            thickness_s_prev = thickness_s_new;
            thickness_i_prev = thickness_i_new;
            prev_surface_err = surface_err;

            // recalculate conductivity of ice and snow
            NumType k_i = Params<NumType>::Conductivity(this->ice_k_param, salinity_i, 0.5*(T_ib + T_is_prev), rho_i);
            NumType k_s = Params<NumType>::Conductivity(this->snow_k_param, 0.0, 0.5*(T_is_prev + T_ss_prev), rho_s);

            // recalculate interface temperature as function of surface temperature
            FuncPtr<NumType> int_temp = [thickness_s_prev, thickness_i_prev, T_ib, k_i, k_s](NumType T)
                                        {
                                            NumType c1 = k_s/thickness_s_prev;
                                            NumType c2 = k_i/thickness_i_prev;
                                            return (NumType)((c1*T + c2*T_ib)/(c1 + c2));
                                        };

            // compute new value of snow surface temperature
            T_ss_new = this->T_from_BC_0D(int_temp,
                                          thickness_s_prev,
                                          k_s,
                                          0.0,
                                          0.0,
                                          rho_s,
                                          false,
                                          true);

            // recalculate interface temperature
            T_is_new  = int_temp(T_ss_new);

            // compute new value of omega at the base
            omega_ib = this->W_from_BC_0D(T_ib,
                                          T_is_new,
                                          thickness_i_prev,
                                          k_i,
                                          salinity_i,
                                          rho_i,
                                          true,
                                          false);
            
            


            // force the convergence of snow surface temperature
            surface_err = std::abs(T_ss_new - T_ss_prev);
            
            //std::abs(T_ss_new - T_ss_prev);///(std::abs(T_ss) + (NumType)0.1);
            
            //if (surface_err <= prev_surface_err)
            //{
            //    T_ss_history.push_back(T_ss_new);
            //}
            //else
            //{
            //    T_ss_new = sum_vec<NumType>(T_ss_history)/T_ss_history.size();
            //}

            // recalculate ice thickness
            thickness_i_new = this->Update_dz_0D(thickness_i_old,
                                                 omega_ib,
                                                 0.0);

            // recalculate snow thickness
            thickness_s_new = this->Update_dz_0D(thickness_s_old,
                                                 0.0,
                                                 0.0);


            // evaluate the error
            std::vector<NumType> prev_temps = {T_is_prev, T_ss_prev};
            std::vector<NumType> new_temps = {T_is_new, T_ss_new};

            current_err = L2_norm(new_temps - prev_temps)/(L2_norm(std::vector<NumType>{T_is, T_ss}) + (NumType)0.1);

            //FuncPtr<NumType> F_lws = this->F_lw;
            //FuncPtr<NumType> F_lwis = this->F_lwi;
            //FuncPtr<NumType> F_sws = this->F_sw;
            //FuncPtr<NumType> F_shs = this->F_sh;
            //FuncPtr<NumType> F_lhs = this->F_lh;
            //FuncPtr<NumType> F_Ps = this->F_P;
            //FuncPtr<NumType> F_total = this->F_up;

            //NumType T_ss_neww = -24.5;

            //std::cout << "snow temp:" << T_ss_neww << std::endl; 
            //std::cout << "T_ib:" << T_ib << std::endl;
            //std::cout << "C1:" << k_s/thickness_s_new << std::endl;
            //std::cout << "C2:" << k_i/thickness_i_new << std::endl;
            //std::cout << "T_is_new:" << T_is_new << std::endl;
            //std::cout << "F_cond:" <<-SnowConsts<NumType>::k0_s*(T_ss_neww - int_temp(T_ss_neww))/thickness_s_new << std::endl;
            //std::cout << "F_lws:" << GenConsts<NumType>::emissivity*F_lws(T_ss_neww) << std::endl;
            //std::cout << "F_lwis:" << -GenConsts<NumType>::emissivity*F_lwis(T_ss_neww) << std::endl;
            //std::cout << "F_sws:" << (1- SnowConsts<NumType>::albedo_dry_s)*F_sws(T_ss_neww) << std::endl;
            //std::cout << "F_shs:" << F_shs(T_ss_neww) << std::endl;
            //std::cout << "F_lhs:" << F_lhs(T_ss_neww) << std::endl;
            //std::cout << "F_Ps:" << F_Ps(T_ss_neww) << std::endl;
            //std::cout << "F_total:" << F_total(T_ss_neww) << std::endl;
            //std::cout << std::endl;
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }
        
        return {{T_is_new, thickness_i_new}, {T_ss_new, thickness_s_new}};
    }

    template<typename NumType>
    std::pair<std::pair<NumType, NumType>, NumType> SeaIce_Solver<NumType>::seaice0d_snow0d_melting
    (NumType T_ib,
     NumType T_is,
     NumType thickness_i,
     NumType salinity_i,
     NumType rho_i,
     NumType T_ss,
     NumType thickness_s,
     NumType rho_s,
     NumType precipitation_rate,
     NumType atm_temperature,
     int max_n_its,
     NumType tol)
    {
        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType thickness_s_old = thickness_s;
        NumType thickness_s_new = thickness_s;
        NumType thickness_s_prev = thickness_s;

        NumType thickness_i_old = thickness_i;
        NumType thickness_i_new = thickness_i;
        NumType thickness_i_prev = thickness_i;

        NumType omega_ib = (NumType)0.0;
        NumType omega_ss = (NumType)0.0;

        NumType current_err = std::numeric_limits<NumType>::max();

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_is_prev = T_is_new;
            thickness_s_prev = thickness_s_new;
            thickness_i_prev = thickness_i_new;
            
            // recalculate conductivity of snow
            NumType k_s = Params<NumType>::Conductivity(this->snow_k_param, 0.0, 0.5*(T_is_prev + T_ss), rho_s);

            // recalculate conductivity of ice
            NumType k_i = Params<NumType>::Conductivity(this->ice_k_param, salinity_i, 0.5*(T_ib + T_is_prev), rho_i);

            // recalculate function for interface temperature
            FuncPtr<NumType> int_temp = [thickness_s_prev, thickness_i_prev, T_ib, k_i, k_s](NumType T)
                                        {
                                            NumType c1 = k_s/thickness_s_prev;
                                            NumType c2 = k_i/thickness_i_prev;
                                            return (NumType)((c1*T + c2*T_ib)/(c1 + c2));
                                        };
            
            T_is_new = int_temp(T_ss);

            // compute new value of omega at the ice base
            omega_ib = this->W_from_BC_0D(T_ib,
                                          T_is_new,
                                          thickness_i_prev,
                                          k_i,
                                          salinity_i,
                                          rho_i,
                                          true,
                                          false);

            

            // compute new value of omega at snow surface
            omega_ss = this->W_from_BC_0D(T_ss,
                                          T_is_new,
                                          thickness_s_prev,
                                          k_s,
                                          0.0,
                                          rho_s,
                                          false,
                                          true);

            // recalculate ice thickness
            thickness_i_new = this->Update_dz_0D(thickness_i_old,
                                                 omega_ib,
                                                 0.0);

            // recalculate snow thickness
            thickness_s_new = this->Update_dz_0D(thickness_s_old,
                                                 0.0,
                                                 omega_ss);

            // evaluate the error
            std::vector<NumType> prev_temps = {T_is_prev};
            std::vector<NumType> new_temps = {T_is_new};

            current_err = L2_norm(new_temps - prev_temps)/(L2_norm(std::vector<NumType>{T_is}) + (NumType)0.1);
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }
        
        return {{T_is_new, thickness_i_new}, thickness_s_new};
    }

    // 1d sea-ice with 0d snow freezing mode
    template<typename NumType>
    std::pair<ThreeVecs<NumType>, TwoVecs<NumType>> SeaIce_Solver<NumType>::seaice1d_snow0d_freezing(NumType T_ib,
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
                                                                                                     int max_n_its,
                                                                                                     NumType tol)
    {
        std::vector<NumType> T_i_cells_new = T_i_cells;
        std::vector<NumType> T_i_cells_prev = T_i_cells;

        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType T_ss_new = T_ss;
        NumType T_ss_prev = T_ss;

        NumType omega_ib;

        std::vector<NumType> dz_i_cells_new = dz_i_cells;
        NumType h_s_new = thickness_s;

        std::vector<NumType> T_ss_history = {T_is};
        NumType current_err = std::numeric_limits<NumType>::max();
        NumType surface_err = std::numeric_limits<NumType>::max();
        NumType prev_surface_err = surface_err;

        std::vector<NumType> full_temp_vec(T_i_cells.size() + 2);
        std::vector<NumType> prev_temp_vec(T_i_cells.size() + 2);
        std::vector<NumType> old_temp_vec = concatenate<NumType>({T_i_cells, std::vector<NumType>{T_is}, std::vector<NumType>{T_ss}});

        NumType r_w = (Configured()) ? GetConfigConsts<NumType>()->WaterConsts.rho_w : WaterConsts<NumType>::rho_w;
        NumType r_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.rho_s :SnowConsts<NumType>::rho_s;
        NumType omega_ss =  (atm_temperature < (NumType)0.0) ? -precipitation_rate*r_w/r_s : (NumType)0.0;

        std::vector<NumType> radiation_nodes_ice(T_i_cells.size() + 1);

        NumType omega_ss_sublim;
        NumType omega_interface = (NumType)0.0;

        NumType kappa_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.kappa_i : IceConsts<NumType>::kappa_i;
        NumType kappa_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.kappa_s :SnowConsts<NumType>::kappa_s;
        NumType i0_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.i0_i : IceConsts<NumType>::i0_i;
        NumType i0_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.i0_s : SnowConsts<NumType>::i0_s;
        

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_i_cells_prev = T_i_cells_new;
            T_is_prev = T_is_new;
            T_ss_prev = T_ss_new;
            prev_surface_err = surface_err;

            // compute new value of omega at the base
            omega_ib = this->W_from_BC(T_ib,
                                       T_i_cells_prev,
                                       dz_i_cells_new,
                                       salinity_i_cells,
                                       rho_i_cells,
                                       true,
                                       false);
            
            // recalculate conductivity of snow
            NumType k_s = Params<NumType>::Conductivity(this->snow_k_param, (NumType)0.5*(T_ss_prev + T_is_prev), (NumType)0.0, rho_s);

            // add sublimation 
            NumType Ls = (Configured()) ? GetConfigConsts<NumType>()->GenConsts.L_sublim : GenConsts<NumType>::L_sublim;
            omega_ss_sublim = (this->is_sublimation) ?
                -this->F_lh(T_ss_new)/(rho_s*Ls) :
                (NumType)0.0;

            NumType k_i = Params<NumType>::Conductivity(this->ice_k_param, salinity_i_cells.back(), T_is_new, rho_i_cells.back());

            // recalculate interface temperature as function of surface temperature
            FuncPtr<NumType> int_temp = [h_s_new, dz_i_cells_new, T_i_cells_new, k_i, k_s](NumType T)
                                        {
                                            NumType c1 = k_s/h_s_new;
                                            NumType c2 = k_i/(0.5*dz_i_cells_new.back());
                                            return (NumType)((c1*T + c2*T_i_cells_new.back())/(c1 + c2));
                                        };

            // compute new value of snow surface temperature
            T_ss_new = this->T_from_BC_0D([int_temp](NumType T){return int_temp(T);},
                                          h_s_new,
                                          k_s,
                                          (NumType)0.0,
                                          0.0,
                                          rho_s,
                                          false,
                                          true);

            // force the convergence of surface temperature
            if ((Configured()) ? GetConfigConsts<NumType>()->SolverRelaxation.force_surf_conv : true)
            {
                surface_err = std::abs(T_ss_new - T_ss_prev)/std::abs(T_ss + (NumType)0.1);

                if (surface_err < prev_surface_err)
                {
                    T_ss_history.push_back(T_ss_new);
                }
                else
                {
                    T_ss_new = sum_vec<NumType>(T_ss_history)/T_ss_history.size();
                }
            }

            // compute ice-snow interface omega value
            if (this->si_transition_mode == SnowIceTransition::SnowAging)
            {
                omega_interface = ((NumType)1.0/GenConsts<NumType>::si_trans_time_scale)*h_s_new;
            }

            // recalculate ice thickness
            dz_i_cells_new = this->Update_dz(dz_i_cells,
                                             omega_ib + (*(this->omega_ib_expl_mass)),
                                             omega_interface);

            // recalculate snow thickness
            h_s_new = this->Update_dz_0D(thickness_s,
                                         omega_interface,
                                         omega_ss + omega_ss_sublim);

            // recalculate ice radiation
            if (this->is_radiation)
            {

                NumType h_i = sum_vec<NumType>(dz_i_cells_new);
                NumType Tfi = Params<NumType>::TempFusion((*(this->Si_cells)).back());
                auto ice_albedo_param = this->ice_albedo_param;
                FuncPtr<NumType> alb_i = [ice_albedo_param, h_i, Tfi] (NumType T)
                {
                    return Params<NumType>::Albedo(ice_albedo_param, T, h_i, Tfi);
                };

                
                auto snow_albedo_param = this->snow_albedo_param;
                FuncPtr<NumType> alb_s = [snow_albedo_param, h_s_new, h_i, Tfi] (NumType T)
                {
                    return Params<NumType>::Albedo(snow_albedo_param, T, h_s_new, 0.0);
                };

                radiation_nodes_ice = this->Compute_radiation_nodes(this->F_sw(T_ss_new),
                                                                    alb_i(T_is_new),
                                                                    i0_i,
                                                                    kappa_i,
                                                                    dz_i_cells_new,
                                                                    alb_s(T_ss_new),
                                                                    i0_s,
                                                                    kappa_s,
                                                                    std::vector<NumType>{h_s_new}).first;
            }
            
            // recalculate ice temperature profile
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_i_cells_prev, T_i_cells,
                                                                T_is_new, T_is_prev, T_is,
                                                                T_ib, T_ib, T_ib, 
                                                                omega_ib, (NumType)0.0, 
                                                                dz_i_cells_new, dz_i_cells,
                                                                salinity_i_cells,
                                                                rho_i_cells,
                                                                radiation_nodes_ice,
                                                                true);

            // solve linear system and update temperatures
            T_i_cells_new = thomas_solver<NumType>(std::get<0>(matrix_rhs),
                                                   std::get<1>(matrix_rhs),
                                                   std::get<2>(matrix_rhs),
                                                   std::get<3>(matrix_rhs));


            // recalculate ice-snow interface temperature 
            NumType top_ice_k = Params<NumType>::Conductivity(this->ice_k_param, T_i_cells_new.back(), salinity_i_cells.back(), rho_i_cells.back());
            NumType ratio = (top_ice_k*h_s_new)/(k_s*dz_i_cells_new.back()*(NumType)0.5);
            T_is_new  = ((NumType)1.0/(ratio + (NumType)1.0))*T_ss_new + ((ratio)/(ratio + (NumType)1.0))*T_i_cells_new.back();

            
            // evaluate the error
            prev_temp_vec = concatenate<NumType>({T_i_cells_prev, std::vector<NumType>{T_is_prev}, std::vector<NumType>{T_ss_prev}});
            
            full_temp_vec = concatenate<NumType>({T_i_cells_new, std::vector<NumType>{T_is_new}, std::vector<NumType>{T_ss_new}});

            current_err = L2_norm(full_temp_vec - prev_temp_vec)/(L2_norm(old_temp_vec) + (NumType)0.1);

            if (current_err < tol)
                break;
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }

        return {{T_i_cells_new, std::vector<NumType>{T_is_new}, dz_i_cells_new}, {std::vector<NumType>{T_ss_new}, std::vector<NumType>{h_s_new}}};
    }

    // 1d sea-ice with 0d snow melting mode
    template<typename NumType>
    std::pair<ThreeVecs<NumType>, NumType> SeaIce_Solver<NumType>::seaice1d_snow0d_melting(NumType T_ib,
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
                                                                                           int max_n_its,
                                                                                           NumType tol)
    {
        std::vector<NumType> T_i_cells_new = T_i_cells;
        std::vector<NumType> T_i_cells_prev = T_i_cells;

        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType omega_ib;
        NumType omega_ss = (NumType)0.0;

        std::vector<NumType> dz_i_cells_new = dz_i_cells;
        NumType h_s_new = thickness_s;

        NumType k_s = Params<NumType>::Conductivity(this->snow_k_param, T_ss_old, (NumType)0.0, rho_s);

        NumType current_err = std::numeric_limits<NumType>::max();

        std::vector<NumType> full_temp_vec(T_i_cells.size() + 2);
        std::vector<NumType> prev_temp_vec(T_i_cells.size() + 2);

        std::vector<NumType> old_temp_vec = concatenate<NumType>({T_i_cells, std::vector<NumType>{T_is}, std::vector<NumType>{T_ss_old}});

        std::vector<NumType> radiation_nodes_ice(T_i_cells.size() + 1);

        NumType omega_ss_sublim;
        NumType omega_interface = (NumType)0.0;

        NumType kappa_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.kappa_i : IceConsts<NumType>::kappa_i;
        NumType kappa_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.kappa_s :SnowConsts<NumType>::kappa_s;
        NumType i0_i = (Configured()) ? GetConfigConsts<NumType>()->IceConsts.i0_i : IceConsts<NumType>::i0_i;
        NumType i0_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.i0_s : SnowConsts<NumType>::i0_s;

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_i_cells_prev = T_i_cells_new;
            T_is_prev = T_is_new;

            // compute new value of omega at the ice base
            omega_ib = this->W_from_BC(T_ib,
                                       T_i_cells_prev,
                                       dz_i_cells_new,
                                       salinity_i_cells,
                                       rho_i_cells,
                                       true,
                                       false);
            
            

            // compute new value of omega at the snow surface
            omega_ss = this->W_from_BC_0D((NumType)0.0,
                                          T_is_prev,
                                          h_s_new,
                                          k_s,
                                          (NumType)0.0,
                                          rho_s,
                                          false,
                                          true);
            

            // add sublimation 
            NumType r_s = (Configured()) ? GetConfigConsts<NumType>()->SnowConsts.rho_s : SnowConsts<NumType>::rho_s;
            NumType Ls = (Configured()) ? GetConfigConsts<NumType>()->GenConsts.L_sublim : GenConsts<NumType>::L_sublim;
            omega_ss_sublim = (this->is_sublimation) ?
                -this->F_lh((NumType)0.0)/(r_s*Ls) :
                (NumType)0.0;
            
            NumType r_w = (Configured()) ? GetConfigConsts<NumType>()->WaterConsts.rho_w : WaterConsts<NumType>::rho_w;

            if (atm_temperature < (NumType)0.0)
                omega_ss -= precipitation_rate*r_w/r_s;

            // compute ice-snow interface omega value
            if (this->si_transition_mode == SnowIceTransition::SnowAging)
            {
                omega_interface = ((NumType)1.0/GenConsts<NumType>::si_trans_time_scale)*h_s_new;
            }

            // recalculate ice thickness
            dz_i_cells_new = this->Update_dz(dz_i_cells,
                                             omega_ib + (*(this->omega_ib_expl_mass)),
                                             omega_interface);

            // recalculate snow thickness
            h_s_new = this->Update_dz_0D(thickness_s,
                                         omega_interface,
                                         omega_ss + omega_ss_sublim);

            // recalculate ice radiation
            if (this->is_radiation)
            {
                NumType h_i = sum_vec<NumType>(dz_i_cells_new);
                NumType Tfi = Params<NumType>::TempFusion((*(this->Si_cells)).back());
                auto ice_albedo_param = this->ice_albedo_param;
                FuncPtr<NumType> alb_i = [ice_albedo_param, h_i, Tfi] (NumType T)
                {
                    return Params<NumType>::Albedo(ice_albedo_param, T, h_i, Tfi);
                };

                
                auto snow_albedo_param = this->snow_albedo_param;
                FuncPtr<NumType> alb_s = [snow_albedo_param, h_s_new, h_i, Tfi] (NumType T)
                {
                    return Params<NumType>::Albedo(snow_albedo_param, T, h_s_new, 0.0);
                };

                radiation_nodes_ice = this->Compute_radiation_nodes(this->F_sw((NumType)0.0),
                                                                    alb_i(T_is_new),
                                                                    i0_i,
                                                                    kappa_i,
                                                                    dz_i_cells_new,
                                                                    alb_s(0.0),
                                                                    i0_s,
                                                                    kappa_s,
                                                                    std::vector<NumType>{h_s_new}).first;
            }
            
            // recalculate ice temperature profile
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_i_cells_prev, T_i_cells,
                                                                T_is_new, T_is_prev, T_is,
                                                                T_ib, T_ib, T_ib, 
                                                                omega_ib, (NumType)0.0, 
                                                                dz_i_cells_new, dz_i_cells,
                                                                salinity_i_cells,
                                                                rho_i_cells,
                                                                radiation_nodes_ice,
                                                                true);

            // solve linear system and update temperatures
            T_i_cells_new = thomas_solver<NumType>(std::get<0>(matrix_rhs),
                                                   std::get<1>(matrix_rhs),
                                                   std::get<2>(matrix_rhs),
                                                   std::get<3>(matrix_rhs));


            // recalculate ice-snow interface temperature 
            NumType top_ice_k = Params<NumType>::Conductivity(this->ice_k_param, T_i_cells_new.back(), salinity_i_cells.back(), rho_i_cells.back());
            NumType ratio = (top_ice_k*h_s_new)/(k_s*dz_i_cells_new.back()*(NumType)0.5);
            T_is_new  = ((NumType)1.0/(ratio + (NumType)1.0))*(NumType)0.0 + ((ratio)/(ratio + (NumType)1.0))*T_i_cells_new.back();

            
            // evaluate the error
            prev_temp_vec = concatenate<NumType>({T_i_cells_prev, std::vector<NumType>{T_is_prev}, std::vector<NumType>{T_ss_old}});
            
            full_temp_vec = concatenate<NumType>({T_i_cells_new, std::vector<NumType>{T_is_new}, std::vector<NumType>{0.0}});

            current_err = L2_norm(full_temp_vec - prev_temp_vec)/(L2_norm(old_temp_vec) + (NumType)0.1);

            //if (current_err < tol)
            //    break;
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }

        return {{T_i_cells_new, std::vector<NumType>{T_is_new}, dz_i_cells_new}, h_s_new};
    }

    // explicit instantiation
    template class SeaIce_Solver<float>;
    template class SeaIce_Solver<double>;
}