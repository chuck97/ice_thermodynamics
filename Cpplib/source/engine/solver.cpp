#include "solver.hpp"

namespace icethermo
{
    template <typename NumType>
    ThermoSolver<NumType>::ThermoSolver(Mesh<NumType>* mesh_ice_,
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
                                        SnowIceTransition si_transition_mode_)
    {
        this->mesh_ice = mesh_ice_;
        this->mesh_snow = mesh_snow_;
        this->time_step = time_step_;
        this->min_ice_thick = min_ice_thick_;
        this->min_snow_thick = min_snow_thick_;
        this->grad_approx_order = grad_approx_order_;
        this->ice_k_param = ice_k_param_;
        this->ice_c_eff_param = ice_c_eff_param_;
        this->ice_E_param = ice_E_param_;
        this->ice_L_param = ice_L_param_;
        this->ice_albedo_param = ice_albedo_param_;
        this->snow_k_param = snow_k_param_;
        this->snow_c_eff_param = snow_c_eff_param_;
        this->snow_E_param = snow_E_param_;
        this->snow_L_param = snow_L_param_;
        this->snow_albedo_param = snow_albedo_param_;
        this->prec_rate = std::make_shared<NumType>(0.0);
        this->atm_temp = std::make_shared<NumType>(0.0);
        this->atm_press = std::make_shared<NumType>(101325.0); // Pa
        this->atm_humid = std::make_shared<NumType>(0.0); // g/kg
        this->abs_wind_speed = std::make_shared<NumType>(5.0);

        NumType C_sh = GenConsts<NumType>::C_sh;
        this->sh_trans_coeff = std::make_shared<NumType>(C_sh);

        NumType C_lh = GenConsts<NumType>::C_lh;
        this->lh_trans_coeff = std::make_shared<NumType>(C_lh);

        NumType atm_dens = AirConsts<NumType>::rho_a;
        this->atm_dens = std::make_shared<NumType>(atm_dens);

        this->is_radiation = is_radiation_;
        this->is_sublimation = is_sublimation_;
        this->si_transition_mode = si_transition_mode_;
        this->is_verbose = is_verbose_;
    }

    template <typename NumType>
    ThermoSolver<NumType>::~ThermoSolver()
    {
        mesh_ice = NULL;
        mesh_snow = NULL;

        Ti_cells = NULL;
        dzi_cells = NULL;
        Si_cells = NULL;
        rhoi_cells = NULL;    
        Ti_s = NULL;
        Ti_b = NULL;

        Ts_cells = NULL;
        dzs_cells = NULL;
        rhos_cells = NULL;
        Ts_s = NULL;
        Ts_b = NULL; 

        F_up = [](NumType T){return 0.0;};
        F_down = [](NumType T){return 0.0;};
        F_sw = [](NumType T){return 0.0;};
        F_lh = [](NumType T){return 0.0;};
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateUpperFlux(FuncPtr<NumType> F_up_)
    {
        this->F_up = F_up_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateUpperFlux()
    {
        if ((this->mesh_snow != NULL) and 
            ((this->mesh_snow)->GetTotalThickness() > this->min_snow_thick))
        {
            // assign common parameterizations
            this->UpdateLatentHeatFlux();
            this->UpdateSensibleHeatFlux();
            this->UpdatePrecipitationHeatFlux();
            this->UpdateEmittingHeatFlux();

            NumType em_s = GenConsts<NumType>::emissivity;
            NumType i0_s = SnowConsts<NumType>::i0_s;
            auto al_param = this->snow_albedo_param;
            NumType hs = sum_vec<NumType>(*(this->dzs_cells));

            FuncPtr<NumType> F_lws = this->F_lw;
            FuncPtr<NumType> F_lwis = this->F_lwi;
            FuncPtr<NumType> F_sws = this->F_sw;
            FuncPtr<NumType> F_shs = this->F_sh;
            FuncPtr<NumType> F_lhs = this->F_lh;
            FuncPtr<NumType> F_Ps = this->F_P;

            //std::cout << "snow temp:" << *(this->Ts_s) << std::endl;
            //std::cout << "F_lws:" << F_lws(*(this->Ts_s)) << std::endl;
            //std::cout << "F_lwis:" << F_lwis(*(this->Ts_s)) << std::endl;
            //std::cout << "F_sws:" << F_sws(*(this->Ts_s)) << std::endl;
            //std::cout << "F_shs:" << F_shs(*(this->Ts_s)) << std::endl;
            //std::cout << "F_lhs:" << F_lhs(*(this->Ts_s)) << std::endl;
            //std::cout << "F_Ps:" << F_Ps(*(this->Ts_s)) << std::endl;

            FuncPtr<NumType> alb_s = [al_param, hs] (NumType T)
            {
                //return SnowConsts<NumType>::albedo_wet_s;
                return Params<NumType>::Albedo(al_param, T, hs, 0.0);
            };

            this->F_up = [em_s,
                          i0_s,
                          alb_s,
                          F_lws,
                          F_lwis,
                          F_sws,
                          F_shs,
                          F_lhs,
                          F_Ps
                          ](NumType T)
            {
                return    em_s*F_lws(T)
                        - em_s*F_lwis(T)
                        + (1.0 - alb_s(T))*(1.0 - i0_s)*F_sws(T)
                        + F_shs(T)  
                        + F_lhs(T) 
                        //+ F_Ps(T)
                        ;
            };
        }
        else if (this->mesh_ice != NULL)
        {
            // assign common parameterizations
            this->UpdateLatentHeatFlux();
            this->UpdateSensibleHeatFlux();
            this->UpdatePrecipitationHeatFlux();
            this->UpdateEmittingHeatFlux();

            NumType em_i = GenConsts<NumType>::emissivity;
            
            IceConsts<NumType>::albedo_i;
            NumType i0_i = IceConsts<NumType>::i0_i;

            FuncPtr<NumType> F_lwi = this->F_lw;
            FuncPtr<NumType> F_lwii = this->F_lwi;
            FuncPtr<NumType> F_swi = this->F_sw;
            FuncPtr<NumType> F_shi = this->F_sh;
            FuncPtr<NumType> F_lhi = this->F_lh;
            FuncPtr<NumType> F_Pi = this->F_P;

            auto al_param = this->ice_albedo_param;
            NumType hi = sum_vec<NumType>(*(this->dzi_cells));
            NumType Tfi = GenConsts<NumType>::TempFusion((*(this->Si_cells)).back());

            FuncPtr<NumType> alb_i = [al_param, hi, Tfi] (NumType T)
            {
                return Params<NumType>::Albedo(al_param, T, hi, Tfi);
            };

            //std::cout << "ice temp:" << *(this->Ti_s) << std::endl;
            //std::cout << "F_lwi:" << F_lwi(*(this->Ti_s)) << std::endl;
            //std::cout << "F_lwii:" << F_lwii(*(this->Ti_s)) << std::endl;
            //std::cout << "F_swi:" << F_swi(*(this->Ti_s)) << std::endl;
            //std::cout << "F_shi:" << F_shi(*(this->Ti_s)) << std::endl;
            //std::cout << "F_lhi:" << F_lhi(*(this->Ti_s)) << std::endl;
            //std::cout << "F_Pi:" << F_Pi(*(this->Ti_s)) << std::endl;

            this->F_up = [em_i,
                          alb_i,
                          i0_i,
                          F_lwi,
                          F_lwii,
                          F_swi,
                          F_shi,
                          F_lhi,
                          F_Pi
                          ](NumType T)
            {
                return    em_i*F_lwi(T)
                        - em_i*F_lwii(T) 
                        + (1.0 - alb_i(T))*(1.0 - i0_i)*F_swi(T) 
                       + F_shi(T)   
                       + F_lhi(T) 
                       //+ F_Pi(T)
                       ;
            };
        }
        else
        {
            THERMO_ERR((std::string) "At least ice mesh should be assigned to assemble default atm flux\n");
        }
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateLowerFlux(FuncPtr<NumType> F_down_)
    {
        this->F_down = F_down_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateShortWaveRadiation(FuncPtr<NumType> F_sw_)
    {
        this->F_sw = F_sw_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateLongWaveRadiation(FuncPtr<NumType> F_lw_)
    {
        this->F_lw = F_lw_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateLatentHeatFlux(FuncPtr<NumType> F_lh_)
    {
        this->F_lh = F_lh_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateLatentHeatFlux()
    {
        auto press = this->atm_press;
        auto ws = this->abs_wind_speed;
        auto ah = this->atm_humid;
        auto Clh = this->lh_trans_coeff;
        auto ra = this->atm_dens;

        NumType c1 = IceConsts<NumType>::c1_i;
        NumType c2 = IceConsts<NumType>::c2_i;
        NumType T_0 = GenConsts<NumType>::T0;
        NumType Ls = GenConsts<NumType>::L_s;

        this->F_lh = [c1, c2, T_0, Ls, press, ws, ah, Clh, ra](NumType T)
        {   
            //NumType es = (NumType)(6.11*exp(c1*T/(T + T_0 - c2)));
            //NumType q_surf = 0.622*es/((*press)*1e-2 - 0.378*es);


            NumType pb = 21.85/(T + GenConsts<NumType>::T0 - 7.65);
            NumType qmax_dup = 3.80042e-3*exp(pb*T)/((*press)*1e-5);
            NumType q_surf = qmax_dup*0.98;
             
            return (*ra)*Ls*(*Clh)*(*ws)*((*ah)*1e-3 - q_surf);
        };
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateSensibleHeatFlux(FuncPtr<NumType> F_sh_)
    {
        this->F_sh = F_sh_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateSensibleHeatFlux()
    {
        NumType cpa = AirConsts<NumType>::cp_a;

        auto ws = this->abs_wind_speed;
        auto at = this->atm_temp;
        auto Csh = this->sh_trans_coeff;
        auto ra = this->atm_dens;

        this->F_sh = [cpa, ws, at, Csh, ra](NumType T)
        {
            return (*ra)*cpa*(*Csh)*(*ws)*((*at) - T);
        };
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdatePrecipitationHeatFlux(FuncPtr<NumType> F_p_)
    {
        this->F_P = F_p_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdatePrecipitationHeatFlux()
    {
        NumType rw = WaterConsts<NumType>::rho_w;
        NumType cpw = WaterConsts<NumType>::c_pw;

        auto pr = this->prec_rate;
        auto at = this->atm_temp;

        this->F_P = [rw, cpw, pr, at](NumType T)
        {
            return 0.0;
            //return (((*at) - T) > 0.0) ? (*pr)*rw*cpw*((*at) - T) : 0.0;
        };
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateEmittingHeatFlux(FuncPtr<NumType> F_lwi_)
    {
        this->F_lwi = F_lwi_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateEmittingHeatFlux()
    {
        NumType sig = GenConsts<NumType>::sigma;
        NumType T_0 = GenConsts<NumType>::T0;

        this->F_lwi = [sig, T_0](NumType T)
        {
            return sig*(T + T_0)*(T + T_0)*(T + T_0)*(T + T_0);
        };
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdatePrecipitationRate(NumType prec_rate_mm_sm1_)
    {
        *(this->prec_rate) = prec_rate_mm_sm1_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateAtmosphereTemperature(NumType atm_temp_)
    {
        *(this->atm_temp) = atm_temp_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateAtmospherePressure(NumType atm_press_)
    {
        *(this->atm_press) = atm_press_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateAirSpecificHumidity(NumType atm_humid_)
    {
        *(this->atm_humid) = atm_humid_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateAbsWindSpeed(NumType abs_wind_speed_)
    {
        *(this->abs_wind_speed) = abs_wind_speed_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateShTransCoeff(NumType sh_coeff_)
    {
        *(this->sh_trans_coeff) = sh_coeff_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateLhTransCoeff(NumType lh_coeff_)
    {
        *(this->lh_trans_coeff) = lh_coeff_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::UpdateAirDensity(NumType atm_dens_)
    {
        *(this->atm_dens) = atm_dens_;
    }

    template<typename NumType>
    void ThermoSolver<NumType>::AddPrecipitation()
    {
        // if there are precipitations with less than zero atm temp add snow
        if (*(this->atm_temp) < (NumType)0.0)
        {
            NumType snow_thick_before = sum_vec(*(this->dzs_cells));

            if ((*(this->dzs_cells)).size() == 1)
            {
                // update snow thickness according to precipitation rate
                (*(this->dzs_cells))[0] = this->Update_dz_0D((*(this->dzs_cells))[0],
                                                              (NumType)0.0, 
                                                              -(*(this->prec_rate))*
                                                              WaterConsts<NumType>::rho_w/
                                                              SnowConsts<NumType>::rho_s);
            }
            else
            {
                // update snow thickness according to precipitation rate
                (*(this->dzs_cells)) = this->Update_dz((*(this->dzs_cells)),
                                                          (NumType)0.0, 
                                                          -(*(this->prec_rate))*
                                                          WaterConsts<NumType>::rho_w/
                                                          SnowConsts<NumType>::rho_s);
            }

                
            // if snow appeared initialize snow temperatures
            if ((sum_vec(*(this->dzs_cells)) > this->min_snow_thick) and 
                 (snow_thick_before < this->min_snow_thick))
            {
                *(this->Ts_s) = *(this->atm_temp);
                *(this->Ts_b) = *(this->Ti_s);
                for (int i = 0; i < (*(this->dzs_cells)).size(); ++i)
                {
                    (*(this->Ts_cells))[i] = (*(this->Ts_b)) + 1.0*(i + 0.5)/((*(this->dzs_cells)).size())*((*(this->Ts_s)) - (*(this->Ts_b)));
                }
            }
        }
    }

    template<typename NumType>
    std::vector<NumType> ThermoSolver<NumType>::Update_dz(const std::vector<NumType>& dz_cells_old,
                                                          NumType omega_down,
                                                          NumType omega_up)
    {
        NumType thickness = sum_vec(dz_cells_old);
        return dz_cells_old - (this->time_step*(omega_up - omega_down)/thickness)*dz_cells_old;
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::Update_dz_0D(NumType dz_old,
                                                NumType omega_down,
                                                NumType omega_up)
    {
        return dz_old - this->time_step*(omega_up - omega_down);
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::T_from_BC(const std::vector<NumType>& T_cells,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             const std::vector<NumType>& rho_cells,
                                             NumType omega_value,
                                             bool is_ice,
                                             bool is_surface)
    {
        Kparam kparam = (is_ice) ? this->ice_k_param: this->snow_k_param;
        Lparam Lparam = (is_ice) ? this->ice_L_param: this->snow_L_param;

        // surface gradient approximation
        if (is_surface)
        {
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T - T_cells.back())/(0.5*dz_cells.back());};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                int dz_size =  dz_cells.size();
                NumType h1 = (NumType)0.5*dz_cells.back();
                NumType h2 = (NumType)dz_cells.back() + (NumType)0.5*dz_cells[dz_size-2];
                grad = [&T_cells, h1, h2, dz_size](NumType T) {return (T*(h2*h2 - h1*h1) - T_cells.back()*h2*h2 + T_cells[dz_size-2]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // physical constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells.back(), rho_cells.back());};
            FuncPtr<NumType> L = [&salinity_cells, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells.back());};

            FuncPtr<NumType> nonlin_func = [&rho_cells, omega_value, this, &k, &grad, &L](NumType T){return rho_cells.back()*L(T)*omega_value + this->F_up(T) - k(T)*grad(T);};
            
            // solve nonlinear 1D equation
            auto secant_res = bisection_solver<NumType>(nonlin_func, (NumType)-50.0, (NumType)1.0);

           return std::get<0>(secant_res);
        }
        else
        {
            // surface gradient approximation
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T_cells[0] - T)/(0.5*dz_cells[0]);};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                NumType h1 = (NumType)0.5*dz_cells[0];
                NumType h2 = dz_cells[0] + (NumType)0.5*dz_cells[1];
                grad = [&T_cells, h1, h2](NumType T) {return (-T*(h2*h2 - h1*h1) + T_cells[0]*h2*h2 - T_cells[1]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // physical constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells[0], rho_cells[0]);};
            FuncPtr<NumType> L = [&salinity_cells, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells[0]);};

            // assemble nonlinear function for 1D solver
            FuncPtr<NumType> nonlin_func = [&rho_cells, omega_value, this, &k, &grad, &L](NumType T){return rho_cells[0]*L(T)*omega_value + this->F_down(T) - k(T)*grad(T);};
            
            // solve nonlinear 1D equation
            auto secant_res = bisection_solver<NumType>(nonlin_func, (NumType)-50.0, (NumType)1.0);
            
            return std::get<0>(secant_res);
        }   
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::T_from_BC_0D(FuncPtr<NumType> T_op_interface,
                                                NumType thickness,
                                                NumType k_value,
                                                NumType sal_value,
                                                NumType omega_value,
                                                NumType rho_value,
                                                bool is_ice,
                                                bool is_surface)
    {
        Lparam Lparam = (is_ice) ? this->ice_L_param: this->snow_L_param;
        FuncPtr<NumType> L = [&sal_value, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, sal_value);};
        FuncPtr<NumType> grad;

        if (is_surface)
        {
            grad = [T_op_interface, &thickness](NumType T) {return (T - T_op_interface(T))/thickness;};

            // assemble nonlinear function for 1D solver
            FuncPtr<NumType> nonlin_func = [&rho_value, omega_value, this, &k_value, &grad, &L](NumType T){return this->F_up(T) - k_value*grad(T);};

            // solve nonlinear 1D equation
            auto res = bisection_solver<NumType>(nonlin_func, (NumType)-50.0, (NumType)1.0);

            return std::get<0>(res);
        }
        else
        {
            grad = [T_op_interface, &thickness](NumType T) {return (T_op_interface(T) - T)/thickness;};

            // assemble nonlinear function for 1D solver
            FuncPtr<NumType> nonlin_func = [&rho_value, omega_value, this, &k_value, &grad, &L](NumType T){return this->F_down(T) - k_value*grad(T);};

            // solve nonlinear 1D equation
            auto res = bisection_solver<NumType>(nonlin_func, (NumType)-50.0, (NumType)1.0);

            return std::get<0>(res);
        }
    }

    template<typename NumType>
    NumType ThermoSolver<NumType>::W_from_BC(NumType T_bnd,
                                             const std::vector<NumType>& T_cells,
                                             const std::vector<NumType>& dz_cells,
                                             const std::vector<NumType>& salinity_cells,
                                             const std::vector<NumType>& rho_cells,
                                             bool is_ice,
                                             bool is_surface)
    {
        Kparam kparam = (is_ice) ? this->ice_k_param: this->snow_k_param;
        Lparam Lparam = (is_ice) ? this->ice_L_param: this->snow_L_param;

        if (is_surface)
        {
            // surface gradient approximation
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T - T_cells.back())/(0.5*dz_cells.back());};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                int dz_size =  dz_cells.size();
                NumType h1 = 0.5*dz_cells.back();
                NumType h2 = dz_cells.back() + 0.5*dz_cells[dz_size-2];
                grad = [&T_cells, h1, h2, dz_size](NumType T) {return (T*(h2*h2 - h1*h1) - T_cells.back()*h2*h2 + T_cells[dz_size-2]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells.back(), rho_cells.back());};
            FuncPtr<NumType> L = [&salinity_cells, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, salinity_cells.back());};

            // calculate omega from boundary conditions
            return (k(T_cells.back())*grad(T_bnd) - this->F_up(T_bnd))/(rho_cells.back()*L(T_cells.back()));
        }
        else
        {
            // surface gradient approximation
            FuncPtr<NumType> grad;

            if (this->grad_approx_order == ApproxOrder::first)
            {
                grad = [&T_cells, &dz_cells](NumType T) {return (T_cells[0] - T)/(0.5*dz_cells[0]);};
            }
            else if (this->grad_approx_order == ApproxOrder::second)
            {
                NumType h1 = 0.5*dz_cells[0];
                NumType h2 = dz_cells[0] + 0.5*dz_cells[1];
                grad = [&T_cells, h1, h2](NumType T) {return (-T*(h2*h2 - h1*h1) + T_cells[0]*h2*h2 - T_cells[1]*h1*h1)/(h1*h2*(h2-h1));};
            }
            else
            {
                THERMO_ERR("Available gradient approximation order: first, second!");
            }

            // physical constants
            FuncPtr<NumType> k = [&salinity_cells, &rho_cells, kparam](NumType T){return Params<NumType>::Conductivity(kparam, T, salinity_cells[0], rho_cells[0]);};
            std::function<NumType(NumType, NumType)> L = [&salinity_cells, Lparam](NumType T, NumType S){return Params<NumType>::FusionHeat(Lparam, T, S);};

            // calculate omega from boundary conditions
            if (kparam == Kparam::BubblyBrine or kparam == Kparam::Untersteiner)
            {
                if ((k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd)) > 0)
                {
                    // basal melt with 4 psu salinity
                    return (k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd))/(rho_cells[0]*L(T_cells[0], 4.0));
                }
                else
                {
                    // basal growth with 10 psu salinity
                    return (k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd))/(rho_cells[0]*L(T_cells[0], 10.0));
                }
            }
            else
            {
                return (k(T_cells[0])*grad(T_bnd) + this->F_down(T_bnd))/(rho_cells[0]*L(T_cells[0], salinity_cells[0]));
            }
        }
    }

    // find omega value consistent with boundary conditions (0D)
    template<typename NumType>
    NumType ThermoSolver<NumType>::W_from_BC_0D(NumType T_bnd,
                                                NumType T_op_interface,
                                                NumType thickness,
                                                NumType k_value,
                                                NumType sal_value,
                                                NumType rho_value,
                                                bool is_ice,
                                                bool is_surface)
    {
        Lparam Lparam = (is_ice) ? this->ice_L_param: this->snow_L_param;
        FuncPtr<NumType> L = [&sal_value, Lparam](NumType T){return Params<NumType>::FusionHeat(Lparam, T, sal_value);};
        FuncPtr<NumType> grad;
        if (is_surface)
        {
            grad = [&T_op_interface, &thickness](NumType T) {return (T - T_op_interface)/thickness;};
            return (k_value*grad(T_bnd) - this->F_up(T_bnd))/(rho_value*L(T_bnd));
        }
        else
        {
            grad = [&T_op_interface, &thickness](NumType T) {return (T_op_interface - T)/thickness;};
            return (k_value*grad(T_bnd) - this->F_down(T_bnd))/(rho_value*L(T_bnd));
        }
    }

    // assembling of tridiagonal advection-diffusion matrix and rhs
    template<typename NumType>
    FourVecs<NumType> ThermoSolver<NumType>::Assemble_advdiff_martix_rhs(const std::vector<NumType>& T_cells_prev,
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
                                                                         bool is_ice)
    {
        // get all parameterizations
        Kparam kparam = (is_ice) ? this->ice_k_param: this->snow_k_param;
        Cparam cparam = (is_ice) ? this->ice_c_eff_param: this->snow_c_eff_param;
        Eparam Eparam = (is_ice) ? this->ice_E_param: this->snow_E_param;

        // get number of cells
        int N = T_cells_old.size();

        // check the size of input arrays
        if (T_cells_prev.size() != N)
        {
            THERMO_ERR("Wrong size of \'T_cells_prev\'!");
        }

        if (dz_cells_new.size() != N)
        {
            THERMO_ERR("Wrong size of \'dz_cells_new\'!")
        }

        if (dz_cells_old.size() != N)
        {
            THERMO_ERR("Wrong size of \'dz_cells_old\'!");
        }

        if (radiation_nodes.size() != N+1)
        {
            THERMO_ERR("Wrong size of \'radiation_nodes\'!");
        }

        // compute effective thermal conductivity at nodes
        std::vector<NumType> eff_k_nodes(N+1);
        
        
        eff_k_nodes[0] = 2.0*Params<NumType>::Conductivity(kparam, T_cells_prev[0], salinity_cells[0], rho_cells[0])/dz_cells_new[0];
        
        for (int i = 1; i < N; ++i)
        {
            NumType k_prev = Params<NumType>::Conductivity(kparam, T_cells_prev[i-1], salinity_cells[i-1], rho_cells[i-1]);
            NumType k_forw = Params<NumType>::Conductivity(kparam, T_cells_prev[i], salinity_cells[i], rho_cells[i]);
            eff_k_nodes[i] = 2.0*k_prev*k_forw/(k_prev*dz_cells_new[i] + k_forw*dz_cells_new[i-1]);
        }

        eff_k_nodes[N] = 2.0*Params<NumType>::Conductivity(kparam, T_cells_prev[N-1], salinity_cells[N-1], rho_cells[N-1])/dz_cells_new[N-1];

        // compute effective heat capacity and enthalpy at the cells and interface nodes
        std::vector<NumType> eff_c_cells(N);
        std::vector<NumType> E_cells(N);


        NumType eff_c_down = Params<NumType>::EffCapacity(cparam, T_down_prev, T_down_old, salinity_cells[0]);
        NumType E_down = Params<NumType>::Enthalpy(Eparam, T_down_old, salinity_cells[0]); 
        
        for (int i = 0; i < N; ++i)
        {
            eff_c_cells[i] = Params<NumType>::EffCapacity(cparam, T_cells_prev[i], T_cells_old[i], salinity_cells[i]);
            E_cells[i] = Params<NumType>::Enthalpy(Eparam, T_cells_old[i],  salinity_cells[i]);
        }

        NumType eff_c_up = Params<NumType>::EffCapacity(cparam, T_up_prev, T_up_old, salinity_cells[N-1]);
        NumType E_up = Params<NumType>::Enthalpy(Eparam, T_up_old, salinity_cells[N-1]);

        // compute nodal values of omega
        std::vector<NumType> omega_nodes(N+1);

        for (int i = 0; i < N+1; ++i)
        {
            omega_nodes[i] = omega_down + (sum_vec(dz_cells_new, 0, i)/sum_vec(dz_cells_new))*(omega_up - omega_down);
        }

        // construct diagonals of matrix and rhs vector
        std::vector<NumType> A(N);
        std::vector<NumType> B(N);
        std::vector<NumType> C(N);
        std::vector<NumType> RHS(N);

        // first row
        B[0] = rho_cells[0]*eff_c_cells[0]*dz_cells_new[0]/this->time_step + eff_k_nodes[1] + eff_k_nodes[0];

        C[0] = -eff_k_nodes[1];

        RHS[0] = rho_cells[0]*eff_c_cells[0]*dz_cells_new[0]*T_cells_old[0]/this->time_step - 
                  rho_cells[0]*E_cells[0]*(dz_cells_new[0] - dz_cells_old[0])/this->time_step +
                  (radiation_nodes[1] - radiation_nodes[0]) + 
                  eff_k_nodes[0]*T_down_new;
        
        if (omega_nodes[0] >= 0)
        {
            RHS[0] += rho_cells[0]*eff_c_down*T_down_new*omega_nodes[0] - 
                      rho_cells[0]*(eff_c_down*T_down_old - E_down)*omega_nodes[0];
        }
        else
        {
            B[0] += -rho_cells[0]*eff_c_cells[0]*omega_nodes[0];

            RHS[0] += -rho_cells[0]*(eff_c_cells[0]*T_cells_old[0] - E_cells[0])*omega_nodes[0];
        }

        if (omega_nodes[1] >= 0)
        {
            B[0] += rho_cells[0]*eff_c_cells[0]*omega_nodes[1];

            RHS[0] += rho_cells[0]*(eff_c_cells[0]*T_cells_old[0] - E_cells[0])*omega_nodes[1];
        }
        else
        {
            C[0] += rho_cells[0]*eff_c_cells[1]*omega_nodes[1];

            RHS[0] += rho_cells[0]*(eff_c_cells[1]*T_cells_old[1] - E_cells[1])*omega_nodes[1];
        }

        // middle rows
        for (int i = 1; i < N-1; ++i)
        {
            A[i] = -eff_k_nodes[i];

            B[i] = rho_cells[i]*eff_c_cells[i]*dz_cells_new[i]/this->time_step + eff_k_nodes[i+1] + eff_k_nodes[i];
            
            C[i] =  -eff_k_nodes[i+1];
            
            RHS[i] = rho_cells[i]*eff_c_cells[i]*dz_cells_new[i]*T_cells_old[i]/this->time_step -
                      rho_cells[i]*E_cells[i]*(dz_cells_new[i] - dz_cells_old[i])/this->time_step +
                      (radiation_nodes[i+1] - radiation_nodes[i]);

            if (omega_nodes[i] >= 0)
            {
                A[i] += -rho_cells[i]*eff_c_cells[i-1]*omega_nodes[i];

                RHS[i] += -rho_cells[i]*(eff_c_cells[i-1]*T_cells_old[i-1] - E_cells[i-1])*omega_nodes[i];
            }
            else
            {
                B[i] += -rho_cells[i]*eff_c_cells[i]*omega_nodes[i];

                RHS[i] += -rho_cells[i]*(eff_c_cells[i]*T_cells_old[i] - E_cells[i])*omega_nodes[i];
            }
            
            if (omega_nodes[i+1] >= 0)
            {
                B[i] += rho_cells[i]*eff_c_cells[i]*omega_nodes[i+1];

                RHS[i] += rho_cells[i]*(eff_c_cells[i]*T_cells_old[i] - E_cells[i])*omega_nodes[i+1];
            }
            else
            {
                C[i] += rho_cells[i]*eff_c_cells[i+1]*omega_nodes[i+1];

                RHS[i] += rho_cells[i]*(eff_c_cells[i+1]*T_cells_old[i+1] - E_cells[i+1])*omega_nodes[i+1];
            }
        }

        // last row
        A[N-1] = -eff_k_nodes[N-1];

        B[N-1] = rho_cells[N-1]*eff_c_cells[N-1]*dz_cells_new[N-1]/this->time_step + eff_k_nodes[N] + eff_k_nodes[N-1];

        RHS[N-1] = rho_cells[N-1]*eff_c_cells[N-1]*T_cells_old[N-1]*dz_cells_new[N-1]/this->time_step - 
                    rho_cells[N-1]*E_cells[N-1]*(dz_cells_new[N-1] - dz_cells_old[N-1])/this->time_step -
                    (radiation_nodes[N] - radiation_nodes[N-1]) +
                    eff_k_nodes[N]*T_up_new; 
        
        if (omega_nodes[N-1] >= 0)
        {
            A[N-1] += -rho_cells[N-1]*eff_c_cells[N-2]*omega_nodes[N-1];
            RHS[N-1] += -rho_cells[N-1]*(eff_c_cells[N-2]*T_cells_old[N-2] - E_cells[N-2])*omega_nodes[N-1];
        }
        else
        {
            B[N-1] += -rho_cells[N-1]*eff_c_cells[N-1]*omega_nodes[N-1];
            RHS[N-1] += -rho_cells[N-1]*(eff_c_cells[N-1]*T_cells_old[N-1] - E_cells[N-1])*omega_nodes[N-1];
        }

        if (omega_nodes[N] >= 0)
        {
            B[N-1] += rho_cells[N-1]*eff_c_cells[N-1]*omega_nodes[N];
            RHS[N-1] += rho_cells[N-1]*(eff_c_cells[N-1]*T_cells_old[N-1] - E_cells[N-1])*omega_nodes[N];
        }
        else
        {
            RHS[N-1] += -rho_cells[N-1]*eff_c_up*T_up_new*omega_nodes[N] + rho_cells[N-1]*(eff_c_up*T_up_old - E_up)*omega_nodes[N];
        }

        return {std::vector<NumType>{A.begin()+1, A.end()},
                B,
                std::vector<NumType>{C.begin(), C.end()-1},
                RHS};
    }

    // assembling radiation on nodes according to Beer's Law
    template<typename NumType>
    std::pair<std::vector<NumType>, std::vector<NumType>> ThermoSolver<NumType>::Compute_radiation_nodes(NumType F_sw_value,
                                                                                                         NumType albedo_ice,
                                                                                                         NumType i0_ice,
                                                                                                         NumType kappa_ice,
                                                                                                         const std::vector<NumType>& dz_cells_ice,
                                                                                                         NumType albedo_snow,
                                                                                                         NumType i0_snow,
                                                                                                         NumType kappa_snow,
                                                                                                         const std::vector<NumType>& dz_cells_snow)
    {
        std::vector<NumType> dzi_extended = concatenate(std::vector<NumType>{0.0}, reverse_vec(dz_cells_ice));
        if (dz_cells_snow.size() == 0)
        {
            std::vector<NumType> radiation_nodes_ice = F_sw_value*((NumType)1.0 - albedo_ice)*i0_ice*exp_vec(-kappa_ice*cumsum(dzi_extended));
            return {reverse_vec(radiation_nodes_ice), std::vector<NumType>{0.0}};
        }
        else
        {
            std::vector<NumType> ones_array(dz_cells_ice.size() + 1);
            for (int i = 0; i < dz_cells_ice.size() + 1; ++i)
            {
                ones_array[i] = 1.0;
            }

            std::vector<NumType> dzs_extended = concatenate(std::vector<NumType>{0.0}, reverse_vec(dz_cells_snow));
            std::vector<NumType> radiation_nodes_snow = F_sw_value*((NumType)1.0 - albedo_snow)*i0_snow*exp_vec(-kappa_snow*cumsum(dzs_extended));
            std::vector<NumType> radiation_nodes_ice = F_sw_value*((NumType)1.0 - albedo_snow)*i0_snow*
                                                       exp_vec(-kappa_ice*cumsum(dzi_extended) - ones_array*sum_vec(kappa_snow*dz_cells_snow));
            return {reverse_vec(radiation_nodes_ice), reverse_vec(radiation_nodes_snow)};
        }
    }

    // 1d sea-ice freezing mode
    template<typename NumType>
    FourVecs<NumType> ThermoSolver<NumType>::glacier1d_freezing(NumType T_ib,
                                                                const std::vector<NumType>& T_cells,
                                                                NumType T_is,
                                                                const std::vector<NumType>& dz_cells,
                                                                const std::vector<NumType>& salinity_cells,
                                                                const std::vector<NumType>& rho_cells,
                                                                int max_n_its,
                                                                NumType tol)
    {
        NumType T_ib_new = T_ib;
        NumType T_ib_prev = T_ib;

        std::vector<NumType> T_cells_new = T_cells;
        std::vector<NumType> T_cells_prev = T_cells;

        NumType T_is_new = T_is;
        NumType T_is_prev = T_is;

        NumType omega_ib = (NumType)0.0;

        std::vector<NumType> dz_cells_new = dz_cells;

        std::vector<NumType> radiation_nodes(T_cells.size() + 1);

        std::vector<NumType> T_is_history = {T_is};
        std::vector<NumType> T_ib_history = {T_ib};

        NumType current_err = std::numeric_limits<NumType>::max();

        NumType surface_err = std::numeric_limits<NumType>::max();
        NumType base_err = std::numeric_limits<NumType>::max();

        NumType prev_surface_err = surface_err;
        NumType prev_base_err = base_err;

        std::vector<NumType> full_temp_vec(T_cells.size() + 2);
        std::vector<NumType> prev_temp_vec(T_cells.size() + 2);
        std::vector<NumType> old_temp_vec = concatenate<NumType>({std::vector<NumType>{T_ib}, T_cells, std::vector<NumType>{T_is}});

        NumType omega_ss;

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;

            T_ib_prev = T_ib_new;
            T_cells_prev = T_cells_new;
            T_is_prev = T_is_new;
            prev_surface_err = surface_err;
            prev_base_err = base_err;

            // compute new value of ice base temperature
            T_ib_new = this->T_from_BC(T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       0.0,
                                       true,
                                       false);
            
            // add sublimation
            omega_ss = (this->is_sublimation) ? 
                -this->F_lh(T_is_new)/(rho_cells.back()*GenConsts<NumType>::L_s) :
                (NumType)0.0;
            
            // compute new value of ice surface temperature
            T_is_new = this->T_from_BC(T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       omega_ss,
                                       true,
                                       true);
            
            // force the convergence of base temperature
            base_err = std::abs(T_ib_new - T_ib_prev)/(std::abs(T_ib) + (NumType)0.1);
            
            if (base_err < prev_base_err)
            {
                T_ib_history.push_back(T_ib_new);
            }
            else
            {
                T_ib_new = sum_vec<NumType>(T_ib_history)/T_ib_history.size();
            }

            // force the convergence of surface temperature
            surface_err = std::abs(T_is_new - T_is_prev)/(std::abs(T_is) + (NumType)0.1);
            
            if (surface_err < prev_surface_err)
            {
                T_is_history.push_back(T_is_new);
            }
            else
            {
                T_is_new = sum_vec<NumType>(T_is_history)/T_is_history.size();
            }

            // recalculate ice thickness
            dz_cells_new = this->Update_dz(dz_cells,
                                           0.0,
                                           omega_ss);
            
            // recalculate radiation
            if (this->is_radiation)
            {
                radiation_nodes = this->Compute_radiation_nodes(F_sw(T_is_new),
                                                                IceConsts<NumType>::albedo_i,
                                                                IceConsts<NumType>::i0_i,
                                                                IceConsts<NumType>::kappa_i,
                                                                dz_cells_new).first;
            }

            // assemble matrix and rhs for ice 
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_cells_prev, T_cells,
                                                                T_is_new, T_is_prev, T_is,
                                                                T_ib_new, T_ib_prev, T_ib, 
                                                                0.0, omega_ss, 
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
            prev_temp_vec = concatenate<NumType>({std::vector<NumType>{T_ib_prev}, T_cells_prev, std::vector<NumType>{T_is_prev}});
            
            full_temp_vec = concatenate<NumType>({std::vector<NumType>{T_ib_new}, T_cells_new, std::vector<NumType>{T_is_new}}); 

            current_err = L2_norm(full_temp_vec - prev_temp_vec)/(L2_norm(old_temp_vec) + (NumType)0.1);

            if (current_err < tol)
                break;
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }

        return {std::vector<NumType>{T_ib_new}, T_cells_new, std::vector<NumType>{T_is_new}, dz_cells_new};
    }

    // 1d glacier melting mode
    template<typename NumType>
    ThreeVecs<NumType> ThermoSolver<NumType>::glacier1d_melting(NumType T_ib,
                                                                NumType T_is,
                                                                const std::vector<NumType>& T_cells,
                                                                NumType T_is_old,
                                                                const std::vector<NumType>& dz_cells,
                                                                const std::vector<NumType>& salinity_cells,
                                                                const std::vector<NumType>& rho_cells,
                                                                int max_n_its,
                                                                NumType tol)
    {
        NumType T_ib_prev = T_ib;
        NumType T_ib_new = T_ib;
        std::vector<NumType> T_cells_new = T_cells;
        std::vector<NumType> T_cells_prev = T_cells;

        NumType omega_ib = 0.0;
        NumType omega_is = 0.0;

        std::vector<NumType> dz_cells_new = dz_cells;

        std::vector<NumType> T_ib_history = {T_ib};

        NumType base_err = std::numeric_limits<NumType>::max();

        NumType prev_base_err = base_err;

        std::vector<NumType> radiation_nodes(T_cells.size() + 1);

        NumType current_err = std::numeric_limits<NumType>::max();

        std::vector<NumType> full_temp_vec(T_cells.size() + 1);
        std::vector<NumType> prev_temp_vec(T_cells.size() + 1);
        std::vector<NumType> old_temp_vec = concatenate<NumType>(std::vector<NumType>{T_ib}, T_cells);

        NumType omega_is_sublim;

        int npseudo = 0;

        for (int pseudoit = 0; pseudoit < max_n_its; ++pseudoit)
        {
            npseudo++;
            
            T_cells_prev = T_cells_new;
            T_ib_prev = T_ib_new;
            
            // compute new value of ice base temperature
            T_ib_new = this->T_from_BC(T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       0.0,
                                       true,
                                       false);
            
            // force the convergence of base temperature
            base_err = std::abs(T_ib_new - T_ib_prev)/(std::abs(T_ib) + (NumType)0.1);
            
            if (base_err < prev_base_err)
            {
                T_ib_history.push_back(T_ib_new);
            }
            else
            {
                T_ib_new = sum_vec<NumType>(T_ib_history)/T_ib_history.size();
            }

            // compute new value of omega at the base
            omega_is = this->W_from_BC(T_is,
                                       T_cells_prev,
                                       dz_cells_new,
                                       salinity_cells,
                                       rho_cells,
                                       true,
                                       true);

            // add sublimation
            omega_is_sublim = (this->is_sublimation) ? 
                -this->F_lh(T_is)/(rho_cells.back()*GenConsts<NumType>::L_s) :
                (NumType)0.0;

            // recalculate ice thickness
            dz_cells_new = this->Update_dz(dz_cells,
                                           (NumType)0.0,
                                           omega_is + omega_is_sublim);
            
            // recalculate radiation
            if (this->is_radiation)
            {
                radiation_nodes = this->Compute_radiation_nodes(this->F_sw(T_is),
                                                                IceConsts<NumType>::albedo_i,
                                                                IceConsts<NumType>::i0_i,
                                                                IceConsts<NumType>::kappa_i,
                                                                dz_cells_new).first;
            }

            // assemble matrix and rhs for ice 
            auto matrix_rhs = this->Assemble_advdiff_martix_rhs(T_cells_prev, T_cells,
                                                                T_is, T_is, T_is_old,
                                                                T_ib_new, T_ib_prev, T_ib, 
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
            prev_temp_vec = concatenate<NumType>(std::vector<NumType>{T_ib_prev}, T_cells_prev);
            full_temp_vec = concatenate<NumType>(std::vector<NumType>{T_ib_prev}, T_cells_prev);
            current_err = L2_norm(full_temp_vec - prev_temp_vec)/(L2_norm(old_temp_vec) + (NumType)0.1);

            if (current_err < tol)
                break;
        }

        if (this->is_verbose)
        {
            std::cout << "nits:" << npseudo << ", err:" << current_err << std::endl;
        }

        return {std::vector<NumType>{T_ib_new}, T_cells_new, dz_cells_new};
    }

    // explicit instantiation
    template class ThermoSolver<float>;
    template class ThermoSolver<double>;
}