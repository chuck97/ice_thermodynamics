#include "ice_thermo.h"

using namespace std;

IceThermo::IceThermo(unsigned int n_cells,
                     double T_b,
                     double T_su,
                     double h):
    n_cells_(n_cells),
    n_nodes_(n_cells + 1)
{
    T_cells_old_.resize(n_cells_);
    T_cells_last_.resize(n_cells_);
    T_cells_new_.resize(n_cells_);
    dz_cells_old_.resize(n_cells_);
    dz_cells_new_.resize(n_cells_);
    c_cells_.resize(n_cells_);
    E_cells_.resize(n_cells_);
    dR_cells_.resize(n_cells_);
    S_cells_.resize(n_cells_);

    T_nodes_old_.resize(n_nodes_);
    T_nodes_last_.resize(n_nodes_);
    T_nodes_new_.resize(n_nodes_);
    w_nodes_.resize(n_nodes_);
    c_nodes_.resize(n_nodes_);
    k_nodes_.resize(n_nodes_);
    E_nodes_.resize(n_nodes_);
    a_nodes_.resize(n_nodes_);
    b_nodes_.resize(n_nodes_);
    S_nodes_.resize(n_nodes_);


    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        dz_cells_old_[i] = h/n_cells_;
        T_cells_old_[i] = T_b + ((i + 0.5)/n_cells_)*(T_su - T_b);
        T_cells_new_[i] = T_cells_old_[i];
        T_cells_last_[i] = T_cells_old_[i];
    }
    dz_cells_new_ = dz_cells_old_;

    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        S_cells_[i] = consts.sal_ice;
    }
    RecalculateSalNodes(S_cells_, S_nodes_);
    RecalculateCoeffs();   
}

IceThermo::IceThermo(unsigned int n_cells,
                     const vector<double>& T_init_cells,
                     const vector<double>& dz_init_cells):
    n_cells_(n_cells),
    n_nodes_(n_cells + 1)
{
    T_cells_old_.resize(n_cells_);
    T_cells_last_.resize(n_cells_);
    T_cells_new_.resize(n_cells_);
    dz_cells_old_.resize(n_cells_);
    dz_cells_new_.resize(n_cells_);
    c_cells_.resize(n_cells_);
    E_cells_.resize(n_cells_);
    dR_cells_.resize(n_cells_);
    S_cells_.resize(n_cells_);

    T_nodes_old_.resize(n_nodes_);
    T_nodes_last_.resize(n_nodes_);
    T_nodes_new_.resize(n_nodes_);
    w_nodes_.resize(n_nodes_);
    c_nodes_.resize(n_nodes_);
    k_nodes_.resize(n_nodes_);
    E_nodes_.resize(n_nodes_);
    a_nodes_.resize(n_nodes_);
    b_nodes_.resize(n_nodes_);
    S_nodes_.resize(n_nodes_);

    if (T_init_cells.size() != n_cells_)
    {
        ERR("wrong size of inial T array");
    }

    if (dz_init_cells.size() != n_cells_)
    {
        ERR("wrong size of inial dz array");
    }
    dz_cells_old_ = dz_init_cells;
    dz_cells_new_ = dz_init_cells;

    T_cells_old_ = T_init_cells;
    T_cells_new_ = T_init_cells;
    T_cells_last_ = T_init_cells;
    

    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        S_cells_[i] = consts.sal_ice;
    }

    RecalculateSalNodes(S_cells_, S_nodes_);
    RecalculateCoeffs();
}

void IceThermo::SetTimeStep(double time_step)
{
    dt_ = time_step;    
}

void IceThermo::UpdateFluxes(double flux_boundary, double flux_surface)
{
    F_b_ = flux_boundary;
    F_su_ = flux_surface;
}

void IceThermo::WriteVTU(const string& path,
                         unsigned int step_num) const
{
    ThermoMesh t_m(n_cells_);
	t_m.SetWidth(0.2);
	t_m.AssignCellThickness(dz_cells_old_);

	t_m.AssignCellData("T cells", T_cells_old_);
    t_m.AssignCellData("dz cells", dz_cells_new_);
    t_m.AssignCellData("c cells", c_cells_);
    t_m.AssignCellData("E cells", E_cells_);
    t_m.AssignCellData("dR cells", dR_cells_);
    t_m.AssignCellData("S cells", S_cells_);

    t_m.AssignNodeData("T nodes", T_nodes_old_);
    t_m.AssignNodeData("W nodes", w_nodes_);
    t_m.AssignNodeData("c nodes", c_nodes_);
    t_m.AssignNodeData("k nodes", k_nodes_);
    t_m.AssignNodeData("E nodes", E_nodes_);
    t_m.AssignNodeData("S nodes", S_nodes_);
    
    stringstream ss;
    ss << setfill('0') << setw(5) << step_num;

	t_m.PlotVertical(path + "vert" + ss.str() + ".vtu");
	t_m.PlotSigma(path + "sigm" + ss.str() + ".vtu");
}

void IceThermo::RecalculateCoeffs()
{
    if (!IS_COEFFS_SIGMA)
    {
        // a0, b0
        a_nodes_[0] = (dz_cells_new_[1] + 2.0*dz_cells_new_[0])/(dz_cells_new_[0] + dz_cells_new_[1]);
        b_nodes_[0] = (-dz_cells_new_[0])/(dz_cells_new_[0] + dz_cells_new_[1]);

        // ai, bi
        for (unsigned int i = 1; i < (n_cells_ - 1); ++i)
        {
            a_nodes_[i] = dz_cells_new_[i]/(dz_cells_new_[i] + dz_cells_new_[i-1]);
            b_nodes_[i] = dz_cells_new_[i-1]/(dz_cells_new_[i] + dz_cells_new_[i-1]);
        }

        // aN, bN
        a_nodes_[n_nodes_ - 1] = (-dz_cells_new_[n_nodes_ - 1])/(dz_cells_new_[n_nodes_ - 1] + dz_cells_new_[n_nodes_ - 2]);
        b_nodes_[n_nodes_ - 1] = (dz_cells_new_[n_nodes_ - 2] + 2.0*dz_cells_new_[n_nodes_ - 1])/(dz_cells_new_[n_nodes_ - 1] + dz_cells_new_[n_nodes_ - 2]);
    }
    else
    {
        // calculate total thickness
        double h = 0.0;
        for (unsigned int i = 0; i < n_cells_; ++i)
        {
            h += dz_cells_new_[i];
        }

        // a0, b0
        a_nodes_[0] = 0.5 + (2.0*h/n_cells_)/(dz_cells_new_[0] + dz_cells_new_[1]);
        b_nodes_[0] = 0.5 - (2.0*h/n_cells_)/(dz_cells_new_[0] + dz_cells_new_[1]);

        // aN, bN
        a_nodes_[n_nodes_ - 1] = 0.5 - (2.0*h/n_cells_)/(dz_cells_new_[n_nodes_ - 1] + dz_cells_new_[n_nodes_ - 2]);
        b_nodes_[n_nodes_ - 1] = 0.5 + (2.0*h/n_cells_)/(dz_cells_new_[n_nodes_ - 1] + dz_cells_new_[n_nodes_ - 2]);;

        // ai, bi
        for (unsigned int i = 0; i < (n_cells_ - 1); ++i)
        {
            a_nodes_[i] = 0.5;
            b_nodes_[i] = 0.5;
        }
    }
}

double IceThermo::CalculateE(double T, double S)
{
    double Tf = -consts.mu_ice*S;
    
    if (abs(Tf) < ABS_T_MIN)
    {
        if (Tf <= 0)
        {
            Tf = -ABS_T_MIN;
        }
        else
        {
            Tf = ABS_T_MIN;
        }
        
    }

    if (abs(T) < ABS_T_MIN)
    {
        if (T <= 0)
        {
            T = -ABS_T_MIN;
        }
        else
        {
            T = ABS_T_MIN;
        }
    }

    return (consts.c0_ice*(T - Tf) - consts.L0_ice*(1.0 - Tf/T) + consts.cw_water*Tf);
}

double IceThermo::CalculateK(double T, double S)
{
    if (abs(T) < ABS_T_MIN)
    {
        if (T <= 0)
        {
            T = -ABS_T_MIN;
        }
        else
        {
            T = ABS_T_MIN;
        }
    }
    return (consts.k0_ice + consts.k1_ice*(S/T));
}

double IceThermo::CalculateC(double T, double T_old, double S)
{
    if (abs(T) < ABS_T_MIN)
    {
        if (T <= 0)
        {
            T = -ABS_T_MIN;
        }
        else
        {
            T = ABS_T_MIN;
        }
    }

    if (abs(T_old) < ABS_T_MIN)
    {
        if (T_old <= 0)
        {
            T_old = -ABS_T_MIN;
        }
        else
        {
            T_old = ABS_T_MIN;
        }
    }

    double Tf = -consts.mu_ice*S;
    return (consts.c0_ice - consts.L0_ice*Tf/(T*T_old));
}

void IceThermo::UpdateW()
{
    // bottom w
    double dT_dz_b = (2.0*(1.0 - a_nodes_[0])/dz_cells_new_[0])*T_cells_old_[0] +
                     (-2.0*b_nodes_[0]/dz_cells_new_[0])*T_cells_old_[1];
    
    double T_b = a_nodes_[0]*T_cells_old_[0] + b_nodes_[0]*T_cells_old_[1];
    double E_b = CalculateE(T_b, S_nodes_[0]);
    double k_b = CalculateK(T_b, S_nodes_[0]);
    double w_b = (F_b_ - k_b*dT_dz_b)/(consts.rho_ice*E_b);

    // surface w
    double dT_dz_s = ((2.0*a_nodes_[n_nodes_-1])/dz_cells_new_[n_cells_-1])*T_cells_old_[n_cells_-2] +
                     ((2.0*b_nodes_[n_nodes_-1] - 1.0)/dz_cells_new_[n_cells_-1])*T_cells_old_[n_cells_-1];
    
    double T_s = a_nodes_[n_nodes_-1]*T_cells_old_[n_cells_-2] + b_nodes_[n_nodes_-1]*T_cells_old_[n_cells_-1];
    double E_s = CalculateE(T_s, S_nodes_[n_nodes_-1]);
    double k_s = CalculateK(T_s, S_nodes_[n_nodes_-1]);
    double w_s = (F_su_ - k_s*dT_dz_s)/(consts.rho_ice*E_s);

    // nodes linear interpolation
    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        w_nodes_[i] = w_b + (i/(n_nodes_-1))*(w_s - w_b);
    }   
}

void IceThermo::UpdateDeltaZ()
{
    vector<double> delz_prev = dz_cells_new_;
    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        dz_cells_new_[i] = dz_cells_new_[i] - dt_*(w_nodes_[i+1] - w_nodes_[i]);
    }
    dz_cells_old_ = delz_prev;
}

void IceThermo::RecalculateTempNodes(const vector<double>& T_cells_array,
                                           vector<double>& T_nodes_array)
{
    // T_0
    T_nodes_array[0] = a_nodes_[0]*T_cells_array[0] + b_nodes_[0]*T_cells_array[1]; 

    // T_i
    for (unsigned int i = 1; i < (n_nodes_-1); ++i)
    {
        T_nodes_array[i] = a_nodes_[i]*T_cells_array[i-1] + b_nodes_[i]*T_cells_array[i]; 
    }

    // T_N
    T_nodes_array[n_nodes_-1] = a_nodes_[n_nodes_-1]*T_cells_array[n_nodes_-2] + b_nodes_[n_nodes_-1]*T_cells_array[n_nodes_-1];
}

void IceThermo::RecalculateSalNodes(const vector<double>& S_cells_array,
                                          vector<double>& S_nodes_array)
{
    // S_0
    S_nodes_array[0] = a_nodes_[0]*S_cells_array[0] + b_nodes_[0]*S_cells_array[1]; 

    // S_i
    for (unsigned int i = 1; i < (n_nodes_-1); ++i)
    {
        S_nodes_array[i] = a_nodes_[i]*S_cells_array[i-1] + b_nodes_[i]*S_cells_array[i]; 
    }

    // S_N
    S_nodes_array[n_nodes_-1] = a_nodes_[n_nodes_-1]*S_cells_array[n_nodes_-2] + b_nodes_[n_nodes_-1]*S_cells_array[n_nodes_-1];
}

void IceThermo::UpdateC_Nodes(const vector<double>& T_nodes_array,
                              const vector<double>& T_old_nodes_array,
                              const vector<double>& Sal_nodes_array,
                                    vector<double>& C_nodes_array)
{
    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        C_nodes_array[i] = CalculateC(T_nodes_array[i], T_old_nodes_array[i], Sal_nodes_array[i]);
    }
}

void IceThermo::UpdateC_Cells(const vector<double>& T_cells_array,
                              const vector<double>& T_old_cells_array,
                              const vector<double>& Sal_cells_array,
                                    vector<double>& C_cells_array)
{
    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        C_cells_array[i] = CalculateC(T_cells_array[i], T_old_cells_array[i], Sal_cells_array[i]);
    }
}

void IceThermo::UpdateE_Nodes(const vector<double>& T_nodes_array,
                              const vector<double>& Sal_nodes_array,
                                    vector<double>& E_nodes_array)
{
    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        E_nodes_array[i] = CalculateE(T_nodes_array[i], Sal_nodes_array[i]);
    }
}

void IceThermo::UpdateE_Cells(const vector<double>& T_cells_array,
                              const vector<double>& Sal_cells_array,
                                    vector<double>& E_cells_array)
{
    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        E_cells_array[i] = CalculateE(T_cells_array[i], Sal_cells_array[i]);
    }
}

void IceThermo::UpdateK_Nodes(const vector<double>& T_nodes_array,
                              const vector<double>& Sal_nodes_array,
                                    vector<double>& K_nodes_array)
{
    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        K_nodes_array[i] = CalculateK(T_nodes_array[i], Sal_nodes_array[i]);
    }
}


void IceThermo::UpdateT()
{
    // update T^n nodes
    RecalculateTempNodes(T_cells_old_, T_nodes_old_);

    // update E^n nodes
    UpdateE_Nodes(T_nodes_old_, S_nodes_, E_nodes_);

    // update E^n cells
    UpdateE_Cells(T_cells_old_, S_cells_, E_cells_);

    for (unsigned int n_pseudoiter = 0; n_pseudoiter < N_PSEUDOITS; ++n_pseudoiter)
    {
        // update T^p nodes
        RecalculateTempNodes(T_cells_last_, T_nodes_last_);

        // update c^p nodes
        UpdateC_Nodes(T_nodes_last_, T_nodes_old_, S_nodes_, c_nodes_);

        // update c^p cells
        UpdateC_Cells(T_cells_last_, T_cells_old_, S_cells_, c_cells_);

        // update k^p cells
        UpdateK_Nodes(T_cells_last_, S_nodes_, k_nodes_);

        // assemble A vector
        vector<double> A_vector(n_cells_ - 1);

        for (unsigned int i = 0; i < (n_cells_-2); ++i)
        {
            A_vector[i] = -consts.rho_ice*c_nodes_[i+1]*w_nodes_[i+1]*a_nodes_[i+1] - 
                          2.0*k_nodes_[i+1]/(dz_cells_new_[i] + dz_cells_new_[i+1]);
        }

        A_vector[n_cells_-2] = consts.rho_ice*c_nodes_[n_nodes_-1]*w_nodes_[n_nodes_-1]*a_nodes_[n_nodes_-1] -
                               consts.rho_ice*c_nodes_[n_nodes_-2]*w_nodes_[n_nodes_-2]*a_nodes_[n_nodes_-2] -
                               2.0*k_nodes_[n_nodes_-1]*a_nodes_[n_nodes_-1]/(dz_cells_new_[n_cells_-1]) -
                               2.0*k_nodes_[n_nodes_-2]/(dz_cells_new_[n_nodes_-2] + dz_cells_new_[n_nodes_-1]);

        // assemble B vector
        vector<double> B_vector(n_cells_);
        B_vector[0] = consts.rho_ice*c_cells_[0]*dz_cells_new_[0]/dt_ + 
                      consts.rho_ice*c_nodes_[1]*w_nodes_[1]*a_nodes_[1] -
                      consts.rho_ice*c_nodes_[0]*w_nodes_[0]*a_nodes_[0] +
                      2.0*k_nodes_[1]/(dz_cells_new_[0] + dz_cells_new_[1]) +
                      2.0*k_nodes_[0]*(1.0 - a_nodes_[0])/(dz_cells_new_[0]);
        
        for (unsigned int i = 1; i < (n_cells_-1); ++i)
        {
            B_vector[i] = consts.rho_ice*c_cells_[i-1]*dz_cells_new_[i-1]/dt_ +
                          consts.rho_ice*c_nodes_[i]*w_nodes_[i]*a_nodes_[i] -
                          consts.rho_ice*c_nodes_[i-1]*w_nodes_[i-1]*b_nodes_[i-1] +
                          2.0*k_nodes_[i]/(dz_cells_new_[i-1] + dz_cells_new_[i]) +
                          2.0*k_nodes_[i-1]/(dz_cells_new_[i-1] + dz_cells_new_[i]);
        }

        B_vector[n_cells_-1] = consts.rho_ice*c_cells_[n_cells_-1]*dz_cells_new_[n_cells_-1]/dt_ + 
                               consts.rho_ice*c_nodes_[n_nodes_-1]*w_nodes_[n_nodes_-1]*b_nodes_[n_nodes_-1] -
                               consts.rho_ice*c_nodes_[n_nodes_-2]*w_nodes_[n_nodes_-2]*b_nodes_[n_nodes_-2] -
                               2.0*k_nodes_[n_nodes_-1]*(b_nodes_[n_nodes_-1])/(dz_cells_new_[n_cells_-1]) +
                               2.0*k_nodes_[n_nodes_-2]/(dz_cells_new_[n_cells_-2] + dz_cells_new_[n_cells_-1]);
        
        // assemble C vector
        vector<double> C_vector(n_cells_ - 1);

        C_vector[0] = consts.rho_ice*c_nodes_[1]*w_nodes_[1]*b_nodes_[1] -
                      consts.rho_ice*c_nodes_[0]*w_nodes_[0]*b_nodes_[0] -
                      2.0*k_nodes_[1]/(dz_cells_new_[0] + dz_cells_new_[1]) -
                      2.0*k_nodes_[0]*b_nodes_[0]/dz_cells_new_[0];
        
        for (unsigned int i = 1; i < (n_cells_ - 1); ++i)
        {
            C_vector[i] = consts.rho_ice*c_nodes_[i]*w_nodes_[i]*b_nodes_[i] - 
                          2.0*k_nodes_[i]/(dz_cells_new_[i-1] + dz_cells_new_[i]);
        }

        // assemble D vector
        vector<double> D_vector(n_cells_);

        for (unsigned int i = 0; i < n_cells_; ++i)
        {
            D_vector[i] = consts.rho_ice*c_cells_[i]*dz_cells_new_[i]*T_cells_old_[i]/dt_ +
                          consts.rho_ice*(c_nodes_[i+1]*w_nodes_[i+1]*T_nodes_old_[i+1] - c_nodes_[i]*w_nodes_[i]*T_nodes_old_[i]) -
                          consts.rho_ice*E_cells_[i]*(dz_cells_new_[i] - dz_cells_old_[i])/dt_ -
                          consts.rho_ice*(E_nodes_[i+1]*w_nodes_[i+1] - E_nodes_[i]*w_nodes_[i]) + 
                          dR_cells_[i]; 
        }

        // solve linear system
        vector<double> new_T_cells = thomas_solver({A_vector, B_vector, C_vector}, D_vector);

        // update T_cells_new_
        for (unsigned int i = 0; i < n_cells_; ++i)
        {
            T_cells_new_[i] = new_T_cells[i];
        }

        // calculate T diff norm
        cout << "preudoiter " << n_pseudoiter 
             << " -> Rel C norm = " << vec_C_norm(T_cells_new_ - T_cells_last_)/vec_C_norm(T_cells_old_) << ";"
             << " Rel L2 norm = " << vec_L2_norm(T_cells_new_ - T_cells_last_)/vec_L2_norm(T_cells_old_) << ";" << endl;

        // update T_cells_last_
        for (unsigned int i = 0; i < n_cells_; ++i)
        {
            T_cells_last_[i] = T_cells_new_[i];
        }
    }

    // update T_cells_old_
    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        T_cells_old_[i] = T_cells_new_[i];
    }

    // update T_nodes_old_
    RecalculateTempNodes(T_cells_old_, T_nodes_old_);

}

void IceThermo::Evaluate()
{
    UpdateW();
    UpdateDeltaZ();
    RecalculateCoeffs();
    UpdateT();
}