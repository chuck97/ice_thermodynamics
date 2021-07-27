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
        dz_cells_new_[i] = dz_cells_old_[i];
        T_cells_old_[i] = T_b + ((i + 0.5)/n_cells_)*(T_su - T_b);
        T_cells_new_[i] = T_cells_old_[i];
        T_cells_last_[i] = T_cells_old_[i];
    }

    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        S_cells_[i] = consts.sal_ice;
    }

    RecalculateCoeffs(dz_cells_new_);
    UpdateW(T_cells_new_, S_cells_, w_nodes_);
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

    RecalculateCoeffs(dz_cells_new_);
    UpdateW(T_cells_new_, S_cells_, w_nodes_);
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
    t_m.AssignCellThickness(dz_cells_new_);

    t_m.AssignCellData("T cells", K_to_C(T_cells_new_));
    t_m.AssignCellData("dz cells", dz_cells_new_);
    t_m.AssignCellData("c cells", c_cells_);
    t_m.AssignCellData("E cells", E_cells_);
    t_m.AssignCellData("dR cells", dR_cells_);
    t_m.AssignCellData("S cells", S_cells_);

    t_m.AssignNodeData("T nodes", K_to_C(T_nodes_new_));
    t_m.AssignNodeData("omega nodes", w_nodes_);
    t_m.AssignNodeData("c nodes", c_nodes_);
    t_m.AssignNodeData("k nodes", k_nodes_);
    t_m.AssignNodeData("E nodes", E_nodes_);
    t_m.AssignNodeData("S nodes", S_nodes_);
    
    stringstream ss;
    ss << setfill('0') << setw(5) << step_num;

	t_m.PlotVertical(path + "vert" + ss.str() + ".vtu");
	t_m.PlotSigma(path + "sigm" + ss.str() + ".vtu");
}

void IceThermo::RecalculateCoeffs(const std::vector<double>& dz_array)
{
    if (!IS_COEFFS_SIGMA)
    {
        // a0, b0
        a_nodes_[0] = (dz_array[1] + 2.0*dz_array[0])/(dz_array[0] + dz_array[1]);
        b_nodes_[0] = (-dz_array[0])/(dz_array[0] + dz_array[1]);

        // ai, bi
        for (unsigned int i = 1; i < (n_nodes_ - 1); ++i)
        {
            a_nodes_[i] = dz_array[i]/(dz_array[i] + dz_array[i-1]);
            b_nodes_[i] = dz_array[i-1]/(dz_array[i] + dz_array[i-1]);
        }

        // aN, bN
        a_nodes_[n_nodes_ - 1] = (-dz_array[n_nodes_ - 1])/(dz_array[n_nodes_ - 1] + dz_array[n_nodes_ - 2]);
        b_nodes_[n_nodes_ - 1] = (dz_array[n_nodes_ - 2] + 2.0*dz_array[n_nodes_ - 1])/(dz_array[n_nodes_ - 1] + dz_array[n_nodes_ - 2]);
    }
    else
    {
        // calculate total thickness
        double h = 0.0;
        for (unsigned int i = 0; i < n_cells_; ++i)
        {
            h += dz_array[i];
        }

        // a0, b0
        a_nodes_[0] = 0.5 + (2.0*h/n_cells_)/(dz_array[0] + dz_array[1]);
        b_nodes_[0] = 0.5 - (2.0*h/n_cells_)/(dz_array[0] + dz_array[1]);

        // ai, bi
        for (unsigned int i = 0; i < (n_nodes_ - 1); ++i)
        {
            a_nodes_[i] = 0.5;
            b_nodes_[i] = 0.5;
        }

        // aN, bN
        a_nodes_[n_nodes_ - 1] = 0.5 - (2.0*h/n_cells_)/(dz_array[n_nodes_ - 1] + dz_array[n_nodes_ - 2]);
        b_nodes_[n_nodes_ - 1] = 0.5 + (2.0*h/n_cells_)/(dz_array[n_nodes_ - 1] + dz_array[n_nodes_ - 2]);;

    }
}

double IceThermo::CalculateE(double T, double S)
{
    double Tf = -consts.mu_ice*S;
    return (consts.c0_ice*(T - Tf) - consts.L0_ice*(1.0 - Tf/T) + consts.cw_water*Tf);
}

double IceThermo::CalculateK(double T, double S)
{
    return (consts.k0_ice + consts.k1_ice*(S/T));
}

double IceThermo::CalculateC(double T_last, double T_old, double S)
{
    double Tf = -consts.mu_ice*S;
    return (consts.c0_ice - consts.L0_ice*Tf/(T_last*T_old));
}

void IceThermo::UpdateW(const vector<double>& T_cells_array,
                        const vector<double>& S_cells_array,
                              vector<double>& w_nodes_array)
{
    // bottom w
    double dT_dz_b = (2.0*(1.0 - a_nodes_[0])/dz_cells_new_[0])*T_cells_array[0] +
                     (-2.0*b_nodes_[0]/dz_cells_new_[0])*T_cells_array[1];
    
    double T_b = a_nodes_[0]*T_cells_array[0] + b_nodes_[0]*T_cells_array[1];
    double S_b = a_nodes_[0]*S_cells_array[0] + b_nodes_[0]*S_cells_array[1];
    double E_b = CalculateE(T_b, S_b);
    double k_b = CalculateK(T_b, S_b);
    double w_b = (F_b_ - k_b*dT_dz_b)/(consts.rho_ice*E_b);

    // surface w
    double dT_dz_s = ((2.0*a_nodes_[n_nodes_-1])/dz_cells_new_[n_cells_-1])*T_cells_array[n_cells_-2] +
                     (2.0*(b_nodes_[n_nodes_-1] - 1.0)/dz_cells_new_[n_cells_-1])*T_cells_array[n_cells_-1];
    
    double T_s = a_nodes_[n_nodes_-1]*T_cells_array[n_cells_-2] + b_nodes_[n_nodes_-1]*T_cells_array[n_cells_-1];
    double S_s = a_nodes_[n_nodes_-1]*S_cells_array[n_cells_-2] + b_nodes_[n_nodes_-1]*S_cells_array[n_cells_-1];
    double E_s = CalculateE(T_s, S_s);
    double k_s = CalculateK(T_s, S_s);
    double w_s = (F_su_ - k_s*dT_dz_s)/(consts.rho_ice*E_s);

    // nodes linear interpolation
    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        w_nodes_array[i] = w_b + (i/(n_nodes_-1))*(w_s - w_b);
    }   
}

void IceThermo::UpdateDeltaZ()
{
    dz_cells_old_ = dz_cells_new_;
    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        dz_cells_new_[i] = dz_cells_old_[i] - dt_*(w_nodes_[i+1] - w_nodes_[i]);
    }
}

void IceThermo::InterpolateScalarNodes(const vector<double>& scalar_cells_array,
                                             vector<double>& scalar_nodes_array)
{
    // scalar_0
    scalar_nodes_array[0] = a_nodes_[0]*scalar_cells_array[0]
                          + b_nodes_[0]*scalar_cells_array[1];

    // scalar_i
    for (unsigned int i = 1; i < (n_nodes_-1); ++i)
    {
        scalar_nodes_array[i] = a_nodes_[i]*scalar_cells_array[i-1]
                              + b_nodes_[i]*scalar_cells_array[i];
    }

    // scalar_N
    scalar_nodes_array[n_nodes_-1] = a_nodes_[n_nodes_-1]*scalar_cells_array[n_nodes_-2]
                                   + b_nodes_[n_nodes_-1]*scalar_cells_array[n_nodes_-1];
}

void IceThermo::UpdateC(const std::vector<double>& T_old_array,
                        const std::vector<double>& T_new_array,
                        const std::vector<double>& S_array,
                              std::vector<double>& C_array)
{
    if (T_old_array.size() != T_new_array.size() ||
        T_new_array.size() != S_array.size() ||
        S_array.size() != C_array.size())
    {
        ERR("Sizes of T_old, T_new, S and C arrays in UpdateC are not equal!");
    }

    for (unsigned int i = 0; i < C_array.size(); ++i)
    {
        C_array[i] = CalculateC(T_new_array[i], T_old_array[i], S_array[i]);
    }
}

void IceThermo::UpdateE(const std::vector<double>& T_array,
                        const std::vector<double>& S_array,
                              std::vector<double>& E_array)
{
    if (T_array.size() != S_array.size() || S_array.size() != E_array.size())
    {
        ERR("Sizes of S, T and E arrays in UpdateE are not equal!");
    }

    for (unsigned int i = 0; i < E_array.size(); ++i)
    {
        E_array[i] = CalculateE(T_array[i], S_array[i]);
    }
}

void IceThermo::UpdateK(const std::vector<double>& T_array,
                        const std::vector<double>& S_array,
                              std::vector<double>& K_array)
{
    if (T_array.size() != S_array.size() || S_array.size() != K_array.size())
    {
        ERR("Sizes of S, T and K arrays in UpdateK are not equal!");
    }

    for (unsigned int i = 0; i < K_array.size(); ++i)
    {
        K_array[i] = CalculateK(T_array[i], S_array[i]);
    }
}


void IceThermo::UpdateT()
{
    // update E^n cells
    UpdateE(T_cells_old_, S_cells_, E_cells_);

    // update T^n nodes
    InterpolateScalarNodes(T_cells_old_, T_nodes_old_);

    // update E^n nodes
    UpdateE(T_nodes_old_, S_nodes_, E_nodes_);

    // initialize T_last (T^p) with T_old (T^n)
    T_cells_last_ = T_cells_old_;

    for (unsigned int n_pseudoiter = 0; n_pseudoiter < N_PSEUDOITS; ++n_pseudoiter)
    {
        // update T^p nodes
        InterpolateScalarNodes(T_cells_last_, T_nodes_last_);

        // update c^p nodes
        UpdateC(T_nodes_old_, T_nodes_last_, S_nodes_, c_nodes_);

        // update c^p cells
        UpdateC(T_cells_old_, T_cells_last_, S_cells_, c_cells_);

        // update k^p cells
        UpdateK(T_nodes_last_, S_nodes_, k_nodes_);

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
                               2.0*k_nodes_[n_nodes_-2]/(dz_cells_new_[n_cells_-2] + dz_cells_new_[n_cells_-1]);

        // assemble B vector
        vector<double> B_vector(n_cells_);

        B_vector[0] = consts.rho_ice*c_cells_[0]*dz_cells_new_[0]/dt_ + 
                      consts.rho_ice*c_nodes_[1]*w_nodes_[1]*a_nodes_[1] -
                      consts.rho_ice*c_nodes_[0]*w_nodes_[0]*a_nodes_[0] +
                      2.0*k_nodes_[1]/(dz_cells_new_[0] + dz_cells_new_[1]) +
                      2.0*k_nodes_[0]*(1.0 - a_nodes_[0])/(dz_cells_new_[0]);
        
        for (unsigned int i = 1; i < (n_cells_-1); ++i)
        {
            B_vector[i] = consts.rho_ice*c_cells_[i]*dz_cells_new_[i]/dt_ +
                          consts.rho_ice*c_nodes_[i+1]*w_nodes_[i+1]*a_nodes_[i+1] -
                          consts.rho_ice*c_nodes_[i]*w_nodes_[i]*b_nodes_[i] +
                          2.0*k_nodes_[i+1]/(dz_cells_new_[i] + dz_cells_new_[i+1]) +
                          2.0*k_nodes_[i]/(dz_cells_new_[i-1] + dz_cells_new_[i]);
        }

        B_vector[n_cells_-1] = consts.rho_ice*c_cells_[n_cells_-1]*dz_cells_new_[n_cells_-1]/dt_ + 
                               consts.rho_ice*c_nodes_[n_nodes_-1]*w_nodes_[n_nodes_-1]*b_nodes_[n_nodes_-1] -
                               consts.rho_ice*c_nodes_[n_nodes_-2]*w_nodes_[n_nodes_-2]*b_nodes_[n_nodes_-2] -
                               2.0*k_nodes_[n_nodes_-1]*(b_nodes_[n_nodes_-1] - 1.0)/(dz_cells_new_[n_cells_-1]) +
                               2.0*k_nodes_[n_nodes_-2]/(dz_cells_new_[n_cells_-2] + dz_cells_new_[n_cells_-1]);
        
        // assemble C vector
        vector<double> C_vector(n_cells_ - 1);

        C_vector[0] = consts.rho_ice*c_nodes_[1]*w_nodes_[1]*b_nodes_[1] -
                      consts.rho_ice*c_nodes_[0]*w_nodes_[0]*b_nodes_[0] -
                      2.0*k_nodes_[1]/(dz_cells_new_[0] + dz_cells_new_[1]) -
                      2.0*k_nodes_[0]*b_nodes_[0]/dz_cells_new_[0];
        
        for (unsigned int i = 1; i < (n_cells_ - 1); ++i)
        {
            C_vector[i] = consts.rho_ice*c_nodes_[i+1]*w_nodes_[i+1]*b_nodes_[i+1] -
                          2.0*k_nodes_[i+1]/(dz_cells_new_[i] + dz_cells_new_[i+1]);
        }

        // assemble D vector
        vector<double> D_vector(n_cells_);

        for (unsigned int i = 0; i < n_cells_; ++i)
        {
            D_vector[i] = consts.rho_ice*c_cells_[i]*dz_cells_new_[i]*T_cells_old_[i]/dt_ +
                          consts.rho_ice*(c_nodes_[i+1]*w_nodes_[i+1]*T_nodes_old_[i+1] -
                                          c_nodes_[i]*w_nodes_[i]*T_nodes_old_[i]) -
                          consts.rho_ice*E_cells_[i]*(dz_cells_new_[i] - dz_cells_old_[i])/dt_ -
                          consts.rho_ice*(E_nodes_[i+1]*w_nodes_[i+1] - E_nodes_[i]*w_nodes_[i]) + 
                          dR_cells_[i]; 
        }

//        cout << "preudoiter " << n_pseudoiter
//             << " -> A C norm = " << vec_C_norm(A_vector) << ";"
//             << " -> B C norm = " << vec_C_norm(B_vector) << ";"
//             << " -> C C norm = " << vec_C_norm(C_vector) << ";"
//             << " -> D C norm = " << vec_C_norm(D_vector) << ";" << endl;

        // solve linear system
        vector<double> T_cells_curr = thomas_solver({A_vector, B_vector, C_vector}, D_vector);

        // calculate T diff norm
        cout << "preudoiter " << n_pseudoiter
             << " -> Rel C norm = " << vec_C_norm(T_cells_curr - T_cells_last_)/vec_C_norm(T_cells_old_) << ";"
             << " Rel L2 norm = " << vec_L2_norm(T_cells_curr - T_cells_last_)/vec_L2_norm(T_cells_old_) << ";" << endl;

        // update T_cells_last_
        T_cells_last_ = T_cells_curr;
    }

    T_cells_new_ = T_cells_last_;
}

void IceThermo::Evaluate()
{
    // new variables becomes old
    dz_cells_old_ = dz_cells_new_;
    T_cells_old_ = T_cells_new_;

    UpdateDeltaZ();
    RecalculateCoeffs(dz_cells_new_);
    UpdateT();
    UpdateW(T_cells_new_, S_cells_, w_nodes_);
}

double C_to_K(double T)
{
    return T + 273.15;
}

double K_to_C(double T)
{
    return T - 273.15;
}

vector<double> C_to_K(const std::vector<double>& T_C_array)
{
    vector<double> T_K_array;
    for (double T_C: T_C_array)
    {
        T_K_array.push_back(C_to_K(T_C));
    }

    return T_K_array;
}

vector<double> K_to_C(const std::vector<double>& T_K_array)
{
    vector<double> T_C_array;
    for (double T_K: T_K_array)
    {
        T_C_array.push_back(K_to_C(T_K));
    }

    return T_C_array;
}
