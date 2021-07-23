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
    }
    dz_cells_new_ = dz_cells_old_;

    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        S_nodes_[i] = consts.sal_ice;
    }

    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        S_cells_[i] = consts.sal_ice;
    }

    UpdateCoeffs();
}

IceThermo::IceThermo(unsigned int n_cells,
                     const std::vector<double>& T_init_cells,
                     const std::vector<double>& dz_init_cells):
    n_cells_(n_cells),
    n_nodes_(n_cells + 1)
{
    T_cells_old_.resize(n_cells_);
    T_cells_last_.resize(n_cells_);
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
    T_cells_old_ = T_init_cells;
    dz_cells_new_ = dz_cells_old_;

    for (unsigned int i = 0; i < n_nodes_; ++i)
    {
        S_nodes_[i] = consts.sal_ice;
    }

    for (unsigned int i = 0; i < n_cells_; ++i)
    {
        S_cells_[i] = consts.sal_ice;
    }

    UpdateCoeffs();
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

void IceThermo::WriteVTU(const std::string& path,
                         unsigned int step_num) const
{
    ThermoMesh t_m(n_cells_);
	t_m.SetWidth(0.2);
	t_m.AssignCellThickness(dz_cells_old_);

	t_m.AssignCellData("T cells", T_cells_old_);
    t_m.AssignCellData("dz cells", dz_cells_old_);
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

void IceThermo::UpdateCoeffs()
{
    if (!IS_COEFFS_SIGMA)
    {
        // a0, b0
        a_nodes_[0] = (dz_cells_new_[1] + 2.0*dz_cells_new_[0])/(dz_cells_new_[0] + dz_cells_new_[1]);
        b_nodes_[0] = (-dz_cells_new_[0])/(dz_cells_new_[0] + dz_cells_new_[1]);

        // aN, bN
        a_nodes_[n_nodes_ - 1] = (-dz_cells_new_[n_nodes_ - 1])/(dz_cells_new_[n_nodes_ - 1] + dz_cells_new_[n_nodes_ - 2]);
        b_nodes_[n_nodes_ - 1] = (dz_cells_new_[n_nodes_ - 2] + 2.0*dz_cells_new_[n_nodes_ - 1])/(dz_cells_new_[n_nodes_ - 1] + dz_cells_new_[n_nodes_ - 2]);

        // ai, bi
        for (unsigned int i = 0; i < (n_cells_ - 1); ++i)
        {
            a_nodes_[i] = dz_cells_new_[i]/(dz_cells_new_[i] + dz_cells_new_[i-1]);
            b_nodes_[i] = dz_cells_new_[i-1]/(dz_cells_new_[i] + dz_cells_new_[i-1]);
        }
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
    return (consts.c0_ice*(T - Tf) - consts.L0_ice*(1.0 - Tf/T) + consts.cw_water*Tf);
}

double IceThermo::CalculateK(double T, double S)
{
    return (consts.k0_ice + consts.k1_ice*(S/T));
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

void IceThermo::UpdateT()
{
    // to do
}

void IceThermo::Evaluate()
{
    UpdateW();
    UpdateDeltaZ();
    UpdateCoeffs();
    UpdateT();
}
