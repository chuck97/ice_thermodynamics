#include "icethermo.hpp"
#include <fstream>
#include <vector>
#include <iostream>
#include <string>

using namespace icethermo;

template<typename NumType>
FourVecs<NumType> Assemble_advdiff_martix_rhs(double time_step,
                                              const std::vector<NumType>& T_cells_prev,
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
                                              const std::vector<NumType>& radiation_nodes)
{
    // get all parameterizations
    Kparam kparam = Kparam::BubblyBrine;
    Cparam cparam = Cparam::SeaIce;
    Eparam Eparam = Eparam::SeaIce;

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
    B[0] = rho_cells[0]*eff_c_cells[0]*dz_cells_new[0]/time_step + eff_k_nodes[1] + eff_k_nodes[0];

    C[0] = -eff_k_nodes[1];

    RHS[0] = rho_cells[0]*eff_c_cells[0]*dz_cells_new[0]*T_cells_old[0]/time_step - 
              rho_cells[0]*E_cells[0]*(dz_cells_new[0] - dz_cells_old[0])/time_step +
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

        B[i] = rho_cells[i]*eff_c_cells[i]*dz_cells_new[i]/time_step + eff_k_nodes[i+1] + eff_k_nodes[i];
            
        C[i] =  -eff_k_nodes[i+1];
            
        RHS[i] = rho_cells[i]*eff_c_cells[i]*dz_cells_new[i]*T_cells_old[i]/time_step -
                 rho_cells[i]*E_cells[i]*(dz_cells_new[i] - dz_cells_old[i])/time_step +
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

    B[N-1] = rho_cells[N-1]*eff_c_cells[N-1]*dz_cells_new[N-1]/time_step + eff_k_nodes[N] + eff_k_nodes[N-1];

    RHS[N-1] = rho_cells[N-1]*eff_c_cells[N-1]*T_cells_old[N-1]*dz_cells_new[N-1]/time_step - 
               rho_cells[N-1]*E_cells[N-1]*(dz_cells_new[N-1] - dz_cells_old[N-1])/time_step -
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


std::vector<double> ReadVecFromFile(const::std::string& filename)
{
    std::fstream s(filename);

    if (!s.is_open()) 
    {
        std::cout << "failed to open " << filename << '\n';
        exit;
    }

    std::vector<double> vec;

    while(true)
    {
        double item;
        s >> item;
        if (s.eof())
            break;
        vec.push_back(item);
    }
    return vec;
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec)
{
    for (int i = 0; i < vec.size(); ++i)
    {
        os << vec[i] << " ";
    }
    return os;
} 

int main()
{
    // read reference matrix
    std::vector<double> left_diag_ref = ReadVecFromFile("../../../tests/matrix_test/matrix_left_diag.txt");
    std::vector<double> cent_diag_ref = ReadVecFromFile("../../../tests/matrix_test/matrix_main_diag.txt");
    std::vector<double> right_diag_ref = ReadVecFromFile("../../../tests/matrix_test/matrix_right_diag.txt");
    std::vector<double> rhs_ref = ReadVecFromFile("../../../tests/matrix_test/matrix_rhs.txt");

    // compute Cpp matrix
    double rho_i = IceConsts<double>::rho_i;
    auto res = Assemble_advdiff_martix_rhs<double>(3600.0,
                                                   std::vector<double> {-1.55, -5.0, -9.0, -13.0, -17.0},
                                                   std::vector<double> {-1.50, -4.95, -8.95, -12.95, -16.95},
                                                   -20.0, 
                                                   -21.0, 
                                                   -22.0,
                                                   -1.5, 
                                                   -1.4, 
                                                   -1.3, 
                                                   -0.00001, 
                                                   -0.0002, 
                                                   std::vector<double> {0.1, 0.1, 0.1, 0.1, 0.1}, 
                                                   std::vector<double> {0.98, 0.98, 0.98, 0.98, 0.98},
                                                   std::vector<double> {0.0, 0.0, 0.0, 0.0, 0.0},
                                                   std::vector<double> {rho_i, rho_i, rho_i, rho_i, rho_i},
                                                   std::vector<double> {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    
    std::cout << "Python vs Cpp left diad:" << std::endl;
    std::cout << left_diag_ref << std::endl;
    std::cout << std::get<0>(res) << std::endl;

    std::cout << "Python vs Cpp centr diad:" << std::endl;
    std::cout << cent_diag_ref << std::endl;
    std::cout << std::get<1>(res) << std::endl;

    std::cout << "Python vs Cpp right diad:" << std::endl;
    std::cout << right_diag_ref << std::endl;
    std::cout << std::get<2>(res) << std::endl;

    std::cout << "Python vs Cpp rhs:" << std::endl;
    std::cout << rhs_ref << std::endl;
    std::cout << std::get<3>(res) << std::endl;
}