#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include "Basics.h"
#include "multithread_loop.h"
#include "Auxillaries.h"
#include "SolveModel.h"

// #include "Interpolation_spline.h"
//
// ////
// //
namespace alias {

    /**************************************************************
    * Value function in the initial period
    **************************************************************/
    //// Equilibrium value/policy functions for period t == 0
    struct EquStateV0 {
        ArrayXd EVal_Lur0_error;  // dim = N_phi * N_K
        ArrayXd EVal_Lur0;  // dim = N_phi * N_K
        ArrayXd EVpart_Lur0; // dim = N_phi * N_K
        ArrayXd OptLur0;  // dim = N_phi * N_K

        ArrayXd EVal0; // dim = N_phi
        ArrayXd OptK0; // dim = N_phi
    };

    /**************************************************************
    * Solve the value/policy function in the first period
    **************************************************************/
    EquStateV0 solveV0(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & Vprime,
        MultiThreads::Threads_Management & threadsManagement);

    /**************************************************************
    * Solve policy function and value function at t = 0
    **************************************************************/
    //// At period 0, choose the optimal Lur
    tuple<ArrayXd,ArrayXd,ArrayXd> solveOptLur0(const ParaEst & para_est, const ParaVec & para_vec,
        const ArrayXd & EVal_PD, const ArrayXXd & RevOpt_Prod_mat,
        MultiThreads::Threads_Management & threadsManagement);
    tuple<double,double,double,int> solve1D_L0_LinearSpline_RevOpt(const double & wstar,
        const double & phi, const double & K, const ParaEst & para_est,
        const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const ArrayXd & RevOpt_Prod_vec);
    tuple<double, double, double> SolveMax_Linear1D_L0_RevOpt(const ArrayXd & coef, const double & dRevOpt_dL_coef,
            const double & x0,const double & x1, const double & TotV0,const double & TotV1,
            const double & wstar, const double & phi, const double & K, const ParaEst & para_est);

    //// At period 0, choose the optimal K
    tuple<ArrayXd,ArrayXd> solveOptK0(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & EVpart_Lur0,
        const ArrayXd & OptLur0, MultiThreads::Threads_Management & threadsManagement);
    tuple<double,double,int> solve1D_K_LinearSpline(const ArrayXd & vec_K, const ArrayXd & EVpart_Lur0_vec,
        const ArrayXd & OptLur0_vec, const ParaEst & para_est, const double & phi);
    tuple<double, double> SolveMax_Linear1D_K0_RevOpt(const ArrayXd & coef, const double & x0,const double & x1,
        const double & TotV0, const double & TotV1, const double & EVpart_Lur0_0, const double & EVpart_Lur0_1,
        const double & OptLur0, const double & OptLur1, const double & phi, const ParaEst & para_est);

    /**************************************************************
    * Solve the entry cost
    **************************************************************/
    //// The equilibrium state of the initial period
    struct EquState0 {
        double F_Entry;
        double FirmMass;
    };
    double SolveFEntry_StatusQuoEqu(const ParaVec &para_vec, const EquStateV0 &EquV0);
    //

}
