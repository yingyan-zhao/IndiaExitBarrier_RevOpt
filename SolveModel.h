#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include "Basics.h"
#include "multithread_loop.h"
#include "Auxillaries.h"

// #include "Interpolation_spline.h"
//
// ////
// //
namespace alias
{
    /**************************************************************
    * Data Structure to store the value function and policy functions
    **************************************************************/
    //// Equilibrium value/policy functions for period t >= 1
    struct EquStateV {
        /// value function
        ArrayXd EVal_PD_Lur_error;
        ArrayXd EVal_PD_OptLurP; // the value function after deciding labor employment given the firm decides to produce
        ArrayXd EVal_PD_OptLurD; // the value function after deciding labor employment given the firm decides to stay dormant
        ArrayXd EVal_PD_ProdDorm; // value function after deciding production or dormancy
        ArrayXd EVal_PD; // the value function at the beginning of each period

        ArrayXd OptLur_P;
        ArrayXd OptLur_D;
        //
        /// policy function
        ArrayXd ProbPD_P;
        ArrayXd ProbPD_D;

        ArrayXd ProbPD_S;
        ArrayXd ProbPD_E;

        //
        // ArrayXXd EVal_PD_Lur_mat;
        // ArrayXXd EVal_PD_SE_mat;
    };

    /**************************************************************
    * solve value function for period >=1
    **************************************************************/
    EquStateV solveV(const ParaEst & para_est, const ParaVec & para_vec, MultiThreads::Threads_Management & threadsManagement);

    /**************************************************************
    * Some useful functions
    **************************************************************/
    ////// Calculate the residual value of firms if exiting
    ArrayXd calResidual_Value_Exit(const ParaEst & para_est, const ParaVec & para_vec);
    ////// define hiring and firing cost function form
    tuple<double,double> CalHFcost_Lbase(const double & c_H_ur, const double & c_high_F_ur, const double & c_low_F_ur,
        const double & Lbase);
    //// Calculate the flow value added for every grids in the space
    ArrayXXd RevOpt_Prod_mat(const ParaEst & para_est, const ParaVec & para_vec);
    ArrayXXd RevOpt_Prod_vec(const ParaEst & para_est, const ParaVec & para_vec);
    tuple<double,double,double> RevOpt_Prod_ELur(const ParaEst & para_est, const double & phi, const double & K, const double & Lur);
    tuple<double,double,double> RevOpt_Prod_K(const ParaEst & para_est, const double & phi, const double & K, const double & Lur);

    /**************************************************************
    * Within each loop: we solve it by backward induction
    **************************************************************/
    //// the expected value function with the productivity transition
    ArrayXd CalEVprime_mat(const int & PD, const ArrayXd & Vprime, const ArrayXXd & tran_phi);

    //// Step 3. solve L_ur
    /* With a chosen target, firms' expected value depends on the shocks on employment */
    ArrayXd CalEVal_Lur_error(const int & PD, const ParaVec & para_vec, const ArrayXd & EVal, const double & sigmaL_ur);
    /* Choose the targeted Lur to maximize the expected value */
    tuple<ArrayXd,ArrayXd,ArrayXd,ArrayXd> solveOptLur(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
        const ArrayXd & EVal_PD, const ArrayXXd & RevOpt_Prod_mat, MultiThreads::Threads_Management & threadsManagement);
    //// data structure to save position
    struct Pos {
        size_t i_phi;
        size_t i_K;
        size_t i_Lur;
    };
    tuple<double,double,int> solve1D_L_LinearSpline(const double & H, const double & F, const double & wstar,
        const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec);
    tuple<double, double> SolveMax_Linear1D_L(const ArrayXd & coef,const double & x0,const double & x1,
        const double & TotV0,const double & TotV1, const double & H, const double & F, const double & wstar,
        const double & Lbase);
    tuple<double, double> SolveMax_Linear1D_L_part(const ArrayXd & coef,const double & Adj,
        const double & wstar, const double & Lbase);

    tuple<double,double,int> solve1D_L_LinearSpline_RevOpt(const double & H, const double & F, const double & wstar,
        const double & Lbase, const double & phi, const double & K, const ParaEst & para_est,
        const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const ArrayXd & RevOpt_Prod_vec);
    tuple<double, double> SolveMax_Linear1D_L_RevOpt(const ArrayXd & coef, const double & dRevOpt_dL_coef,
        const double & x0,const double & x1, const double & TotV0,const double & TotV1, const double & H, const double & F,
        const double & wstar, const double & Lbase, const double & phi, const double & K, const ParaEst & para_est);

    //// Step 2.  Firms decide if produce or be dormant
    tuple<ArrayXd,ArrayXd,ArrayXd> calEV_ProbStatus_PD(const ParaEst & para_est,
        const ArrayXd & Vprime_LurP, const ArrayXd & Vprime_LurD, MultiThreads::Threads_Management & threadsManagement);

    //// Entry and exit (Step 1)
    tuple<ArrayXd,ArrayXd,ArrayXd> calEV_ProbStatus_SE_logit(const ParaEst & para_est,
        const ArrayXd & ResVal_E, const ArrayXd & EVal, MultiThreads::Threads_Management & threadsManagement);

}
// ////    using namespace std;
// ////    using namespace alias;
// //////
// ////    using ArrayXd = Eigen::ArrayXd;
// ////    using ArrayXXd = Eigen::ArrayXXd;
// ////    using ArrayXi = Eigen::ArrayXi;
// ////    using namespace Ipopt_Wrapper;
// ////

//

//     ////// Calculate the flow value added for every grids in the space
//     struct Rev {
//         ArrayXd RevOpt_Prod;
//         ArrayXXd RevOpt_Prod_mat;
//     };
//     Rev RevOpt_Prod_vec(const ParaEst & para_est, const ParaVec & para_vec);
//

//
//     /**************************************************************
//     * Per loop, update the value function
//     **************************************************************/
//     EquStateV solveV_OneLoop(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & Vprime_ini,
//         const ArrayXd & ResVal_E, MultiThreads::Threads_Management & threadsManagement);
//


//     //// Find the optimal solution (Lur) in Step 3
//     tuple<double,double,int> solve1D_L_LinearSpline(const double & H, const double & F, const double & wstar,
//         const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec);
//     double solve1D_L_LinearSpline_Diff(const double & H, const double & F, const double & wstar,
//         const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const int & L_index_max);
//

//
//
//
//

//     // //// Find the optimal solution (K)
//     // tuple<ArrayXd,ArrayXd> alias::solveOptK(const int & PD, const ParaEst & para_est,const ParaVec & para_vec,
//     //     const ArrayXd & EVal_PD_K, MultiThreads::Threads_Management & threadsManagement);
//
// //
// //
// ////    /**************************************************************
// ////     * Calculate the value function of the end period:
// ////     * Approximation: Assume that firms stop change capital or employment after 40 period. This is for computation only
// ////    **************************************************************/
// ////    ArrayXd calVEnd(const ParaEst & para_est, const ParaVec & para_vec, const Rev & RevOpt);
// ////
// //
// ////
// //
// ////    //// Firms decide if produce or be dormant
// ////    tuple<ArrayXd,ArrayXd,ArrayXd> calEV_ProbStatus_PD(const ParaEst & para_est, const ParaVec & para_vec,
// ////        const ArrayXd & Vprime, const Rev & RevOpt, MultiThreads::Threads_Management & threadsManagement);
// ////
// ////    //// solve Luc (Step 4)
// ////    /* With a chosen target, firms' expected value depends on the shocks on employment */
// ////    ArrayXd CalEVal_Luc_error(const int & PD, const ParaVec & para_vec, const ArrayXd & EVal, const double & sigmaL_uc,
// ////        MultiThreads::Threads_Management & threadsManagement);
// ////    /* Choose the targeted Luc to maximize the expected value */
// ////    tuple<ArrayXd,ArrayXd> solveOptLuc(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
// ////        const ArrayXd & EVal_Luc, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement);
// ////    /* solve the optimal Luc by iteration */
// ////    tuple<ArrayXd,ArrayXd> solveOptLuc_iteration(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
// ////        const ArrayXd & EVal_Luc, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement);
// ////
// ////    //// solve L_ur (Step 3)
// ////    /* With a chosen target, firms' expected value depends on the shocks on employment */
// ////    ArrayXd CalEVal_Lur_error(const int & PD, const ParaVec & para_vec, const ArrayXd & EVal, const double & sigmaL_ur,
// ////        MultiThreads::Threads_Management & threadsManagement);
// ////    /* Choose the targeted Lur to maximize the expected value */
// ////    tuple<ArrayXd,ArrayXd> solveOptLur(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
// ////        const ArrayXd & EVal_Lur, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement);
// ////    /* solve the optimal Luc by iteration */
// ////    tuple<ArrayXd,ArrayXd> solveOptLur_iteration(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
// ////        const ArrayXd & EVal_Lur, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement);
// ////
// ////    //// Find the optimal solution (Luc & Lur) in Step 3 and Step 4
// ////    tuple<double,double,int> solve1D_L_LinearSpline(const double & H, const double & F, const double & wstar,
// ////        const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec);
// ////
// ////    //// solve K (Step 2)
// ////    /* With a chosen target, firms' expected value depends on the shocks on capital */
// ////    ArrayXd CalEVal_K_error(const int & PD, const ParaVec & para_vec, const ArrayXd & EVal, const double & sigmaK,
// ////        MultiThreads::Threads_Management & threadsManagement);
// ////    /* Choose the targeted Luc to maximize the expected value */
// ////    tuple<ArrayXd,ArrayXd> solveOptK(const int & PD, const ParaEst & para_est,const ParaVec & para_vec,
// ////        const ArrayXd & EVal_K, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement);
// ////    tuple<ArrayXd,ArrayXd> solveOptK_iteration(const int & PD, const ParaEst & para_est,const ParaVec & para_vec,
// ////        const ArrayXd & EVal_K, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement);
// ////
// ////    //// Find the optimal solution (K) in Step 2
// ////    tuple<double,double,int> solve1D_K_LinearSpline(const double & H, const double & F, const double & c_K,
// ////        const double & rstar, const double & p_K, const double & Kbase, const ArrayXd & vec_K, const ArrayXd & EVal_K_vec);
// ////
// //
// ////    /**************************************************************
// ////    * Testing the concacvity of value functions
// ////    **************************************************************/
// ////    void testConcavity_L(const ArrayXd & Ymat, const ArrayXd & vec_Luc, const ArrayXd & vec_Lur);
// ////    void testConcavity_K(const ArrayXXd & Y, const ArrayXd & vec_K);
// ////    void testConcavity_L_loop(const ArrayXd & Eval, const ArrayXd & vec_Luc, const ArrayXd & vec_Lur);
// ////    void testConcavity_K_loop(const ArrayXd & Eval, const ArrayXd & vec_K);
// ////
// ////    void testMonotonicity_by_K_loop(const ArrayXd & OptKL);
// ////    void testMonotonicity_by_Lur_loop(const ArrayXd & OptKL);
// ////    void testMonotonicity_by_Luc_loop(const ArrayXd & OptKL);
// ////
// //////    tuple<ArrayXd,ArrayXd,ArrayXd> calEV_ProbStatus_SE_lognormal(const ParaEst & para_est, const ParaVec & para_vec,
// //////        const ArrayXd & ResVal_E, const ArrayXd & EVal, MultiThreads::Threads_Management & threadsManagement);
// //////
// }
