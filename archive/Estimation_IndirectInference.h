#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <map>
#include "multithread_loop.h"

#include "Basics.h"
#include "SimulationData.h"
#include "EstimationAuxiliary.h"
//
//#include "DFS.h"
//#include "direct_search.h"
//#include "soo.h"
//
//
//
//
////#include "Simulation.h"
//////
//
////#include "SolveModel_V0_FEntry.h"
////#include "SolveModel.h"
////
//
//
namespace alias {
    using namespace std;
    using namespace alias;
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;

    /**************************************************************************
    * In the second stage of estimate, we first simulate data with a guessed parameter.
    * If we guessed it right, theta_a should be the maximizer of the auxilliary model with the guessed parameter.
    * That is, the derivative of the maximum likelihood function of the auxiliary model with simulated data
    * at theta_a should be as close to zero as possible.
    **************************************************************************/

    /**************************************************************
    * Estimation Part 2: Second Stage model
    **************************************************************/
    ArrayXd EstimationIndirectInference_Output(const ArrayXd & theta1_a, const ArrayXd & theta2_a, const ArrayXd & theta3_a,
        MultiThreads::Threads_Management & threadsManagement);

    struct EquV_Auxiliary {
        ParaEst para_est_a_point;
        std::vector<ParaEst> para_est_a_plus = std::vector<ParaEst>(para.dim3);

        ParaVec para_vec_a_point;
        std::vector<ParaVec> para_vec_a_plus = std::vector<ParaVec>(para.dim3);
        std::vector<ParaVec> para_vec_a_minus = std::vector<ParaVec>(para.dim3);

        EquStateV EquV_a_point;
        std::vector<EquStateV> EquV_a_plus = std::vector<EquStateV>(para.dim3);
        std::vector<EquStateV> EquV_a_minus = std::vector<EquStateV>(para.dim3);

        EquStateVmat Evalmat_a_point;
        std::vector<EquStateVmat> Evalmat_a_plus = std::vector<EquStateVmat>(para.dim3);
        std::vector<EquStateVmat> Evalmat_a_minus = std::vector<EquStateVmat>(para.dim3);

        ArrayXd eps_diff = ArrayXd::Zero(para.dim3);
    };

    EquV_Auxiliary CalEquV_plus_minus(const ArrayXd &theta_a, const int &GoodState, const int &LaborIntensive,
        MultiThreads::Threads_Management &threadsManagement);

    /**************************************************************************
    * Estimation in the second stage
    **************************************************************************/
    ArrayXd EstimationIndirectInference(const SimVar & sim_var_2step_GoodState_LaborInt,
        const SimVar & sim_var_2step_BadState_LaborInt,const SimVar & sim_var_2step_GoodState_CapitalInt,
        const SimVar & sim_var_2step_BadState_CapitalInt,
        const ArrayXd & theta1_Step2_ini, const ArrayXd & theta2_Step2_ini, const ArrayXd & theta3_Step2_ini,
        const ArrayXd & theta_a,
        const EquV_Auxiliary & EquV_a_GoodState_LaborInt, const EquV_Auxiliary & EquV_a_BadState_LaborInt,
        const EquV_Auxiliary & EquV_a_GoodState_CapitalInt, const EquV_Auxiliary & EquV_a_BadState_CapitalInt,
        const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
        const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
        MultiThreads::Threads_Management & threadsManagement);

    /**************************************************************************
    * Estimation in the second stage
    * Objective function calculation
    **************************************************************************/
    double EstimationIndirectInference_FullVersion_Obj(std::vector<double> x_vec,
        const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
        const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
        const ArrayXd theta_a,
        const EquV_Auxiliary & EquV_a_GoodState_LaborInt, const EquV_Auxiliary & EquV_a_BadState_LaborInt,
        const EquV_Auxiliary & EquV_a_GoodState_CapitalInt, const EquV_Auxiliary & EquV_a_BadState_CapitalInt,
        const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
        const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
        MultiThreads::Threads_Management & threadsManagement);

    /**************************************************************************
    * Simulate the data for the second stage
    **************************************************************************/
    // the function
    SimData SimulationData_Step2Estimation_FullVersion(const ArrayXd & thetaData, const SimVar & sim_var,
        const double & p_K, const SimData & SimDataPanel, const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);
    // within the function: simulate the full data
    SimData Simulation_Step2Estimation_FullVersion(const ParaEst & para_est, const ParaVec & para_vec,
        const EquStateV & EquV, const EquStateVmat & Evalmat, const SimVar & sim_var, const SimData & SimDataPanel,
        const int & SimN, const int & TBar, const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);
    // within the function: generate the missing data
    SimData Simulation_Missing_Step2Estimation_FullVersion(const SimVar & sim_var, const SimData & sim_data,
        const int & SimN, const int & TBar, const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);
    // initialize the simulated variables
    SimData InitializeSimData(const int & SimN, const int & TBar);
    // calculate revenue given productivity,capital,employment
    double RevOpt1_sim_Realized(const ParaEst & para_est, const double & phi, const double & K, const double & L_ur,
        const double & L_uc);
    // calculate price
    double calPrice(const ParaEst & para_est, const double & alpha_M, const double & phi, const double & K,
        const double & L_ur, const double & L_uc, const double & PI);
    // Backout phi based on theta1 : need to differentiate whether a firm is produciton/dormant in the previous period
    ArrayXXd Backout_phi_a_PD(const ArrayXXd & Revenue, const ArrayXXd & Capital, const ArrayXXd & Employ_ur,
        const ArrayXXd & Employ_uc, const double & alpha_tilde_K,const double & alpha_tilde_L, const double & alpha_Lr,
        const double & alpha_Lc);
//
    // Simulate capital/employment given the previous year
    tuple<double,double,double,double,double,double> calOptKLurLuc(const ParaEst & para_est,
        const ParaVec & para_vec, const SimVar & sim_var, const PhiW & phi_w,
        const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state, const int & n, const int & t,
        const ArrayXXd & EVal_PD_K_mat, const ArrayXXd & EVal_PD_Lur_mat, const ArrayXXd & EVal_PD_Luc_mat,
        const int & PD_lag, const double & p_K, const int & InitialPeriod);
    //
    /**************************************************************************
    * Part1 of objection function in the second stage estimation
    **************************************************************************/
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
        const std::vector<double> & theta1, const SimData & sim_data,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);
    /**************************************************************************
    * Part2 of objection function in the second stage estimation
    **************************************************************************/
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(
        const std::vector<double> & theta2, const SimData & sim_data, MultiThreads::Threads_Management & threadsManagement);
    /**************************************************************************
    * Part3 of objection function in the second stage estimation
    **************************************************************************/
    // original value
    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(
        const SimData & sim_data, const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV,
        const EquStateVmat & Evalmat, const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);
    // take difference
    tuple<double,ArrayXXd> EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(const SimData & sim_data,
        const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV, const EquStateVmat & Evalmat,
        const ArrayXXi & K_index_max_P_mat, const ArrayXXi & Lur_index_max_P_mat, const ArrayXXi & Luc_index_max_P_mat,
        const ArrayXXd & K_max_P_mat, const ArrayXXd & Lur_max_P_mat, const ArrayXXd & Luc_max_P_mat,
        const ArrayXXi & K_index_max_D_mat, const ArrayXXi & Lur_index_max_D_mat, const ArrayXXi & Luc_index_max_D_mat,
        const ArrayXXd & K_max_D_mat, const ArrayXXd & Lur_max_D_mat, const ArrayXXd & Luc_max_D_mat,
        const int & GoodState, const int & LaborIntensive,MultiThreads::Threads_Management & threadsManagement);

    tuple<ArrayXXd,double> EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2_likelihood(
        const SimData & sim_data, const ParaEst & para_est, const ParaVec & para_vec, const EquStateVmat & Evalmat,
        const ArrayXXd & lnphi_a_P, const ArrayXXd & lnphi_a_D,
        const ArrayXXd & ProbP_S_mat, const ArrayXXd & ProbP_E_mat, const ArrayXXd & lnProbP_S_mat, const ArrayXXd & lnProbP_E_mat,
        const ArrayXXd & ProbD_S_mat, const ArrayXXd & ProbD_E_mat, const ArrayXXd & lnProbD_S_mat, const ArrayXXd & lnProbD_E_mat,
        const ArrayXXd & K_opt_P_mat, const ArrayXXd & Lur_opt_P_mat, const ArrayXXd & Luc_opt_P_mat,
        const ArrayXXd & K_opt_D_mat, const ArrayXXd & Lur_opt_D_mat, const ArrayXXd & Luc_opt_D_mat,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);

    // Calculate the probability of Lur/Luc/K for each data point; take the original optimal value as a reference
    tuple<double,double,double,int,int,int> cal_Opt_KLurLuc_Data( const ParaEst & para_est, const ParaVec & para_vec,
        const EquStateVmat & Evalmat, const PhiW & phi_w, const LKState & K_state, const LKState & Lur_state,
        const LKState & Luc_state, const int & PD_lag);
    tuple<double,double,double> cal_Opt_KLurLuc_Data_Diff( const ParaEst & para_est, const ParaVec & para_vec,
        const EquStateVmat & Evalmat, const PhiW & phi_w, const LKState & K_state, const LKState & Lur_state,
        const LKState & Luc_state, const int & K_index_max, const int & Lur_index_max, const int & Luc_index_max,
        const int & PD_lag);

    // solve the optimal K; take the original optimal value as a reference
    double sol_opt_K_Diff_simulation(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd EVal_K_mat, const LKState & Lur_state, const LKState & Luc_state, const double & K_state_val,
        const int & K_index_max);
    // solve the optimal Lur; take the original optimal value as a reference
    double sol_opt_Lur_Diff_simulation(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd & EVal_Lur_mat, const LKState & K_state, const LKState & Luc_state, const double & Lur_state_val,
        const int & Lur_index_max);
    // solve the optimal Luc; take the original optimal value as a reference
    double sol_opt_Luc_Diff_simulation(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd & EVal_Luc_mat, const LKState & K_state, const LKState & Lur_state, const double & Luc_state_val,
        const int & Luc_index_max);
    // the solver for sol_opt_K_Diff_simulation/sol_opt_Lur_Diff_simulation/sol_opt_Luc_Diff_simulation
    tuple<double,double> solve1D_K_LinearSpline_Diff_simulation(const double & H, const double & F, const double & c_K,
        const double & rstar, const double & p_K, const double & Kbase, const ArrayXd & vec_K, const ArrayXd & EVal_K_vec,
        const int & K_index_max);
    tuple<double,double> solve1D_L_LinearSpline_Diff_simulation(const double & H, const double & F, const double & wstar,
        const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const int & L_index_max);
}
