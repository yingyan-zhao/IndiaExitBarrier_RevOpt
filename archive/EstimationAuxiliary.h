#pragma once
#include <iostream>
#include <cmath>
#include <map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Ipopt_wrapper.h"
#include "DFS.h"
#include "direct_search.h"
#include "soo.h"

#include "SimulationData.h"
#include "Basics.h"
#include "SolveModel.h"
//
////#include "SolveModel_V0_FEntry.h"
//
namespace alias {
    using namespace std;
    using namespace alias;
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;

    /**************************************************************
    * Version 2: Randomly choose inital value and see the estimation
    **************************************************************/
    tuple<ArrayXd,ArrayXd,ArrayXd> EstimationAuxiliaryModel_Output_Version2(const SimData & RealData_GoodState_LaborInt,
        const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
        const SimData & RealData_BadState_CapitalInt, MultiThreads::Threads_Management & threadsManagement);

    /**************************************************************
    * Initial Parameters for the Auxilliary model
    **************************************************************/
    tuple<ArrayXd,ArrayXd,ArrayXd> SetupInitialParaGuess_forEstimation(const ArrayXd & theta_RandomInitial_select);

    /**************************************************************
    * Estimating the auxilliary model
    **************************************************************/
    tuple<ArrayXd,ArrayXd,ArrayXd,double> EstimationAuxiliaryModel(const SimData & RealData_GoodState_LaborInt,
        const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
        const SimData & RealData_BadState_CapitalInt,
        const ArrayXd & theta1_a_ini, const ArrayXd & theta2_a_ini, const ArrayXd & theta3_a_ini,
        MultiThreads::Threads_Management & threadsManagement);

    /******************************************************************
     * Auxilliary Model Estimation: Part 1
     *******************************************************************/
    //// Given the objective function, find the estimates by maximize the objective function
    ArrayXd EstimationAuxiliaryModel_part1_theta1(const SimData & RealData_GoodState_LaborInt,
        const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
        const SimData & RealData_BadState_CapitalInt, const ArrayXd & theta1_a_ini,
        MultiThreads::Threads_Management & threadsManagement);
    //// Compute the objective function
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> EstimationAuxiliaryModel_part1_ObjFun(
        const std::vector<double> & theta1, const SimData & RealData, const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);

    //// Backout phi based on theta1_a_ini
    ArrayXXd Backout_phi_a(const SimData & RealData,const double & alpha_tilde_K,const double & alpha_tilde_L,
                           const double & alpha_Lr,const double & alpha_Lc);

    /********************************************************************************************************************
    * Auxilliary Model Estimation: Part 2
    ********************************************************************************************************************/
    //// Given the objective function, find the estimates by maximize the objective function
    ArrayXd EstimationAuxiliaryModel_part2_theta2(const SimData & RealData_GoodState_LaborInt,
        const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
        const SimData & RealData_BadState_CapitalInt, const ArrayXd & theta1_a, const ArrayXd & theta2_a_ini,
        MultiThreads::Threads_Management & threadsManagement);

    //// Compute the objective function
    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> EstimationAuxiliaryModel_part2_ObjFun(
        const std::vector<double> & theta2, const ParaEst1 & para_est1, const SimData & RealData, const ArrayXXd & lnphi_a);

    /********************************************************************************************************************
    * Auxilliary Model Estimation: Part 3
    ********************************************************************************************************************/
    //// Given the objective function, find the estimates by maximize the objective function
    tuple<ArrayXd,double> EstimationAuxiliaryModel_part3_theta3(const SimData & RealData_GoodState_LaborInt,
        const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
        const SimData & RealData_BadState_CapitalInt,
        const ArrayXd & theta1_a, const ArrayXd & theta2_a, const ArrayXd & theta3_a_ini,
        MultiThreads::Threads_Management & threadsManagement);

    //// Compute the objective function
    tuple<double,int,EquStateV> EstimationAuxiliaryModel_part3_ObjFun(const std::vector<double> & theta3,
        const SimData & RealData, const ArrayXd & theta1_a, const ArrayXd & theta2_a,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);

    //// solve the value function/policy function for a given set of parameters
    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        const std::vector<double> & theta3, const ArrayXd & theta1_a, const ArrayXd & theta2_a,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);

    //// Calculate the likelihood for the data
    tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi> EstimationAuxiliaryModel_part3_likelihood(
        const SimData & RealData, const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV,
        const EquStateVmat & Evalmat, const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);

    tuple<double,int,ArrayXXd> EstimationAuxiliaryModel_part3_likelihood_Diff(const SimData & RealData,
        const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV, const EquStateVmat & Evalmat,
        const ArrayXXi & K_index_max_mat, const ArrayXXi & Lur_index_max_mat, const ArrayXXi & Luc_index_max_mat,
        const int & GoodState, const int & LaborIntensive,
        MultiThreads::Threads_Management & threadsManagement);

    ////**** Define the state of productivity ****////
    ////**** Define the state of labor and capital ****////
    struct PhiW {
        double w1;
        double w2;
        int phi_index;
        int phi_grid_choice;
        double phi;
    };
    struct LKState {
        double val_state;
        int i_val_state;

        double w_1;
        double w_2;
    };
    PhiW definePhiWeights(const ArrayXd & vec_phi, const double & phi_state);
    LKState defineLKState(const ArrayXd & vec_Dim,const double & Dim_state);

    ////**** Calculate the probability of Lur/Luc/K for each data point  ****////
    //// the main function to calculate the probability of Lur/Luc/K
    tuple<double,double,double,double,double,double,double,double,double,double,int,int,int> cal_prob_KLurLuc_Data(
        const ParaEst & para_est, const ParaVec & para_vec, const EquStateVmat & Evalmat, const PhiW & phi_w,
        const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state, const int & PD_lag,
        const double & K_state_val, const double & Lur_state_val, const double & Luc_state_val,
        const int & K_missing, const int & Lur_missing, const int & Luc_missing);

    tuple<double, double, double, double, double, double, double> cal_prob_KLurLuc_Data_Diff(
        const ParaEst & para_est, const ParaVec & para_vec, const EquStateVmat & Evalmat, const PhiW & phi_w,
        const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state, const int & PD_lag,
        const double & K_state_val, const double & Lur_state_val, const double & Luc_state_val,
        const int & K_missing, const int & Lur_missing, const int & Luc_missing,
        const int & K_index_max, const int & Lur_index_max, const int & Luc_index_max);

    //// solve the optimal Capital target given the states in the data
    tuple<double,int> sol_opt_K(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd EVal_K_mat, const LKState & Lur_state, const LKState & Luc_state, const double & K_state_val,
        const int & InitialPeriod);
    double sol_opt_K_Diff(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd EVal_K_mat, const LKState & Lur_state, const LKState & Luc_state, const double & K_state_val,
        const int & InitialPeriod, const int & K_index_max);

    //// solve the optimal Lur target given the states in the data
    tuple<double,int> sol_opt_Lur(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd & EVal_Lur_mat, const LKState & K_state, const LKState & Luc_state, const double & Lur_state_val,
        const int & InitialPeriod);
    double sol_opt_Lur_Diff(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd & EVal_Lur_mat, const LKState & K_state, const LKState & Luc_state, const double & Lur_state_val,
        const int & InitialPeriod, const int & Lur_index_max);
    //// solve the optimal Luc target given the states in the data
    tuple<double,int> sol_opt_Luc(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd & EVal_Luc_mat, const LKState & K_state, const LKState & Lur_state, const double & Luc_state_val,
        const int & InitialPeriod);
    double sol_opt_Luc_Diff(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
        const ArrayXXd & EVal_Luc_mat, const LKState & K_state, const LKState & Lur_state, const double & Luc_state_val,
        const int & InitialPeriod, const int & Luc_index_max);

    ////**** Calculate the exit/stay probability for each data point  ****////
    //// with logit distribution assumption
    tuple<double,double,double,double> calStayExit_Prob_logit(const ParaEst & para_est, const ParaVec & para_vec,
        const EquStateVmat & Evalmat, const PhiW & phi_w, const LKState & K_state, const LKState & Lur_state,
        const LKState & Luc_state, const int & PD_state);
    ////**** calculation for a specific capital / labor employment  ****////
    double calResidual_Value_Exit_double(const ParaEst & para_est, const double & p_K, const double & K,
        const double & Lur, const double & Luc);

    ////**** Interpolation the value function (EVal_K_mat; EVal_Lur_mat; EVal_Luc_mat)  ****////
    ArrayXd cal_EVal_K_phi_interpolation(const ArrayXXd & EVal_K_mat, const PhiW & phi_w,
        const LKState & Lur_state, const LKState & Luc_state, const int & NK);
    ArrayXd cal_EVal_Lur_phi_interpolation(const ArrayXXd & EVal_Lur_mat, const PhiW & phi_w,
        const LKState & K_state, const LKState & Luc_state, const int & NLur);
    ArrayXd cal_EVal_Luc_phi_interpolation(const ArrayXXd & EVal_Luc_mat, const PhiW & phi_w,
        const LKState & K_state, const LKState & Lur_state, const int & NLuc);

    tuple<double,double,double,double,double,double> calProductionDormancy_Prob_lognormal(const ParaEst & para_est,
        const ParaVec & para_vec, const EquStateVmat & Evalmat, const int & PD_lag, const double & phi,
        const double & K_val_state, const double & Lur_val_state, const double & Luc_val_state);

    double RevOpt1_sim(const ParaEst & para_est, const double & phi, const double & K, const double & L_ur,
        const double & L_uc);



    void PrintResultEstimationAuxiliaryModel(const std::vector<double> & x,const EquStateV & EquV);
//
//    //// Store the model predicted optimal capital/labor
//    struct EstArray {
//        ArrayXXd Est_KOpt;
//        ArrayXXd Est_LurOpt;
//        ArrayXXd Est_LucOpt;
//
//        ArrayXXd Est_lnpdf_KOpt;
//        ArrayXXd Est_lnpdf_LurOpt;
//        ArrayXXd Est_lnpdf_LucOpt;
//
//        ArrayXXd Est_lnProbP;
//        ArrayXXd Est_lnProbD;
//
//        ArrayXXd Est_lnProbS;
//        ArrayXXd Est_lnProbE;
//
//        ArrayXXd Est_ProbP;
//        ArrayXXd Est_ProbD;
//
//        ArrayXXd Est_ProbS;
//        ArrayXXd Est_ProbE;
//
//        ArrayXXd Est_Prob_E_Miss;
//        ArrayXXd Est_Prob_DP_Miss;
//        ArrayXXd Prob_Miss_DP;
//    };
//    EstArray initializeEstArray(const size_t & N);
//

//






//

//


//    /**************************************************************
//    * Version 1: Randomly choose inital value and see the estimation
//    **************************************************************/
//    tuple<ArrayXd,ArrayXd,ArrayXd> EstimationAuxiliaryModel_Output_Version1(const SimData & RealData_GoodState_LaborInt,
//        const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
//        const SimData & RealData_BadState_CapitalInt, MultiThreads::Threads_Management & threadsManagement);
//

////
////
////
////
////
////
//////
////
//////    tuple<double,ArrayXXd,EstArray> EstimationAuxiliaryModel_part3_Simulation_Diff(const SimData & RealData,
//////        const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV, const DiffEVal & diff_Eval,
//////        const ArrayXXi & K_index_max_mat, const ArrayXXi & Lur_index_max_mat, const ArrayXXi & Luc_index_max_mat,
//////        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);
//////
////
//////
//////    tuple<double,ArrayXXd,EstArray,EquStateV,ArrayXd> EstimationAuxiliaryModel_part3_ObjFun_Diff(
//////        const std::vector<double> & theta3, const SimData & RealData, const ArrayXd & theta1_a,
//////        const ArrayXd & theta2_a, const EquStateV & EquV_point, const DiffEVal & diff_Eval_point,
//////        const ArrayXXi & K_index_max_mat, const ArrayXXi & Lur_index_max_mat, const ArrayXXi & Luc_index_max_mat,
//////        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);
////
//
//////
//////
//////    void EstimationAuxiliaryModel_part1_theta1_checking(const std::vector<double> & theta1,const SimData & sim_data,
//////        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);
//////    void EstimationAuxiliaryModel_part2_theta2_checking(const std::vector<double> & theta2,const ArrayXd & theta1,
//////        const SimData & sim_data, const int & GoodState, const int & LaborIntensive);
//////    void EstimationAuxiliaryModel_part3_theta3_checking(const ArrayXd & theta3, const SimData & RealData,
//////        const ArrayXd & theta1_a, const ArrayXd & theta2_a, const int & GoodState, const int & LaborIntensive,
//////        MultiThreads::Threads_Management & threadsManagement);
//
//
}