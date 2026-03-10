#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <map>
#include "multithread_loop.h"

#include "Basics.h"
#include "SolveModel.h"
#include "SimulationData.h"

//
//
//
//#include "Estimation_IndirectInference.h"
//
////#include "DFS.h"
////#include "direct_search.h"
////#include "soo.h"
////
////
//
////
//
////
namespace alias {
    using namespace std;
    using namespace alias;
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;

    /**************************************************************
    * **************************************************************
    * Solve the Status quo equilibrium and simulate firm distribution
    * **************************************************************
    **************************************************************/
    tuple<CounterFactVariable, CounterFactVariable, CounterFactVariable, CounterFactVariable>
        main_SimulateStatusQuo(const ArrayXd &theta_Est_S, const SimVar &sim_var, const int &SimTbar,
        MultiThreads::Threads_Management &threadsManagement);

    tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat, EquState0, double, SimData>
        SolveSimulate_StatusQuo_FEntry_V0(
        const ArrayXd &theta_Est, const SimVar &sim_var, const int &GoodState, const int &LaborIntensive,
        MultiThreads::Threads_Management &threadsManagement);

}
//
//    /** **************************************************************
//    * Variables that we care about in counterfactural
//    * **************************************************************/
//    struct CounterFactVariable {
//        double F_Entry;
//        double p_K;
//        double p_K_Delta;
//        double NumFirm;
//        double NumFirm_Delta;
//        double FirmMass;
//        double FirmMass_Delta;
//        double CapitalDemand;
//        double CapitalDemand_Delta;
//        /*** Average Productivity ***/
//        double AvgProd;
//        double AvgProd_Delta;
//        double WAvgProd;
//        double WAvgProd_Delta;
//        double EntrantProd;
//        double EntrantProd_Delta;
//        double ExitProd;
//        double ExitProd_Delta;
//        /*** Total Labor Employment ***/
//        double TotEmploy_ur;
//        double TotEmploy_ur_Delta;
//        double TotEmploy_uc;
//        double TotEmploy_uc_Delta;
//        double TotEmploy;
//        double TotEmploy_Delta;
//        /*** Total Output ***/
//        double TotOutput;
//        double TotOutput_Delta;
//        double Output_per_Capital;
//        double Output_per_Capital_Delta;
//        double Output_per_Employ;
//        double Output_per_Employ_Delta;
//        /*** Average length of dormancy ***/
//        double DormLength;
//        /*** Average length of production ***/
//        double ProdLength;
//        /*** Average length of Survival ***/
//        double SurvivalLength;
//        /*** Price Index ***/
//        double PriceIndex;
//        double PriceIndex_Delta;
//        /*** Average Exit Rate ***/
//        double AvgExitRate;
//        /*** Capital Labor Ratio ***/
//        double K_L_ratio;
//        double K_L_ratio_Delta;
//
//        double Welfare_Manu;
//        double Welfare_Manu_Delta;
//
//
//        double AggValueAdded;
//    };
//

//
////
//    ///* Solve the entry cost */
//    double SolveFEntry_StatusQuoEqu(const ParaVec &para_vec, const EquStateV0 &EquV0);
//
//    /**********************************************************************************************
//    * Simulation of firm distribution for Status quo / counterfactual
//    **********************************************************************************************/
//    SimData Simulation_StatusQuo_Counterfactual(const ParaEst &para_est, const ParaVec &para_vec,
//        const EquStateV &EquV, const EquStateVmat &Evalmat,
//        const EquStateV0 &EquV0, const EquStateV0mat &Eval0mat,
//        const SimVar &sim_var, const int &SimN, const int &TBar,
//        const int &GoodState, const int &LaborIntensive,
//        MultiThreads::Threads_Management &threadsManagement);
//
//
//    /**********************************************************************************************
//    * Calculate Aggregates for counterfactuals
//    **********************************************************************************************/
//    CounterFactVariable CalEconomicVariable(const double &p_K, const double &FirmMass, const double &F_Entry,
//        const double &CapitalDemand, SimData &sim_data);
//
//    /*** Calculate Average Productivity ***/
//    double calAvgProductivity(SimData &sim_data_C, const double FirmMass);
//
//    double calWAvgProductivity(SimData &sim_data_C, const double FirmMass);
//
//    double calEntrantProductivity(SimData &sim_data_C, const double FirmMass);
//
//    double calExitProductivity(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Total Labor Employment ***/
//    tuple<double, double> calTotalEmployment(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Total Output ***/
//    double calTotalOutput(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Average length of dormancy ***/
//    double calLengthDormancy(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Average length of production ***/
//    double calLengthProduction(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Average length of Survival ***/
//    double calLengthSurvival(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Price Index ***/
//    double calPriceIndex(SimData &sim_data_C, const double FirmMass);
//
//    /*** Calculate Average Exit Rate ***/
//    double calAvgExitRate(SimData &sim_data_C, const double FirmMass);
//
//    double calCapitalLaborRatio(SimData &sim_data_C, const double FirmMass);
//
//    /******** Auxiliary function for counterfactual ********/
//    struct CFDeltaVec {
//        ArrayXd FirmMass_Delta;
//        ArrayXd NumFirm_Delta;
//        ArrayXd CapitalDemand_Delta;
//        ArrayXd p_K_Delta;
//        ArrayXd AvgProd_Delta;
//        ArrayXd WAvgProd_Delta;
//        ArrayXd EntrantProd_Delta;
//        ArrayXd ExitProd_Delta;
//        ArrayXd TotEmploy_ur_Delta;
//        ArrayXd TotEmploy_uc_Delta;
//        ArrayXd TotEmploy_Delta;
//        ArrayXd TotOutput_Delta;
//        ArrayXd Output_per_Capital_Delta;
//        ArrayXd Output_per_Employ_Delta;
//        ArrayXd DormLength;
//        ArrayXd ProdLength;
//        ArrayXd SurvivalLength;
//        ArrayXd PriceIndex_Delta;
//        ArrayXd AvgExitRate;
//        ArrayXd K_L_ratio_Delta;
//        ArrayXd firing_share;
//        ArrayXd ResidualValue;
//        ArrayXd EntryCost;
//        ArrayXd ExitSubsidy;
//        ArrayXd AggValueAdded;
//
//        ArrayXd Welfare_Manu_Delta;
//    };
//
//    CFDeltaVec InitializeCFDeltaVec(const int &N_C);
//
//    ArrayXd Assign_thetaC(const ArrayXd &theta_Est_S, const ArrayXd &ResidualValue_row, const ArrayXd &firing_share_row,
//                          const int &GoodState, const int &LaborIntensive);
//
//    CounterFactVariable CalEconomicVariableDelta(const CounterFactVariable &Value_C, const CounterFactVariable &Value_StatusQuo);
//
//    int writeToCSVfileValueDelta(const string &filename, const CFDeltaVec &CFDeltaVec_Status);
//
//
//    /**************************************************************
//    * **************************************************************
//    * Solve the Partial Counterfactural equilibrium
//    * **************************************************************
//    **************************************************************/
//    /**************************************************************
//    *** Counterfactual Exercise 1: Partial Equilibrium: Varying the exit barrier ***
//    **************************************************************/
//    int SolveSimulate_CounterfacturalResult_1_v1(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//        const string & FKL, const ArrayXd & TargetExitRate_vec, MultiThreads::Threads_Management & threadsManagement);
//
//    int SolveSimulate_CounterfacturalResult_1_v2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//        const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas,
//        const string &FKL, const ArrayXd &ResVal_vec, MultiThreads::Threads_Management &threadsManagement);
//
//    /**************************************************************
//    *** Counterfactual Exercise 2: Partial Equilibrium: Varying the exit barrier ***
//    **************************************************************/
//    int SolveSimulate_CounterfacturalResult_2_v1(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//        const string & FKL, const ArrayXd & TargetExitRate_vec, MultiThreads::Threads_Management & threadsManagement);
//
//    int SolveSimulate_CounterfacturalResult_2_v2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//        const string & FKL, const ArrayXd & firingshare_vec, MultiThreads::Threads_Management & threadsManagement);
//
//    /**************************************************************
//    * Counterfactual Exercise 3: Partial Equilibrium: Change KsupplyElas and Change the exit barrier to match the US exit rate
//    **************************************************************/
//    int SolveSimulate_CounterfacturalResult_3_Weighted(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN,
//        const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    int SolveSimulate_CounterfacturalResult_3_Aggregate(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN,
//        const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    /**************************************************************
//    * Counterfactual Exercise 4: Partial Equilibrium: Change KsupplyElas and Change the firing cost to match the US exit rate
//    **************************************************************/
//    int SolveSimulate_CounterfacturalResult_4_Weighted(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN,
//        const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    int SolveSimulate_CounterfacturalResult_4_Aggregate(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN,
//        const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    /**************************************************************
//    * Counterfactual Exercise 5: Partial Equilibrium: Varying both the exit barrier and firing costs at the same time
//    **************************************************************/
//    int SolveSimulate_CounterfacturalResult_5(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//        const double & Residual_up, const double & firingshare_low, const int & N_C,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    int SolveSimulate_CounterfacturalResult_6(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    /**************************************************************
//    * **************************************************************
//    * Partial Equilibrium: functions.
//    * **************************************************************
//    **************************************************************/
//    double SolveResidualValueMatchingTargetExitRate(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo,
//        const SimVar &sim_var, const int &GoodState, const int &LaborIntensive,
//        const double &TargetExitRate,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    double SolvefiringshareMatchingTargetExitRate(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo,
//        const SimVar &sim_var, const int &GoodState, const int &LaborIntensive,
//        const double &TargetExitRate,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat, EquState0, double, double, SimData>
//        SolveSimulate_Counterfactural_FEntry_V0_EquPrice(const ArrayXd &theta_C, const double &F_Entry,
//        const double &CapitalStock, const double &PriceIndex_S, const double &FirmMass_S, const SimVar &sim_var_C,
//        const int &GoodState, const int &LaborIntensive, const double &KsupplyElas,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    double Solve_Counterfactural_FEntry_V0_Solve_EquPrice(const ParaEst &para_est_C, const double &F_Entry,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    tuple<ArrayXXd,ArrayXXd> CalSimulation_CounterfactualResult(const ArrayXd &theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//        const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas,
//        const ArrayXXd &ResidualValue, const ArrayXXd &firing_share, const ArrayXXd &EntryCost, const string &FKL,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    double CalEntrySubsidy_per_Entrants(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_Segment,
//        const SimVar &sim_var, const int SimTbar, const int SimN, const int &GoodState, const int &LaborIntensive,
//        const double &KsupplyElas, const double &ExitSubsidy, const double & NumEntrants,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    /**************************************************************
//    * **************************************************************
//    * Solve the General Counterfactural equilibrium
//    * **************************************************************
//    **************************************************************/
//    /**************************************************************
//    *** Counterfactual Exercise 1: General Equilibrium: Varying the exit barrier ***
//    **************************************************************/
//    int SolveSimulate_GECounterfacturalResult_1(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//        const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas, const ArrayXd & ResVal_vec,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    int SolveSimulate_GECounterfacturalResult_2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//        const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas, const ArrayXd &firingshare_vec,
//        MultiThreads::Threads_Management &threadsManagement);
//
//    tuple<ArrayXXd,ArrayXXd> CalSimulation_GECounterfactualResult(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
//        const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
//        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//        const ArrayXXd & ResidualValue, const ArrayXXd & firing_share, const ArrayXXd & EntryCost, const string & FKL,
//        const double & p_K_C_low, const double & p_K_C_high, MultiThreads::Threads_Management & threadsManagement);
//
//    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat> Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(
//        const ArrayXd &theta_C, const double & F_Entry, const double & p_K, const int & GoodState, const int & LaborIntensive,
//        MultiThreads::Threads_Management & threadsManagement);
//
//    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,double> Solve_GECounterfactural_FEntry_V0_Cal_VEntry(
//        const ArrayXd &theta_C, const double & PI, const double & F_Entry, const double & p_K,
//        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);
//
//    tuple<double,double,double,double> CalFirmMass_Given_sim_data(const ParaEst &para_est_C_GoodLaborInt,
//        const ParaEst &para_est_C_BadLaborInt, const ParaEst &para_est_C_GoodCapitalInt, const ParaEst &para_est_C_BadCapitalInt,
//        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//        SimData & sim_data_C_GoodLaborInt,SimData & sim_data_C_BadLaborInt,
//        SimData & sim_data_C_GoodCapitalInt,SimData & sim_data_C_BadCapitalInt);
//
//    tuple<double,double,double,double> CalFirmMass_Given_sim_data_old(const ParaEst &para_est_C_GoodLaborInt,
//        const ParaEst &para_est_C_BadLaborInt, const ParaEst &para_est_C_GoodCapitalInt, const ParaEst &para_est_C_BadCapitalInt,
//        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//        SimData & sim_data_C_GoodLaborInt,SimData & sim_data_C_BadLaborInt,
//        SimData & sim_data_C_GoodCapitalInt,SimData & sim_data_C_BadCapitalInt);
//
////
//
////
//
//
////    int SolveSimulate_GECounterfacturalResult_2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
////        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////        const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas,
////        const int &N_C, MultiThreads::Threads_Management &threadsManagement);
////        //
////    int SolveSimulate_GECounterfacturalResult_3(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
////        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////        const SimVar &sim_var, const int SimTbar, const int SimN,
////        const int &N_C, MultiThreads::Threads_Management &threadsManagement);
////
////
////
//////    ArrayXd Assign_thetaC_SingleSegment(const ArrayXd &theta_Est_S, const double & ResidualValue,
//////        const double & firing_share, const int &GoodState, const int &LaborIntensive);
//////
//////    tuple<ArrayXXd,ArrayXXd> CalSimulation_CounterfactualResult_Test_SingleSegment(const ArrayXd & theta_Est_S,
//////        const CounterFactVariable & Value_StatusQuo, const int &GoodState, const int &LaborIntensive,
//////        const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//////        const ArrayXd & ResidualValue, const ArrayXd & firing_share, const ArrayXd & EntryCost, const string & FKL,
//////        MultiThreads::Threads_Management & threadsManagement);
//////
//////    tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat, EquState0, double, double, SimData>
//////        SolveSimulate_Counterfactural_FEntry_V0_EquPrice_SingleSegment(const ArrayXd &theta_C, const double &F_Entry,
//////        const double &CapitalStock, const double &PriceIndex_S, const double &FirmMass_S,
//////        const SimVar &sim_var_C, const int &GoodState, const int &LaborIntensive, const double &KsupplyElas,
//////        MultiThreads::Threads_Management &threadsManagement);
//////
//////    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat> Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(
//////        const ArrayXd &theta_C, const double & F_Entry, const double & p_K, const int & GoodState, const int & LaborIntensive,
//////        MultiThreads::Threads_Management & threadsManagement);
//////
//////    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,double> Solve_Counterfactural_FEntry_V0_Cal_VEntry(
//////        const ArrayXd &theta_C, const double & PI, const double & F_Entry, const double & p_K,
//////        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement);
//////
//////    tuple<double,double,double,double> CalFirmMass_Given_sim_data(const ParaEst &para_est_C_GoodLaborInt,
//////        const ParaEst &para_est_C_BadLaborInt, const ParaEst &para_est_C_GoodCapitalInt, const ParaEst &para_est_C_BadCapitalInt,
//////        const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
//////        const CounterFactVariable &Value_StatusQuo_BadLaborInt, const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
//////        const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//////        SimData & sim_data_C_GoodLaborInt, SimData & sim_data_C_BadLaborInt,
//////        SimData & sim_data_C_GoodCapitalInt, SimData & sim_data_C_BadCapitalInt);
//////
//////    double CalFirmMass_Given_sim_data_SingleSegment(const ParaEst &para_est_C, const CounterFactVariable &Value_StatusQuo,
//////        SimData & sim_data_C_GoodLaborInt, const int &GoodState, const int &LaborIntensive);
//////
//////    /**************************************************************
//////    * Counterfactual test -- Single Segment: Varying the exit barrier
//////    **************************************************************/
//////    int SolveSimulate_CounterfacturalResult_Test_SingleSegment(const ArrayXd &theta_Est_S,
//////        const CounterFactVariable & Value_StatusQuo, const int &GoodState, const int &LaborIntensive,
//////        const SimVar &sim_var, const int SimTbar, const int SimN,
//////        const double &KsupplyElas, const int &N_C, MultiThreads::Threads_Management &threadsManagement);
//////
//////    /**************************************************************
//////    * Counterfactual Exercise 1: Varying the exit barrier
//////    **************************************************************/
//////    int SolveSimulate_CounterfacturalResult_1(const ArrayXd &theta_Est_S, const SegmentShr &shr_RealData,
//////        const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
//////        const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
//////        const SimVar &sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
//////        const int &N_C, MultiThreads::Threads_Management &threadsManagement);
//////
//////
////
////////
////////    /**************************************************************
////////    * Counterfactual Exercise 2: Varying the exit barrier
////////    **************************************************************/
////////    int
////////    SolveSimulate_CounterfacturalResult_2_v1(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_All,
////////                                             const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                             const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                             const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                             const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                             const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                             const double &KsupplyElas,
////////                                             const int &N_C, MultiThreads::Threads_Management &threadsManagement);
////////
////////    int
////////    SolveSimulate_CounterfacturalResult_2_v2(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_All,
////////                                             const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                             const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                             const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                             const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                             const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                             const double &KsupplyElas,
////////                                             const int &N_C, MultiThreads::Threads_Management &threadsManagement);
////////
////////    /**************************************************************
////////    * Counterfactual Exercise 3: Change KsupplyElas and Change the exit barrier to match the US exit rate
////////    **************************************************************/
////////    int
////////    SolveSimulate_CounterfacturalResult_3(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_All,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                          const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                          const ArrayXd &TargetExitRate_LaborInt_vec,
////////                                          const ArrayXd &TargetExitRate_CapitalInt_vec,
////////                                          MultiThreads::Threads_Management &threadsManagement);
////////
////////    /**************************************************************
////////    * Counterfactual Exercise 4: Change KsupplyElas and Change firing costs to match the US exit rate
////////    **************************************************************/
////////    int
////////    SolveSimulate_CounterfacturalResult_4(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_All,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                          const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                          const ArrayXd &TargetExitRate_LaborInt_vec,
////////                                          const ArrayXd &TargetExitRate_CapitalInt_vec,
////////                                          MultiThreads::Threads_Management &threadsManagement);
////////
////////    /**************************************************************
////////    * Counterfactual Exercise 5: Varying both the exit barrier and firing costs at the same time
////////    **************************************************************/
////////    int
////////    SolveSimulate_CounterfacturalResult_5(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_All,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                          const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                          const double &KsupplyElas, const double &Residual_up,
////////                                          const double &firingshare_low,
////////                                          const int &N_C, MultiThreads::Threads_Management &threadsManagement);
////////
////////    /**************************************************************
////////    * Counterfactual Exercise 8: Entry Cost vs Exit Cost
////////    **************************************************************/
////////    int
////////    SolveSimulate_CounterfacturalResult_8(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_All,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                          const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                          const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                          const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                          const double &KsupplyElas,
////////                                          MultiThreads::Threads_Management &threadsManagement);
////////
////////    double CalExitSubsidy(const ArrayXd &theta_Est_S, const CounterFactVariable &Value_StatusQuo_Segment,
////////        const SimVar &sim_var, const int SimTbar, const int SimN, const int GoodState, const int LaborIntensive,
////////        const double & KsupplyElas, const ArrayXd & ResidualValue, const ArrayXd & firing_share,
////////        MultiThreads::Threads_Management &threadsManagement);
////////
//
////////
////////    double CalEntrySubsidy_per_Entrants_ForAll(const ArrayXd &theta_Est_S,
////////                                               const CounterFactVariable &Value_StatusQuo_GoodLaborInt,
////////                                               const CounterFactVariable &Value_StatusQuo_BadLaborInt,
////////                                               const CounterFactVariable &Value_StatusQuo_GoodCapitalInt,
////////                                               const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
////////                                               const SimVar &sim_var, const int SimTbar, const int SimN,
////////                                               const double &KsupplyElas, const double &ExitSubsidy_Tot,
////////                                               const double & NumEntrants_Tot,
////////                                               MultiThreads::Threads_Management &threadsManagement);
//////////
//}
//////////tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
//////////SolveSimulate_Counterfactural_FEntry_V0(const ArrayXd & theta_C, const double & F_Entry,
//////////                                        const double & CapitalStock, const double & PriceIndex_S, const double & FirmMass_S, const SimVar & sim_var_C,
//////////                                        const int & GoodState, const int & LaborIntensive, const IndustryPara & para_ind, const double & KsupplyElas,
//////////                                        MultiThreads::Threads_Management & threadsManagement);
//////////
///////////**************************************************************
//////////* Counterfactual Exercise 6: With a target exit rate, Draw the PPF of exit barrier and firing costs, with employment loss as the cost curve
//////////**************************************************************/
//////////int SolveSimulate_CounterfacturalResult_6(const ArrayXd & theta_Est_S, const CounterFactVariable & Val_StatusQuo,
//////////                                          const SimVar & sim_var, const int SimTbar, const int SimN, const int & GoodState, const int & LaborIntensive,
//////////                                          const IndustryPara & para_ind, const double & KsupplyElas, const string & State, const double & TargetExitRate,
//////////                                          const int & N_C, MultiThreads::Threads_Management & threadsManagement);
//////////
///////////**************************************************************
//////////* Counterfactual Exercise 7: With a target exit rate, Draw the PPF of exit barrier and firing costs, with employment loss as the cost curve
//////////**************************************************************/
//////////int SolveSimulate_CounterfacturalResult_7(const ArrayXd & theta_Est_S, const CounterFactVariable & Val_StatusQuo,
//////////                                          const SimVar & sim_var, const int SimTbar, const int SimN, const int & GoodState, const int & LaborIntensive,
//////////                                          const IndustryPara & para_ind, const double & KsupplyElas, const string & State, const double & TargetOutput_Delta,
//////////                                          const int & N_C, MultiThreads::Threads_Management & threadsManagement);
