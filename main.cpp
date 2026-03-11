#include <iostream>
#include <Eigen/Dense>
#include <thread>
#include "multithread_loop.h"

#include "Basics.h"
#include "SolveModel.h"
#include "SimulationData.h"
#include "CounterfactualExercise.h"

//// #include "EstimationAuxiliary.h"
// //#include "Estimation_IndirectInference.h"
// //#include "EstimationIndirectInference_StandardDeviation.h"
// //
// //
// ////#include "Ipopt_wrapper.h"
// ////
// //////
using namespace std;
using namespace Eigen;
using namespace alias;
// ////using namespace Ipopt_Wrapper;
// //
// //
int main() {
    /**************************************************************
    * Set up the environment
    **************************************************************/
    /* Prepare Multithread for parallel computation  */
    const int n_threads = std::thread::hardware_concurrency();
    cout << "n_threads = " << n_threads << endl;
    MultiThreads::Threads_Management threadsManagement;
    //
    // int GoodState = 1;
    // int LaborIntensive = 1;
    // ArrayXd theta_Est_S(para.dim1+para.dim2+para.dim3);
    // theta_Est_S << 103.318, 218.532, 178.482, 256.496, 200.002, 0.381529, -413.768,-626.425, 0.5, 528.122,
    //      0.198573, 0.180748, 0.892402, 0.713,
    //      -12.5651, 11.973, 9.1589;
    // cout << para.vec_norm_phi << endl;
    // ParaEst para_est = constructParaEst(theta_Est_S,GoodState,LaborIntensive);
    // para_est.PI = 1.0;
    // ParaVec para_vec = constructParaFull(para_est);
    //
    //  //// solve value function for t > 1
    //  EquStateV EquV = solveV(para_est,para_vec,VEnd,threadsManagement);
//
// //
// //    /**************************************************************
// //    * Import the data
// //    **************************************************************/
// //    int GoodState = 1;
// //    int LaborIntensive = 1;
// //    SimData data_panel_GoodState_LaborInt = createRealData("", GoodState, LaborIntensive);
// //    int N_RealData_GoodState_LaborInt = data_panel_GoodState_LaborInt.good.rows(); // number of observations
// //
// //    GoodState = 0;
// //    LaborIntensive = 1;
// //    SimData data_panel_BadState_LaborInt = createRealData("", GoodState, LaborIntensive);
// //    int N_RealData_BadState_LaborInt = data_panel_BadState_LaborInt.good.rows(); // number of observations
// //
// //    GoodState = 1;
// //    LaborIntensive = 0;
// //    SimData data_panel_GoodState_CapitalInt = createRealData("", GoodState, LaborIntensive);
// //    int N_RealData_GoodState_CapitalInt = data_panel_GoodState_CapitalInt.good.rows(); // number of observations
// //
// //    GoodState = 0;
// //    LaborIntensive = 0;
// //    SimData data_panel_BadState_CapitalInt = createRealData("", GoodState, LaborIntensive);
// //    int N_RealData_BadState_CapitalInt = data_panel_BadState_CapitalInt.good.rows(); // number of observations
// //
// //    cout << "/*** Load the actual data done.  ***/" << endl;
// //    cout << "Number of observation - N_RealData_GoodState = " << N_RealData_GoodState_LaborInt << endl;
// //    cout << "Number of observation - N_RealData_BadState = " << N_RealData_BadState_LaborInt << endl;
// //    cout << "Number of observation - N_RealData_GoodState = " << N_RealData_GoodState_CapitalInt << endl;
// //    cout << "Number of observation - N_RealData_BadState = " << N_RealData_BadState_CapitalInt << endl;
// ////    throw runtime_error("Import the Data: 39");
// //
// //    SegmentShr shr_RealData;
// //    shr_RealData.shr_GoodState_LaborInt = double(N_RealData_GoodState_LaborInt)
// //        / double(N_RealData_GoodState_LaborInt+N_RealData_BadState_LaborInt+N_RealData_GoodState_CapitalInt+N_RealData_BadState_CapitalInt);
// //    shr_RealData.shr_BadState_LaborInt = double(N_RealData_BadState_LaborInt)
// //        / double(N_RealData_GoodState_LaborInt+N_RealData_BadState_LaborInt+N_RealData_GoodState_CapitalInt+N_RealData_BadState_CapitalInt);
// //    shr_RealData.shr_GoodState_CapitalInt = double(N_RealData_GoodState_CapitalInt)
// //        / double(N_RealData_GoodState_LaborInt+N_RealData_BadState_LaborInt+N_RealData_GoodState_CapitalInt+N_RealData_BadState_CapitalInt);
// //    shr_RealData.shr_BadState_CapitalInt = double(N_RealData_BadState_CapitalInt)
// //        / double(N_RealData_GoodState_LaborInt+N_RealData_BadState_LaborInt+N_RealData_GoodState_CapitalInt+N_RealData_BadState_CapitalInt);
// //    cout << "shr_RealData = " << shr_RealData.shr_GoodState_LaborInt << "; " << shr_RealData.shr_BadState_LaborInt
// //         << "; " << shr_RealData.shr_GoodState_CapitalInt << "; " << shr_RealData.shr_BadState_CapitalInt << endl;
// //
// //    /**************************************************************
// //    * Estimation Part 1: Auxilliary model
// //    **************************************************************/
// ////    /*** Version 1: Randomly choose inital value and see the estimation ***/
// ////    // Generate 10 initial guess and esimtation
// //////    tuple<ArrayXd,ArrayXd,ArrayXd> t_EstAux_V1 = EstimationAuxiliaryModel_Output_Version1(
// //////            data_panel_GoodState_LaborInt,data_panel_BadState_LaborInt,
// //////            data_panel_GoodState_CapitalInt,data_panel_BadState_CapitalInt,
// //////            threadsManagement);
// //////    ArrayXd theta1_a = get<0>(t_EstAux_V1);
// //////    ArrayXd theta2_a = get<1>(t_EstAux_V1);
// //////    ArrayXd theta3_a = get<2>(t_EstAux_V1);
// //////////
// /////
// //    //  /*** Version 2: Choose a specific initial value and see the estimation ***/
// //    //  tuple<ArrayXd, ArrayXd, ArrayXd> t_EstAux_V2 = EstimationAuxiliaryModel_Output_Version2(
// //    //          data_panel_GoodState_LaborInt, data_panel_BadState_LaborInt,
// //    //          data_panel_GoodState_CapitalInt, data_panel_BadState_CapitalInt,
// //    //          threadsManagement);
// //    //  ArrayXd theta1_a = get<0>(t_EstAux_V2);
// //    //  ArrayXd theta2_a = get<1>(t_EstAux_V2);
// //    //  ArrayXd theta3_a = get<2>(t_EstAux_V2);
// //    //  cout << "theta1_a = " << theta1_a.transpose() << endl;
// //    //  cout << "theta2_a = " << theta2_a.transpose() << endl;
// //    //  cout << "theta3_a = " << theta3_a.transpose() << endl;
// //    //  cout << "/*** Auxilliary model Estimation done: 65 ***/" << endl;
// //    // throw runtime_error("Auxilliary model Estimation done: 78");
// //// ////////////
// ////
// ////       /**************************************************************
// ////       * Estimation Part 2: Second Stage model
// ////       **************************************************************/
// ////        // ArrayXd theta1_a = ArrayXd::Zero(para.dim1);
// ////        theta1_a << 0.198573, 0.180748, 0.892402, 0.713;
// ////        // ArrayXd theta2_a = ArrayXd::Zero(para.dim2);
// ////        theta2_a << -12.5651, 11.973, 9.1589;
// ////       // ArrayXd theta3_a = ArrayXd::Zero(para.dim3);
// ////       //  theta3_a << 61.6331, 179.402, 1662.86,103.318,218.532,178.482,256.496,200.002,
// ////       //          5.19303,38.2064,36.8367,0.381529,0.3766,0.680907,-413.768,-626.425,528.122; // full sample
// ////       theta3_a << 64.3592, 168.407, 1662.58, 115.746, 247.061, 178.294, 273.716, 195.673, 5.36315, 31.725, 30.3631,
// ////           0.380091, 0.377282, 0.68308, -414.518, -621.988, 508.867; // only single unit firms
// //// // //
// //// // //
// ////      theta1_a << 0.19857250857,0.180748214758,0.892402465991,0.712999509591;
// ////      theta2_a << -12.5650976354,11.9729819971,9.15889702071;
// ////      theta3_a << 64.3591932411,168.407327531,1662.57579102,115.745612122,247.061467265,178.294113248,273.716457993,195.673380953,
// ////          5.36315306194,31.7249969942,30.3631367811,0.380090876926,0.377281880635,0.683080451129,
// ////          -414.517906307,-621.987707083,508.86704831;
// ////      ArrayXd theta_a(para.dim3+para.dim2+para.dim1);
// ////      theta_a << theta3_a, theta1_a, theta2_a;
// //// // //
// ////       // ArrayXd theta_Est = EstimationIndirectInference_Output(theta1_a, theta2_a, theta3_a, threadsManagement);
// ////       cout << "/*** The Second stage of estimation done: 85 ***/" << endl;
// ////       // throw runtime_error("The Second stage of estimation done: 109");
// ////     ArrayXd theta_Est(para.dim3+para.dim2+para.dim1);
// ////     // x = 36.3564290951, 174.906301424, 1663.01504136, 120.745553674, 246.71777459, 176.793953137, 272.19883312, 199.267361457, 1.61319803097, 39.9749920247, 34.6326721631, 0.377231346324, 0.376776459801, 0.683424640327, -409.955406269, -622.205480495, 521.117048314, 0.191979283045, 0.174697166717, 0.892388848844, 0.71317358153, -13.0025976353, 17.4729819971, 10.4870220206, ;
// //
// ////     // //
// ////      /**************************************************************
// ////      * Standard Deviation
// ////      **************************************************************/
// ////      /** estimates from the second stage **/
// ////      ArrayXd theta1_Est = ArrayXd::Zero(para.dim1);
// ////      theta1_Est << 0.191979283045, 0.174697166717, 0.892388848844, 0.71317358153;
// ////      ArrayXd theta2_Est = ArrayXd::Zero(para.dim2);
// ////      theta2_Est << -13.0025976353, 17.4729819971, 10.4870220206;
// ////      ArrayXd theta3_Est = ArrayXd::Zero(para.dim3);
// ////      theta3_Est << 36.3564290951, 174.906301424, 1663.01504136, 120.745553674, 246.71777459, 176.793953137, 272.19883312, 199.267361457, 1.61319803097, 39.9749920247, 34.6326721631, 0.377231346324, 0.376776459801, 0.683424640327, -409.955406269, -622.205480495, 521.117048314;
// ////      // ArrayXd theta_Est(para.dim3+para.dim2+para.dim1);
// ////      theta_Est << theta3_Est, theta1_Est, theta2_Est;
// //// //
// ////     cout << "144" << endl;
// ////      //
// ////      ArrayXd theta_SD = EstimationIndirectInference_StandardDeviation(theta_Est,theta_a,
// ////          data_panel_GoodState_LaborInt, data_panel_BadState_LaborInt,
// ////          data_panel_GoodState_CapitalInt,data_panel_BadState_CapitalInt,
// ////          threadsManagement);
// //// //     //
// ////      cout << "/*** The standard deviation of the Second stage of estimation done: 129 ***/" << endl;
// ////      throw runtime_error("The standard deviation of the Second stage of estimation done: 139");
// //// //
// //// //
    /**************************************************************
    * Counterfactual Exercise: Current Status
    **************************************************************/
    int SimTbar = 50;
    int SimN = 100000;
    // simulate random numbers for simulation
    SimVar sim_var = SimulationBase(SimN, SimTbar);
    // Value of Estimates
    ArrayXd theta_Est_S(para.dim1 + para.dim2 + para.dim3);
    theta_Est_S << 103.318, 218.532, 178.482, 256.496, 200.002, 0.381529, -413.768,-626.425, 0.5, 528.122,
    0.198573, 0.180748, 0.892402, 0.713,
    -12.5651, 11.973, 9.1589;

    /**************************************************************************************************************************/
    /*** Current Status:: Simulate the data and Get the value for the Status Quo as a base for comparison later ***/
    tuple<CounterFactVariable, CounterFactVariable, CounterFactVariable, CounterFactVariable> t_StatusQuo =
        main_SimulateStatusQuo(theta_Est_S, sim_var, SimTbar, threadsManagement);
    CounterFactVariable Value_StatusQuo_GoodLaborInt = get<0>(t_StatusQuo);
    CounterFactVariable Value_StatusQuo_BadLaborInt = get<1>(t_StatusQuo);
    CounterFactVariable Value_StatusQuo_GoodCapitalInt = get<2>(t_StatusQuo);
    CounterFactVariable Value_StatusQuo_BadCapitalInt = get<3>(t_StatusQuo);
// //
// //     cout << "EconVal.K_L_ratio: GoodLabor = " << Value_StatusQuo_GoodLaborInt.K_L_ratio <<
// //         "; BadLabor " << Value_StatusQuo_BadLaborInt.K_L_ratio <<
// //         "; GoodCapital " << Value_StatusQuo_GoodCapitalInt.K_L_ratio <<
// //         "; BadCapital " << Value_StatusQuo_BadCapitalInt.K_L_ratio << endl;
// //
// //     cout << "Dormlength: GoodLabor = " << Value_StatusQuo_GoodLaborInt.DormLength
// //         << "; BadLabor = " << Value_StatusQuo_BadLaborInt.DormLength
// //         << "; GoodCapital = " << Value_StatusQuo_GoodCapitalInt.DormLength
// //         << "; BadCapital = " << Value_StatusQuo_BadCapitalInt.DormLength << endl;
// //
// //     cout << "SurvivaLength: GoodLabor = " << Value_StatusQuo_GoodLaborInt.SurvivalLength
// //         << "; BadLabor = " << Value_StatusQuo_BadLaborInt.SurvivalLength
// //         << "; GoodCapital = " << Value_StatusQuo_GoodCapitalInt.SurvivalLength
// //         << "; BadCapital = " << Value_StatusQuo_BadCapitalInt.SurvivalLength << endl;
// //
// //     cout << "Average Exit Rate: GoodLabor = " << Value_StatusQuo_GoodLaborInt.AvgExitRate
// //          << "; BadLabor = " << Value_StatusQuo_BadLaborInt.AvgExitRate
// //          << "; GoodCapital = " << Value_StatusQuo_GoodCapitalInt.AvgExitRate
// //          << "; BadCapital = " << Value_StatusQuo_BadCapitalInt.AvgExitRate << endl;
// //
// //     cout << "Value-added: GoodLaborInt = " << Value_StatusQuo_GoodLaborInt.TotOutput
// //         << "; BadLaborInt" << Value_StatusQuo_BadLaborInt.TotOutput
// //         << "; GoodCapitalInt" << Value_StatusQuo_GoodCapitalInt.TotOutput
// //         << "; BadCapitalInt" << Value_StatusQuo_BadCapitalInt.TotOutput << endl;
// //
// //     cout << "Price Index: GoodLaborInt = " << Value_StatusQuo_GoodLaborInt.PriceIndex
// //         << "; BadLaborInt" << Value_StatusQuo_BadLaborInt.PriceIndex
// //         << "; GoodCapitalInt" << Value_StatusQuo_GoodCapitalInt.PriceIndex
// //         << "; BadCapitalInt" << Value_StatusQuo_BadCapitalInt.PriceIndex << endl;
// //     cout << "/*** Current Status:: Simulate the data and Get the value for the Status Quo as a base for comparison later ***/" << endl;
// //     // throw runtime_error("117");
// ////
// ////      /**************************************************************
// ////      * Counterfactual Exercise: Partial Equilibrium
// ////      **************************************************************/
// ////       /*** Counterfactual 1.1 (Table): Varying the exit barrier ***/
// ////       double KsupplyElas = 0.0;
// ////       int N_C = 8;
// ////       ArrayXd TargetExitRate_vec(N_C);
// ////       TargetExitRate_vec << Value_StatusQuo_GoodLaborInt.AvgExitRate,Value_StatusQuo_BadLaborInt.AvgExitRate,
// ////           Value_StatusQuo_GoodCapitalInt.AvgExitRate, Value_StatusQuo_BadCapitalInt.AvgExitRate,
// ////           0.045,0.054,0.063,0.072;
// //// //
// ////       string FKL = "_Partial_C1_ResVal_v1";
// //// //      int PartialCounterfactual1_v1_Output = SolveSimulate_CounterfacturalResult_1_v1(theta_Est_S, shr_RealData,
// //// //          Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //// //          Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //// //          sim_var, SimTbar, SimN, KsupplyElas, FKL, TargetExitRate_vec, threadsManagement);
// //// //      cout << "Counterfactual 1 v1: targeting exit rate for each segment: Varying ResidualVal done" << endl;
// //// //       // throw runtime_error("181");
// //// // //
// ////        /*** Counterfactual 1.2 (Graph): Varying the exit barrier ***/
// ////       KsupplyElas = 0.0;
// ////       N_C = 9;
// ////       ArrayXd ResVal1_vec = ArrayXd::LinSpaced(N_C, -1000, 0);
// ////       ArrayXd ResVal2_vec = ArrayXd::LinSpaced(N_C, 50, 1000);
// ////       ArrayXd ResVal_vec(N_C+N_C);
// ////       ResVal_vec << ResVal1_vec,ResVal2_vec;
// ////       // ResVal_vec << -300,-150,0,150,250,350;
// ////       //
// ////       FKL = "_Partial_C1_ResVal_v2";
// ////       int PartialCounterfactual1_v2_Output = SolveSimulate_CounterfacturalResult_1_v2(theta_Est_S, shr_RealData,
// ////           Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////           Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////           sim_var, SimTbar, SimN, KsupplyElas, FKL, ResVal_vec, threadsManagement);
// ////       cout << "Counterfactual 1 v2: targeting exit rate for each segment: Varying ResidualVal done" << endl;
// ////       // throw runtime_error("239");
// //// // // // ////
// ////       // /*** Counterfactual 2 (Table): Varying firing costs ***/
// ////       // KsupplyElas = 0.0;
// ////       // FKL = "_Partial_C2_FiringCost_v1";
// ////       // N_C = 8;
// ////       // TargetExitRate_vec << Value_StatusQuo_GoodLaborInt.AvgExitRate,Value_StatusQuo_BadLaborInt.AvgExitRate,
// ////       //     Value_StatusQuo_GoodCapitalInt.AvgExitRate, Value_StatusQuo_BadCapitalInt.AvgExitRate,
// ////       //     0.045,0.054,0.063,0.072;
// ////       // cout << "TargetExitRate_vec = " << TargetExitRate_vec.transpose() << endl;
// ////       // int PartialCounterfactual2_v1_Output = SolveSimulate_CounterfacturalResult_2_v1(theta_Est_S, shr_RealData,
// ////       //     Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////       //     Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////       //     sim_var, SimTbar, SimN, KsupplyElas, FKL, TargetExitRate_vec, threadsManagement);
// ////       // cout << "Counterfactual 2 v1: Varying Firing Cost done" << endl;
// //// // // //
// //// //     // /*** Counterfactual 2 (Graph): Varying firing costs ***/
// //// //      KsupplyElas = 0.0;
// //// //      FKL = "_Partial_C2_FiringCost_v2";
// //// //      N_C = 9;
// //// //      ArrayXd firingshare1_vec = ArrayXd::LinSpaced(N_C, 0.01, 1);
// //// //      ArrayXd firingshare2_vec = ArrayXd::LinSpaced(N_C, 1.05, 5);
// //// //      ArrayXd firingshare_vec(N_C+N_C);
// //// //      firingshare_vec << firingshare1_vec,firingshare2_vec;
// //// // //
// //// //      int PartialCounterfactual2_v2_Output = SolveSimulate_CounterfacturalResult_2_v2(theta_Est_S, shr_RealData,
// //// //          Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //// //          Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //// //          sim_var, SimTbar, SimN, KsupplyElas, FKL, firingshare_vec, threadsManagement);
// //// //      cout << "Counterfactual 2 v2: Varying Firing Cost done" << endl;
// //// //      throw runtime_error("239");
// //// //
// ////     /*** Counterfactual 3: Change KsupplyElas and Change the exit barrier to match the US exit rate ***/
// ////     ArrayXd TargetExitRate_LaborInt_vec(1);
// ////     TargetExitRate_LaborInt_vec << 0.045;
// ////     ArrayXd TargetExitRate_CapitalInt_vec(1);
// ////     TargetExitRate_CapitalInt_vec << 0.045;
// //// // //
// //// //     int Counterfactual3_Output_Weighted = SolveSimulate_CounterfacturalResult_3_Weighted(theta_Est_S, shr_RealData,
// //// //         Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //// //         Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //// //         sim_var, SimTbar, SimN, TargetExitRate_LaborInt_vec, TargetExitRate_CapitalInt_vec, threadsManagement);
// ////     //
// ////     int Counterfactual3_Output_Aggregate = SolveSimulate_CounterfacturalResult_3_Aggregate(theta_Est_S, shr_RealData,
// ////         Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////         Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////         sim_var, SimTbar, SimN, TargetExitRate_LaborInt_vec, TargetExitRate_CapitalInt_vec, threadsManagement);
// ////      cout << "Counterfactual 3: Change KsupplyElas and Change the exit barrier to match the US exit rate done" << endl;
// ////      // throw runtime_error("265");
// //// //
// ////     // ArrayXd TargetExitRate_LaborInt_vec(1);
// ////     TargetExitRate_LaborInt_vec << 0.045;
// //// //    ArrayXd TargetExitRate_CapitalInt_vec(1);
// ////     TargetExitRate_CapitalInt_vec << 0.045;
// ////     /*** Counterfactual 4: Change KsupplyElas and Change firing costs to match the US exit rate ***/
// ////     // int Counterfactual4_Output_Weighted = SolveSimulate_CounterfacturalResult_4_Weighted(theta_Est_S, shr_RealData,
// ////     //     Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////     //     Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////     //     sim_var, SimTbar, SimN, TargetExitRate_LaborInt_vec, TargetExitRate_CapitalInt_vec, threadsManagement);
// ////     // cout << "Counterfactual 4: Change KsupplyElas and Change the firing costs to match the US exit rate done" << endl;
// ////
// ////     int Counterfactual4_Output_Aggregate = SolveSimulate_CounterfacturalResult_4_Aggregate(theta_Est_S, shr_RealData,
// ////     Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////     Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////     sim_var, SimTbar, SimN, TargetExitRate_LaborInt_vec, TargetExitRate_CapitalInt_vec, threadsManagement);
// ////     cout << "Counterfactual 4: Change KsupplyElas and Change the firing costs to match the US exit rate done" << endl;
// ////     throw runtime_error("273");
// //// // //
// //// //      // /*** Counterfactual 5: Varying both the exit barrier and firing costs at the same time ***/
// //// //      // KsupplyElas = 0.0;
// //// //      // double Residual_up = 500; double firingshare_low = 0.01;
// //// //      // N_C = 21;
// //// //      // int Counterfactual5_Output = SolveSimulate_CounterfacturalResult_5(theta_Est_S,shr_RealData,
// //// //      //     Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //// //      //     Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //// //      //     sim_var, SimTbar, SimN, KsupplyElas,
// //// //      //     Residual_up, firingshare_low, N_C, threadsManagement);
// //// //      // cout << "Counterfactual 5: Varying both the exit barrier and firing costs at the same time done" << endl;
// //// // //     // throw runtime_error("128");
// //// //     //
// //// //     // /*** Counterfactual 6: Entry Cost vs Exit Cost ***/
// //// //     // int Counterfactual6_Output = SolveSimulate_CounterfacturalResult_6(theta_Est_S,shr_RealData,
// //// //     //     Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //// //     //     Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //// //     //     sim_var, SimTbar, SimN, KsupplyElas, threadsManagement);
// //// //     // throw runtime_error("Counterfactual 6: Entry Cost vs Exit Cost");
// //// // // ////
// //     /**************************************************************
// //     * Counterfactual Exercise: General Equilibrium
// //     **************************************************************/
// //     /*** Counterfactual Graph: Varying the exit barrier ***/
// //     int N_C = 7;
// //     double KsupplyElas = 0.75;
// //     ArrayXd ResVal_vec(1);
// //     ResVal_vec << 300;
// //     //
// //     // int GECounterfactual1_Output = SolveSimulate_GECounterfacturalResult_1(theta_Est_S, shr_RealData,
// //     //     Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //     //     Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //     //     sim_var, SimTbar, SimN, KsupplyElas, ResVal_vec, threadsManagement);
// //     // cout << "Counterfactual 1: targeting exit rate for each segment: Varying ResidualVal done" << endl;
// //     // throw runtime_error("239");
// ////
// ////     FKL = "_GEPartial_C1_ResVal";
// ////     int GEPartialCounterfactual1_Output = SolveSimulate_CounterfacturalResult_1_v2(theta_Est_S, shr_RealData,
// ////         Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////         Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////         sim_var, SimTbar, SimN, KsupplyElas, FKL, ResVal_vec, threadsManagement);
// ////     cout << "Counterfactual 1 GE (partial version): targeting exit rate for each segment: Varying ResidualVal done" << endl;
// ////
// //     N_C = 5;
// //     KsupplyElas = 0.75;
// //     ArrayXd firingshare_vec(1);
// //     firingshare_vec << 0.5;
// //
// //     int GECounterfactual2_Output = SolveSimulate_GECounterfacturalResult_2(theta_Est_S, shr_RealData,
// //         Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //         Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //         sim_var, SimTbar, SimN, KsupplyElas, firingshare_vec, threadsManagement);
// //     cout << "Counterfactual 2: targeting exit rate for each segment: Varying firing cost done" << endl;
// //     // throw runtime_error("239");
// ////
// ////     FKL = "_GEPartial_C2_FiringCost";
// ////     int GEPartialCounterfactual2_Output = SolveSimulate_CounterfacturalResult_2_v2(theta_Est_S, shr_RealData,
// ////         Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// ////         Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////         sim_var, SimTbar, SimN, KsupplyElas, FKL, firingshare_vec, threadsManagement);
// ////     cout << "Counterfactual 2 GE (Partial version): targeting exit rate for each segment: Varying firing cost done" << endl;
// //
    return 0;
}
