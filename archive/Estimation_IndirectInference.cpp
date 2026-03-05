#include "Estimation_IndirectInference.h"

using namespace Eigen;
using namespace std;
using namespace alias;
////
/**************************************************************************
* In the second stage of estimate, we first simulate data with a guessed parameter.
* If we guessed it right, theta_a should be the maximizer of the auxilliary model with the guessed parameter.
* That is, the derivative of the maximum likelihood function of the auxiliary model with simulated data
* at theta_a should be as close to zero as possible.
**************************************************************************/
/**************************************************************
* Estimation Part 2: Second Stage model
**************************************************************/
ArrayXd alias::EstimationIndirectInference_Output(const ArrayXd & theta1_a, const ArrayXd & theta2_a,
    const ArrayXd & theta3_a, MultiThreads::Threads_Management & threadsManagement) {

    /*** Estimates from The Second stage of estimation: deflated ***/
//     set the estimates from the auxiliary model as the initial states.
    ArrayXd theta1_Para_ini(para.dim1);
    theta1_Para_ini << theta1_a;
    ArrayXd theta2_Para_ini(para.dim2);
    theta2_Para_ini << theta2_a;
    ArrayXd theta3_Para_ini(para.dim3);
    theta3_Para_ini << theta3_a;

    ArrayXd theta_a(para.dim1 + para.dim2 + para.dim3);
    theta_a << theta3_a,theta1_a,theta2_a;
    cout << "theta_a = " << theta_a.transpose() << endl;

    ArrayXd theta_Est(para.dim1+para.dim2+para.dim3);
    theta_Est << theta3_Para_ini,theta1_Para_ini,theta2_Para_ini;
    cout << "theta_Est = " << theta_Est.transpose() << endl;
//    throw runtime_error("147");

    /** In the second stage of estimate, we first simulate data with a guessed parameter.
     * If we guessed it right, theta_a should be the maximizer of the auxilliary model with the guessed parameter.
     * That is, the derivative of the maximum likelihood function of the auxiliary model with simulated data
     * at theta_a should be as close to zero as possible. **/

    /** Simulated Data: we take the initial period of each firm in the data as given, then simulate the panel **/
    // we can duplicate the firms in the real data to increase the size of simulated data */
    int GoodState = 1; int LaborIntensive = 1;
    SimData simdata_panel_GoodState_LaborInt = createRealData("_ForSim",GoodState,LaborIntensive);
    GoodState = 0; LaborIntensive = 1;
    SimData simdata_panel_BadState_LaborInt = createRealData("_ForSim",GoodState,LaborIntensive);
    GoodState = 1; LaborIntensive = 0;
    SimData simdata_panel_GoodState_CapitalInt = createRealData("_ForSim",GoodState,LaborIntensive);
    GoodState = 0; LaborIntensive = 0;
    SimData simdata_panel_BadState_CapitalInt = createRealData("_ForSim",GoodState,LaborIntensive);
    cout << "Data imported for simulation" << endl;
//    throw runtime_error("147");

    // to prepare for the calculation of the derivative at theta_a, we calculate the value/policy function at theta_a and theta_a + epsilon
    GoodState = 1; LaborIntensive = 1;
    EquV_Auxiliary EquV_a_GoodState_LaborInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 1;
    EquV_Auxiliary EquV_a_BadState_LaborInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    GoodState = 1; LaborIntensive = 0;
    EquV_Auxiliary EquV_a_GoodState_CapitalInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 0;
    EquV_Auxiliary EquV_a_BadState_CapitalInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    cout << "EquV_a done" << endl;

    /** simulate a series of random number **/
    int SimN_2step_GoodState_LaborInt = simdata_panel_GoodState_LaborInt.good.rows();
    SimVar sim_var_2step_GoodState_LaborInt = SimulationBase(SimN_2step_GoodState_LaborInt,para.SimTbar);
    int SimN_2step_BadState_LaborInt = simdata_panel_BadState_LaborInt.good.rows();
    SimVar sim_var_2step_BadState_LaborInt = SimulationBase(SimN_2step_BadState_LaborInt,para.SimTbar);

    int SimN_2step_GoodState_CapitalInt = simdata_panel_GoodState_CapitalInt.good.rows();
    SimVar sim_var_2step_GoodState_CapitalInt = SimulationBase(SimN_2step_GoodState_CapitalInt,para.SimTbar);
    int SimN_2step_BadState_CapitalInt = simdata_panel_BadState_CapitalInt.good.rows();
    SimVar sim_var_2step_BadState_CapitalInt = SimulationBase(SimN_2step_BadState_CapitalInt,para.SimTbar);
    cout << "Simulate a series of random number" << endl;

    theta_Est = EstimationIndirectInference(sim_var_2step_GoodState_LaborInt,
        sim_var_2step_BadState_LaborInt,sim_var_2step_GoodState_CapitalInt, sim_var_2step_BadState_CapitalInt,
        theta1_Para_ini,theta2_Para_ini,theta3_Para_ini,theta_a,
        EquV_a_GoodState_LaborInt,EquV_a_BadState_LaborInt,EquV_a_GoodState_CapitalInt,EquV_a_BadState_CapitalInt,
        simdata_panel_GoodState_LaborInt,simdata_panel_BadState_LaborInt,
        simdata_panel_GoodState_CapitalInt,simdata_panel_BadState_CapitalInt,
        threadsManagement);
    cout << "Estimates from The Second stage of estimation:" << endl;
    cout << "Group 1 parameters: theta1_Est = " << theta_Est(seqN(0,para.dim1)).transpose() << endl;
    cout << "Group 2 parameters: theta2_Est = " << theta_Est(seqN(para.dim1,para.dim2)).transpose() << endl;
    cout << "Group 3 parameters: theta2_Est = " << theta_Est(seqN(para.dim1+para.dim2,para.dim3)).transpose() << endl;

    return theta_Est;
}

EquV_Auxiliary alias::CalEquV_plus_minus(const ArrayXd & theta_a, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    // cout << "theta_a = " << theta_a.transpose() << endl;
    EquV_Auxiliary EquV_a;

    ArrayXd theta1_a(para.dim1);
    for (size_t i = 0; i < para.dim1; ++i) {theta1_a(i) = theta_a(para.dim3+i);}
    // cout << "theta1_a = " << theta1_a.transpose() << endl;

    ArrayXd theta2_a(para.dim2);
    for (size_t i = 0; i < para.dim2; ++i) {theta2_a(i) = theta_a(para.dim3+para.dim1+i);}
    // cout << "theta2_a = " << theta2_a.transpose() << endl;

    ArrayXd theta_a_point = theta_a;
    std::vector<double> theta3_a_point(para.dim3);
    for (size_t i = 0; i < para.dim3; ++i) {
        theta3_a_point[i] = theta_a_point(i);
    }
    // cout << "theta3_a = " << theta_a_point.transpose() << endl;

    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_point = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
    EquV_a.para_est_a_point = get<0>(t_Equ_point);
    EquV_a.para_vec_a_point = get<1>(t_Equ_point);
    EquV_a.EquV_a_point = get<2>(t_Equ_point);
    EquV_a.Evalmat_a_point = get<3>(t_Equ_point);
    // cout << "******* Prob_E = ";
    // for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV_a.EquV_a_point.ProbPD_E(i_phi*para.N_KL) << ", "; }
    // cout << ";" << endl;

    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
        //// value function of the auxilliary model
        ArrayXd theta_a_plus = theta_a;
        EquV_a.eps_diff(i_a) = abs(theta_a(i_a)) * 1e-8;
        theta_a_plus(i_a) = theta_a(i_a) + EquV_a.eps_diff(i_a);

        std::vector<double> theta3_a_plus(para.dim3);
        for (size_t i = 0; i < para.dim3; ++i) {
            theta3_a_plus[i] = theta_a_plus(i);
        }

        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_plus = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
            theta3_a_plus, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
        EquV_a.para_est_a_plus[i_a] = get<0>(t_Equ_plus);
        EquV_a.para_vec_a_plus[i_a] = get<1>(t_Equ_plus);
        EquV_a.EquV_a_plus[i_a] = get<2>(t_Equ_plus);
        EquV_a.Evalmat_a_plus[i_a] = get<3>(t_Equ_plus);

        // cout << "******* Prob_E = ";
        // for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV_a.EquV_a_plus[i_a].ProbPD_E(i_phi*para.N_KL) << ", "; }
        // cout << ";" << endl;
        // cout << "******* diff EquV_a.EquV_a_plus[i_a] = " << abs(EquV_a.Evalmat_a_plus[i_a].EVal_PD_K_mat-EquV_a.Evalmat_a_point.EVal_PD_K_mat).maxCoeff() << endl;
    }
    return EquV_a;
}

/**************************************************************************
* Estimation in the second stage
**************************************************************************/
ArrayXd alias::EstimationIndirectInference(const SimVar & sim_var_2step_GoodState_LaborInt,
    const SimVar & sim_var_2step_BadState_LaborInt,const SimVar & sim_var_2step_GoodState_CapitalInt,
    const SimVar & sim_var_2step_BadState_CapitalInt,
    const ArrayXd & theta1_Step2_ini, const ArrayXd & theta2_Step2_ini, const ArrayXd & theta3_Step2_ini,
    const ArrayXd & theta_a,
    const EquV_Auxiliary & EquV_a_GoodState_LaborInt, const EquV_Auxiliary & EquV_a_BadState_LaborInt,
    const EquV_Auxiliary & EquV_a_GoodState_CapitalInt, const EquV_Auxiliary & EquV_a_BadState_CapitalInt,
    const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    MultiThreads::Threads_Management & threadsManagement) {

    int dim = para.dim1+para.dim2+para.dim3;
    //// objective function
    auto f_derivative_free
        = [&sim_var_2step_GoodState_LaborInt,&sim_var_2step_BadState_LaborInt,
        &sim_var_2step_GoodState_CapitalInt,&sim_var_2step_BadState_CapitalInt,
        &theta_a,
        &EquV_a_GoodState_LaborInt,&EquV_a_BadState_LaborInt,&EquV_a_GoodState_CapitalInt,&EquV_a_BadState_CapitalInt,
        &simdata_panel_GoodState_LaborInt,&simdata_panel_BadState_LaborInt,
        &simdata_panel_GoodState_CapitalInt,&simdata_panel_BadState_CapitalInt,
        &threadsManagement] (std::vector<double> x){
            cout << "x = " << std::setprecision(12);
            for (size_t n = 0; n < para.dim1+para.dim2+para.dim3; ++n) {cout << x[n] << ", ";}
            cout << "; " << endl;

            double obj = EstimationIndirectInference_FullVersion_Obj(x,
                sim_var_2step_GoodState_LaborInt,sim_var_2step_BadState_LaborInt,
                sim_var_2step_GoodState_CapitalInt,sim_var_2step_BadState_CapitalInt,
                theta_a,
                EquV_a_GoodState_LaborInt, EquV_a_BadState_LaborInt,
                EquV_a_GoodState_CapitalInt, EquV_a_BadState_CapitalInt,
                simdata_panel_GoodState_LaborInt, simdata_panel_BadState_LaborInt,
                simdata_panel_GoodState_CapitalInt, simdata_panel_BadState_CapitalInt,
                threadsManagement);

            cout << "Iteration : obj = " << obj << endl;
            return obj;
        };

    std::vector<double> x0_vec(dim);
    for (size_t n = 0; n < para.dim3; ++n) { x0_vec[n] = theta3_Step2_ini(n);}
    for (size_t n = 0; n < para.dim1; ++n) { x0_vec[para.dim3+n] = theta1_Step2_ini(n);}
    for (size_t n = 0; n < para.dim2; ++n) { x0_vec[para.dim3+para.dim1+n] = theta2_Step2_ini(n);}

    int nn = 1;
    std::vector<double> lb(dim);
    lb[0] = -100;
    for (size_t k = nn+0; k < nn+para.dim3-6; ++k) { lb[k] = -100; }
    lb[para.dim3-6] = 0.01; lb[para.dim3-5] = 0.01; lb[para.dim3-4] = 0.01;
    lb[para.dim3-3] = -5000; lb[para.dim3-2] = -5000;
    lb[para.dim3-1] = 0.1;

    lb[para.dim3+0] = -5; lb[para.dim3+1] = -5; lb[para.dim3+2] = -0.98; lb[para.dim3+3] = 0.01;
    lb[para.dim3+para.dim1+0] = -100; lb[para.dim3+para.dim1+1] = -100; lb[para.dim3+para.dim1+2] = 0.1;


    std::vector<double> ub(dim);
    ub[0] = 5000;
    for (size_t k = nn+0; k < nn+para.dim3-6; ++k) { ub[k] = 5000; }
    ub[para.dim3-6] = 2; ub[para.dim3-5] = 2; ub[para.dim3-4] = 2;
    ub[para.dim3-3] = 5000; ub[para.dim3-2] = 5000;
    ub[para.dim3-1] = 5000;

    ub[para.dim3+0] = 5; ub[para.dim3+1] = 5; ub[para.dim3+2] = 0.98; ub[para.dim3+3] = 5;
//    double alpha_M = std::max(para.alpha_M_G_L,para.alpha_M_B_L);
//    ub[3] = (1-alpha_M) * (para.sigma-1)/para.sigma / (1 - (para.sigma-1)/para.sigma * alpha_M);
//    ub[4] = 1;
    ub[para.dim3+para.dim1+0] = 100; ub[para.dim3+para.dim1+1] = 100; ub[para.dim3+para.dim1+2] = 50;


    DFS::Optimization_Problem problem(dim, f_derivative_free);
    problem.set_x_lb(lb);
    problem.set_x_ub(ub);

    DFS::Direct_Search_Option option;
    option.display_level = 5;
    auto result = DFS::Coordinate_Direct_Search(problem, x0_vec, option);
    std::cout << "optimal_x = " ;
    for (size_t n = 0; n < dim; ++n) { cout << result.optimal_x[n] << ", "; }
    cout << ";" << endl;
    std::cout << "optimal f = " << result.optimal_f << std::endl;

    ArrayXd thetaPara(dim);
    for (size_t n = 0; n < dim; ++n) { thetaPara(n) = result.optimal_x[n]; }
    cout << ";" << endl;
    return thetaPara;
}

/**************************************************************************
* Estimation in the second stage
* Objective function calculation
**************************************************************************/
double alias::EstimationIndirectInference_FullVersion_Obj(std::vector<double> x_vec,
    const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
    const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
    const ArrayXd theta_a,
    const EquV_Auxiliary & EquV_a_GoodState_LaborInt, const EquV_Auxiliary & EquV_a_BadState_LaborInt,
    const EquV_Auxiliary & EquV_a_GoodState_CapitalInt, const EquV_Auxiliary & EquV_a_BadState_CapitalInt,
    const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd thetaPara(para.dim1 + para.dim2 + para.dim3);
    for (size_t i = 0; i < para.dim1 + para.dim2 + para.dim3; ++i) { thetaPara(i) = x_vec[i]; }

    cout << "thetaPara = " << thetaPara.transpose() << endl;
    /** simulate data for the second stage **/
    int GoodState = 1; int LaborIntensive = 1;
    SimData sim_data_Para_GoodState_LaborInt = SimulationData_Step2Estimation_FullVersion(thetaPara,
        sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);

    GoodState = 0; LaborIntensive = 1;
    SimData sim_data_Para_BadState_LaborInt = SimulationData_Step2Estimation_FullVersion(thetaPara,
        sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);

    GoodState = 1; LaborIntensive = 0;
    SimData sim_data_Para_GoodState_CapitalInt = SimulationData_Step2Estimation_FullVersion(thetaPara,
        sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);

    GoodState = 0; LaborIntensive = 0;
    SimData sim_data_Para_BadState_CapitalInt = SimulationData_Step2Estimation_FullVersion(thetaPara,
        sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    // throw runtime_error("279");

    /** Part 1 of the full model estimation **/
    std::vector<double> theta_a_part1(para.dim1);
    for (size_t n = 0; n < para.dim1; ++n) { theta_a_part1[n] = theta_a(para.dim3+n);}

    GoodState = 1; LaborIntensive = 1;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Good_LaborInt = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
        theta_a_part1,sim_data_Para_GoodState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Good_LaborInt = get<0>(t_Part1_Good_LaborInt);
    std::vector<double> grad_a_Part1_sum_Good_LaborInt = get<1>(t_Part1_Good_LaborInt);
    int Ncount_part1_Good_LaborInt = get<2>(t_Part1_Good_LaborInt);

    GoodState = 0; LaborIntensive = 1;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Bad_LaborInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(theta_a_part1,sim_data_Para_BadState_LaborInt,
        GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Bad_LaborInt = get<0>(t_Part1_Bad_LaborInt);
    std::vector<double> grad_a_Part1_sum_Bad_LaborInt = get<1>(t_Part1_Bad_LaborInt);
    int Ncount_part1_Bad_LaborInt = get<2>(t_Part1_Bad_LaborInt);
////
////    double Nobs = double(Ncount_part1_Good_LaborInt+Ncount_part1_Bad_LaborInt);
////    double y_part1 = (y_part1_sum_Good_LaborInt+y_part1_sum_Bad_LaborInt)/Nobs;
////    ArrayXd grad_a1(para.dim1);
////    grad_a1(0) = (grad_a_Part1_sum_Good_LaborInt[0]) / Nobs;
////    grad_a1(1) = (grad_a_Part1_sum_Bad_LaborInt[0]) / Nobs;
////    grad_a1(2) = (grad_a_Part1_sum_Good_LaborInt[1]+grad_a_Part1_sum_Bad_LaborInt[1]) / Nobs;
////    grad_a1(3) = (grad_a_Part1_sum_Good_LaborInt[2]+grad_a_Part1_sum_Bad_LaborInt[2]) / Nobs;
//
    GoodState = 1; LaborIntensive = 0;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Good_CapitalInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(theta_a_part1,sim_data_Para_GoodState_CapitalInt,
        GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Good_CapitalInt = get<0>(t_Part1_Good_CapitalInt);
    std::vector<double> grad_a_Part1_sum_Good_CapitalInt = get<1>(t_Part1_Good_CapitalInt);
    int Ncount_part1_Good_CapitalInt = get<2>(t_Part1_Good_CapitalInt);

    GoodState = 0; LaborIntensive = 0;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Bad_CapitalInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(theta_a_part1,sim_data_Para_BadState_CapitalInt,
        GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Bad_CapitalInt = get<0>(t_Part1_Bad_CapitalInt);
    std::vector<double> grad_a_Part1_sum_Bad_CapitalInt = get<1>(t_Part1_Bad_CapitalInt);
    int Ncount_part1_Bad_CapitalInt = get<2>(t_Part1_Bad_CapitalInt);

    double Nobs = double(Ncount_part1_Good_LaborInt+Ncount_part1_Bad_LaborInt+Ncount_part1_Good_CapitalInt+Ncount_part1_Bad_CapitalInt);
    double y_part1 = (y_part1_sum_Good_LaborInt+y_part1_sum_Bad_LaborInt+y_part1_sum_Good_CapitalInt+y_part1_sum_Bad_CapitalInt)
            /Nobs;
    ArrayXd grad_a1(para.dim1);
    grad_a1(0) = (grad_a_Part1_sum_Good_LaborInt[0]+grad_a_Part1_sum_Good_CapitalInt[0]) / Nobs;
    grad_a1(1) = (grad_a_Part1_sum_Bad_LaborInt[0]+grad_a_Part1_sum_Bad_CapitalInt[0])  / Nobs;
    grad_a1(2) = (grad_a_Part1_sum_Good_LaborInt[1]+grad_a_Part1_sum_Bad_LaborInt[1]
        +grad_a_Part1_sum_Good_CapitalInt[1]+grad_a_Part1_sum_Bad_CapitalInt[1]) / Nobs;
    grad_a1(3) = (grad_a_Part1_sum_Good_LaborInt[2]+grad_a_Part1_sum_Bad_LaborInt[2]
        +grad_a_Part1_sum_Good_CapitalInt[2]+grad_a_Part1_sum_Bad_CapitalInt[2]) / Nobs;
    cout << "y_part1 = " << y_part1 << "; grad_a1 = " << grad_a1.transpose() *  Nobs << "; sum = " << (grad_a1*grad_a1).sum() << endl;
//    throw runtime_error("Part 1 of the full model estimation done");
////    cout << "186 : Part 1 of the full model estimation done" << endl;
//
    /** Part 2 of the full model estimation **/
    std::vector<double> theta1(para.dim1);
    for (size_t k = 0; k < para.dim1; ++k) { theta1[k] = theta_a(para.dim3+k);}
    std::vector<double> theta_a_part2(para.dim2);
    for (size_t n = 0; n < para.dim2; ++n) { theta_a_part2[n] = theta_a(para.dim3+para.dim1+n);}

    GoodState = 1; LaborIntensive = 1;
    ParaEst1 para_est1_Good_LaborInt = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_LaborInt = Backout_phi_a_PD(sim_data_Para_GoodState_LaborInt.Revenue_P,
        sim_data_Para_GoodState_LaborInt.Capital_P, sim_data_Para_GoodState_LaborInt.Employ_ur_P,
        sim_data_Para_GoodState_LaborInt.Employ_uc_P,
        para_est1_Good_LaborInt.alpha_tilde_K, para_est1_Good_LaborInt.alpha_tilde_L,
        para_est1_Good_LaborInt.alpha_Lr, para_est1_Good_LaborInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Good_LaborInt
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_a_part2,
        sim_data_Para_GoodState_LaborInt,threadsManagement);
    double y_part2_sum_Good_LaborInt = get<0>(t_Part2_Good_LaborInt);
    std::vector<double> grad_a_Part2_sum_Good_LaborInt = get<1>(t_Part2_Good_LaborInt);
    int Ncount_part2_Good_LaborInt = get<2>(t_Part2_Good_LaborInt);

    GoodState = 0; LaborIntensive = 1;
    ParaEst1 para_est1_Bad_LaborInt = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_BadState_LaborInt = Backout_phi_a_PD(sim_data_Para_BadState_LaborInt.Revenue_P,
        sim_data_Para_BadState_LaborInt.Capital_P, sim_data_Para_BadState_LaborInt.Employ_ur_P,
        sim_data_Para_BadState_LaborInt.Employ_uc_P,
        para_est1_Bad_LaborInt.alpha_tilde_K, para_est1_Bad_LaborInt.alpha_tilde_L,
        para_est1_Bad_LaborInt.alpha_Lr, para_est1_Bad_LaborInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Bad_LaborInd
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_a_part2,
        sim_data_Para_BadState_LaborInt,threadsManagement);
    double y_part2_sum_Bad_LaborInt = get<0>(t_Part2_Bad_LaborInd);
    std::vector<double> grad_a_Part2_sum_Bad_LaborInt = get<1>(t_Part2_Bad_LaborInd);
    int Ncount_part2_Bad_LaborInt = get<2>(t_Part2_Bad_LaborInd);

////    double y_part2 = (y_part2_sum_Good_LaborInt + y_part2_sum_Bad_LaborInt)/double(Ncount_part2_Good_LaborInt+Ncount_part2_Bad_LaborInt);
////    ArrayXd grad_a2(para.dim2);
////    for (size_t n = 0; n < para.dim2; ++n) {
////    grad_a2(n) = (grad_a_Part2_sum_Good_LaborInt[n]+grad_a_Part2_sum_Bad_LaborInt[n])
////                 / double(Ncount_part2_Good_LaborInt+Ncount_part2_Bad_LaborInt);
////    }

    GoodState = 1; LaborIntensive = 0;
    ParaEst1 para_est1_Good_CapitalInt = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_CapitalInt = Backout_phi_a_PD(sim_data_Para_GoodState_CapitalInt.Revenue_P,
        sim_data_Para_GoodState_CapitalInt.Capital_P, sim_data_Para_GoodState_CapitalInt.Employ_ur_P,
        sim_data_Para_GoodState_CapitalInt.Employ_uc_P,
        para_est1_Good_CapitalInt.alpha_tilde_K, para_est1_Good_CapitalInt.alpha_tilde_L,
        para_est1_Good_CapitalInt.alpha_Lr, para_est1_Good_CapitalInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Good_CapitalInd
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_a_part2,
        sim_data_Para_GoodState_CapitalInt,threadsManagement);
    double y_part2_sum_Good_CapitalInt = get<0>(t_Part2_Good_CapitalInd);
    std::vector<double> grad_a_Part2_sum_Good_CapitalInt = get<1>(t_Part2_Good_CapitalInd);
    int Ncount_part2_Good_CapitalInt = get<2>(t_Part2_Good_CapitalInd);

    GoodState = 0; LaborIntensive = 0;
    ParaEst1 para_est1_Bad_CapitalInt = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_BadState_CapitalInt = Backout_phi_a_PD(sim_data_Para_BadState_CapitalInt.Revenue_P,
        sim_data_Para_BadState_CapitalInt.Capital_P, sim_data_Para_BadState_CapitalInt.Employ_ur_P,
        sim_data_Para_BadState_CapitalInt.Employ_uc_P,
        para_est1_Bad_CapitalInt.alpha_tilde_K, para_est1_Bad_CapitalInt.alpha_tilde_L,
        para_est1_Bad_CapitalInt.alpha_Lr, para_est1_Bad_CapitalInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Bad_CapitalInd
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_a_part2,
        sim_data_Para_BadState_CapitalInt,threadsManagement);
    double y_part2_sum_Bad_CapitalInt = get<0>(t_Part2_Bad_CapitalInd);
    std::vector<double> grad_a_Part2_sum_Bad_CapitalInt = get<1>(t_Part2_Bad_CapitalInd);
    int Ncount_part2_Bad_CapitalInt = get<2>(t_Part2_Bad_CapitalInd);

    double y_part2 = (y_part2_sum_Good_LaborInt + y_part2_sum_Bad_LaborInt + y_part2_sum_Good_CapitalInt + y_part2_sum_Bad_CapitalInt)
        / double(Ncount_part2_Good_LaborInt+Ncount_part2_Bad_LaborInt+Ncount_part2_Good_CapitalInt+Ncount_part2_Bad_CapitalInt);
    ArrayXd grad_a2(para.dim2);
    for (size_t n = 0; n < para.dim2; ++n) {
        grad_a2(n) = (grad_a_Part2_sum_Good_LaborInt[n]+grad_a_Part2_sum_Bad_LaborInt[n]
                +grad_a_Part2_sum_Good_CapitalInt[n]+grad_a_Part2_sum_Bad_CapitalInt[n])
            / double(Ncount_part2_Good_LaborInt+Ncount_part2_Bad_LaborInt+Ncount_part2_Good_CapitalInt+Ncount_part2_Bad_CapitalInt);
    }
    cout << "y_part2 = " << y_part2 << "; grad_a2 = "
        << grad_a2.transpose() * double(Ncount_part2_Good_LaborInt+Ncount_part2_Bad_LaborInt+Ncount_part2_Good_CapitalInt+Ncount_part2_Bad_CapitalInt) << "; sum = " << (grad_a2*grad_a2).sum() << endl;
//    throw runtime_error("Part 2 of the full model estimation done");

    /** Part 3 of the full model estimation **/
    GoodState = 1; LaborIntensive = 1;
    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_point_Good_LaborInt
        = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_GoodState_LaborInt,
        EquV_a_GoodState_LaborInt.para_est_a_point,EquV_a_GoodState_LaborInt.para_vec_a_point,
        EquV_a_GoodState_LaborInt.EquV_a_point,EquV_a_GoodState_LaborInt.Evalmat_a_point,
        GoodState, LaborIntensive,threadsManagement);
    double obj_point_sum_Good_LaborInt = get<0>(t_point_Good_LaborInt);
    int Ncount_part3_Good_LaborInt = get<1>(t_point_Good_LaborInt);
    ArrayXXi K_index_max_P_mat_Good_LaborInt = get<4>(t_point_Good_LaborInt);
    ArrayXXi Lur_index_max_P_mat_Good_LaborInt = get<5>(t_point_Good_LaborInt);
    ArrayXXi Luc_index_max_P_mat_Good_LaborInt = get<6>(t_point_Good_LaborInt);
    ArrayXXd K_max_P_mat_Good_LaborInt = get<7>(t_point_Good_LaborInt);
    ArrayXXd Lur_max_P_mat_Good_LaborInt = get<8>(t_point_Good_LaborInt);
    ArrayXXd Luc_max_P_mat_Good_LaborInt = get<9>(t_point_Good_LaborInt);
    ArrayXXi K_index_max_D_mat_Good_LaborInt = get<10>(t_point_Good_LaborInt);
    ArrayXXi Lur_index_max_D_mat_Good_LaborInt = get<11>(t_point_Good_LaborInt);
    ArrayXXi Luc_index_max_D_mat_Good_LaborInt = get<12>(t_point_Good_LaborInt);
    ArrayXXd K_max_D_mat_Good_LaborInt = get<13>(t_point_Good_LaborInt);
    ArrayXXd Lur_max_D_mat_Good_LaborInt = get<14>(t_point_Good_LaborInt);
    ArrayXXd Luc_max_D_mat_Good_LaborInt = get<15>(t_point_Good_LaborInt);
//    throw runtime_error(" Part 3 of the full model estimation ");

    ArrayXd grad_a3_sum_Good_LaborInt = ArrayXd::Zero(para.dim3);
    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
//        cout << "i_a = " << i_a << endl;
        double obj_plus_sum_Good_LaborInt = 0.0;
        if (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15) {
            obj_plus_sum_Good_LaborInt = 0.0;
            grad_a3_sum_Good_LaborInt(i_a) = 0.0;
        }
        else {
            tuple<double, ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(
                sim_data_Para_GoodState_LaborInt,
                EquV_a_GoodState_LaborInt.para_est_a_plus[i_a], EquV_a_GoodState_LaborInt.para_vec_a_plus[i_a],
                EquV_a_GoodState_LaborInt.EquV_a_plus[i_a], EquV_a_GoodState_LaborInt.Evalmat_a_plus[i_a],
                K_index_max_P_mat_Good_LaborInt, Lur_index_max_P_mat_Good_LaborInt, Luc_index_max_P_mat_Good_LaborInt,
                K_max_P_mat_Good_LaborInt, Lur_max_P_mat_Good_LaborInt, Luc_max_P_mat_Good_LaborInt,
                K_index_max_D_mat_Good_LaborInt, Lur_index_max_D_mat_Good_LaborInt, Luc_index_max_D_mat_Good_LaborInt,
                K_max_D_mat_Good_LaborInt, Lur_max_D_mat_Good_LaborInt, Luc_max_D_mat_Good_LaborInt,
                GoodState, LaborIntensive, threadsManagement);
            obj_plus_sum_Good_LaborInt = get<0>(t_plus);
            grad_a3_sum_Good_LaborInt(i_a) = (obj_plus_sum_Good_LaborInt - obj_point_sum_Good_LaborInt)
                    / EquV_a_GoodState_LaborInt.eps_diff(i_a);
        }
    }
    cout << "grad_a3_sum_Good_LaborInt = " << grad_a3_sum_Good_LaborInt.transpose() << endl;


    GoodState = 0; LaborIntensive = 1;
    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_point_Bad_LaborInt
        = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(
        sim_data_Para_BadState_LaborInt,EquV_a_BadState_LaborInt.para_est_a_point,
        EquV_a_BadState_LaborInt.para_vec_a_point,EquV_a_BadState_LaborInt.EquV_a_point,
        EquV_a_BadState_LaborInt.Evalmat_a_point,GoodState, LaborIntensive,threadsManagement);
    double obj_point_sum_Bad_LaborInt = get<0>(t_point_Bad_LaborInt);
    int Ncount_part3_Bad_LaborInt = get<1>(t_point_Bad_LaborInt);
    ArrayXXi K_index_max_P_mat_Bad_LaborInt = get<4>(t_point_Bad_LaborInt);
    ArrayXXi Lur_index_max_P_mat_Bad_LaborInt = get<5>(t_point_Bad_LaborInt);
    ArrayXXi Luc_index_max_P_mat_Bad_LaborInt = get<6>(t_point_Bad_LaborInt);
    ArrayXXd K_max_P_mat_Bad_LaborInt = get<7>(t_point_Bad_LaborInt);
    ArrayXXd Lur_max_P_mat_Bad_LaborInt = get<8>(t_point_Bad_LaborInt);
    ArrayXXd Luc_max_P_mat_Bad_LaborInt = get<9>(t_point_Bad_LaborInt);
    ArrayXXi K_index_max_D_mat_Bad_LaborInt = get<10>(t_point_Bad_LaborInt);
    ArrayXXi Lur_index_max_D_mat_Bad_LaborInt = get<11>(t_point_Bad_LaborInt);
    ArrayXXi Luc_index_max_D_mat_Bad_LaborInt = get<12>(t_point_Bad_LaborInt);
    ArrayXXd K_max_D_mat_Bad_LaborInt = get<13>(t_point_Bad_LaborInt);
    ArrayXXd Lur_max_D_mat_Bad_LaborInt = get<14>(t_point_Bad_LaborInt);
    ArrayXXd Luc_max_D_mat_Bad_LaborInt = get<15>(t_point_Bad_LaborInt);
//    throw runtime_error(" Part 3 of the full model estimation ");

    ArrayXd grad_a3_sum_Bad_LaborInt = ArrayXd::Zero(para.dim3);
    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
//        cout << "i_a = " << i_a << endl;
        double obj_plus_sum_Bad_LaborInt;
        if (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14) {
//            cout << "thetaPara = " << thetaPara(i_a) << endl;
            obj_plus_sum_Bad_LaborInt = 0.0;
            grad_a3_sum_Bad_LaborInt(i_a) = 0.0;
        }
        else {
            tuple<double, ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(
                sim_data_Para_BadState_LaborInt, EquV_a_BadState_LaborInt.para_est_a_plus[i_a],
                EquV_a_BadState_LaborInt.para_vec_a_plus[i_a], EquV_a_BadState_LaborInt.EquV_a_plus[i_a],
                EquV_a_BadState_LaborInt.Evalmat_a_plus[i_a],
                K_index_max_P_mat_Bad_LaborInt, Lur_index_max_P_mat_Bad_LaborInt, Luc_index_max_P_mat_Bad_LaborInt,
                K_max_P_mat_Bad_LaborInt, Lur_max_P_mat_Bad_LaborInt, Luc_max_P_mat_Bad_LaborInt,
                K_index_max_D_mat_Bad_LaborInt, Lur_index_max_D_mat_Bad_LaborInt, Luc_index_max_D_mat_Bad_LaborInt,
                K_max_D_mat_Bad_LaborInt, Lur_max_D_mat_Bad_LaborInt, Luc_max_D_mat_Bad_LaborInt,
                GoodState, LaborIntensive, threadsManagement);
            obj_plus_sum_Bad_LaborInt = get<0>(t_plus);
            grad_a3_sum_Bad_LaborInt(i_a) = (obj_plus_sum_Bad_LaborInt - obj_point_sum_Bad_LaborInt)
                    / EquV_a_BadState_LaborInt.eps_diff(i_a);
        }
    }
    cout << "grad_a3_sum_Bad_LaborInt = " << grad_a3_sum_Bad_LaborInt.transpose() << endl;

//    double obj_point = (obj_point_sum_Good_LaborInt + obj_point_sum_Bad_LaborInt)
//            / double(Ncount_part3_Good_LaborInt + Ncount_part3_Bad_LaborInt);
//    ArrayXd grad_a3(para.dim2+para.dim3);
//    for (size_t n = 0; n < para.dim2+para.dim3; ++n) {
//        grad_a3(n) = (grad_a3_sum_Good_LaborInt[n]+grad_a3_sum_Bad_LaborInt[n])
//            / double(Ncount_part3_Good_LaborInt + Ncount_part3_Bad_LaborInt);
//    }

    GoodState = 1; LaborIntensive = 0;
    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_point_Good_CapitalInt
        = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_GoodState_CapitalInt,
        EquV_a_GoodState_CapitalInt.para_est_a_point,EquV_a_GoodState_CapitalInt.para_vec_a_point,
        EquV_a_GoodState_CapitalInt.EquV_a_point,EquV_a_GoodState_CapitalInt.Evalmat_a_point,
        GoodState, LaborIntensive,threadsManagement);
    double obj_point_sum_Good_CapitalInt = get<0>(t_point_Good_CapitalInt);
    int Ncount_part3_Good_CapitalInt = get<1>(t_point_Good_CapitalInt);
    ArrayXXi K_index_max_P_mat_Good_CapitalInt = get<4>(t_point_Good_CapitalInt);
    ArrayXXi Lur_index_max_P_mat_Good_CapitalInt = get<5>(t_point_Good_CapitalInt);
    ArrayXXi Luc_index_max_P_mat_Good_CapitalInt = get<6>(t_point_Good_CapitalInt);
    ArrayXXd K_max_P_mat_Good_CapitalInt = get<7>(t_point_Good_CapitalInt);
    ArrayXXd Lur_max_P_mat_Good_CapitalInt = get<8>(t_point_Good_CapitalInt);
    ArrayXXd Luc_max_P_mat_Good_CapitalInt = get<9>(t_point_Good_CapitalInt);
    ArrayXXi K_index_max_D_mat_Good_CapitalInt = get<10>(t_point_Good_CapitalInt);
    ArrayXXi Lur_index_max_D_mat_Good_CapitalInt = get<11>(t_point_Good_CapitalInt);
    ArrayXXi Luc_index_max_D_mat_Good_CapitalInt = get<12>(t_point_Good_CapitalInt);
    ArrayXXd K_max_D_mat_Good_CapitalInt = get<13>(t_point_Good_CapitalInt);
    ArrayXXd Lur_max_D_mat_Good_CapitalInt = get<14>(t_point_Good_CapitalInt);
    ArrayXXd Luc_max_D_mat_Good_CapitalInt = get<15>(t_point_Good_CapitalInt);
//    throw runtime_error(" Part 3 of the full model estimation ");

    ArrayXd grad_a3_sum_Good_CapitalInt = ArrayXd::Zero(para.dim3);
    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
//        cout << "i_a = " << i_a << endl;
        double obj_plus_sum_Good_CapitalInt;
        if (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15) {
            obj_plus_sum_Good_CapitalInt = 0.0;
            grad_a3_sum_Good_CapitalInt(i_a) = 0.0;
        }
        else {
            tuple<double, ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(
                sim_data_Para_GoodState_CapitalInt, EquV_a_GoodState_CapitalInt.para_est_a_plus[i_a],
                EquV_a_GoodState_CapitalInt.para_vec_a_plus[i_a], EquV_a_GoodState_CapitalInt.EquV_a_plus[i_a],
                EquV_a_GoodState_CapitalInt.Evalmat_a_plus[i_a],
                K_index_max_P_mat_Good_CapitalInt, Lur_index_max_P_mat_Good_CapitalInt, Luc_index_max_P_mat_Good_CapitalInt,
                K_max_P_mat_Good_CapitalInt, Lur_max_P_mat_Good_CapitalInt, Luc_max_P_mat_Good_CapitalInt,
                K_index_max_D_mat_Good_CapitalInt, Lur_index_max_D_mat_Good_CapitalInt, Luc_index_max_D_mat_Good_CapitalInt,
                K_max_D_mat_Good_CapitalInt, Lur_max_D_mat_Good_CapitalInt, Luc_max_D_mat_Good_CapitalInt,
                GoodState, LaborIntensive, threadsManagement);
            obj_plus_sum_Good_CapitalInt = get<0>(t_plus);
            grad_a3_sum_Good_CapitalInt(i_a) = (obj_plus_sum_Good_CapitalInt - obj_point_sum_Good_CapitalInt)
                    / EquV_a_GoodState_CapitalInt.eps_diff(i_a);
        }
    }
    cout << "grad_a3_sum_Good_CapitalInt = " << grad_a3_sum_Good_CapitalInt.transpose() << endl;

    GoodState = 0; LaborIntensive = 0;
    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_point_Bad_CapitalInt
        = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_BadState_CapitalInt,
        EquV_a_BadState_CapitalInt.para_est_a_point,EquV_a_BadState_CapitalInt.para_vec_a_point,
        EquV_a_BadState_CapitalInt.EquV_a_point,EquV_a_BadState_CapitalInt.Evalmat_a_point,
        GoodState, LaborIntensive,threadsManagement);
    double obj_point_sum_Bad_CapitalInt = get<0>(t_point_Bad_CapitalInt);
    int Ncount_part3_Bad_CapitalInt = get<1>(t_point_Bad_CapitalInt);
    ArrayXXi K_index_max_P_mat_Bad_CapitalInt = get<4>(t_point_Bad_CapitalInt);
    ArrayXXi Lur_index_max_P_mat_Bad_CapitalInt = get<5>(t_point_Bad_CapitalInt);
    ArrayXXi Luc_index_max_P_mat_Bad_CapitalInt = get<6>(t_point_Bad_CapitalInt);
    ArrayXXd K_max_P_mat_Bad_CapitalInt = get<7>(t_point_Bad_CapitalInt);
    ArrayXXd Lur_max_P_mat_Bad_CapitalInt = get<8>(t_point_Bad_CapitalInt);
    ArrayXXd Luc_max_P_mat_Bad_CapitalInt = get<9>(t_point_Bad_CapitalInt);
    ArrayXXi K_index_max_D_mat_Bad_CapitalInt = get<10>(t_point_Bad_CapitalInt);
    ArrayXXi Lur_index_max_D_mat_Bad_CapitalInt = get<11>(t_point_Bad_CapitalInt);
    ArrayXXi Luc_index_max_D_mat_Bad_CapitalInt = get<12>(t_point_Bad_CapitalInt);
    ArrayXXd K_max_D_mat_Bad_CapitalInt = get<13>(t_point_Bad_CapitalInt);
    ArrayXXd Lur_max_D_mat_Bad_CapitalInt = get<14>(t_point_Bad_CapitalInt);
    ArrayXXd Luc_max_D_mat_Bad_CapitalInt = get<15>(t_point_Bad_CapitalInt);
//    throw runtime_error(" Part 3 of the full model estimation ");

    ArrayXd grad_a3_sum_Bad_CapitalInt = ArrayXd::Zero(para.dim3);
    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
//        cout << "i_a = " << i_a << endl;
        double obj_plus_sum_Bad_CapitalInt;
        if (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14) {
            obj_plus_sum_Bad_CapitalInt = 0.0;
            grad_a3_sum_Bad_CapitalInt(i_a) = 0.0;
        }
        else {
            tuple<double, ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(
                sim_data_Para_BadState_CapitalInt, EquV_a_BadState_CapitalInt.para_est_a_plus[i_a],
                EquV_a_BadState_CapitalInt.para_vec_a_plus[i_a], EquV_a_BadState_CapitalInt.EquV_a_plus[i_a],
                EquV_a_BadState_CapitalInt.Evalmat_a_plus[i_a],
                K_index_max_P_mat_Bad_CapitalInt, Lur_index_max_P_mat_Bad_CapitalInt, Luc_index_max_P_mat_Bad_CapitalInt,
                K_max_P_mat_Bad_CapitalInt, Lur_max_P_mat_Bad_CapitalInt, Luc_max_P_mat_Bad_CapitalInt,
                K_index_max_D_mat_Bad_CapitalInt, Lur_index_max_D_mat_Bad_CapitalInt, Luc_index_max_D_mat_Bad_CapitalInt,
                K_max_D_mat_Bad_CapitalInt, Lur_max_D_mat_Bad_CapitalInt, Luc_max_D_mat_Bad_CapitalInt,
                GoodState, LaborIntensive, threadsManagement);
            obj_plus_sum_Bad_CapitalInt = get<0>(t_plus);
            grad_a3_sum_Bad_CapitalInt(i_a) = (obj_plus_sum_Bad_CapitalInt - obj_point_sum_Bad_CapitalInt)
                / EquV_a_BadState_CapitalInt.eps_diff(i_a);
        }
    }
    cout << "grad_a3_sum_Bad_CapitalInt = " << grad_a3_sum_Bad_CapitalInt.transpose() << endl;


    double obj_point = (obj_point_sum_Good_LaborInt + obj_point_sum_Bad_LaborInt
            + obj_point_sum_Good_CapitalInt + obj_point_sum_Bad_CapitalInt)
            / double(Ncount_part3_Good_LaborInt + Ncount_part3_Bad_LaborInt
            + Ncount_part3_Good_CapitalInt + Ncount_part3_Bad_CapitalInt);
    ArrayXd grad_a3(para.dim3);
    for (size_t n = 0; n < para.dim3; ++n) {
        grad_a3(n) = (grad_a3_sum_Good_LaborInt[n]+grad_a3_sum_Bad_LaborInt[n]
            +grad_a3_sum_Good_CapitalInt[n]+grad_a3_sum_Bad_CapitalInt[n])
            / double(Ncount_part3_Good_LaborInt + Ncount_part3_Bad_LaborInt
            + Ncount_part3_Good_CapitalInt + Ncount_part3_Bad_CapitalInt);
    }
    cout << "y_part3 = " << obj_point << "; grad_a3 = " << grad_a3.transpose() * double(Ncount_part3_Good_LaborInt + Ncount_part3_Bad_LaborInt
            + Ncount_part3_Good_CapitalInt + Ncount_part3_Bad_CapitalInt) << "; sum = " << (grad_a3*grad_a3).sum() << endl;

    double Obj = (grad_a1*grad_a1).sum() + (grad_a2*grad_a2).sum() + (grad_a3*grad_a3).sum();
    cout << "############ Obj = " << Obj << endl;
    // throw runtime_error("453");

    return Obj;
}

/**************************************************************************
* Simulate the data for the second stage
**************************************************************************/
/** the function **/
SimData alias::SimulationData_Step2Estimation_FullVersion(const ArrayXd & thetaData, const SimVar & sim_var,
    const double & p_K, const SimData & SimDataPanel, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    ParaEst para_est_Data = constructParaEst(thetaData,GoodState,LaborIntensive);
    para_est_Data.PI = 1.0;
    ParaVec para_vec_Data = constructParaFull(para_est_Data, p_K);
//    cout << "para_est = " << "F_P = " << para_est_Data.F_P_PP << ", " << para_est_Data.F_P_DP << ", " << para_est_Data.sigma_PD << ", "
//        << "Capital adj: " << para_est_Data.H_K_Cons << ", " << para_est_Data.H_K << ", "
//        << para_est_Data.c_HK << "," << para_est_Data.F_K_Cons << "," << para_est_Data.F_K << "," << para_est_Data.c_FK << ","
//        << "; RegWorker adj: " << para_est_Data.H_ur_Cons << "," << para_est_Data.H_ur << "," << para_est_Data.c_H_ur << ","
//        << para_est_Data.F_high_ur_Cons << "," << para_est_Data.F_high_ur << "," << para_est_Data.c_high_F_ur << ","
//        << para_est_Data.F_low_ur_Cons << "," << para_est_Data.F_low_ur << "," << para_est_Data.c_low_F_ur << ","
//        << "; NonRegWorker adj: " << para_est_Data.H_uc_Cons << "," << para_est_Data.H_uc << "," << para_est_Data.c_H_uc << ","
//        << para_est_Data.F_uc_Cons << "," << para_est_Data.F_uc << "," << para_est_Data.c_F_uc << ","
//        << "; Sigma adj: " << para_est_Data.sigma_Kerror << "," << para_est_Data.sigma_Lerror_ur << "," << para_est_Data.sigma_Lerror_uc << ","
//        << "; Exit Cost: " << para_est_Data.F_E << "," << para_est_Data.F_E_c_FK << "," << para_est_Data.F_E_c_F_uc << ","
//        << para_est_Data.F_E_c_F_ur << "," << para_est_Data.sigma_SE << endl;

//    cout << "start to solve the value function" << endl;
    ArrayXd VEnd = ArrayXd::Zero(para.N*2);
    tuple<EquStateV,EquStateVmat> t_Data = solveV(para_est_Data,para_vec_Data,VEnd,threadsManagement);
    EquStateV EquV_Data = get<0>(t_Data);
    EquStateVmat Evalmat = get<1>(t_Data);
//    cout << "Good States:: the value function solved" << endl;

    cout << "******* Eval = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV_Data.EVal_PD(i_phi*para.N_KL) << "; "; }
    cout << ";" << endl;
    cout << "******* Good States: Prob_E = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV_Data.ProbPD_E(i_phi*para.N_KL) << "; "; }
    cout << ";" << endl;
    cout << "******* Good States Production: Prob_D = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV_Data.ProbPD_D(i_phi*para.N_KL + 5*para.N_Lur*para.N_Luc + 5*para.N_Lur + 5) << "; "; }
    cout << ";" << endl;
    cout << "******* Good States Dormancy: Prob_D = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV_Data.ProbPD_D(para.N + i_phi*para.N_KL + 5*para.N_Lur*para.N_Luc + 5*para.N_Lur + 5) << "; "; }
    cout << ";" << endl;
//    throw runtime_error("419");
//
    //// simulate full data
//    cout << "Start simulation without missing 421" << endl;
    int SimN = SimDataPanel.good.rows();
    SimData data_panel_full = Simulation_Step2Estimation_FullVersion(para_est_Data,para_vec_Data,
        EquV_Data,Evalmat,sim_var,SimDataPanel,SimN,para.SimTbar,GoodState,LaborIntensive,
        threadsManagement);
//
//    ArrayXd avgExitProb = data_panel_full.Prob_E.cast<double>().colwise().sum() / data_panel_full.Prob_S.cast<double>().colwise().sum();
//    cout << "avgExitProb = " << avgExitProb.transpose() << endl;
////
//    ArrayXXd tempP = (data_panel_full.ProbP_P(all,seqN(0,para.SimTbar-1)) > 0).cast<double>();
//    ArrayXd avgProd_preProd =
//        ( data_panel_full.ProbP_P(all,seqN(0,para.SimTbar-1)) * tempP ).cast<double>().colwise().sum()
//        / tempP.cast<double>().colwise().sum();
//    cout << "avgProd_preProd = " << avgProd_preProd.transpose() << endl;
//
//    ArrayXXd tempD = (data_panel_full.ProbD_P(all,seqN(0,para.SimTbar-1)) > 0.00001).cast<double>();
//    ArrayXd avgProd_preDorm =
//        ( data_panel_full.ProbD_P(all,seqN(0,para.SimTbar-1)) * tempD ).cast<double>().colwise().sum()
//        / tempD.cast<double>().colwise().sum();
//    cout << "avgProd_preDorm = " << avgProd_preDorm.transpose() << endl;
//    cout << "Done simulation without missing 433" << endl;
////    throw runtime_error("431");
//
    //// simulate the missing pattern
//    cout << "Start simulation with missing 435" << endl;
    SimData data_panel = Simulation_Missing_Step2Estimation_FullVersion(sim_var,data_panel_full,
        SimN,para.SimTbar,GoodState,LaborIntensive,threadsManagement);
//    ArrayXd avgProb_P_nonMiss = data_panel.Prob_P_nonMiss.cast<double>().colwise().mean();
//    cout << "avgProb_P_nonMiss = " << avgProb_P_nonMiss.transpose() << endl;
//    cout << "Done simulation with missing 440" << endl;
    throw runtime_error("456");

    return data_panel;
}

/** within the function: simulate the full data **/
SimData alias::Simulation_Step2Estimation_FullVersion(const ParaEst & para_est, const ParaVec & para_vec,
    const EquStateV & EquV, const EquStateVmat & Evalmat, const SimVar & sim_var, const SimData & SimDataPanel,
    const int & SimN, const int & TBar, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    double alpha_M;
    if (GoodState == 1){alpha_M = para.alpha_M_G_L;}
    else {alpha_M = para.alpha_M_B_L;}

    SimData sim_data = InitializeSimData(SimN,TBar);
    sim_data.CapitalMiss = ArrayXXi::Zero(SimN,TBar);
    sim_data.Employ_urMiss = ArrayXXi::Zero(SimN,TBar);
    sim_data.Employ_ucMiss = ArrayXXi::Zero(SimN,TBar);

    sim_data.phi.col(0) = exp(sim_var.lnphi_sim.col(0) * para_vec.prod_sigma0 + para_vec.prod_mu0);
//    auto worker = [&](size_t n, unsigned thread_id) {
//    size_t thread_id = 0;
    for (size_t n = 0; n < SimN; ++n) {
//        cout << "n = " << n << endl;
        sim_data.good(n) = SimDataPanel.good(n);
        sim_data.labour_intensive_ind(n) = SimDataPanel.labour_intensive_ind(n);
        sim_data.first_year_in_Data(n) = SimDataPanel.first_year_in_Data(n);
        sim_data.last_year_in_Data(n) = SimDataPanel.last_year_in_Data(n);
        sim_data.state_code(n) = SimDataPanel.state_code(n);
        sim_data.industry(n) = SimDataPanel.industry(n);

        sim_data.last_year_in_Data(n) = 0;
        for (size_t t = sim_data.first_year_in_Data(n); t <= TBar-1; ++t) {
//            cout << "n = " << n << "; t = " << t << endl;
            if (t == sim_data.first_year_in_Data(n)) {
                sim_data.phi(n, t) = exp(sim_var.lnphi_sim(n, t) * para_vec.prod_sigma0 + para_vec.prod_mu0);

                sim_data.Capital_P(n, t) = SimDataPanel.Capital(n, t);
                sim_data.Employ_ur_P(n, t) = SimDataPanel.Employ_ur(n, t);
                sim_data.Employ_uc_P(n, t) = SimDataPanel.Employ_uc(n, t);

                sim_data.Capital_D(n, t) = SimDataPanel.Capital(n, t);
                sim_data.Employ_ur_D(n, t) = SimDataPanel.Employ_ur(n, t);
                sim_data.Employ_uc_D(n, t) = SimDataPanel.Employ_uc(n, t);

                sim_data.phi(n, t+1) = exp(para_est.gamma0 + para_est.gamma1 * log(sim_data.phi(n, t))
                    + sim_var.lnphi_sim(n, t+1) * para_est.sigma_phi_eps);

                sim_data.Revenue_P(n, t) = RevOpt1_sim_Realized(para_est,sim_data.phi(n, t),
                    sim_data.Capital_P(n, t),sim_data.Employ_ur_P(n, t),
                    sim_data.Employ_uc_P(n, t));

                sim_data.Revenue_D(n, t) = RevOpt1_sim_Realized(para_est,sim_data.phi(n, t),
                    sim_data.Capital_D(n, t),sim_data.Employ_ur_D(n, t),
                    sim_data.Employ_uc_D(n, t));

                sim_data.Price_P(n, t) = calPrice(para_est,alpha_M,sim_data.phi(n, t),
                    sim_data.Capital_P(n, t),sim_data.Employ_ur_P(n, t),
                    sim_data.Employ_uc_P(n, t), 1.0);

                sim_data.Price_D(n, t) = calPrice(para_est,alpha_M,sim_data.phi(n, t),
                    sim_data.Capital_D(n, t),sim_data.Employ_ur_D(n, t),
                    sim_data.Employ_uc_D(n, t), 1.0);

                sim_data.Prob_P(n, t) = 1.0;
                sim_data.Prob_D(n,t) = 0.0;

                sim_data.Prob_S(n, t) = 1.0;
                sim_data.Prob_E(n,t) = 0.0;

                sim_data.ProbP_P(n, t) = 1.0;
                sim_data.ProbP_D(n, t) = 0.0;
                sim_data.ProbD_P(n, t) = 1.0;
                sim_data.ProbD_D(n, t) = 0.0;

                sim_data.ProbP_S(n,t) = 1.0;
                sim_data.ProbP_E(n,t) = 0.0;
                sim_data.ProbD_S(n,t) = 1.0;
                sim_data.ProbD_E(n,t) = 0.0;

                sim_data.ProbP_P_nonMiss(n,t) = 1.0;
                sim_data.ProbP_D_nonMiss(n,t) = 0.0;
                sim_data.ProbP_DP_Miss_UnObs(n,t) = 0.0;
                sim_data.ProbP_E_Miss_UnObs(n,t) = 0.0;

                sim_data.ProbD_P_nonMiss(n,t) = 1.0;
                sim_data.ProbD_D_nonMiss(n,t) = 0.0;
                sim_data.ProbD_DP_Miss_UnObs(n,t) = 0.0;
                sim_data.ProbD_E_Miss_UnObs(n,t) = 0.0;

                sim_data.Prob_P_nonMiss(n, t) = 1.0;
                sim_data.Prob_D_nonMiss(n, t) = 0.0;
                sim_data.Prob_DP_Miss_UnObs(n,t) = 0.0;
                sim_data.Prob_E_Miss_UnObs(n,t) = 0.0;
            }
            else {
                sim_data.phi(n, t+1) = exp(para_est.gamma0 + para_est.gamma1 * log(sim_data.phi(n, t))
                    + sim_var.lnphi_sim(n, t+1) * para_est.sigma_phi_eps);
                PhiW phi_w = definePhiWeights(para_vec.vec_phi, sim_data.phi(n, t));

                //// production in the previous period
//                cout << "548 production in the previous period" << endl;
                LKState K_state_P = defineLKState(para_vec.vec_K, sim_data.Capital_P(n, t - 1));
                LKState Lur_state_P = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur_P(n, t - 1));
                LKState Luc_state_P = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc_P(n, t - 1));

                tuple<double,double,double,double> t_SE_P = calStayExit_Prob_logit(para_est,para_vec,Evalmat,phi_w,
                    K_state_P, Lur_state_P, Luc_state_P, 0);
                sim_data.ProbP_S(n,t) = get<0>(t_SE_P);
                sim_data.ProbP_E(n,t) = get<1>(t_SE_P);
                double lnProbP_S = get<2>(t_SE_P);
                double lnProbP_E = get<3>(t_SE_P);

                tuple<double,double,double,double,double,double> t_LK_P = calOptKLurLuc(para_est,
                    para_vec,sim_var,phi_w,K_state_P,Lur_state_P,Luc_state_P,n,t,
                    Evalmat.EVal_PD_K_mat,Evalmat.EVal_PD_Lur_mat,Evalmat.EVal_PD_Luc_mat,0,para_vec.p_K,0);
                sim_data.Capital_P(n, t) = get<0>(t_LK_P);
                sim_data.Employ_ur_P(n, t) = get<1>(t_LK_P);
                sim_data.Employ_uc_P(n, t) = get<2>(t_LK_P);
                sim_data.K_sol_P(n,t) = get<3>(t_LK_P);
                sim_data.Lur_sol_P(n,t) = get<4>(t_LK_P);
                sim_data.Luc_sol_P(n,t) = get<5>(t_LK_P);
//                cout << "569 production in the previous period" << endl;

                //// dormancy in the previous period
//                cout << "572 dormancy in the previous period" << endl;
                LKState K_state_D = defineLKState(para_vec.vec_K, sim_data.Capital_D(n, t - 1));
                LKState Lur_state_D = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur_D(n, t - 1));
                LKState Luc_state_D = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc_D(n, t - 1));

                tuple<double,double,double,double> t_SE_D = calStayExit_Prob_logit(para_est,para_vec,Evalmat,phi_w,
                    K_state_D, Lur_state_D, Luc_state_D, 1);
                sim_data.ProbD_S(n,t) = get<0>(t_SE_D);
                sim_data.ProbD_E(n,t) = get<1>(t_SE_D);
                double lnProbD_S = get<2>(t_SE_D);
                double lnProbD_E = get<3>(t_SE_D);

                tuple<double,double,double,double,double,double> t_LK_D = calOptKLurLuc(para_est,
                    para_vec,sim_var,phi_w,K_state_D,Lur_state_D,Luc_state_D,n,t,
                    Evalmat.EVal_PD_K_mat,Evalmat.EVal_PD_Lur_mat,Evalmat.EVal_PD_Luc_mat,1,para_vec.p_K,0);
                sim_data.Capital_D(n, t) = get<0>(t_LK_D);
                sim_data.Employ_ur_D(n, t) = get<1>(t_LK_D);
                sim_data.Employ_uc_D(n, t) = get<2>(t_LK_D);
                sim_data.K_sol_D(n,t) = get<3>(t_LK_D);
                sim_data.Lur_sol_D(n,t) = get<4>(t_LK_D);
                sim_data.Luc_sol_D(n,t) = get<5>(t_LK_D);
//                cout << "593 dormancy in the previous period" << endl;

                //// production in the previous period
                tuple<double,double,double,double,double,double> t_PD_P = calProductionDormancy_Prob_lognormal(para_est,
                    para_vec,Evalmat,0,sim_data.phi(n, t),
                    sim_data.Capital_P(n, t),sim_data.Employ_ur_P(n, t),
                    sim_data.Employ_uc_P(n, t));
                sim_data.ProbP_P(n, t) = get<0>(t_PD_P); sim_data.ProbP_D(n, t) = get<1>(t_PD_P);

                //// dormancy in the previous period
                tuple<double,double,double,double,double,double> t_PD_D = calProductionDormancy_Prob_lognormal(para_est,
                    para_vec,Evalmat,1,sim_data.phi(n, t),
                    sim_data.Capital_D(n, t),sim_data.Employ_ur_D(n, t),
                    sim_data.Employ_uc_D(n, t));
                sim_data.ProbD_P(n, t) = get<0>(t_PD_D); sim_data.ProbD_D(n, t) = get<1>(t_PD_D);

                //// the probability evolution given the conditional prod/dorm probability
                sim_data.Prob_P(n,t) =
                    sim_data.Prob_P(n,t-1)*sim_data.ProbP_S(n, t)*sim_data.ProbP_P(n, t)
                    + sim_data.Prob_D(n,t-1)*sim_data.ProbD_S(n, t)*sim_data.ProbD_P(n, t);
                sim_data.Prob_D(n,t) =
                    sim_data.Prob_P(n,t-1)*sim_data.ProbP_S(n, t)*sim_data.ProbP_D(n, t)
                    + sim_data.Prob_D(n,t-1)*sim_data.ProbD_S(n, t)*sim_data.ProbD_D(n, t);
                sim_data.Prob_S(n,t) = sim_data.Prob_P(n,t) + sim_data.Prob_D(n,t);
                sim_data.Prob_E(n,t) = sim_data.Prob_E(n,t-1)
                    + sim_data.Prob_P(n,t-1)*sim_data.ProbP_E(n, t)
                    + sim_data.Prob_D(n,t-1)*sim_data.ProbD_E(n, t);

                //// Revenue and Price
                sim_data.Revenue_P(n, t) = RevOpt1_sim_Realized(para_est,sim_data.phi(n, t),
                    sim_data.Capital_P(n, t),sim_data.Employ_ur_P(n, t),
                    sim_data.Employ_uc_P(n, t));
                sim_data.Revenue_D(n, t) = RevOpt1_sim_Realized(para_est,sim_data.phi(n, t),
                    sim_data.Capital_D(n, t),sim_data.Employ_ur_D(n, t),
                    sim_data.Employ_uc_D(n, t));

                sim_data.Price_P(n, t) = calPrice(para_est,alpha_M,sim_data.phi(n, t),
                    sim_data.Capital_P(n, t),sim_data.Employ_ur_P(n, t),
                    sim_data.Employ_uc_P(n, t), 1.0);
                sim_data.Price_D(n, t) = calPrice(para_est,alpha_M,sim_data.phi(n, t),
                    sim_data.Capital_D(n, t),sim_data.Employ_ur_D(n, t),
                    sim_data.Employ_uc_D(n, t), 1.0);

                // probability of Missing (at period t):
                // Case 1: sim_data_missing.State_P_nonMiss(n,t-1) = 1
                sim_data.Prob_Missing_P1(n,t) = ProbMissing(sim_data.Capital_P(n, t-1),
                     sim_data.Employ_ur_P(n, t-1),
                     sim_data.Employ_uc_P(n, t-1),
                     sim_data.Revenue_P(n, t-1),
                     0,0,0,
                     0,0,1,0,
                     sim_data.state_code(n), sim_data.industry(n), t,
                     0, GoodState, LaborIntensive);
                sim_data.Prob_Missing_D1(n,t) = ProbMissing(sim_data.Capital_D(n, t-1),
                     sim_data.Employ_ur_D(n, t-1),
                     sim_data.Employ_uc_D(n, t-1),
                     sim_data.Revenue_D(n, t-1),
                     0,0,0,
                     1,0,0,0,
                     sim_data.state_code(n), sim_data.industry(n),t,
                     0, GoodState, LaborIntensive);

                // Case 3: sim_data_missing.State_DP_nonMiss(n,t-1) = 1
                sim_data.Prob_Missing_DP1(n,t) = ProbMissing(sim_data.Capital_D(n, t-1),
                                                      sim_data.Employ_ur_D(n, t-1),
                                                      sim_data.Employ_uc_D(n, t-1),
                                                      sim_data.Revenue_D(n, t-1),
                                                      0,0,0,
                                                      0,1,0,0,
                                                      sim_data.state_code(n), sim_data.industry(n),t,
                                                      0, GoodState, LaborIntensive);

                //// Transition probability
                sim_data.ProbP_P_nonMiss(n,t) = sim_data.ProbP_S(n,t) * sim_data.ProbP_P(n,t) * (1.0-sim_data.Prob_Missing_P1(n,t));
                sim_data.ProbP_D_nonMiss(n,t) = sim_data.ProbP_S(n,t) * sim_data.ProbP_D(n,t) * (1.0-sim_data.Prob_Missing_P1(n,t));
                sim_data.ProbP_DP_Miss_UnObs(n,t) = sim_data.ProbP_S(n,t) * sim_data.ProbP_P(n,t) * sim_data.Prob_Missing_P1(n,t)
                    + sim_data.ProbP_S(n,t) * sim_data.ProbP_D(n,t) * sim_data.Prob_Missing_P1(n,t);
                sim_data.ProbP_E_Miss_UnObs(n,t) = sim_data.ProbP_E(n,t);

                sim_data.ProbD_P_nonMiss(n,t) = sim_data.ProbD_S(n,t) * sim_data.ProbD_P(n,t) * (1-sim_data.Prob_Missing_D1(n,t));
                sim_data.ProbD_D_nonMiss(n,t) = sim_data.ProbD_S(n,t) * sim_data.ProbD_D(n,t) * (1-sim_data.Prob_Missing_D1(n,t));
                sim_data.ProbD_DP_Miss_UnObs(n,t) = sim_data.ProbD_S(n,t) * sim_data.ProbD_P(n,t) * sim_data.Prob_Missing_D1(n,t)
                    + sim_data.ProbD_S(n,t) * sim_data.ProbD_D(n,t) * sim_data.Prob_Missing_D1(n,t);
                sim_data.ProbD_E_Miss_UnObs(n,t) = sim_data.ProbD_E(n,t);

                sim_data.Prob_P_nonMiss(n, t) =
                    sim_data.Prob_P_nonMiss(n,t-1) * sim_data.ProbP_P_nonMiss(n,t)
                    + sim_data.Prob_D_nonMiss(n,t-1) * sim_data.ProbD_P_nonMiss(n,t)
                    + (sim_data.Prob_P(n,t-1) - sim_data.Prob_P_nonMiss(n,t-1)) * (1.0-sim_data.Prob_Missing_DP1(n,t)) * (1.0 - sim_data.ProbP_E(n,t)) * sim_data.ProbP_P(n,t)
                    + (sim_data.Prob_D(n,t-1) - sim_data.Prob_D_nonMiss(n,t-1)) * (1.0-sim_data.Prob_Missing_DP1(n,t)) * (1.0 - sim_data.ProbD_E(n,t)) * sim_data.ProbD_P(n,t);
                sim_data.Prob_D_nonMiss(n, t) =
                    sim_data.Prob_P_nonMiss(n,t-1) * sim_data.ProbP_D_nonMiss(n,t)
                    + sim_data.Prob_D_nonMiss(n,t-1) * sim_data.ProbD_D_nonMiss(n,t)
                    + (sim_data.Prob_P(n,t-1) - sim_data.Prob_P_nonMiss(n,t-1)) * (1.0-sim_data.Prob_Missing_DP1(n,t)) * (1.0 - sim_data.ProbP_E(n,t)) * sim_data.ProbP_D(n,t)
                    + (sim_data.Prob_D(n,t-1) - sim_data.Prob_D_nonMiss(n,t-1)) * (1.0-sim_data.Prob_Missing_DP1(n,t)) * (1.0 - sim_data.ProbD_E(n,t)) * sim_data.ProbD_D(n,t);
                sim_data.Prob_DP_Miss_UnObs(n, t) =
                    sim_data.Prob_P_nonMiss(n,t-1) * sim_data.ProbP_DP_Miss_UnObs(n,t)
                    + sim_data.Prob_D_nonMiss(n,t-1) * sim_data.ProbD_DP_Miss_UnObs(n,t)
                    + (sim_data.Prob_P(n,t-1) - sim_data.Prob_P_nonMiss(n,t-1)) * (1.0 - sim_data.ProbP_E(n,t)) * sim_data.Prob_Missing_DP1(n,t)
                    + (sim_data.Prob_D(n,t-1) - sim_data.Prob_D_nonMiss(n,t-1)) * (1.0 - sim_data.ProbD_E(n,t)) * sim_data.Prob_Missing_DP1(n,t);
                sim_data.Prob_E_Miss_UnObs(n, t) = sim_data.Prob_E(n,t-1)
                    + sim_data.Prob_P_nonMiss(n,t-1) * sim_data.ProbP_E_Miss_UnObs(n,t)
                    + sim_data.Prob_D_nonMiss(n,t-1) * sim_data.ProbD_E_Miss_UnObs(n,t)
                    + (sim_data.Prob_P(n,t-1) - sim_data.Prob_P_nonMiss(n,t-1)) * sim_data.ProbP_E(n,t)
                    + (sim_data.Prob_D(n,t-1) - sim_data.Prob_D_nonMiss(n,t-1)) * sim_data.ProbD_E(n,t);

                if ( abs(sim_data.Prob_P_nonMiss(n,t)+sim_data.Prob_D_nonMiss(n,t)
                         +sim_data.Prob_DP_Miss_UnObs(n,t)+sim_data.Prob_E_Miss_UnObs(n,t) - 1.0) > 0.001) {
                    cout << "sim_data.Prob_Missing_P1(n,t) = " << sim_data.Prob_Missing_P1(n,t)
                         << "; sim_data.Prob_Missing_D1(n,t) = " << sim_data.Prob_Missing_D1(n,t)
                         << "; sim_data.Prob_Missing_DP1(n,t) = " << sim_data.Prob_Missing_DP1(n,t) << endl;

                    cout << "Prob_P_nonMiss(n,t-1) = " << sim_data.Prob_P_nonMiss(n,t-1)
                         << "; Prob_D_nonMiss(n,t-1) = " << sim_data.Prob_D_nonMiss(n,t-1)
                         << "; Prob_DP_Miss_UnObs(n,t-1) = " << sim_data.Prob_DP_Miss_UnObs(n,t-1)
                         << "; Prob_E_Miss_UnObs(n,t-1) = " << sim_data.Prob_E_Miss_UnObs(n,t-1)
                         << "; sum = " << sim_data.Prob_P_nonMiss(n,t-1) + sim_data.Prob_D_nonMiss(n,t-1)
                         + sim_data.Prob_DP_Miss_UnObs(n,t-1)+ sim_data.Prob_E_Miss_UnObs(n,t-1) << endl;

                    cout << "*******  Transition matrix *********" << endl;
                    cout << "sim_data.ProbP_P_nonMiss(n,t) = " << sim_data.ProbP_P_nonMiss(n,t)
                        << "; sim_data.ProbP_D_nonMiss(n,t) = " << sim_data.ProbP_D_nonMiss(n,t)
                        << "; sim_data.ProbP_DP_Miss_UnObs(n,t) = " << sim_data.ProbP_DP_Miss_UnObs(n,t)
                        << "; sim_data.ProbP_E_Miss_UnObs(n,t) = " << sim_data.ProbP_E_Miss_UnObs(n,t) << endl;

                    cout << "sim_data.ProbD_P_nonMiss(n,t) = " << sim_data.ProbP_P_nonMiss(n,t)
                         << "; sim_data.ProbD_D_nonMiss(n,t) = " << sim_data.ProbP_D_nonMiss(n,t)
                         << "; sim_data.ProbD_DP_Miss_UnObs(n,t) = " << sim_data.ProbP_DP_Miss_UnObs(n,t)
                         << "; sim_data.ProbD_E_Miss_UnObs(n,t) = " << sim_data.ProbP_E_Miss_UnObs(n,t) << endl;

                    cout << "sim_data.ProbP_S(n,t) = " << sim_data.ProbP_S(n,t)
                         << "; sim_data.ProbP_E(n,t) = " << sim_data.ProbP_E(n,t)
                         << "; sim_data.ProbP_S(n,t) + sim_data.ProbP_E(n,t) = " << sim_data.ProbP_S(n,t) + sim_data.ProbP_E(n,t) << endl;

                    cout << "sim_data.ProbD_S(n,t) = " << sim_data.ProbD_S(n,t)
                         << "; sim_data.ProbD_E(n,t) = " << sim_data.ProbD_E(n,t)
                         << "; sim_data.ProbD_S(n,t) + sim_data.ProbD_E(n,t) = " << sim_data.ProbD_S(n,t) + sim_data.ProbD_E(n,t) << endl;

                    cout << "sim_data.ProbP_P(n,t) = " << sim_data.ProbP_P(n,t)
                         << "; sim_data.ProbP_D(n,t) = " << sim_data.ProbP_D(n,t)
                         << "; sim_data.ProbP_P(n,t) + sim_data.ProbP_D(n,t) = " << sim_data.ProbP_P(n,t) + sim_data.ProbP_D(n,t) << endl;

                    cout << "sim_data.ProbD_P(n,t) = " << sim_data.ProbD_P(n,t)
                         << "; sim_data.ProbD_D(n,t) = " << sim_data.ProbD_D(n,t)
                         << "; sim_data.ProbD_P(n,t) + sim_data.ProbD_D(n,t) = " << sim_data.ProbD_P(n,t) + sim_data.ProbD_D(n,t) << endl;

                    cout << "; sim_data.Prob_P(n,t) = " << sim_data.Prob_P(n,t)
                        << "; sim_data.Prob_D(n,t) = " << sim_data.Prob_D(n,t)
                        << "; sim_data.Prob_S(n,t) = " << sim_data.Prob_S(n,t)
                        << "; sim_data.Prob_E(n,t) = " << sim_data.Prob_E(n,t) << endl;

                    cout << "Prob_P_nonMiss(n,t) = " << sim_data.Prob_P_nonMiss(n,t)
                         << "; Prob_D_nonMiss(n,t) = " << sim_data.Prob_D_nonMiss(n,t)
                         << "; Prob_DP_Miss_UnObs(n,t) = " << sim_data.Prob_DP_Miss_UnObs(n,t)
                         << "; Prob_E_Miss_UnObs(n,t) = " << sim_data.Prob_E_Miss_UnObs(n,t)
                         << "; sum = " << sim_data.Prob_P_nonMiss(n,t) + sim_data.Prob_D_nonMiss(n,t)
                         + sim_data.Prob_DP_Miss_UnObs(n,t)+ sim_data.Prob_E_Miss_UnObs(n,t) << endl;

                    cout << "sim_data.Prob_E(n,t-1) + sim_data.Prob_P(n,t-1) + sim_data.Prob_D(n,t-1) = "
                         << sim_data.Prob_E(n,t-1) + sim_data.Prob_P(n,t-1) + sim_data.Prob_D(n,t-1) << endl;
                    cout << "sim_data.Prob_E(n,t) + sim_data.Prob_D(n,t) + sim_data.Prob_P(n,t) = "
                         << sim_data.ProbP_E(n,t) + sim_data.ProbP_D(n,t) + sim_data.ProbP_P(n,t) << endl;
                    cout << "sim_data.ProbP_D(n,t) + sim_data.ProbP_P(n,t) = "
                         << sim_data.ProbP_D(n,t) + sim_data.ProbP_P(n,t) << endl;

                    throw runtime_error("871");
                }
            }
        }
    }
//    };
//    MultiThreads::simple_parallel_for(worker,SimN, threadsManagement);

    sim_data.cost_P = para_est.w_uc * sim_data.Employ_uc_P + para_est.w_ur * sim_data.Employ_ur_P
        + para_est.r_K * sim_data.Capital_P;

    sim_data.cost_D = para_est.w_uc * sim_data.Employ_uc_D + para_est.w_ur * sim_data.Employ_ur_D
        + para_est.r_K * sim_data.Capital_D;


    // sim_data.Prob_E; sim_data.Prob_S; sim_data.Prob_P; sim_data.Prob_D
    // are the probability tracking states evolution without considering missing probability
    return sim_data;
}

/** within the function: generate the missing data **/
SimData alias::Simulation_Missing_Step2Estimation_FullVersion(const SimVar & sim_var, const SimData & sim_data,
    const int & SimN, const int & TBar, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    SimData sim_data_missing = sim_data;
//    cout << "808 Simulate sim_data_missing.Prob_E_Miss, sim_data_missing.Prob_DP_Miss, ..." << endl;
    //// calculate the probability of disappearing forever
    //// sim_data_missing.Prob_E_Miss; sim_data_missing.Prob_DP_Miss; sim_data_missing.Prob_P_nonMiss;
    //// sim_data_missing.Prob_D_nonMiss are the probabilities of states taking the disappearing probability into account;
    auto worker2 = [&](size_t n, unsigned thread_id) {
//    size_t thread_id = 0;
//    for (size_t n = 0; n < SimN; ++n) {
        for (int t = TBar-1; t >= sim_data_missing.first_year_in_Data(n); --t) {
//                    cout << "n = " << n << "; t = " << t << endl;
            if (t == TBar-1) {
                sim_data_missing.Prob_E_Miss(n,t) = sim_data_missing.Prob_DP_Miss_UnObs(n,t)
                    + sim_data_missing.Prob_E_Miss_UnObs(n,t);
                sim_data_missing.Prob_DP_Miss(n, t) = 0.0;

                sim_data_missing.ProbP_E_Miss = sim_data_missing.ProbP_DP_Miss_UnObs(n,t)
                    + sim_data_missing.ProbP_E_Miss_UnObs(n,t);
                sim_data_missing.ProbP_DP_Miss(n, t) = 0.0;

                sim_data_missing.ProbD_E_Miss = sim_data_missing.ProbD_DP_Miss_UnObs(n,t)
                    + sim_data_missing.ProbD_E_Miss_UnObs(n,t);
                sim_data_missing.ProbD_DP_Miss(n, t) = 0.0;

                if ( abs(sim_data_missing.Prob_E_Miss(n,t)+sim_data_missing.Prob_DP_Miss(n, t)
                         +sim_data_missing.Prob_P_nonMiss(n,t)+sim_data_missing.Prob_D_nonMiss(n,t) - 1.0) > 0.0001) {
                    cout << "Prob_E_Miss(n,t) = " << sim_data_missing.Prob_E_Miss(n,t)
                         << "; Prob_DP_Miss(n,t) = " << sim_data_missing.Prob_DP_Miss(n, t)
                         << "; Prob_P_nonMiss(n,t) = " << sim_data_missing.Prob_P_nonMiss(n,t)
                         << "; Prob_D_nonMiss(n,t) = " << sim_data_missing.Prob_D_nonMiss(n,t)
                         << "; sum = " << sim_data_missing.Prob_E_Miss(n,t)+sim_data_missing.Prob_DP_Miss(n, t)
                         +sim_data_missing.Prob_P_nonMiss(n,t)+sim_data_missing.Prob_D_nonMiss(n,t) << endl;
                    throw runtime_error("817");
                }
            }
            else {
                sim_data_missing.Prob_E_Miss(n,t) =
                    sim_data_missing.Prob_DP_Miss_UnObs(n,t) * sim_data_missing.Prob_E_Miss(n,t+1)
                    + sim_data_missing.Prob_E_Miss_UnObs(n,t);
                sim_data_missing.Prob_DP_Miss(n, t) =
                    sim_data_missing.Prob_DP_Miss_UnObs(n,t) * (1 - sim_data_missing.Prob_E_Miss(n,t+1));

                sim_data_missing.ProbP_E_Miss =
                    sim_data_missing.ProbP_DP_Miss_UnObs(n,t) * sim_data_missing.Prob_E_Miss(n,t+1)
                    + sim_data_missing.ProbP_E_Miss_UnObs(n,t);
                sim_data_missing.ProbP_DP_Miss(n, t) =
                    sim_data_missing.ProbP_DP_Miss_UnObs(n,t) * (1 - sim_data_missing.Prob_E_Miss(n,t+1));

                sim_data_missing.ProbD_E_Miss =
                    sim_data_missing.ProbD_DP_Miss_UnObs(n,t) * sim_data_missing.Prob_E_Miss(n,t+1)
                    + sim_data_missing.ProbD_E_Miss_UnObs(n,t);
                sim_data_missing.ProbD_DP_Miss(n, t) =
                    sim_data_missing.ProbD_DP_Miss_UnObs(n,t) * (1 - sim_data_missing.Prob_E_Miss(n,t+1));

                if ( abs(sim_data_missing.Prob_E_Miss(n,t)+sim_data_missing.Prob_DP_Miss(n, t)
                         +sim_data_missing.Prob_P_nonMiss(n,t)+sim_data_missing.Prob_D_nonMiss(n,t) - 1.0) > 0.0001) {
                    cout << "Prob_E_Miss(n,t) = " << sim_data_missing.Prob_E_Miss(n,t)
                         << "; Prob_DP_Miss(n,t) = " << sim_data_missing.Prob_DP_Miss(n, t)
                         << "; Prob_P_nonMiss(n,t) = " << sim_data_missing.Prob_P_nonMiss(n,t)
                         << "; Prob_D_nonMiss(n,t) = " << sim_data_missing.Prob_D_nonMiss(n,t)
                         << "; sum = " << sim_data_missing.Prob_E_Miss(n,t)+sim_data_missing.Prob_DP_Miss(n, t)
                                          +sim_data_missing.Prob_P_nonMiss(n,t)+sim_data_missing.Prob_D_nonMiss(n,t) << endl;
                    throw runtime_error("837");
                }
            }
        }
//    }
    };
    MultiThreads::simple_parallel_for(worker2,SimN, threadsManagement);
//
////
//    ArrayXXd ProbAll = sim_data_missing.Prob_P_nonMiss + sim_data_missing.Prob_D_nonMiss + sim_data_missing.Prob_DP_Miss
//        + sim_data_missing.Prob_E_Miss;
//    cout << "min and max of sim_data_missing.Prob_P_nonMiss = " << sim_data_missing.Prob_P_nonMiss.minCoeff()
//        << "; " << sim_data_missing.Prob_P_nonMiss.maxCoeff() << endl;
//    cout << "min and max of sim_data_missing.Prob_D_nonMiss = " << sim_data_missing.Prob_D_nonMiss.minCoeff()
//         << "; " << sim_data_missing.Prob_D_nonMiss.maxCoeff() << endl;
//    cout << "min and max of sim_data_missing.Prob_DP_Miss = " << sim_data_missing.Prob_DP_Miss.minCoeff()
//         << "; " << sim_data_missing.Prob_DP_Miss.maxCoeff() << endl;
//    cout << "min and max of sim_data_missing.Prob_E_Miss = " << sim_data_missing.Prob_E_Miss.minCoeff()
//         << "; " << sim_data_missing.Prob_E_Miss.maxCoeff() << endl;
////    throw runtime_error("867");
////
//     writeToCSVfile("phi_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.phi.cast<double>().matrix());
//     writeToCSVfile("Employ_uc_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Employ_uc_P.cast<double>().matrix());
//     writeToCSVfile("Employ_uc_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Employ_uc_D.cast<double>().matrix());
//
//     writeToCSVfile("Employ_ur_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Employ_ur_P.cast<double>().matrix());
//     writeToCSVfile("Employ_ur_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Employ_ur_D.cast<double>().matrix());
//
//     writeToCSVfile("Revenue_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Revenue_P.cast<double>().matrix());
//     writeToCSVfile("Revenue_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Revenue_D.cast<double>().matrix());
//
//     writeToCSVfile("Capital_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Capital_P.cast<double>().matrix());
//     writeToCSVfile("Capital_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Capital_D.cast<double>().matrix());
//
//     writeToCSVfile("ProbP_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_P.cast<double>().matrix());
//     writeToCSVfile("ProbP_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_D.cast<double>().matrix());
//     writeToCSVfile("ProbP_E_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_E.cast<double>().matrix());
//
//     writeToCSVfile("ProbD_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_P.cast<double>().matrix());
//     writeToCSVfile("ProbD_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_D.cast<double>().matrix());
//     writeToCSVfile("ProbD_E_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_E.cast<double>().matrix());
//
//     writeToCSVfile("cost_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.cost_P.cast<double>().matrix());
//     writeToCSVfile("cost_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.cost_D.cast<double>().matrix());
//
//     writeToCSVfile("Price_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Price_P.cast<double>().matrix());
//     writeToCSVfile("Price_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Price_D.cast<double>().matrix());
//
//     writeToCSVfile("Prob_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_P.cast<double>().matrix());
//     writeToCSVfile("Prob_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_D.cast<double>().matrix());
//     writeToCSVfile("Prob_S_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_S.cast<double>().matrix());
//     writeToCSVfile("Prob_E_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_E.cast<double>().matrix());
//
//     writeToCSVfile("ProbP_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_P.cast<double>().matrix());
//     writeToCSVfile("ProbP_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_D.cast<double>().matrix());
//     writeToCSVfile("ProbD_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_P.cast<double>().matrix());
//     writeToCSVfile("ProbD_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_P.cast<double>().matrix());
//
//     writeToCSVfile("Prob_P_nonMiss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_P_nonMiss.cast<double>().matrix());
//     writeToCSVfile("Prob_D_nonMiss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_D_nonMiss.cast<double>().matrix());
//     writeToCSVfile("Prob_DP_Miss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_DP_Miss.cast<double>().matrix());
//     writeToCSVfile("Prob_E_Miss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_E_Miss.cast<double>().matrix());
//
//     writeToCSVfile("ProbP_P_nonMiss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_P_nonMiss.cast<double>().matrix());
//     writeToCSVfile("ProbP_D_nonMiss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_D_nonMiss.cast<double>().matrix());
//     writeToCSVfile("ProbP_DP_Miss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_DP_Miss.cast<double>().matrix());
//     writeToCSVfile("ProbP_E_Miss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_E_Miss.cast<double>().matrix());
//
//     writeToCSVfile("ProbD_P_nonMiss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_P_nonMiss.cast<double>().matrix());
//     writeToCSVfile("ProbD_D_nonMiss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_D_nonMiss.cast<double>().matrix());
//     writeToCSVfile("ProbD_DP_Miss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_DP_Miss.cast<double>().matrix());
//     writeToCSVfile("ProbD_E_Miss_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_E_Miss.cast<double>().matrix());
//
//     writeToCSVfile("Prob_Missing_P1_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_Missing_P1.cast<double>().matrix());
//     writeToCSVfile("Prob_Missing_D1_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_Missing_D1.cast<double>().matrix());
//     writeToCSVfile("Prob_Missing_DP1_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Prob_Missing_DP1.cast<double>().matrix());
//
//     writeToCSVfile("ProbD_DP_Miss_UnObs_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_DP_Miss_UnObs.cast<double>().matrix());
//     writeToCSVfile("ProbD_E_Miss_UnObs_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbD_E_Miss_UnObs.cast<double>().matrix());
//
//     writeToCSVfile("ProbP_DP_Miss_UnObs_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_DP_Miss_UnObs.cast<double>().matrix());
//     writeToCSVfile("ProbP_E_Miss_UnObs_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.ProbP_E_Miss_UnObs.cast<double>().matrix());
//
//     writeToCSVfile("K_sol_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.K_sol_D.cast<double>().matrix());
//     writeToCSVfile("Lur_sol_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Lur_sol_D.cast<double>().matrix());
//     writeToCSVfile("Luc_sol_D_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Luc_sol_D.cast<double>().matrix());
//
//     writeToCSVfile("K_sol_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.K_sol_P.cast<double>().matrix());
//     writeToCSVfile("Lur_sol_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Lur_sol_P.cast<double>().matrix());
//     writeToCSVfile("Luc_sol_P_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_data_missing.Luc_sol_P.cast<double>().matrix());
//
//     writeToCSVfile("StateCategory_Sim_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_var.StateCategory_Sim.cast<double>().matrix());
//     writeToCSVfile("Misssim_Empirical_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_var.Misssim_Empirical.cast<double>().matrix());
//
//     writeToCSVfile("PDsim_Empirical_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_var.PDsim_Empirical.cast<double>().matrix());
//     writeToCSVfile("SEsim_Empirical_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                    sim_var.SEsim_Empirical.cast<double>().matrix());
// //
//     // throw runtime_error("951");

    return sim_data_missing;
}
//
/** initialize the simulated variables **/
SimData alias::InitializeSimData(const int & SimN, const int & TBar) {
    SimData sim_data;

    sim_data.phi = ArrayXXd::Zero(SimN,TBar+1);
    sim_data.Revenue = ArrayXXd::Zero(SimN,TBar);
    sim_data.Revenue_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Revenue_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.Revenue_Miss = ArrayXXi::Zero(SimN,TBar);

    sim_data.Capital = ArrayXXd::Zero(SimN,TBar);
    sim_data.Capital_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Capital_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.CapitalMiss = ArrayXXi::Zero(SimN,TBar);

    sim_data.Employ_ur = ArrayXXd::Zero(SimN,TBar);
    sim_data.Employ_ur_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Employ_ur_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.Employ_urMiss = ArrayXXi::Zero(SimN,TBar);

    sim_data.Employ_uc = ArrayXXd::Zero(SimN,TBar);
    sim_data.Employ_uc_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Employ_uc_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.Employ_ucMiss = ArrayXXi::Zero(SimN,TBar);

    sim_data.K_sol_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Lur_sol_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Luc_sol_P = ArrayXXd::Zero(SimN,TBar);

    sim_data.K_sol_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.Lur_sol_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.Luc_sol_D = ArrayXXd::Zero(SimN,TBar);

    sim_data.K_sol = ArrayXXd::Zero(SimN,TBar);
    sim_data.Lur_sol = ArrayXXd::Zero(SimN,TBar);
    sim_data.Luc_sol = ArrayXXd::Zero(SimN,TBar);

    sim_data.Price = ArrayXXd::Zero(SimN,TBar);
    sim_data.Price_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Price_D = ArrayXXd::Zero(SimN,TBar);

    sim_data.cost = ArrayXXd::Zero(SimN,TBar);
    sim_data.cost_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.cost_D = ArrayXXd::Zero(SimN,TBar);

    sim_data.Prob_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_S = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_E = ArrayXXd::Zero(SimN,TBar);

    sim_data.ProbP_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_S = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_E = ArrayXXd::Zero(SimN,TBar);

    sim_data.ProbD_P = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_D = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_S = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_E = ArrayXXd::Zero(SimN,TBar);

    sim_data.good = ArrayXi::Zero(SimN);
    sim_data.labour_intensive_ind = ArrayXi::Zero(SimN);
    sim_data.state_code = ArrayXi::Zero(SimN);
    sim_data.industry = ArrayXi::Zero(SimN);
    sim_data.first_year_in_Data = ArrayXi::Zero(SimN);
    sim_data.last_year_in_Data = ArrayXi::Zero(SimN);

    sim_data.State_P = ArrayXXi::Zero(SimN,TBar);
    sim_data.State_D = ArrayXXi::Zero(SimN,TBar);
    sim_data.State_S = ArrayXXi::Zero(SimN,TBar);
    sim_data.State_E = ArrayXXi::Zero(SimN,TBar);

    sim_data.State_P_nonMiss = ArrayXXi::Zero(SimN,TBar);
    sim_data.State_D_nonMiss = ArrayXXi::Zero(SimN,TBar);
    sim_data.State_DP_Miss = ArrayXXi::Zero(SimN,TBar);
    sim_data.State_E_Miss = ArrayXXi::Zero(SimN,TBar);

    sim_data.ProbP_P_nonMiss = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_D_nonMiss = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_DP_Miss_UnObs = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_E_Miss_UnObs = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_DP_Miss = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbP_E_Miss = ArrayXXd::Zero(SimN,TBar);

    sim_data.ProbD_P_nonMiss = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_D_nonMiss = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_DP_Miss_UnObs = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_E_Miss_UnObs = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_DP_Miss = ArrayXXd::Zero(SimN,TBar);
    sim_data.ProbD_E_Miss = ArrayXXd::Zero(SimN,TBar);

    sim_data.Prob_P_nonMiss = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_D_nonMiss = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_DP_Miss = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_E_Miss = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_DP_Miss_UnObs = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_E_Miss_UnObs = ArrayXXd::Zero(SimN,TBar);

    sim_data.Prob_Missing_P1 = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_Missing_D1 = ArrayXXd::Zero(SimN,TBar);
    sim_data.Prob_Missing_DP1 = ArrayXXd::Zero(SimN,TBar);

    sim_data.CutoffVal = ArrayXXd::Zero(SimN,TBar);
    return sim_data;
}

/** calculate revenue given productivity,capital,employment **/
double alias::RevOpt1_sim_Realized(const ParaEst & para_est, const double & phi, const double & K, const double & L_ur,
    const double & L_uc) {

    /*** CES production function: Cobb Douglas ***/
    double L = pow(L_uc, para_est.alpha_Lc) * pow(L_ur, para_est.alpha_Lr);
    double KL = pow(K, para_est.alpha_K) * pow(L, para_est.alpha_L);
    double Rev = phi * pow(KL, (para.sigma-1.0)/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) )
        * pow( para_est.PI,1.0/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) );

    // /*** CES production function: CES ***/
    // double L = pow( para_est.alpha_Lc * pow(L_uc, (para.sigma_Lrc-1)/para.sigma_Lrc)
    //     + para_est.alpha_Lr * pow(L_uc, (para.sigma_Lrc-1)/para.sigma_Lrc) , para.sigma_Lrc/(para.sigma_Lrc-1) );
    // double KL = pow( para_est.alpha_K/(para_est.alpha_K+para_est.alpha_L)* pow(K, (para.sigma_LK-1)/para.sigma_LK)
    //     + para_est.alpha_L/(para_est.alpha_K+para_est.alpha_L) * pow(L, (para.sigma_LK-1)/para.sigma_LK) , para.sigma_LK/(para.sigma_LK-1) );
    // KL = pow(KL,1.0-para_est.alpha_M);
    //
    // double Rev = phi * pow(KL, (para.sigma-1.0)/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) )
    //     * pow( para_est.PI,1.0/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) );

    return Rev;
}

/** calculate price charged by firms **/
double alias::calPrice(const ParaEst & para_est, const double & alpha_M, const double & phi, const double & K,
    const double & L_ur, const double & L_uc, const double & PI) {

    /*** Cobb Douglas Production function: Cobb Douglas ***/
    double L = pow(L_uc, para_est.alpha_Lc) * pow(L_ur, para_est.alpha_Lr);
    double KL = pow(K, para_est.alpha_K) * pow(L, para_est.alpha_L);
    double Price = pow(phi,-1.0/(para.sigma-1.0)) * pow( KL, - 1.0/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) )
        * pow( para_est.PI,(1.0-para_est.alpha_M)/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) );

    // /*** Cobb Douglas Production function: CES ***/
    // double L = pow( para_est.alpha_Lr*pow(L_ur,(para.sigma_Lrc-1)/para.sigma_Lrc)
    //     + para_est.alpha_Lc*pow(L_uc,(para.sigma_Lrc-1)/para.sigma_Lrc), para.sigma_Lrc/(para.sigma_Lrc-1) );
    // double KL = pow( para_est.alpha_L/(para_est.alpha_K+para_est.alpha_L) * pow(L,(para.sigma_LK-1)/para.sigma_LK)
    //     + para_est.alpha_K/(para_est.alpha_K+para_est.alpha_L) * pow(L,(para.sigma_LK-1)/para.sigma_LK), para.sigma_LK/(para.sigma_LK-1) );
    // KL = pow(KL,1.0-para_est.alpha_M);
    //
    // double Price = pow(phi,-1.0/(para.sigma-1.0)) * pow( KL, - 1.0/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) )
    //     * pow( para_est.PI,(1.0-para_est.alpha_M)/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) );

    return Price;
}

/** Backout phi based on theta1 : need to differentiate whether a firm is produciton/dormant in the previous period **/
ArrayXXd alias::Backout_phi_a_PD(const ArrayXXd & Revenue, const ArrayXXd & Capital, const ArrayXXd & Employ_ur,
    const ArrayXXd & Employ_uc, const double & alpha_tilde_K,const double & alpha_tilde_L, const double & alpha_Lr,
    const double & alpha_Lc) {

//    ArrayXXd L = pow( para_est.alpha_Lc*pow(sim_data.Employ_uc,(para_est.mu-1)/para_est.mu)
//        + para_est.alpha_Lr*pow(sim_data.Employ_ur,(para_est.mu-1)/para_est.mu), para_est.mu/(para_est.mu-1));

    ArrayXXd logL = alpha_Lc * log(Employ_uc) + alpha_Lr * log(Employ_ur);
    ArrayXXd lnphi_a = log(Revenue) - alpha_tilde_K*log(Capital)-alpha_tilde_L*logL;

    return lnphi_a;
}
//
/** Simulate capital/employment given the previous year **/
tuple<double,double,double,double,double,double> alias::calOptKLurLuc(const ParaEst & para_est,
    const ParaVec & para_vec, const SimVar & sim_var, const PhiW & phi_w,
    const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state, const int & n, const int & t,
    const ArrayXXd & EVal_PD_K_mat, const ArrayXXd & EVal_PD_Lur_mat, const ArrayXXd & EVal_PD_Luc_mat,
    const int & PD_lag, const double & p_K, const int & InitialPeriod) {

    double Capital; double Employ_ur; double Employ_uc;

    ArrayXXd EVal_K_mat = EVal_PD_K_mat(seqN(PD_lag*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
    tuple<double,int> t_K = sol_opt_K(para_est,para_vec,phi_w, EVal_K_mat, Lur_state,
        Luc_state,K_state.val_state, InitialPeriod);
    double K = get<0>(t_K);
    /*********************************/
    Capital = exp(log(K) + sim_var.lnK_sim(n, t) * para_est.sigma_Kerror);
    Capital = max(0.001, Capital);
    /*********************************/

    LKState K_state_up = defineLKState(para_vec.vec_K,Capital);
    ArrayXXd EVal_Lur_mat = EVal_PD_Lur_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
    tuple<double,int> t_Lur = sol_opt_Lur(para_est,para_vec,phi_w,EVal_Lur_mat,K_state_up,Luc_state,
        Lur_state.val_state, InitialPeriod);
    double Lur = get<0>(t_Lur);
    /*********************************/
    Employ_ur = exp(log(Lur) + sim_var.lnLur_sim(n, t) * para_est.sigma_Lerror_ur);
    Employ_ur = max(0.5, Employ_ur);
    /*********************************/

    LKState Lur_state_up = defineLKState(para_vec.vec_Lur,Employ_ur);
    ArrayXXd EVal_Luc_mat = EVal_PD_Luc_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
    tuple<double,int> t_Luc = sol_opt_Luc(para_est,para_vec,phi_w,EVal_Luc_mat,K_state_up,Lur_state_up,
        Luc_state.val_state, InitialPeriod);
    double Luc = get<0>(t_Luc);
    /*********************************/
    Employ_uc = exp(log(Luc) + sim_var.lnLuc_sim(n, t) * para_est.sigma_Lerror_uc);
    Employ_uc = max(0.5, Employ_uc);
    /*********************************/
    return tuple<double,double,double,double,double,double>(Capital,Employ_ur,Employ_uc,K,Lur,Luc);
}

/********************************************************************************************************************
* Part1 of objection function in the second stage estimation
********************************************************************************************************************/
tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> alias::EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
    const std::vector<double> & theta1, const SimData & sim_data,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {

    double pi = 3.1415926;
    int N = sim_data.good.rows();

    ArrayXi Period_vec = ArrayXi::LinSpaced(para.EstTbar,0,para.EstTbar-1);
    ArrayXXi Period_mat(N,para.EstTbar);
    Period_mat.rowwise() = Period_vec.transpose();

    ArrayXd Ones = ArrayXd::Ones(N);
    ParaEst1 para_est = constructParaEst_Part1(theta1,GoodState,LaborIntensive);

    ArrayXXd logResR_P = ArrayXXd::Zero(N,para.EstTbar);
    logResR_P = log(sim_data.Revenue_P)
        - para_est.alpha_tilde_L*para_est.alpha_Lc*log(sim_data.Employ_uc_P)
        - para_est.alpha_tilde_L*para_est.alpha_Lr*log(sim_data.Employ_ur_P)
        - para_est.alpha_tilde_K*log(sim_data.Capital_P);

//    ArrayXXd logResR_D = ArrayXXd::Zero(N,para.EstTbar);
//    logResR_D = log(sim_data.Revenue_D)
//        - para_est.alpha_tilde_L*para_est.alpha_Lc*log(sim_data.Employ_uc_D)
//        - para_est.alpha_tilde_L*para_est.alpha_Lr*log(sim_data.Employ_ur_D)
//        - para_est.alpha_tilde_K*log(sim_data.Capital_D);

    ArrayXXd lnLtemp_P = para_est.alpha_Lc * log(sim_data.Employ_uc_P) + para_est.alpha_Lr * log(sim_data.Employ_ur_P);
//    ArrayXXd lnLtemp_D = para_est.alpha_Lc * log(sim_data.Employ_uc_D) + para_est.alpha_Lr * log(sim_data.Employ_ur_D);

    ArrayXXd logL = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd dlogL_d0 = ArrayXXd::Zero(N,para.EstTbar);  // gamma0
    ArrayXXd dlogL_d1 = ArrayXXd::Zero(N,para.EstTbar); // gamma1
    ArrayXXd dlogL_d2 = ArrayXXd::Zero(N,para.EstTbar);  // sigma_phi_eps

    ArrayXXi count = ArrayXXi::Zero(N,para.EstTbar);
    ArrayXXd count_prob = ArrayXXd::Zero(N,para.EstTbar);

    //// Other period
    auto worker = [&](size_t n, unsigned thread_id) {
////    size_t thread_id = 0;
    // for (size_t n = 0; n < N; ++n) {
        for (int t = sim_data.first_year_in_Data(n) + 1; t <= para.EstTbar-1; ++t) {
            // cout << "n = " << n << "; t = " << t << "; 770" << endl;
            int diff_period = 1;
            int max_period = t-1;
            double mu1 = -para_est.gamma0 * (1.0 - pow(para_est.gamma1, diff_period)) / (1.0 - para_est.gamma1);
            double mu2_P = logResR_P(n, t) - pow(para_est.gamma1, diff_period) * logResR_P(n, max_period);
            double mu_P = mu1 + mu2_P;

            double var = pow(para_est.sigma_phi_eps, 2.0) * (1.0 - pow(para_est.gamma1, 2.0 * diff_period))
                / (1.0 - pow(para_est.gamma1, 2.0));
            double sigma = sqrt(var);

            double dmu_dgamma0 = -(1.0 - pow(para_est.gamma1, diff_period)) / (1.0 - para_est.gamma1);
            double dmu_dgamma1_P = -diff_period * pow(para_est.gamma1, diff_period - 1) * logResR_P(n, max_period)
                + para_est.gamma0 * diff_period * pow(para_est.gamma1, diff_period - 1) / (1.0 - para_est.gamma1)
                - para_est.gamma0 * (1 - pow(para_est.gamma1, diff_period)) / pow(1.0 - para_est.gamma1, 2);
//            double dmu_dgamma1_D = -diff_period * pow(para_est.gamma1, diff_period - 1) * logResR_D(n, max_period)
//                + para_est.gamma0 * diff_period * pow(para_est.gamma1, diff_period - 1) / (1.0 - para_est.gamma1)
//                - para_est.gamma0 * (1 - pow(para_est.gamma1, diff_period)) / pow(1.0 - para_est.gamma1, 2);

            double dsigma_dgamma1 = -para_est.sigma_phi_eps * diff_period * pow(para_est.gamma1, 2.0 * diff_period - 1)
                * pow(1 - pow(para_est.gamma1, 2 * diff_period), -0.5) * pow(1 - pow(para_est.gamma1, 2), -0.5)
                + para_est.gamma1 * para_est.sigma_phi_eps * pow(1 - pow(para_est.gamma1, 2 * diff_period), 0.5)
                / pow(1 - pow(para_est.gamma1, 2.0), 1.5);
            double dsigma_dsigmaphi = sqrt( (1 - pow(para_est.gamma1, 2.0 * diff_period)) / (1 - pow(para_est.gamma1, 2.0)) );

            double f_P = -0.5 * log(2.0 * pi) - log(sigma) - 0.5 * pow(mu_P / sigma, 2);
//            double f_D = -0.5 * log(2.0 * pi) - log(sigma) - 0.5 * pow(mu_D / sigma, 2);

            double df_d0_P = -mu_P / pow(sigma, 2) * dmu_dgamma0;
//            double df_d0_D = -mu_D / pow(sigma, 2) * dmu_dgamma0;

            double df_d1_P = -1.0 / sigma * dsigma_dgamma1 - mu_P / pow(sigma, 2) * dmu_dgamma1_P
                + pow(mu_P, 2) / pow(sigma, 3) * dsigma_dgamma1; // gamma1
//            double df_d1_D = -1.0 / sigma * dsigma_dgamma1 - mu_D / pow(sigma, 2) * dmu_dgamma1_D
//                + pow(mu_D, 2) / pow(sigma, 3) * dsigma_dgamma1; // gamma1

            double df_d2_P = -1.0 / sigma * dsigma_dsigmaphi + pow(mu_P, 2) / pow(sigma, 3) * dsigma_dsigmaphi; // sigma_phi_eps
//            double df_d2_D = -1.0 / sigma * dsigma_dsigmaphi + pow(mu_D, 2) / pow(sigma, 3) * dsigma_dsigmaphi; // sigma_phi_eps

            logL(n, t) = f_P * sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t);
//            cout << " logL(n, t) = " <<  logL(n, t) << "; f = " << f << "; sim_data.Prob_P_nonMiss(n, t-1) = "
//                << sim_data.Prob_P_nonMiss(n, t-1) << "; sim_data.Prob_P_nonMiss(n, t)" << sim_data.Prob_P_nonMiss(n, t)
//                << "; Prob_S_cum = " << Prob_S_cum << endl;
            if (isfinite(logL(n, t)) == 0) {
                cout << " logResR_P(n, max_period) = " <<  logResR_P(n, max_period)
                    << "; max_period = " << max_period << endl;
                throw runtime_error("832");
            }
            dlogL_d0(n, t) = df_d0_P * sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t);

            dlogL_d1(n, t) = df_d1_P * sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t);

            dlogL_d2(n, t) = df_d2_P * sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t);

            count(n, t) = 1;
            count_prob(n,t) = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t);
        }
    // }
    };
    MultiThreads::simple_parallel_for(worker, N, threadsManagement);

    int Ncount = N;
    double Ncount_prob = count_prob.sum();

//    logL = (logL.isFinite()).select(logL, 0.0);
    double y_sum = -logL.sum();

    std::vector<double> grad_sum(para.dim1);
    dlogL_d0 = (dlogL_d0.isFinite()).select(dlogL_d0, 0.0);
    grad_sum[0] = -dlogL_d0.sum();

    dlogL_d1 = (dlogL_d1.isFinite()).select(dlogL_d1, 0.0);
    grad_sum[1] = -dlogL_d1.sum();

    dlogL_d2 = (dlogL_d2.isFinite()).select(dlogL_d2, 0.0);
    grad_sum[2] = -dlogL_d2.sum();

    return tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd>(
        y_sum,grad_sum,Ncount,Ncount_prob,logL,dlogL_d0,dlogL_d1,dlogL_d2);
}

/********************************************************************************************************************
* Part2 of objection function in the second stage estimation
********************************************************************************************************************/
tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> alias::EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(
    const std::vector<double> & theta2, const SimData & sim_data, MultiThreads::Threads_Management & threadsManagement) {

    int N = sim_data.good.rows();
    ArrayXi Period_vec = ArrayXi::LinSpaced(para.EstTbar,0,para.EstTbar-1);
    ArrayXXi Period_mat(N,para.EstTbar);
    Period_mat.rowwise() = Period_vec.transpose();

    ArrayXd Ones = ArrayXd::Ones(N);
    ArrayXXd ValueAdd = sim_data.Revenue_P;
    ParaEst2 para_est2 = constructParaEst_Part2(theta2);

    ArrayXXd logL = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd dlogL_d0 = ArrayXXd::Zero(N,para.EstTbar);  // F_P_P
    ArrayXXd dlogL_d1 = ArrayXXd::Zero(N,para.EstTbar);  // F_D_P
    ArrayXXd dlogL_d2 = ArrayXXd::Zero(N,para.EstTbar);  // sigma_PD

    ArrayXXd f_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd df_d0_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd df_d1_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd df_d2_P_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXd f_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd df_d0_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd df_d1_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd df_d2_D_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXi count = ArrayXXi::Zero(N,para.EstTbar);
    ArrayXXd Nprob = ArrayXXd::Zero(N,para.EstTbar);
    //// Other period
    auto worker = [&](size_t n, unsigned thread_id) {
        //    for (size_t n = 0; n < N; ++n) {
        //        cout << "n = " << n << endl;
        for (size_t t = sim_data.first_year_in_Data(n)+1; t < para.EstTbar; ++t) {

            //// maximum likelihood function when the previous period is production
            double f_P; double df_d0_P; double df_d1_P; double df_d2_P;
            double VV_PD_P;
            double VV_PD_P_temp = (log(ValueAdd(n, t)) - para_est2.F_P_PP) / para_est2.sigma_PD;
            VV_PD_P = min(VV_PD_P_temp, 7.5);
            VV_PD_P = max(VV_PD_P, -7.5);

            double lnProbP_P = log(normCDF(VV_PD_P, 0, 1));
            double dlnProbP_P_d0 = normPDF(VV_PD_P, 0, 1) / normCDF(VV_PD_P, 0, 1)
                * (-1.0 / para_est2.sigma_PD);
            double dlnProbP_P_d2 = normPDF(VV_PD_P, 0, 1) / normCDF(VV_PD_P, 0, 1)
                * VV_PD_P * (-1.0 / para_est2.sigma_PD);

            double lnProbP_D = log(normCDF(-VV_PD_P, 0, 1));
            double dlnProbP_D_d0 = normPDF(-VV_PD_P, 0, 1) / normCDF(-VV_PD_P, 0, 1)
                * (1.0 / para_est2.sigma_PD);
            double dlnProbP_D_d2 = normPDF(-VV_PD_P, 0, 1) / normCDF(-VV_PD_P, 0, 1)
                * (-VV_PD_P) * (-1.0 / para_est2.sigma_PD);

//            f_P = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * lnProbP_P
//                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * lnProbP_D;
//            df_d0_P = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d0
//                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d0;
//            df_d1_P = 0.0;
//            df_d2_P = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d2
//                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d2;
//
//            f_P = sim_data.Prob_P(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * lnProbP_P
//                  + sim_data.Prob_P(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * lnProbP_D;
//            df_d0_P = sim_data.Prob_P(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d0
//                      + sim_data.Prob_P(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d0;
//            df_d1_P = 0.0;
//            df_d2_P = sim_data.Prob_P(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d2
//                      + sim_data.Prob_P(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d2;

//            f_P_mat(n,t) = sim_data.ProbP_P_nonMiss(n, t) * lnProbP_P
//                + sim_data.ProbP_D_nonMiss(n, t) * lnProbP_D;
//            df_d0_P_mat(n,t) = sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d0
//                + sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d0;
//            df_d1_P_mat(n,t) = 0.0;
//            df_d2_P_mat(n,t) = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d2
//                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d2;

            f_P = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * lnProbP_P
                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * lnProbP_D;
            df_d0_P = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d0
                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d0;
            df_d1_P = 0.0;
            df_d2_P = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t) * dlnProbP_P_d2
                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t) * dlnProbP_D_d2;

            //// maximum likelihood function when the previous period is dormancy
            double f_D; double df_d0_D; double df_d1_D; double df_d2_D;
            double VV_PD_D;
            double VV_PD_D_temp = (log(ValueAdd(n, t)) - para_est2.F_P_DP) / para_est2.sigma_PD;
            VV_PD_D = min(VV_PD_D_temp, 7.5);
            VV_PD_D = max(VV_PD_D, -7.5);

            double lnProbD_P = log(normCDF(VV_PD_D, 0, 1));
            double dlnProbD_P_d1 = normPDF(VV_PD_D, 0, 1) / normCDF(VV_PD_D, 0, 1)
                * (-1.0 / para_est2.sigma_PD);
            double dlnProbD_P_d2 = normPDF(VV_PD_D, 0, 1) / normCDF(VV_PD_D, 0, 1)
                * VV_PD_D * (-1.0 / para_est2.sigma_PD);

            double lnProbD_D = log(normCDF(-VV_PD_D, 0, 1));
            double dlnProbD_D_d1 = normPDF(-VV_PD_D, 0, 1) / normCDF(-VV_PD_D, 0, 1)
                * (1.0 / para_est2.sigma_PD);
            double dlnProbD_D_d2 = normPDF(-VV_PD_D, 0, 1) / normCDF(-VV_PD_D, 0, 1)
                * (-VV_PD_D) * (-1.0 / para_est2.sigma_PD);

//            f_D = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * lnProbD_P
//                  + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * lnProbD_D;
//            df_d0_D = 0.0;
//            df_d1_D = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d1
//                      + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d1;
//            df_d2_D = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d2
//                      + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d2;
//
//            f_D = sim_data.Prob_D(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * lnProbD_P
//                  + sim_data.Prob_D(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * lnProbD_D;
//            df_d0_D = 0.0;
//            df_d1_D = sim_data.Prob_D(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d1
//                      + sim_data.Prob_D(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d1;
//            df_d2_D = sim_data.Prob_D(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d2
//                      + sim_data.Prob_D(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d2;

//            f_D_mat(n,t) = sim_data.ProbD_P_nonMiss(n, t) * lnProbD_P
//                  + sim_data.ProbD_D_nonMiss(n, t) * lnProbD_D;
//            df_d0_D_mat(n,t) = 0.0;
//            df_d1_D_mat(n,t) = sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d1
//                      + sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d1;
//            df_d2_D_mat(n,t) = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d2
//                      + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d2;

            f_D = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * lnProbD_P
                + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * lnProbD_D;
            df_d0_D = 0.0;
            df_d1_D = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d1
                + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d1;
            df_d2_D = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t) * dlnProbD_P_d2
                + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t) * dlnProbD_D_d2;

//            //// maximum likelihood function when the previous period is dormant
//            ArrayXd f_D_prob_tt(para.EstTbar);
//            ArrayXd df_d0_D_prob_tt(para.EstTbar);
//            ArrayXd df_d1_D_prob_tt(para.EstTbar);
//            ArrayXd df_d2_D_prob_tt(para.EstTbar);
//            double prob_P_tt; double prob_D_tt;
//            for (int tt = t-2; tt >= sim_data.first_year_in_Data(n); --tt) {
//                int max_period = tt;
//                int t_lag = max_period;
//                int diff_period = t - max_period;
//
//                double VV_PD_D;
//                double VV_PD_D_temp = (log(ValueAdd(n, tt)) - para_est2.F_P_DP) / para_est2.sigma_PD;
//                VV_PD_D = min(VV_PD_D_temp, 7.5);
//                VV_PD_D = max(VV_PD_D, -7.5);
//                double lnProbD_P = log(normCDF(VV_PD_D, 0, 1));
//                double dlnProbD_P_d1 = normPDF(VV_PD_D, 0, 1) / normCDF(VV_PD_D, 0, 1)
//                    * (-1.0 / para_est2.sigma_PD);
//                double dlnProbD_P_d2 = normPDF(VV_PD_D, 0, 1) / normCDF(VV_PD_D, 0, 1)
//                    * VV_PD_D * (-1.0 / para_est2.sigma_PD);
//
//                double lnProbD_D = log(normCDF(-VV_PD_D, 0, 1));
//                double dlnProbD_D_d1 = normPDF(-VV_PD_D, 0, 1) / normCDF(-VV_PD_D, 0, 1)
//                    * (1.0 / para_est2.sigma_PD);
//                double dlnProbD_D_d2 = normPDF(-VV_PD_D, 0, 1) / normCDF(-VV_PD_D, 0, 1)
//                    * (-VV_PD_D) * (-1.0 / para_est2.sigma_PD);
//
//                if (tt == t - 2) {
//                    prob_P_tt = sim_data.Prob_P_nonMiss(n, t-2) * sim_data.ProbP_D_nonMiss(n, t-1)
//                            * sim_data.ProbD_P_nonMiss(n, t);
//                    prob_D_tt = sim_data.Prob_P_nonMiss(n, t-2) * sim_data.ProbP_D_nonMiss(n, t-1)
//                            * sim_data.ProbD_D_nonMiss(n, t);
//                }
//                else {
//                    double miss_tt = sim_data.ProbP_D_nonMiss(n, tt+1) + sim_data.ProbP_DP_Miss(n, tt+1);
//                    double prob_miss_P = sim_data.ProbP_P(n,tt);
//                    double prob_miss_D = sim_data.ProbP_D(n,tt);
//
//                    for (int ttt = tt+2; ttt < t-1; ++ttt) {
//                        miss_tt = miss_tt
//                            * (prob_miss_P * (sim_data.ProbP_D_nonMiss(n, ttt) + sim_data.ProbP_DP_Miss(n, ttt))
//                            + prob_miss_D * (sim_data.ProbD_D_nonMiss(n, ttt) + sim_data.ProbD_DP_Miss(n, ttt)));
//                        double prob_miss_P_update = prob_miss_P * sim_data.ProbP_P(n,ttt) * sim_data.ProbP_S(n,ttt)
//                            + prob_miss_D * sim_data.ProbD_P(n,ttt) * sim_data.ProbD_S(n,ttt);
//                        double prob_miss_D_update = prob_miss_P * sim_data.ProbP_D(n,ttt) * sim_data.ProbP_S(n,ttt)
//                            + prob_miss_D * sim_data.ProbD_D(n,ttt) * sim_data.ProbD_S(n,ttt);
//                        prob_miss_P = prob_miss_P_update;
//                        prob_miss_D = prob_miss_D_update;
//                    }
//
//                    double prob_pre_P = prob_miss_P;
//                    double prob_pre_D = prob_miss_D;
//                    prob_P_tt = sim_data.Prob_P_nonMiss(n, tt) * miss_tt
//                        * (prob_pre_P * sim_data.ProbP_D_nonMiss(n, t-1) + prob_pre_D * sim_data.ProbD_D_nonMiss(n, t-1))
//                        * sim_data.ProbD_P_nonMiss(n, t);
//                    prob_D_tt = sim_data.Prob_P_nonMiss(n, tt) * miss_tt
//                        * (prob_pre_P * sim_data.ProbP_D_nonMiss(n, t-1) + prob_pre_D * sim_data.ProbD_D_nonMiss(n, t-1))
//                        * sim_data.ProbD_D_nonMiss(n, t);
//                }
//
//                f_D_prob_tt(tt) = prob_P_tt * lnProbD_P + prob_D_tt * lnProbD_D;
//                df_d0_D_prob_tt(tt) = 0.0;
//                df_d1_D_prob_tt(tt) = prob_P_tt * dlnProbD_P_d1 + prob_D_tt * dlnProbD_D_d1;
//                df_d2_D_prob_tt(tt) = prob_P_tt * dlnProbD_P_d2 + prob_D_tt * dlnProbD_D_d2;
//            }
//
//            double f_D = f_D_prob_tt.sum();
//            double df_d0_D = df_d0_D_prob_tt.sum();
//            double df_d1_D = df_d1_D_prob_tt.sum();
//            double df_d2_D = df_d2_D_prob_tt.sum();

            logL(n, t) = f_P + f_D;
            dlogL_d0(n, t) = df_d0_P + df_d0_D;
            dlogL_d1(n, t) = df_d1_P + df_d1_D;
            dlogL_d2(n, t) = df_d2_P + df_d2_D;

            Nprob(n,t) = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t)
                + sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t)
                + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t)
                + sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t);
        }
        //    }
    };
    MultiThreads::simple_parallel_for(worker, N, threadsManagement);
//
//    writeToCSVfile("dlogL_d0_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//        dlogL_d0.cast<double>().matrix());
//    writeToCSVfile("dlogL_d1_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//        dlogL_d1.cast<double>().matrix());
//    writeToCSVfile("dlogL_d2_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//        dlogL_d2.cast<double>().matrix());
//
//    writeToCSVfile("f_P_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   f_P_mat.cast<double>().matrix());
//    writeToCSVfile("df_d0_P_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   df_d0_P_mat.cast<double>().matrix());
//    writeToCSVfile("df_d1_P_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   df_d1_P_mat.cast<double>().matrix());
//    writeToCSVfile("df_d2_P_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   df_d2_P_mat.cast<double>().matrix());
//
//    writeToCSVfile("f_D_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   f_D_mat.cast<double>().matrix());
//    writeToCSVfile("df_d0_D_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   df_d0_D_mat.cast<double>().matrix());
//    writeToCSVfile("df_d1_D_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   df_d1_D_mat.cast<double>().matrix());
//    writeToCSVfile("df_d2_D_mat_State" + to_string(GoodState) + "_Sec" + to_string(LaborIntensive) + ".csv",
//                   df_d2_D_mat.cast<double>().matrix());
//
//    cout << "Ncount = " << N << endl;
////    throw runtime_error("312");

    int Ncount = N;
    double Nprob_sum = Nprob.sum();

    //    logL = (logL.isFinite()).select(logL, 0.0);
    double y_sum = -logL.sum();
    //    cout << "logL.sum() = " << logL.sum() << endl;

    std::vector<double> grad_sum(para.dim2);
    grad_sum[0] = -dlogL_d0.sum();
    //    cout << "dlogL_d0.sum() = " << dlogL_d0.sum() << endl;

    grad_sum[1] = -dlogL_d1.sum();
    //    cout << "dlogL_d1.sum() = " << dlogL_d1.sum() << endl;

    grad_sum[2] = -dlogL_d2.sum();
    //    cout << "dlogL_d2.sum() = " << dlogL_d2.sum() << endl;
    //    cout << "y = " << y << "; grad = " << grad[0] << "; " << grad[1] << endl;
    //    throw runtime_error("312");
    return tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd>(y_sum,grad_sum,Ncount,
        Nprob_sum,logL,dlogL_d0,dlogL_d1,dlogL_d2);
}

/********************************************************************************************************************
* Part3 of objection function in the second stage estimation
********************************************************************************************************************/
tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
    ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd>
    alias::EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(
    const SimData & sim_data, const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV,
    const EquStateVmat & Evalmat, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    int N = sim_data.good.rows();

    //// estimated productivity
    ArrayXXd lnphi_a_P = Backout_phi_a_PD(sim_data.Revenue_P,sim_data.Capital_P,
        sim_data.Employ_ur_P,sim_data.Employ_uc_P,
        para_est.alpha_tilde_K,para_est.alpha_tilde_L,para_est.alpha_Lr,para_est.alpha_Lc);
    ArrayXXd lnphi_a_D = Backout_phi_a_PD(sim_data.Revenue_D,sim_data.Capital_D,
        sim_data.Employ_ur_D,sim_data.Employ_uc_D,
        para_est.alpha_tilde_K,para_est.alpha_tilde_L,para_est.alpha_Lr,para_est.alpha_Lc);

    ArrayXXd K_opt_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Lur_opt_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Luc_opt_P_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXi K_index_max_P_mat = ArrayXXi::Zero(N,para.EstTbar);
    ArrayXXi Lur_index_max_P_mat = ArrayXXi::Zero(N,para.EstTbar);
    ArrayXXi Luc_index_max_P_mat = ArrayXXi::Zero(N,para.EstTbar);
//
    ArrayXXd K_opt_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Lur_opt_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Luc_opt_D_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXi K_index_max_D_mat = ArrayXXi::Zero(N,para.EstTbar);
    ArrayXXi Lur_index_max_D_mat = ArrayXXi::Zero(N,para.EstTbar);
    ArrayXXi Luc_index_max_D_mat = ArrayXXi::Zero(N,para.EstTbar);

    ArrayXXd ProbP_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd ProbP_E_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbP_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbP_E_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXd ProbD_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd ProbD_E_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbD_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbD_E_mat = ArrayXXd::Zero(N,para.EstTbar);

    auto worker = [&](size_t n, unsigned thread_id) {
    // for (size_t n = 0; n < N; ++n) {
        for (size_t t = sim_data.first_year_in_Data(n); t <= para.EstTbar - 1; ++t) {
            // cout << "n = " << n << "; t = " << t << endl;
            //// case 1 Production in the previous period
            PhiW phi_w_P = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_P(n,t)));
            LKState K_state_P = defineLKState(para_vec.vec_K, sim_data.Capital_P(n,t));
            LKState Lur_state_P = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur_P(n,t));
            LKState Luc_state_P = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc_P(n,t));

            tuple<double, double, double, double> t_SE_P = calStayExit_Prob_logit(para_est, para_vec, Evalmat,
                phi_w_P, K_state_P, Lur_state_P,Luc_state_P,0);
            ProbP_S_mat(n,t) = get<0>(t_SE_P);
            ProbP_E_mat(n,t) = get<1>(t_SE_P);
            lnProbP_S_mat(n,t) = get<2>(t_SE_P);
            lnProbP_E_mat(n,t) = get<3>(t_SE_P);

            tuple<double,double,double,int,int,int> t_Opt_P = cal_Opt_KLurLuc_Data(para_est,para_vec,Evalmat,
                phi_w_P,K_state_P, Lur_state_P,Luc_state_P,0);
            K_opt_P_mat(n,t) = get<0>(t_Opt_P);
            Lur_opt_P_mat(n,t) = get<1>(t_Opt_P);
            Luc_opt_P_mat(n,t) = get<2>(t_Opt_P);
            K_index_max_P_mat(n,t) = get<3>(t_Opt_P);
            Lur_index_max_P_mat(n,t) = get<4>(t_Opt_P);
            Luc_index_max_P_mat(n,t) = get<5>(t_Opt_P);

            //// case 2 Dormancy in the previous period
            PhiW phi_w_D = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_D(n,t)));
            LKState K_state_D = defineLKState(para_vec.vec_K, sim_data.Capital_D(n,t));
            LKState Lur_state_D = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur_D(n,t));
            LKState Luc_state_D = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc_D(n,t));

            tuple<double, double, double, double> t_SE_D = calStayExit_Prob_logit(para_est, para_vec, Evalmat,
                phi_w_D, K_state_D, Lur_state_D,Luc_state_D,1);
            ProbD_S_mat(n,t) = get<0>(t_SE_D);
            ProbD_E_mat(n,t) = get<1>(t_SE_D);
            lnProbD_S_mat(n,t) = get<2>(t_SE_D);
            lnProbD_E_mat(n,t) = get<3>(t_SE_D);

            tuple<double,double,double,int,int,int> t_Opt_D = cal_Opt_KLurLuc_Data(para_est,para_vec,Evalmat,
                phi_w_P,K_state_D, Lur_state_D,Luc_state_D,1);
            K_opt_D_mat(n,t) = get<0>(t_Opt_D);
            Lur_opt_D_mat(n,t) = get<1>(t_Opt_D);
            Luc_opt_D_mat(n,t) = get<2>(t_Opt_D);
            K_index_max_D_mat(n,t) = get<3>(t_Opt_D);
            Lur_index_max_D_mat(n,t) = get<4>(t_Opt_D);
            Luc_index_max_D_mat(n,t) = get<5>(t_Opt_D);
        }
    // }
    };
    MultiThreads::simple_parallel_for(worker,N, threadsManagement);
    // throw runtime_error("1962");
    // cout << "EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2_likelihood" << endl;
    tuple<ArrayXXd,double> t_logLikelihood = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2_likelihood(
        sim_data,para_est,para_vec,Evalmat,lnphi_a_P,lnphi_a_D,
        ProbP_S_mat,ProbP_E_mat,lnProbP_S_mat,lnProbP_E_mat,
        ProbD_S_mat,ProbD_E_mat,lnProbD_S_mat,lnProbD_E_mat,
        K_opt_P_mat,Lur_opt_P_mat,Luc_opt_P_mat,
        K_opt_D_mat,Lur_opt_D_mat,Luc_opt_D_mat,
        GoodState,LaborIntensive,threadsManagement);
    ArrayXXd logLikelihood = get<0>(t_logLikelihood);
    double Nprob = get<1>(t_logLikelihood);

    double y_sum = logLikelihood.sum();
    int Ncount = N;

    return tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd>(y_sum,Ncount,Nprob,logLikelihood,
        K_index_max_P_mat,Lur_index_max_P_mat,Luc_index_max_P_mat,K_opt_P_mat,Lur_opt_P_mat,Luc_opt_P_mat,
        K_index_max_D_mat,Lur_index_max_D_mat,Luc_index_max_D_mat,K_opt_D_mat,Lur_opt_D_mat,Luc_opt_D_mat);
}



tuple<double,ArrayXXd> alias::EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(const SimData & sim_data,
    const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV, const EquStateVmat & Evalmat,
    const ArrayXXi & K_index_max_P_mat, const ArrayXXi & Lur_index_max_P_mat, const ArrayXXi & Luc_index_max_P_mat,
    const ArrayXXd & K_max_P_mat, const ArrayXXd & Lur_max_P_mat, const ArrayXXd & Luc_max_P_mat,
    const ArrayXXi & K_index_max_D_mat, const ArrayXXi & Lur_index_max_D_mat, const ArrayXXi & Luc_index_max_D_mat,
    const ArrayXXd & K_max_D_mat, const ArrayXXd & Lur_max_D_mat, const ArrayXXd & Luc_max_D_mat,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {

    int N = sim_data.good.rows();

    //// estimated productivity
    ArrayXXd lnphi_a_P = Backout_phi_a_PD(sim_data.Revenue_P,sim_data.Capital_P,
        sim_data.Employ_ur_P,sim_data.Employ_uc_P,
        para_est.alpha_tilde_K,para_est.alpha_tilde_L,para_est.alpha_Lr,para_est.alpha_Lc);
    ArrayXXd lnphi_a_D = Backout_phi_a_PD(sim_data.Revenue_D,sim_data.Capital_D,
        sim_data.Employ_ur_D,sim_data.Employ_uc_D,
        para_est.alpha_tilde_K,para_est.alpha_tilde_L,para_est.alpha_Lr,para_est.alpha_Lc);

    ArrayXXd K_opt_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Lur_opt_P_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Luc_opt_P_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXd K_opt_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Lur_opt_D_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd Luc_opt_D_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXd ProbP_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd ProbP_E_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbP_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbP_E_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXd ProbD_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd ProbD_E_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbD_S_mat = ArrayXXd::Zero(N,para.EstTbar);
    ArrayXXd lnProbD_E_mat = ArrayXXd::Zero(N,para.EstTbar);

    ArrayXXi count = ArrayXXi::Zero(N,para.EstTbar);

    auto worker = [&](size_t n, unsigned thread_id) {
//    for (size_t n = 0; n < N; ++n) {
        for (size_t t = sim_data.first_year_in_Data(n); t <= para.EstTbar - 1; ++t) {
//            cout << "n = " << n << "; t = " << t << endl;
            //// case 1 Production in the previous period
            PhiW phi_w_P = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_P(n, t)));
            LKState K_state_P = defineLKState(para_vec.vec_K, sim_data.Capital_P(n, t));
            LKState Lur_state_P = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur_P(n, t));
            LKState Luc_state_P = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc_P(n, t));

            tuple<double, double, double, double> t_SE_P = calStayExit_Prob_logit(para_est, para_vec, Evalmat,
                phi_w_P, K_state_P, Lur_state_P,Luc_state_P, 0);
            ProbP_S_mat(n, t) = get<0>(t_SE_P);
            ProbP_E_mat(n, t) = get<1>(t_SE_P);
            lnProbP_S_mat(n, t) = get<2>(t_SE_P);
            lnProbP_E_mat(n, t) = get<3>(t_SE_P);

            tuple<double, double, double> t_Opt_P = cal_Opt_KLurLuc_Data_Diff(para_est, para_vec, Evalmat, phi_w_P,
                K_state_P, Lur_state_P, Luc_state_P,
                K_index_max_P_mat(n, t),
                Lur_index_max_P_mat(n, t),
                Luc_index_max_P_mat(n, t), 0);
            K_opt_P_mat(n, t) = get<0>(t_Opt_P);
            Lur_opt_P_mat(n, t) = get<1>(t_Opt_P);
            Luc_opt_P_mat(n, t) = get<2>(t_Opt_P);

            //// case 2 Dormancy in the previous period
            PhiW phi_w_D = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_D(n, t)));
            LKState K_state_D = defineLKState(para_vec.vec_K, sim_data.Capital_D(n, t));
            LKState Lur_state_D = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur_D(n, t));
            LKState Luc_state_D = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc_D(n, t));

            tuple<double, double, double, double> t_SE_D = calStayExit_Prob_logit(para_est, para_vec, Evalmat,
                phi_w_D, K_state_D, Lur_state_D,
                Luc_state_D, 1);
            ProbD_S_mat(n, t) = get<0>(t_SE_D);
            ProbD_E_mat(n, t) = get<1>(t_SE_D);
            lnProbD_S_mat(n, t) = get<2>(t_SE_D);
            lnProbD_E_mat(n, t) = get<3>(t_SE_D);

            tuple<double, double, double> t_Opt_D = cal_Opt_KLurLuc_Data_Diff(para_est, para_vec, Evalmat, phi_w_D,
                K_state_D, Lur_state_D, Luc_state_D,
                K_index_max_D_mat(n, t),
                Lur_index_max_D_mat(n, t),
                Luc_index_max_D_mat(n, t), 1);
            K_opt_D_mat(n, t) = get<0>(t_Opt_D);
            Lur_opt_D_mat(n, t) = get<1>(t_Opt_D);
            Luc_opt_D_mat(n, t) = get<2>(t_Opt_D);
        }
//    }
    };
    MultiThreads::simple_parallel_for(worker,N, threadsManagement);

    tuple<ArrayXXd,double> t_logLikelihood = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2_likelihood(
        sim_data,para_est,para_vec,Evalmat,lnphi_a_P,lnphi_a_D,
        ProbP_S_mat,ProbP_E_mat,lnProbP_S_mat,lnProbP_E_mat,
        ProbD_S_mat,ProbD_E_mat,lnProbD_S_mat,lnProbD_E_mat,
        K_opt_P_mat,Lur_opt_P_mat,Luc_opt_P_mat,
        K_opt_D_mat,Lur_opt_D_mat,Luc_opt_D_mat,
        GoodState,LaborIntensive,threadsManagement);
    ArrayXXd logLikelihood =  get<0>(t_logLikelihood);
    double y_sum = logLikelihood.sum();
    int Ncount = N;
    return tuple<double,ArrayXXd>(y_sum,logLikelihood);
}

tuple<ArrayXXd,double> alias::EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2_likelihood(
    const SimData & sim_data, const ParaEst & para_est, const ParaVec & para_vec, const EquStateVmat & Evalmat,
    const ArrayXXd & lnphi_a_P, const ArrayXXd & lnphi_a_D,
    const ArrayXXd & ProbP_S_mat, const ArrayXXd & ProbP_E_mat, const ArrayXXd & lnProbP_S_mat, const ArrayXXd & lnProbP_E_mat,
    const ArrayXXd & ProbD_S_mat, const ArrayXXd & ProbD_E_mat, const ArrayXXd & lnProbD_S_mat, const ArrayXXd & lnProbD_E_mat,
    const ArrayXXd & K_opt_P_mat, const ArrayXXd & Lur_opt_P_mat, const ArrayXXd & Luc_opt_P_mat,
    const ArrayXXd & K_opt_D_mat, const ArrayXXd & Lur_opt_D_mat, const ArrayXXd & Luc_opt_D_mat,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {

    int N = sim_data.good.rows();
    ArrayXXd logLikelihood = ArrayXXd::Zero(N, para.EstTbar);

    ArrayXXd Nprob = ArrayXXd::Zero(N, para.EstTbar);

    auto worker = [&](size_t n, unsigned thread_id) {
    //    size_t thread_id = 0;
    // for (size_t n = 0; n < N; ++n) {
        for (size_t t = sim_data.first_year_in_Data(n) + 1; t <= para.EstTbar - 1; ++t) {
            // cout << "n = " << n << "; t = " << t << endl;
            //// Production in the previous period
            PhiW phi_w_P = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_P(n, t-1)));
            double ProbP_S = ProbP_S_mat(n, t-1);
            double ProbP_E = ProbP_E_mat(n, t-1);
            double lnProbP_S = lnProbP_S_mat(n, t-1);
            double lnProbP_E = lnProbP_E_mat(n, t-1);

            tuple<double,double,double,double,double,double> t_PD_P = calProductionDormancy_Prob_lognormal(para_est, para_vec,
                Evalmat, 0, phi_w_P.phi,
                sim_data.Capital_P(n, t),sim_data.Employ_ur_P(n, t),
                sim_data.Employ_uc_P(n, t));
            double ProbP_P = get<0>(t_PD_P);
            double ProbP_D = get<1>(t_PD_P);
            double lnProbP_P = get<2>(t_PD_P);
            double lnProbP_D = get<3>(t_PD_P);

            double lnpdf_State_K_P = -log(para_est.sigma_Kerror) - 0.5 * log(2.0 * 3.1415926)
                - 0.5 * pow((log(sim_data.Capital_P(n, t)) - log(K_opt_P_mat(n, t-1))) / para_est.sigma_Kerror, 2);
            double lnpdf_State_Lur_P = -log(para_est.sigma_Lerror_ur) - 0.5 * log(2.0 * 3.1415926)
                - 0.5 * pow((log(sim_data.Employ_ur_P(n, t)) - log(Lur_opt_P_mat(n, t-1))) / para_est.sigma_Lerror_ur, 2);
            double lnpdf_State_Luc_P = -log(para_est.sigma_Lerror_uc) - 0.5 * log(2.0 * 3.1415926)
                - 0.5 * pow((log(sim_data.Employ_uc_P(n, t)) - log(Luc_opt_P_mat(n, t-1))) / para_est.sigma_Lerror_uc, 2);

            //// Dormancy in the previous period
            PhiW phi_w_D = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_D(n, t-1)));
            double ProbD_S = ProbD_S_mat(n, t-1);
            double ProbD_E = ProbD_E_mat(n, t-1);
            double lnProbD_S = lnProbD_S_mat(n, t-1);
            double lnProbD_E = lnProbD_E_mat(n, t-1);

            tuple<double,double,double,double,double,double> t_PD_D = calProductionDormancy_Prob_lognormal(para_est, para_vec,
                Evalmat, 1, phi_w_D.phi,
                sim_data.Capital_D(n, t),sim_data.Employ_ur_D(n, t),
                sim_data.Employ_uc_D(n, t));
            double ProbD_P = get<0>(t_PD_P);
            double ProbD_D = get<1>(t_PD_P);
            double lnProbD_P = get<2>(t_PD_P);
            double lnProbD_D = get<3>(t_PD_P);

            double lnpdf_State_K_D = -log(para_est.sigma_Kerror) - 0.5 * log(2.0 * 3.1415926)
                - 0.5 * pow((log(sim_data.Capital_D(n, t)) - log(K_opt_D_mat(n, t-1))) / para_est.sigma_Kerror, 2);
            double lnpdf_State_Lur_D = -log(para_est.sigma_Lerror_ur) - 0.5 * log(2.0 * 3.1415926)
                - 0.5 * pow((log(sim_data.Employ_ur_D(n, t)) - log(Lur_opt_D_mat(n, t-1))) / para_est.sigma_Lerror_ur, 2);
            double lnpdf_State_Luc_D = -log(para_est.sigma_Lerror_uc) - 0.5 * log(2.0 * 3.1415926)
                - 0.5 * pow((log(sim_data.Employ_uc_D(n, t)) - log(Luc_opt_D_mat(n, t-1))) / para_est.sigma_Lerror_uc, 2);

            /////// Case 1: sim_data.State_P_nonMiss(n,t-1) = 1
            //// 1.1 production in period t
            double tempP_P_nonMiss = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_P_nonMiss(n, t-1)
                * (lnpdf_State_K_P + lnpdf_State_Lur_P + lnpdf_State_Luc_P + lnProbP_S + lnProbP_P);
            //// 1.2 dormancy in period t
            double tempP_D_nonMiss = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_D_nonMiss(n, t-1)
                * (lnpdf_State_K_P + lnpdf_State_Lur_P + lnpdf_State_Luc_P + lnProbP_S + lnProbP_D);
            //// 1.3 PD missing in period t
            // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
            double Prob_PS_Missing = ProbP_S * sim_data.Prob_Missing_P1(n,t);
            double Prob_P_Missing_cum = ProbP_E;
            for (size_t tt = t; tt < para.SimTbar; ++tt) {
                Prob_P_Missing_cum = Prob_P_Missing_cum + Prob_PS_Missing * ProbP_E;
                Prob_PS_Missing = Prob_PS_Missing * ProbP_S * sim_data.Prob_Missing_DP1(n,t);
            }
            Prob_P_Missing_cum = Prob_P_Missing_cum + Prob_PS_Missing;
            double tempP_E_Miss = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_E_Miss(n, t-1)
                * log(Prob_P_Missing_cum);
            //// 1.4 PD missing in period t
            double tempP_DP_Miss;
            double Prob_S_Missing_P = ProbP_S * sim_data.Prob_Missing_P1(n,t);
            double Prob_S_Missing_next_P = ProbP_S * sim_data.Prob_Missing_DP1(n,t);
            double Prob_Missing_next_cum_P = ProbP_E;
            // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
            for (size_t tt = t+1; tt < para.SimTbar; ++tt) {
                Prob_Missing_next_cum_P = Prob_Missing_next_cum_P + Prob_S_Missing_next_P * ProbP_E;
                Prob_S_Missing_next_P = Prob_S_Missing_next_P * ProbP_S * sim_data.Prob_Missing_DP1(n,t);
            }
            Prob_Missing_next_cum_P = Prob_Missing_next_cum_P + Prob_S_Missing_next_P;
            tempP_DP_Miss = sim_data.Prob_P_nonMiss(n, t-1) * sim_data.ProbP_DP_Miss(n, t-1)
                    * log(Prob_S_Missing_P*(1-Prob_Missing_next_cum_P));


            /////// Case 2: sim_data.State_D_nonMiss(n,t-1) = 1
            //// 1.1 production in period t
            double tempD_P_nonMiss = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_P_nonMiss(n, t-1)
                * (lnpdf_State_K_D + lnpdf_State_Lur_D + lnpdf_State_Luc_D + lnProbD_S + lnProbD_P);
            //// 1.2 dormancy in period t
            double tempD_D_nonMiss = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_D_nonMiss(n, t-1)
                * (lnpdf_State_K_D + lnpdf_State_Lur_D + lnpdf_State_Luc_D + lnProbD_S + lnProbD_D);
            //// 1.3 PD missing in period t
            // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
            double Prob_DS_Missing = ProbD_S * sim_data.Prob_Missing_D1(n,t);
            double Prob_D_Missing_cum = ProbP_E;
            for (size_t tt = t; tt < para.SimTbar; ++tt) {
                Prob_D_Missing_cum = Prob_D_Missing_cum + Prob_DS_Missing * ProbD_E;
                Prob_DS_Missing = Prob_DS_Missing * ProbD_S * sim_data.Prob_Missing_DP1(n,t);
            }
            Prob_D_Missing_cum = Prob_D_Missing_cum + Prob_DS_Missing;
            double tempD_E_Miss = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbD_E_Miss(n, t)
                                  * log(Prob_D_Missing_cum);
            //// 1.4 PD missing in period t
            double tempD_DP_Miss;
            double Prob_S_Missing_D = ProbD_S * sim_data.Prob_Missing_D1(n,t);
            double Prob_S_Missing_next_D = ProbD_S * sim_data.Prob_Missing_DP1(n,t);
            double Prob_Missing_next_cum_D = ProbD_E;
            // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
            for (size_t tt = t+1; tt < para.SimTbar; ++tt) {
                Prob_Missing_next_cum_D = Prob_Missing_next_cum_D + Prob_S_Missing_next_D * ProbD_E;
                Prob_S_Missing_next_D = Prob_S_Missing_next_D * ProbD_S * sim_data.Prob_Missing_DP1(n,t);
            }
            Prob_Missing_next_cum_D = Prob_Missing_next_cum_D + Prob_S_Missing_next_D;
            tempD_DP_Miss = sim_data.Prob_D_nonMiss(n, t-1) * sim_data.ProbP_DP_Miss(n, t)
                    * log(Prob_S_Missing_D*(1-Prob_Missing_next_cum_D));

            double f = tempP_D_nonMiss + tempP_P_nonMiss + tempP_DP_Miss + tempP_E_Miss
                    + tempD_D_nonMiss + tempD_P_nonMiss + tempD_DP_Miss + tempD_E_Miss;

            logLikelihood(n,t) = f;
            Nprob(n,t) = sim_data.Prob_P_nonMiss(n, t-1) + sim_data.Prob_D_nonMiss(n, t-1);
        }
    // }
    };
    MultiThreads::simple_parallel_for(worker,N, threadsManagement);

    double Nprob_sum = Nprob.sum();
    return tuple<ArrayXXd,double>(logLikelihood,Nprob_sum);
}


////**** Calculate the probability of Lur/Luc/K for each data point  ****////
//// the main function to calculate the probability of Lur/Luc/K
tuple<double,double,double,int,int,int> alias::cal_Opt_KLurLuc_Data( const ParaEst & para_est, const ParaVec & para_vec,
    const EquStateVmat & Evalmat, const PhiW & phi_w, const LKState & K_state, const LKState & Lur_state,
    const LKState & Luc_state, const int & PD_lag) {

    ArrayXXd EVal_K_mat = Evalmat.EVal_PD_K_mat(seqN(PD_lag*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
    tuple<double,int> t_K = sol_opt_K(para_est, para_vec, phi_w, EVal_K_mat, Lur_state,Luc_state,K_state.val_state,0);
    double K_opt = get<0>(t_K);
    double K_index_max = get<1>(t_K);

    ArrayXXd EVal_Lur_mat = Evalmat.EVal_PD_Lur_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
    tuple<double,int> t_Lur = sol_opt_Lur(para_est, para_vec, phi_w, EVal_Lur_mat, K_state,Luc_state,Lur_state.val_state,0);
    double Lur_opt = get<0>(t_Lur);
    double Lur_index_max = get<1>(t_Lur);

    ArrayXXd EVal_Luc_mat = Evalmat.EVal_PD_Luc_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
    tuple<double,int> t_Luc = sol_opt_Luc(para_est, para_vec, phi_w, EVal_Luc_mat, K_state,Lur_state,Luc_state.val_state,0);
    double Luc_opt = get<0>(t_Luc);
    double Luc_index_max = get<1>(t_Luc);

    return tuple<double,double,double,int,int,int> (K_opt,Lur_opt,Luc_opt,K_index_max,Lur_index_max,Luc_index_max);
}

tuple<double,double,double> alias::cal_Opt_KLurLuc_Data_Diff( const ParaEst & para_est, const ParaVec & para_vec,
    const EquStateVmat & Evalmat, const PhiW & phi_w,
    const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state,
    const int & K_index_max, const int & Lur_index_max, const int & Luc_index_max, const int & PD_lag) {

    ArrayXXd EVal_K_mat = Evalmat.EVal_PD_K_mat(seqN(PD_lag*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
    double K_opt = sol_opt_K_Diff_simulation(para_est, para_vec, phi_w, EVal_K_mat,
        Lur_state,Luc_state,K_state.val_state,K_index_max);

    ArrayXXd EVal_Lur_mat = Evalmat.EVal_PD_Lur_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
    double Lur_opt = sol_opt_Lur_Diff_simulation(para_est, para_vec, phi_w, EVal_Lur_mat,
        K_state,Luc_state,Lur_state.val_state,Lur_index_max);

    ArrayXXd EVal_Luc_mat = Evalmat.EVal_PD_Luc_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
    double Luc_opt = sol_opt_Luc_Diff_simulation(para_est, para_vec, phi_w, EVal_Luc_mat,
        K_state,Lur_state,Luc_state.val_state,Luc_index_max);

    return tuple<double,double,double> (K_opt,Lur_opt,Luc_opt);
}

////**** solve the optimal K; take the original optimal value as a reference ****////
double alias::sol_opt_K_Diff_simulation(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd EVal_K_mat, const LKState & Lur_state, const LKState & Luc_state, const double & K_state_val,
    const int & K_index_max) {

    ArrayXd EVal_K_vec = cal_EVal_K_phi_interpolation(EVal_K_mat, phi_w, Lur_state, Luc_state, para.N_K);
//    cout << "EVal_K_vec = " << EVal_K_vec.transpose() << endl;

    double rstar = exp(0.5*pow(para_est.sigma_Kerror,2)) * para_vec.p_K * para_est.r_K;
    double Kbase = para.deltaK * K_state_val;

    double H_Kbase = para_est.H_K / pow(Kbase,2) + para_est.c_HK * log(Kbase+1.0) / pow(Kbase,2)
        + para_est.H_K_Cons * para_vec.p_K * para_vec.p_K;
    double F_Kbase = para_est.F_K / pow(Kbase,2) + para_est.c_FK * log(Kbase+1.0) / pow(Kbase,2)
        + para_est.F_K_Cons * para_vec.p_K * para_vec.p_K;

    double c_K = 1.0;

    tuple<double,double> t_K = solve1D_K_LinearSpline_Diff_simulation(H_Kbase,F_Kbase,
        c_K,rstar,para_vec.p_K,Kbase,para_vec.vec_K,EVal_K_vec,K_index_max);
    double K = get<0>(t_K);
    double OptV = get<1>(t_K);

//    cout << "vec_K = " << para_vec.vec_K.transpose() << endl;
//    cout << "369 K_index_max = " << K_index_max << "; K = " << K << endl;
    return K;
}

////**** solve the optimal Lur; take the original optimal value as a reference ****////
double alias::sol_opt_Lur_Diff_simulation(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd & EVal_Lur_mat, const LKState & K_state, const LKState & Luc_state, const double & Lur_state_val,
    const int & Lur_index_max) {

    ArrayXd EVal_Lur_vec;
    EVal_Lur_vec = cal_EVal_Lur_phi_interpolation(EVal_Lur_mat, phi_w, K_state, Luc_state, para.N_Lur);

    //// Value of production
    double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
    double Lbase = Lur_state_val;

    double H_Lbase_ur = para_est.H_ur / pow(Lbase,2) + para_est.c_H_ur*log(Lbase+1.0) / pow(Lbase,2) + para_est.H_ur_Cons;
    double F_Lbase_ur = 0;
    if (Lbase < para.Lur_cutoff1) {
        F_Lbase_ur = para_est.F_low_ur / pow(Lbase,2) + para_est.c_low_F_ur * log(Lbase+1.0) / pow(Lbase,2)
            + para_est.F_low_ur_Cons;
    }
    else {
        F_Lbase_ur = para_est.F_high_ur / pow(Lbase,2) + para_est.c_high_F_ur * log(Lbase+1.0) / pow(Lbase,2)
            + para_est.F_high_ur_Cons;
    }

    tuple<double,double> t_Lur = solve1D_L_LinearSpline_Diff_simulation(H_Lbase_ur, F_Lbase_ur, wstar, Lbase,
        para_vec.vec_Lur, EVal_Lur_vec, Lur_index_max);
    double Lur = get<0>(t_Lur);
    double TotVOpt = get<1>(t_Lur);

    return Lur;
}

////**** solve the optimal Luc; take the original optimal value as a reference ****////
double alias::sol_opt_Luc_Diff_simulation(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd & EVal_Luc_mat, const LKState & K_state, const LKState & Lur_state, const double & Luc_state_val,
    const int & Luc_index_max) {

    ArrayXd EVal_Luc_vec= cal_EVal_Luc_phi_interpolation(EVal_Luc_mat, phi_w, K_state, Lur_state, para.N_Luc);

    double wstar = para_est.w_uc * exp(0.5 * pow(para_est.sigma_Lerror_uc, 2));
    double Lbase = Luc_state_val;

    double H_Lbase_uc = para_est.H_uc / pow(Lbase,2) + para_est.c_H_uc * log(Lbase+1.0) / pow(Lbase,2) + para_est.H_uc_Cons;
    double F_Lbase_uc = para_est.F_uc / pow(Lbase,2) + para_est.c_F_uc * log(Lbase+1.0) / pow(Lbase,2) + para_est.F_uc_Cons;

    tuple<double,double> t_Luc = solve1D_L_LinearSpline_Diff_simulation(H_Lbase_uc, F_Lbase_uc, wstar, Lbase,
        para_vec.vec_Luc, EVal_Luc_vec, Luc_index_max);
    double Luc = get<0>(t_Luc);
    double TotVOpt = get<1>(t_Luc);

    return Luc;
}

/** the solver for sol_opt_K_Diff_simulation/sol_opt_Lur_Diff_simulation/sol_opt_Luc_Diff_simulation **/
tuple<double,double> alias::solve1D_K_LinearSpline_Diff_simulation(const double & H, const double & F, const double & c_K,
    const double & rstar, const double & p_K, const double & Kbase, const ArrayXd & vec_K, const ArrayXd & EVal_K_vec,
    const int & K_index_max) {

    ArrayXd diff_vec_K = vec_K - Kbase;
    ArrayXd TotV = EVal_K_vec - H * diff_vec_K.pow(2) * (diff_vec_K > 0).cast<double>()
        - F * diff_vec_K.pow(2) * (diff_vec_K <= 0).cast<double>()
        - rstar * vec_K
        - p_K * (vec_K - Kbase) * (diff_vec_K > 0).cast<double>()
        - p_K * c_K * (vec_K - Kbase) * (diff_vec_K <= 0).cast<double>();

    double Vdefault = TotV(K_index_max);
    double Kdefault = vec_K(K_index_max);

    double K_opt = vec_K(K_index_max);
    double v_max = Vdefault;

    int K_index = K_index_max;
    if (K_index_max == para.N_K-1) { K_index = K_index_max - 1; }

    ArrayXd Kvec = vec_K.segment(K_index, 2);
    ArrayXd coef_EVal = CalLinearspline1D_Coeff(vec_K(K_index),vec_K(K_index+1),
        EVal_K_vec(K_index),EVal_K_vec(K_index+1));

    tuple<double, double> t_max = SolveMax_Linear1D_K(coef_EVal, vec_K(K_index),
        vec_K(K_index+1),TotV(K_index),TotV(K_index+1),
        H, F, c_K, rstar, p_K, Kbase);
    double K_opt_temp = get<0>(t_max);
    double v_max_temp = get<1>(t_max);

    if (v_max_temp > v_max) {
        K_opt = K_opt_temp; v_max = v_max_temp;
    }
    return tuple<double,double>(K_opt,v_max);
}

tuple<double,double> alias::solve1D_L_LinearSpline_Diff_simulation(const double & H, const double & F, const double & wstar,
    const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const int & L_index_max) {

    int N_L = vec_L.size();

    ArrayXd diff_vec_L = vec_L - Lbase;
    ArrayXd TotV = EVal_L_vec - H * diff_vec_L.pow(2) * (diff_vec_L > 0).cast<double>()
        - F * diff_vec_L.pow(2) * (diff_vec_L <= 0).cast<double>()
        - wstar * vec_L;

    double Vdefault = TotV(L_index_max);
    double Kdefault = vec_L(L_index_max);

    double L_opt = vec_L(L_index_max);
    double v_max = Vdefault;

    int L_index = L_index_max;
    if (L_index_max == N_L-1) { L_index = L_index_max - 1; }

    ArrayXd Lvec = vec_L.segment(L_index, 2);
    ArrayXd coef_EVal = CalLinearspline1D_Coeff(vec_L(L_index),vec_L(L_index+1),
        EVal_L_vec(L_index),EVal_L_vec(L_index+1));

    tuple<double, double> t_max = SolveMax_Linear1D_L(coef_EVal, vec_L(L_index),
        vec_L(L_index+1), TotV(L_index),TotV(L_index+1), H, F,
        wstar, Lbase);
    double L_opt_temp = get<0>(t_max);
    double v_max_temp = get<1>(t_max);

    if (v_max_temp > v_max) {
        L_opt = L_opt_temp; v_max = v_max_temp;
    }
    return tuple<double,double>(L_opt,v_max);
}


////
////////void alias::EstimationIndirectInference_checking(const SimVar & sim_var,
////////    const ArrayXd & theta1_Para_ini, const ArrayXd & theta2_Para_ini, const ArrayXd & theta3_Para_ini,
////////    const ArrayXd & theta_a, const EquV_Auxiliary & EquV_a, const SimData & RealData, const int & GoodState,
////////    const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {
////////
////////    //// checking objective function
////////    std::vector<double> x0_vec(para.dim1 + para.dim2 + para.dim3);
////////    x0_vec[0] = theta1_Para_ini(0); x0_vec[1] = theta1_Para_ini(1); x0_vec[2] = theta1_Para_ini(2);
////////    x0_vec[3] = theta1_Para_ini(3); x0_vec[4] = theta1_Para_ini(4);
////////
////////    x0_vec[5] = theta2_Para_ini(0); x0_vec[6] = theta2_Para_ini(1);
////////
////////    x0_vec[7] = theta3_Para_ini(0); x0_vec[8] = theta3_Para_ini(1); x0_vec[9] = theta3_Para_ini(2);
////////    x0_vec[10] = theta3_Para_ini(3); x0_vec[11] = theta3_Para_ini(4); x0_vec[12] = theta3_Para_ini(5);
////////    x0_vec[13] = theta3_Para_ini(6); x0_vec[14] = theta3_Para_ini(7); x0_vec[15] = theta3_Para_ini(8);
////////    x0_vec[16] = theta3_Para_ini(9); x0_vec[17] = theta3_Para_ini(10); x0_vec[18] = theta3_Para_ini(11);
////////    x0_vec[19] = theta3_Para_ini(12); x0_vec[20] = theta3_Para_ini(13); x0_vec[21] = theta3_Para_ini(14);
////////    x0_vec[22] = theta3_Para_ini(15); x0_vec[23] = theta3_Para_ini(16); x0_vec[24] = theta3_Para_ini(17);
////////
////////    int Ntemp = 5;
////////    ArrayXXd Obj_vec(Ntemp,para.dim1 + para.dim2 + para.dim3);
////////    ArrayXXd ParaSpace_vec(Ntemp,para.dim1 + para.dim2 + para.dim3);
////////    for (size_t i = 0; i < para.dim1 + para.dim2 + para.dim3; ++i) {
////////        ArrayXd tempvec0 = ArrayXd::LinSpaced(Ntemp, x0_vec[i]*0.5, x0_vec[i]*5);
////////        for (size_t n = 0; n < Ntemp; ++n) {
////////            std::vector<double> x1_vec = x0_vec;
////////            x1_vec[i] = tempvec0(n);
////////            //// simulate data with the actual parameters
////////            Obj_vec(n,i) = EstimationIndirectInference_FullVersion_Obj(x1_vec,sim_var,
////////                theta_a, EquV_a, RealData, GoodState, LaborIntensive, threadsManagement);
////////
////////            cout << "x1_vec = " << x1_vec[0] << "; " << x1_vec[1] << "; " << x1_vec[2] << "; " << x1_vec[3] << "; "
////////                 << x1_vec[4] << "; " << x1_vec[5] << "; " << x1_vec[6] << "; " << x1_vec[7] << "; " << x1_vec[8] << "; "
////////                 << x1_vec[9] << "; " << x1_vec[10] << "; " << x1_vec[11] << "; " << x1_vec[12] << "; " << x1_vec[13] << "; "
////////                 << x1_vec[14] << "; " << x1_vec[15] << "; " << x1_vec[16] << "; " << x1_vec[17] << "; " << x1_vec[18] << "; "
////////                 << x1_vec[19] << "; " << x1_vec[20] << "; " << x1_vec[21] << "; " << x1_vec[22] << "; " << x1_vec[23] << "; "
////////                 << x1_vec[24] << endl;
////////            cout << "i = " << i << "; n = " << n << "; Obj_vec(n, i) = " << Obj_vec(n, i) << endl;
////////        }
////////        ParaSpace_vec.col(i) = tempvec0;
////////    }
////////    writeToCSVfile("Obj_II.csv", Obj_vec.matrix());
////////    writeToCSVfile("ParaSpace_II.csv", ParaSpace_vec.matrix());
////////    throw runtime_error("88");}
//
////
//
////
//
