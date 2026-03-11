#include "CounterfactualExercise.h"

using namespace Eigen;
using namespace std;
using namespace alias;

/**************************************************************
* **************************************************************
* Solve the Status quo equilibrium and simulate firm distribution
* **************************************************************
**************************************************************/
tuple<CounterFactVariable,CounterFactVariable,CounterFactVariable,CounterFactVariable>
    alias::main_SimulateStatusQuo(const ArrayXd & theta_Est_S, const SimVar & sim_var, const int & SimTbar,
    MultiThreads::Threads_Management & threadsManagement) {

    int LaborIntensive; int GoodState;
    /*** solve the equilibrium and simulate the data ***/
    LaborIntensive = 1; GoodState = 1;
    tuple<ParaEst,ParaVec,EquStateV,EquStateV0,EquState0,SimData> t_StatusQuo_GoodLaborInt
            = SolveSimulate_StatusQuo_FEntry_V0(theta_Est_S,sim_var,GoodState,LaborIntensive,threadsManagement);
    EquStateV EquV_GoodLaborInt_S = get<2>(t_StatusQuo_GoodLaborInt);
    EquState0 Equ0_GoodLaborInt_S = get<4>(t_StatusQuo_GoodLaborInt);
    SimData sim_data_GoodLaborInt_S = get<5>(t_StatusQuo_GoodLaborInt);
    // CounterFactVariable Value_StatusQuo_GoodLaborInt = CalEconomicVariable(Equ0_GoodLaborInt_S.FirmMass,
    //     Equ0_GoodLaborInt_S.F_Entry,CapitalDemand_GoodLaborInt_S,sim_data_GoodLaborInt_S);
    //
    // /*** Average Exit Rate by Age ***/
    // ArrayXd AvgExitRate_byAge_GoodLaborInt(SimTbar-1);
    // for (size_t t = 0; t < SimTbar-1; ++t) {
    //     AvgExitRate_byAge_GoodLaborInt(t) = sim_data_GoodLaborInt_S.State_E.col(t+1).cast<double>().sum()
    //         / sim_data_GoodLaborInt_S.State_S.col(t).cast<double>().sum();
    // }
    //
    // /*** solve the equilibrium and simulate the data ***/
    // LaborIntensive = 1; GoodState = 0;
    // tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,SimData> t_StatusQuo_BadLaborInt
    //     = SolveSimulate_StatusQuo_FEntry_V0(theta_Est_S,sim_var,GoodState,LaborIntensive,threadsManagement);
    // EquStateV EquV_BadLaborInt_S = get<2>(t_StatusQuo_BadLaborInt);
    // EquState0 Equ0_BadLaborInt_S = get<6>(t_StatusQuo_BadLaborInt);
    // double CapitalDemand_BadLaborInt_S = get<7>(t_StatusQuo_BadLaborInt);
    // SimData sim_data_BadLaborInt_S = get<8>(t_StatusQuo_BadLaborInt);
    // CounterFactVariable Value_StatusQuo_BadLaborInt = CalEconomicVariable(p_K_BadLaborInt_S,
    //     Equ0_BadLaborInt_S.FirmMass,Equ0_BadLaborInt_S.F_Entry,CapitalDemand_BadLaborInt_S,
    //     sim_data_BadLaborInt_S);
    //
    // /*** Average Exit Rate by Age ***/
    // ArrayXd AvgExitRate_byAge_BadLaborInt(SimTbar-1);
    // for (size_t t = 0; t < SimTbar-1; ++t) {
    //     AvgExitRate_byAge_BadLaborInt(t) = sim_data_BadLaborInt_S.State_E.col(t+1).cast<double>().sum()
    //         / sim_data_BadLaborInt_S.State_S.col(t).cast<double>().sum();
    // }
    //
    // /*** solve the equilibrium and simulate the data ***/
    // LaborIntensive = 0; GoodState = 1;
    // tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,SimData> t_StatusQuo_GoodCapitalInt
    //     = SolveSimulate_StatusQuo_FEntry_V0(theta_Est_S,sim_var,GoodState,LaborIntensive,threadsManagement);
    // double p_K_GoodCapitalInt_S = para.p_K_StatusQuo;
    // EquStateV EquV_GoodCapitalInt_S = get<2>(t_StatusQuo_GoodCapitalInt);
    // EquState0 Equ0_GoodCapitalInt_S = get<6>(t_StatusQuo_GoodCapitalInt);
    // double CapitalDemand_GoodCapitalInt_S = get<7>(t_StatusQuo_GoodCapitalInt);
    // SimData sim_data_GoodCapitalInt_S = get<8>(t_StatusQuo_GoodCapitalInt);
    // CounterFactVariable Value_StatusQuo_GoodCapitalInt = CalEconomicVariable(p_K_GoodCapitalInt_S,
    //     Equ0_GoodCapitalInt_S.FirmMass,Equ0_GoodCapitalInt_S.F_Entry,CapitalDemand_GoodCapitalInt_S,
    //     sim_data_GoodCapitalInt_S);
    //
    // /*** Average Exit Rate by Age ***/
    // ArrayXd AvgExitRate_byAge_GoodCapitalInt(SimTbar-1);
    // for (size_t t = 0; t < SimTbar-1; ++t) {
    //     AvgExitRate_byAge_GoodCapitalInt(t) = sim_data_GoodCapitalInt_S.State_E.col(t+1).cast<double>().sum()
    //         / sim_data_GoodCapitalInt_S.State_S.col(t).cast<double>().sum();
    // }
    //
    // /*** solve the equilibrium and simulate the data ***/
    // LaborIntensive = 0; GoodState = 0;
    // tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,SimData> t_StatusQuo_BadCapitalInt
    //         = SolveSimulate_StatusQuo_FEntry_V0(theta_Est_S,sim_var,GoodState,LaborIntensive,threadsManagement);
    // double p_K_BadCapitalInt_S = para.p_K_StatusQuo;
    // EquStateV EquV_BadCapitalInt_S = get<2>(t_StatusQuo_BadCapitalInt);
    // EquState0 Equ0_BadCapitalInt_S = get<6>(t_StatusQuo_BadCapitalInt);
    // double CapitalDemand_BadCapitalInt_S = get<7>(t_StatusQuo_BadCapitalInt);
    // SimData sim_data_BadCapitalInt_S = get<8>(t_StatusQuo_BadCapitalInt);
    // CounterFactVariable Value_StatusQuo_BadCapitalInt = CalEconomicVariable(p_K_BadCapitalInt_S,
    //     Equ0_BadCapitalInt_S.FirmMass,Equ0_BadCapitalInt_S.F_Entry,CapitalDemand_BadCapitalInt_S,
    //     sim_data_BadCapitalInt_S);
    //
    // /*** Average Exit Rate by Age ***/
    // ArrayXd AvgExitRate_byAge_BadCapitalInt(SimTbar-1);
    // for (size_t t = 0; t < SimTbar-1; ++t) {
    //     AvgExitRate_byAge_BadCapitalInt(t) = sim_data_BadCapitalInt_S.State_E.col(t+1).cast<double>().sum()
    //         / sim_data_BadCapitalInt_S.State_S.col(t).cast<double>().sum();
    // }

    // throw runtime_error("97");

    CounterFactVariable Value_StatusQuo_GoodLaborInt;
    CounterFactVariable Value_StatusQuo_BadLaborInt;
    CounterFactVariable Value_StatusQuo_GoodCapitalInt;
    CounterFactVariable Value_StatusQuo_BadCapitalInt;

    return tuple<CounterFactVariable,CounterFactVariable,CounterFactVariable,CounterFactVariable>(
        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt, Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt);
}

tuple<ParaEst,ParaVec,EquStateV,EquStateV0,EquState0,SimData>
    alias::SolveSimulate_StatusQuo_FEntry_V0(const ArrayXd & theta_Est, const SimVar & sim_var,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {

    ParaEst para_est = constructParaEst(theta_Est,GoodState,LaborIntensive);
    para_est.PI = 1.0;
    ParaVec para_vec = constructParaFull(para_est);
    //// solve value function for t > 1
    EquStateV EquV = solveV(para_est,para_vec,threadsManagement);

    Eigen::Map<const Eigen::ArrayXXd> ProbPD_E_mat(EquV.ProbPD_E.data(),para.N_Lur,para.N_phi*para.N_K);
    cout << "ProbPD_E_mat at i_phi 0 = " << ProbPD_E_mat.middleCols(0,para.N_K) << endl;
    cout << "ProbPD_E_mat at i_phi 1 = " << ProbPD_E_mat.middleCols(para.N_K,para.N_K) << endl;

    Eigen::Map<const Eigen::ArrayXXd> ProbPD_P_mat(EquV.ProbPD_P.data(),para.N_Lur,para.N_phi*para.N_K);
    cout << "ProbPD_P_mat at i_phi 0 = " << ProbPD_P_mat.middleCols(0,para.N_K) << endl;
    cout << "ProbPD_P_mat at i_phi 1 = " << ProbPD_P_mat.middleCols(para.N_K,para.N_K) << endl;

    Eigen::Map<const Eigen::ArrayXXd> EVal_PD_mat(EquV.EVal_PD.data(),para.N_Lur,para.N_phi*para.N_K);
    cout << "EVal_PD_mat at i_phi 0 = " << EVal_PD_mat.middleCols(0,para.N_K) << endl;
    cout << "EVal_PD_mat at i_phi 1 = " << EVal_PD_mat.middleCols(para.N_K,para.N_K) << endl;
    // throw runtime_error("135");

    //// With free entry condition, solve value/policy function at the initial t = 1
    //// not done yet.
    EquStateV0 EquV0 = solveV0(para_est,para_vec,EquV.EVal_PD,threadsManagement);


    EquState0 Equ0;
    Equ0.F_Entry = SolveFEntry_StatusQuoEqu(para_vec, EquV0);
    Equ0.FirmMass = 1.0;
    cout << "F_Entry = " << Equ0.F_Entry << endl;
    throw runtime_error("135");

    //// Simulate the data to get the stationary distribution of capital employment
    /** simulate a series of random number **/
    int SimN = sim_var.lnphi_sim.rows();
    int SimTbar = sim_var.lnphi_sim.cols()-1;
    //
    // /** simulate the data under the counterfactual data **/
    SimData sim_data;
    // SimData sim_data = Simulation_StatusQuo_Counterfactual(para_est, para_vec, EquV, EValmat,
    //     EquV0, EVal0mat, sim_var, SimN, SimTbar, GoodState, LaborIntensive, threadsManagement);
    // double CapitalStock = (sim_data.Capital * sim_data.State_S.cast<double>()).sum() * Equ0.FirmMass;
    // cout << "CapitalStock = " << CapitalStock << endl;
// //
// /    ArrayXd avgExitProb = sim_data.Prob_E.cast<double>().colwise().sum() / sim_data.Prob_S.cast<double>().colwise().sum();
// /    cout << "avgExitProb = " << avgExitProb.transpose() << endl;
// ///
// /    ArrayXXd tempP = (sim_data.Prob_P(all,seqN(0,para.SimTbar-1)) > 0).cast<double>();
// /    ArrayXd avgProd_preProd =
// /            ( sim_data.Prob_P(all,seqN(0,para.SimTbar-1)) * tempP ).cast<double>().colwise().sum()
// /            / tempP.cast<double>().colwise().sum();
// /    cout << "avgProd_preProd = " << avgProd_preProd.transpose() << endl;
// ///
// /    ArrayXXd tempD = (sim_data.Prob_P(all,seqN(0,para.SimTbar-1)) > 0.00001).cast<double>();
// /    ArrayXd avgProd_preDorm =
// /            ( sim_data.Prob_P(all,seqN(0,para.SimTbar-1)) * tempD ).cast<double>().colwise().sum()
// /            / tempD.cast<double>().colwise().sum();
// /    cout << "avgProd_preDorm = " << avgProd_preDorm.transpose() << endl;

    return tuple<ParaEst,ParaVec,EquStateV,EquStateV0,EquState0,SimData> (para_est,
        para_vec,EquV,EquV0,Equ0,sim_data);
}
// //
// //
// ///**********************************************************************************************
// //* Simulation of firm distribution for Status quo / counterfactual
// //**********************************************************************************************/
// //SimData alias::Simulation_StatusQuo_Counterfactual(const ParaEst & para_est, const ParaVec & para_vec,
// //    const EquStateV & EquV, const EquStateVmat & Evalmat, const EquStateV0 & EquV0, const EquStateV0mat & Eval0mat,
// //    const SimVar & sim_var, const int & SimN, const int & TBar, const int & GoodState, const int & LaborIntensive,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    double alpha_M;
// //    if (GoodState == 1 & LaborIntensive == 1){alpha_M = para.alpha_M_G_L;}
// //    else if (GoodState == 0 & LaborIntensive == 1){alpha_M = para.alpha_M_B_L;}
// //    else if (GoodState == 1 & LaborIntensive == 0){alpha_M = para.alpha_M_G_K;}
// //    else {alpha_M = para.alpha_M_B_K;}
// //
// //    SimData sim_data = InitializeSimData(SimN,TBar);
// //    sim_data.State_E = ArrayXXi::Zero(SimN,TBar);
// //    sim_data.phi.col(0) = exp(sim_var.lnphi_sim.col(0) * para_vec.prod_sigma0 + para_vec.prod_mu0);
// //    sim_data.CapitalMiss = ArrayXXi::Zero(SimN,TBar);
// //    sim_data.Employ_urMiss = ArrayXXi::Zero(SimN,TBar);
// //    sim_data.Employ_ucMiss = ArrayXXi::Zero(SimN,TBar);
// //
// //    auto worker = [&](size_t n, unsigned thread_id) {
// //        //    size_t thread_id = 0;
// //        //    for (size_t n = 0; n < SimN; ++n) {
// //        sim_data.first_year_in_Data(n) = 0;
// //        for (size_t t = 0; t <= TBar - 1; ++t) {
// //            //            cout << "n = " << n << "; t = " << t << endl;
// //            if (t == 0) {
// //                sim_data.phi(n, t) = exp(sim_var.lnphi_sim(n, t) * para_vec.prod_sigma0 + para_vec.prod_mu0);
// //                PhiW phi_w = definePhiWeights(para_vec.vec_phi, sim_data.phi(n, t));
// //
// //                LKState K_state = defineLKState(para_vec.vec_K, 0.0);
// //                LKState Lur_state = defineLKState(para_vec.vec_Lur, 0.0);
// //                LKState Luc_state = defineLKState(para_vec.vec_Luc, 0.0);
// //                tuple<double, double, double, double, double, double> t_LK = calOptKLurLuc(para_est,
// //                    para_vec, sim_var, phi_w,
// //                    K_state, Lur_state,
// //                    Luc_state, n, t,
// //                    Eval0mat.EVal_PD_K0_mat,
// //                    Eval0mat.EVal_PD_Lur0_mat,
// //                    Eval0mat.EVal_PD_Luc0_mat, 0,
// //                    para_vec.p_K, 1);
// //                double Capital = get<0>(t_LK);
// //                double Employ_ur = get<1>(t_LK);
// //                double Employ_uc = get<2>(t_LK);
// //                double K_P = get<3>(t_LK);
// //                double Lur_P = get<4>(t_LK);
// //                double Luc_P = get<5>(t_LK);
// //
// //                sim_data.K_sol(n, t) = K_P;
// //                sim_data.Lur_sol(n, t) = Lur_P;
// //                sim_data.Luc_sol(n, t) = Luc_P;
// //                sim_data.Capital(n, t) = Capital;
// //                sim_data.Employ_ur(n, t) = Employ_ur;
// //                sim_data.Employ_uc(n, t) = Employ_uc;
// //
// //                sim_data.Revenue(n, t) = RevOpt1_sim_Realized(para_est, sim_data.phi(n, t),
// //                    sim_data.Capital(n, t), sim_data.Employ_ur(n, t),
// //                    sim_data.Employ_uc(n, t));
// //
// //                sim_data.Price(n, t) = calPrice(para_est, alpha_M, sim_data.phi(n, t),
// //                    sim_data.Capital(n, t), sim_data.Employ_ur(n, t),
// //                    sim_data.Employ_uc(n, t), para_est.PI);
// //
// //                //// choose Stay or Exit State
// //                ArrayXd EVal_K_Interp = cal_EVal_K_phi_interpolation(Eval0mat.EVal_PD0_mat, phi_w,
// //                                                                     Lur_state, Luc_state, para.N_K);
// //
// //                double beta = (EVal_K_Interp(K_state.i_val_state + 1) - EVal_K_Interp(K_state.i_val_state))
// //                              / (para_vec.vec_K(K_state.i_val_state + 1) - para_vec.vec_K(K_state.i_val_state));
// //                double EVal_S = beta * (K_state.val_state - para_vec.vec_K(K_state.i_val_state))
// //                                + EVal_K_Interp(K_state.i_val_state);
// //                double ResV_E = calResidual_Value_Exit_double(para_est, para_vec.p_K, K_state.val_state,
// //                                                              Lur_state.val_state, Luc_state.val_state);
// //
// //                double dEVal_SE_sigma = EVal_S / para_est.sigma_SE;
// //                double Prob_S;
// //                double Prob_E;
// //                double lnProb_S;
// //                double lnProb_E;
// //
// //                double dEVal_SE_sigma_temp = max(dEVal_SE_sigma, -700.0);
// //                dEVal_SE_sigma_temp = min(dEVal_SE_sigma_temp, 700.0);
// //                Prob_S = exp(dEVal_SE_sigma_temp - log(exp(dEVal_SE_sigma_temp) + 1.0));
// //                Prob_E = exp(-dEVal_SE_sigma_temp - log(exp(-dEVal_SE_sigma_temp) + 1.0));
// //
// //                sim_data.Prob_S(n, t) = Prob_S;
// //                sim_data.Prob_E(n, t) = Prob_E;
// //                if (sim_var.SEsim_Empirical(n, t) < Prob_S) {
// //                    sim_data.State_S(n, t) = 1;
// //                    sim_data.State_P(n, t) = 1;
// //                    sim_data.State_D(n, t) = 0;
// //                } else {
// //                    sim_data.State_E(n, t) = 1;
// //                    sim_data.State_P(n, t) = 0;
// //                    sim_data.State_D(n, t) = 0;
// //                    sim_data.last_year_in_Data(n) = t;
// //                    break;
// //                }
// //
// //                sim_data.phi(n, t + 1) = exp(para_est.gamma0 + para_est.gamma1 * log(sim_data.phi(n, t))
// //                                             + sim_var.lnphi_sim(n, t + 1) * para_est.sigma_phi_eps);
// //            }
// //            else {
// //                //// for t > 0
// //                PhiW phi_w = definePhiWeights(para_vec.vec_phi, sim_data.phi(n, t));
// //
// //                LKState K_state = defineLKState(para_vec.vec_K, sim_data.Capital(n, t - 1));
// //                LKState Lur_state = defineLKState(para_vec.vec_Lur, sim_data.Employ_ur(n, t - 1));
// //                LKState Luc_state = defineLKState(para_vec.vec_Luc, sim_data.Employ_uc(n, t - 1));
// //
// //                int PD_State;
// //                if (sim_data.State_P(n, t-1) == 1) {PD_State = 0;}
// //                else {PD_State = 1;}
// //
// //                tuple<double, double, double, double> t_SE = calStayExit_Prob_logit(para_est, para_vec, Evalmat, phi_w,
// //                                                                                    K_state, Lur_state, Luc_state, PD_State);
// //                double Prob_S = get<0>(t_SE);
// //                double Prob_E = get<1>(t_SE);
// //
// //                sim_data.Prob_S(n, t) = Prob_S;
// //                sim_data.Prob_E(n, t) = Prob_E;
// //                if (sim_var.SEsim_Empirical(n, t) < Prob_S) {
// //                    sim_data.State_S(n, t) = 1;
// //                } else {
// //                    sim_data.State_E(n, t) = 1;
// //                    sim_data.last_year_in_Data(n) = t;
// //                    break;
// //                }
// //                tuple<double, double, double, double, double, double> t_LK = calOptKLurLuc(para_est,
// //                    para_vec, sim_var, phi_w, K_state, Lur_state, Luc_state, n, t,
// //                    Evalmat.EVal_PD_K_mat, Evalmat.EVal_PD_Lur_mat,
// //                    Evalmat.EVal_PD_Luc_mat,PD_State,para_vec.p_K,0);
// //                double Capital = get<0>(t_LK);
// //                double Employ_ur = get<1>(t_LK);
// //                double Employ_uc = get<2>(t_LK);
// //                double K_P = get<3>(t_LK);
// //                double Lur_P = get<4>(t_LK);
// //                double Luc_P = get<5>(t_LK);
// //
// //                sim_data.K_sol(n, t) = K_P;
// //                sim_data.Lur_sol(n, t) = Lur_P;
// //                sim_data.Luc_sol(n, t) = Luc_P;
// //                sim_data.Capital(n, t) = Capital;
// //                sim_data.Employ_ur(n, t) = Employ_ur;
// //                sim_data.Employ_uc(n, t) = Employ_uc;
// //
// //                tuple<double, double, double, double, double, double> t_PD = calProductionDormancy_Prob_lognormal(para_est,
// //                    para_vec,Evalmat,PD_State,sim_data.phi(n, t),
// //                    Capital,Employ_ur,Employ_uc);
// //                double Prob_P = get<0>(t_PD);
// //                double Prob_D = get<1>(t_PD);
// //                sim_data.CutoffVal(n,t) = get<5>(t_PD);
// //
// //                sim_data.Prob_P(n, t) = Prob_P;
// //                sim_data.Prob_D(n, t) = Prob_D;
// //                if (sim_var.PDsim_Empirical(n, t) < Prob_P) {
// //                    sim_data.State_P(n, t) = 1;
// //                    sim_data.State_D(n, t) = 0;
// //                } else {
// //                    sim_data.State_P(n, t) = 0;
// //                    sim_data.State_D(n, t) = 1;
// //                }
// //
// //                sim_data.Revenue(n, t) = RevOpt1_sim_Realized(para_est, sim_data.phi(n, t),
// //                    sim_data.Capital(n, t), sim_data.Employ_ur(n, t),
// //                    sim_data.Employ_uc(n, t));
// //
// //                sim_data.Price(n, t) = calPrice(para_est, alpha_M, sim_data.phi(n, t),
// //                                                sim_data.Capital(n, t), sim_data.Employ_ur(n, t),
// //                                                sim_data.Employ_uc(n, t), para_est.PI);
// //
// //                sim_data.phi(n, t + 1) = exp(para_est.gamma0 + para_est.gamma1 * log(sim_data.phi(n, t))
// //                                             + sim_var.lnphi_sim(n, t + 1) * para_est.sigma_phi_eps);
// //
// ////                if (n==15497) {
// ////                    cout << "PD_State = " << PD_State << "; sim_data.Revenue(n, t) = " << sim_data.Revenue(n, t)
// ////                        << "; sim_data.Revenue(n, t-1) = " << sim_data.Revenue(n, t-1)
// ////                            << "; sim_data.Capital(n, t) = " << sim_data.Capital(n, t)
// ////                            << "; sim_data.Capital(n, t-1) = " << sim_data.Capital(n, t-1)
// ////                            << "; sim_data.Employ_ur(n, t) = " << sim_data.Employ_ur(n, t)
// ////                            << "; sim_data.Employ_ur(n, t-1) = " << sim_data.Employ_ur(n, t-1)
// ////                            << "; sim_data.Employ_uc(n, t) = " << sim_data.Employ_uc(n, t)
// ////                            << "; sim_data.Employ_uc(n, t-1) = " << sim_data.Employ_uc(n, t-1)
// ////                            << "; Prob_S = " << Prob_S << "; sim_var.PDsim_Empirical(n, t) = " << sim_var.PDsim_Empirical(n, t)
// ////                            << "; Prob_P = " << Prob_P << "; sim_var.PDsim_Empirical(n, t) = " << sim_var.PDsim_Empirical(n, t)
// ////                            << "; sim_data.phi(n, t) = " << sim_data.phi(n, t)
// ////                            << endl;
// ////                }
// //            }
// //        }
// //        //    }
// //    };
// //    MultiThreads::simple_parallel_for(worker, SimN, threadsManagement);
// ////
// ////  for exploring production cost
// //    writeToCSVfile("phi_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.phi.cast<double>().matrix());
// //    writeToCSVfile("Revenue_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.Revenue.cast<double>().matrix());
// //    writeToCSVfile("Prob_P_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.Prob_P.cast<double>().matrix());
// //    writeToCSVfile("Prob_D_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.Prob_D.cast<double>().matrix());
// //    writeToCSVfile("PDsim_Empirical_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_var.PDsim_Empirical.cast<double>().matrix());
// //    writeToCSVfile("SEsim_Empirical" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //               sim_var.SEsim_Empirical.cast<double>().matrix());
// //    writeToCSVfile("State_P_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.State_P.cast<double>().matrix());
// //    writeToCSVfile("State_D_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.State_D.cast<double>().matrix());
// //    writeToCSVfile("State_S_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.State_S.cast<double>().matrix());
// //    writeToCSVfile("State_E_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //               sim_data.State_E.cast<double>().matrix());
// //    writeToCSVfile("Employ_uc_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.Employ_uc.cast<double>().matrix());
// //    writeToCSVfile("Employ_ur_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.Employ_ur.cast<double>().matrix());
// //    writeToCSVfile("Capital_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                   sim_data.Capital.cast<double>().matrix());
// //
// //
// //   writeToCSVfile("CutoffVal_State" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + ".csv",
// //                  sim_data.CutoffVal.cast<double>().matrix());
// //////    throw runtime_error("417");
// //    return sim_data;
// //}
// //
// ///**********************************************************************************************
// //* Calculate Aggregates for counterfactuals
// //**********************************************************************************************/
// /////*** Calculate Economic Variables of interests ***/
// //CounterFactVariable alias::CalEconomicVariable(const double & p_K, const double & FirmMass, const double & F_Entry,
// //    const double & CapitalDemand, SimData & sim_data) {
// //
// //    CounterFactVariable EconVal;
// //
// //    EconVal.FirmMass = FirmMass;
// //    EconVal.F_Entry = F_Entry;
// //    EconVal.p_K = p_K;
// //    EconVal.CapitalDemand = CapitalDemand;
// //    EconVal.NumFirm = EconVal.FirmMass * sim_data.State_S.cast<double>().sum();
// //
// //    /*** Calculate Average Productivity ***/
// //    EconVal.AvgProd = calAvgProductivity(sim_data, EconVal.FirmMass);
// //    EconVal.WAvgProd = calWAvgProductivity(sim_data, EconVal.FirmMass);
// //    EconVal.EntrantProd = calEntrantProductivity(sim_data, EconVal.FirmMass);
// //    EconVal.ExitProd = calExitProductivity(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Total Labor Employment ***/
// //    tuple<double,double> t_TotEmp = calTotalEmployment(sim_data, EconVal.FirmMass);
// //    EconVal.TotEmploy_ur = get<0>(t_TotEmp);
// //    EconVal.TotEmploy_uc = get<1>(t_TotEmp);
// //    EconVal.TotEmploy = EconVal.TotEmploy_ur + EconVal.TotEmploy_uc;
// //
// //    /*** Calculate Total Output ***/
// //    EconVal.TotOutput = calTotalOutput(sim_data, EconVal.FirmMass);
// //    EconVal.Output_per_Capital = EconVal.TotOutput / EconVal.CapitalDemand;
// //    EconVal.Output_per_Employ = EconVal.TotOutput / (EconVal.TotEmploy_ur + EconVal.TotEmploy_uc);
// //
// //    /*** Calculate Average length of dormancy ***/
// //    EconVal.DormLength = calLengthDormancy(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Average length of production ***/
// //    EconVal.ProdLength = calLengthProduction(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Average length of Survival ***/
// //    EconVal.SurvivalLength = calLengthSurvival(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Price Index ***/
// //    EconVal.PriceIndex = calPriceIndex(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Average Exit Rate ***/
// //    EconVal.AvgExitRate = calAvgExitRate(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Capital Labor Ratio ***/
// //    EconVal.K_L_ratio = calCapitalLaborRatio(sim_data, EconVal.FirmMass);
// //
// //    /*** Calculate Welfare ***/
// //    EconVal.Welfare_Manu = EconVal.TotOutput / EconVal.PriceIndex;
// //
// //    return EconVal;
// //}
// //
// ///*** Calculate Average Productivity ***/
// //double alias::calAvgProductivity(SimData & sim_data_C, const double FirmMass) {
// //    int Ncol = sim_data_C.State_S.cols();
// //    double AvgProd = ( sim_data_C.phi(all,seqN(0,Ncol)) * sim_data_C.State_S.cast<double>() ).sum()
// //                     / sim_data_C.State_S.cast<double>().sum();
// //    return AvgProd;
// //}
// //
// //double alias::calWAvgProductivity(SimData & sim_data_C, const double FirmMass) {
// //    int Ncol = sim_data_C.State_S.cols();
// //    double AvgProd = ( sim_data_C.phi(all,seqN(0,Ncol)) * sim_data_C.Revenue ).sum()
// //                     / sim_data_C.Revenue.sum();
// //    return AvgProd;
// //}
// //
// //double alias::calEntrantProductivity(SimData & sim_data_C, const double FirmMass) {
// //    double EntrantProd = ( sim_data_C.phi(all,0) * sim_data_C.State_S(all,0).cast<double>() ).sum()
// //                         / sim_data_C.State_S(all,0).cast<double>().sum();
// //    return EntrantProd;
// //}
// //
// //double alias::calExitProductivity(SimData & sim_data_C, const double FirmMass) {
// //    int Ncol = sim_data_C.State_S.cols();
// //    double ExitProd = ( sim_data_C.phi(all,seqN(0,Ncol)) * sim_data_C.State_E.cast<double>() ).sum()
// //                      / sim_data_C.State_E.cast<double>().sum();
// //    return ExitProd;
// //}
// //
// ///*** Calculate Total Labor Employment ***/
// //tuple<double,double> alias::calTotalEmployment(SimData & sim_data_C, const double FirmMass) {
// //    double TotEmp_ur = ( sim_data_C.Employ_ur * sim_data_C.State_S.cast<double>() ).sum() * FirmMass;
// //    double TotEmp_uc = ( sim_data_C.Employ_uc * sim_data_C.State_S.cast<double>() ).sum() * FirmMass;
// //    cout << "TotEmp_ur = " << ( sim_data_C.Employ_ur * sim_data_C.State_S.cast<double>() ).sum()
// //        << "; TotEmp_uc = " << ( sim_data_C.Employ_uc * sim_data_C.State_S.cast<double>() ).sum()
// //        << "; FirmMass = " << FirmMass << endl;
// //
// //    return tuple<double,double>(TotEmp_ur,TotEmp_uc);
// //}
// //
// ///*** Calculate Total Output ***/
// //double alias::calTotalOutput(SimData & sim_data_C, const double FirmMass) {
// //    double TotOutput = ( sim_data_C.Revenue * sim_data_C.State_P.cast<double>() ).sum() * FirmMass;
// //    return TotOutput;
// //}
// //
// ///*** Calculate Average length of dormancy ***/
// //double alias::calLengthDormancy(SimData & sim_data_C, const double FirmMass) {
// //    ArrayXi DormancyYr = sim_data_C.State_D.rowwise().sum();
// //    double LengthDorm = DormancyYr.cast<double>().mean();
// //    return LengthDorm;
// //}
// //
// ///*** Calculate Average length of production ***/
// //double alias::calLengthProduction(SimData & sim_data_C, const double FirmMass) {
// //    ArrayXi ProductionYr = sim_data_C.State_P.rowwise().sum();
// //    double LengthProd= ProductionYr.cast<double>().mean();
// //    return LengthProd;
// //}
// //
// ///*** Calculate Average length of Survival ***/
// //double alias::calLengthSurvival(SimData & sim_data_C, const double FirmMass) {
// //    ArrayXi SurvivalYr = sim_data_C.State_S.rowwise().sum();
// //    double LengthSurv = SurvivalYr.cast<double>().mean();
// //    return LengthSurv;
// //}
// //
// ///*** Calculate Price Index ***/
// //double alias::calPriceIndex(SimData & sim_data_C, const double FirmMass) {
// //    ArrayXXd Price_temp = ( sim_data_C.Price * sim_data_C.State_P.cast<double>() ).pow(1.0-para.sigma);
// //    Price_temp = (Price_temp.isFinite()).select(Price_temp,0.0);
// //    double PriceIndex = pow(Price_temp.sum()*FirmMass, 1.0/(1.0-para.sigma));
// //    return PriceIndex;
// //}
// //
// ///*** Calculate Average Exit Rate ***/
// //double alias::calAvgExitRate(SimData & sim_data_C, const double FirmMass) {
// //
// //    int SimTbar = sim_data_C.phi.cols()-1;
// //    int SimN = sim_data_C.phi.rows();
// //
// //    double AvgExitRate = sim_data_C.State_E(all,seqN(1,SimTbar-1)).cast<double>().sum()
// //                         / sim_data_C.State_S(all,seqN(0,SimTbar-1)).cast<double>().sum();
// //
// //    return AvgExitRate;
// //}
// //
// ///*** Calculate Capital Labor Ratio ***/
// //double alias::calCapitalLaborRatio(SimData & sim_data_C, const double FirmMass) {
// //
// //    int SimTbar = sim_data_C.phi.cols()-1;
// //    int SimN = sim_data_C.phi.rows();
// //
// //    ArrayXXd Employ = sim_data_C.Employ_ur(all,seqN(1,SimTbar-1)).cast<double>()
// //        + sim_data_C.Employ_uc(all,seqN(1,SimTbar-1)).cast<double>()
// //        + (1.0 - sim_data_C.State_S(all,seqN(1,SimTbar-1)).cast<double>());
// //
// //    ArrayXXd K_L_ratio_mat = sim_data_C.Capital(all,seqN(1,SimTbar-1)).cast<double>() / Employ;
// //
// //    double K_L_ratio = ( K_L_ratio_mat * sim_data_C.State_S(all,seqN(1,SimTbar-1)).cast<double>() ).sum()
// //        / sim_data_C.State_S(all,seqN(1,SimTbar-1)).cast<double>().sum();
// //    return K_L_ratio;
// //}
// //
// /////******** Auxiliary function for counterfactual ********/
// //CFDeltaVec alias::InitializeCFDeltaVec(const int & N_C) {
// //
// //    CFDeltaVec CounterFactDeltaVec;
// //
// //    CounterFactDeltaVec.FirmMass_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.NumFirm_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.CapitalDemand_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.p_K_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.AvgProd_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.WAvgProd_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.EntrantProd_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.ExitProd_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.TotEmploy_ur_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.TotEmploy_uc_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.TotEmploy_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.TotOutput_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.Output_per_Capital_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.Output_per_Employ_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.DormLength = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.ProdLength = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.SurvivalLength = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.PriceIndex_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.AvgExitRate = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.K_L_ratio_Delta = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.AggValueAdded = ArrayXd::Zero(N_C);
// //    CounterFactDeltaVec.Welfare_Manu_Delta = ArrayXd::Zero(N_C);
// //
// //    return CounterFactDeltaVec;
// //}
// //
// //ArrayXd alias::Assign_thetaC(const ArrayXd & theta_Est_S, const ArrayXd & ResidualValue_row, const ArrayXd & firing_share_row,
// //    const int & GoodState, const int & LaborIntensive) {
// //
// //    int n_Est;
// //    if (GoodState == 1) {n_Est = 14;}
// //    else {n_Est = 15;}
// //
// //    int n_Est_high_ur; int n_Est_low_ur; int n_Est_uc; int n_Est_Exit_ur; int n_Est_Exit_uc;
// //    if (GoodState == 1) {
// //        n_Est_high_ur = 4; n_Est_low_ur = 5; n_Est_uc = 9;
// //    }
// //    else {
// //        n_Est_high_ur = 6; n_Est_low_ur = 7; n_Est_uc = 10;
// //    }
// //
// //    int k_ResVFiringShr;
// //    if (GoodState == 1 & LaborIntensive == 1) {k_ResVFiringShr = 0;}
// //    else if (GoodState == 0 & LaborIntensive == 1) {k_ResVFiringShr = 1;}
// //    else if (GoodState == 1 & LaborIntensive == 0) {k_ResVFiringShr = 2;}
// //    else {k_ResVFiringShr = 3;}
// //
// //    ArrayXd theta_C_ResFC = theta_Est_S;
// //    theta_C_ResFC(n_Est) = ResidualValue_row(k_ResVFiringShr) + theta_Est_S(n_Est);
// //
// //    theta_C_ResFC(n_Est_high_ur) = firing_share_row(k_ResVFiringShr) * theta_Est_S(n_Est_high_ur);
// //    theta_C_ResFC(n_Est_low_ur) = firing_share_row(k_ResVFiringShr) * theta_Est_S(n_Est_low_ur);
// //    theta_C_ResFC(n_Est_uc) = firing_share_row(k_ResVFiringShr) * theta_Est_S(n_Est_uc);
// //
// //    return theta_C_ResFC;
// //}
// //
// //CounterFactVariable alias::CalEconomicVariableDelta(const CounterFactVariable & Value_C,
// //    const CounterFactVariable & Value_StatusQuo) {
// //
// //    CounterFactVariable Value_Delta;
// //    Value_Delta.p_K_Delta = Value_C.p_K / Value_StatusQuo.p_K;
// //    Value_Delta.FirmMass_Delta = Value_C.FirmMass / Value_StatusQuo.FirmMass;
// //    Value_Delta.CapitalDemand_Delta = Value_C.CapitalDemand / Value_StatusQuo.CapitalDemand;
// //    Value_Delta.NumFirm_Delta = Value_C.NumFirm / Value_StatusQuo.NumFirm;
// //
// //    Value_Delta.AvgProd_Delta = Value_C.AvgProd / Value_StatusQuo.AvgProd;
// //    Value_Delta.WAvgProd_Delta = Value_C.WAvgProd / Value_StatusQuo.WAvgProd;
// //    Value_Delta.EntrantProd_Delta = Value_C.EntrantProd / Value_StatusQuo.EntrantProd;
// //    Value_Delta.ExitProd_Delta = Value_C.ExitProd / Value_StatusQuo.ExitProd;
// //
// //    Value_Delta.TotEmploy_ur_Delta = Value_C.TotEmploy_ur / Value_StatusQuo.TotEmploy_ur;
// //    Value_Delta.TotEmploy_uc_Delta = Value_C.TotEmploy_uc / Value_StatusQuo.TotEmploy_uc;
// //    Value_Delta.TotEmploy_Delta = Value_C.TotEmploy / Value_StatusQuo.TotEmploy;
// //
// //    Value_Delta.TotOutput_Delta = Value_C.TotOutput / Value_StatusQuo.TotOutput;
// //    Value_Delta.Output_per_Capital_Delta = Value_C.Output_per_Capital / Value_StatusQuo.Output_per_Capital;
// //    Value_Delta.Output_per_Employ_Delta = Value_C.Output_per_Employ / Value_StatusQuo.Output_per_Employ;
// //
// //    Value_Delta.DormLength = Value_C.DormLength;
// //    Value_Delta.ProdLength = Value_C.ProdLength;
// //    Value_Delta.SurvivalLength = Value_C.SurvivalLength;
// //
// //    Value_Delta.PriceIndex_Delta = Value_C.PriceIndex / Value_StatusQuo.PriceIndex;
// //
// //    Value_Delta.AvgExitRate = Value_C.AvgExitRate;
// //    Value_Delta.K_L_ratio_Delta = Value_C.K_L_ratio / Value_StatusQuo.K_L_ratio;
// //
// //    Value_Delta.AggValueAdded = Value_C.TotOutput;
// //
// //    Value_Delta.Welfare_Manu_Delta = 1.0 / (Value_Delta.PriceIndex_Delta);
// //    return Value_Delta;
// //}
// //
// //int alias::writeToCSVfileValueDelta(const string & filename, const CFDeltaVec & CFDeltaVec_Status) {
// //
// //    writeToCSVfile("FirmMass_Delta_" + filename + ".csv", CFDeltaVec_Status.FirmMass_Delta.cast<double>().matrix());
// //    writeToCSVfile("NumFirm_Delta_" + filename + ".csv", CFDeltaVec_Status.NumFirm_Delta.cast<double>().matrix());
// //    writeToCSVfile("CapitalDemand_Delta_" + filename + ".csv", CFDeltaVec_Status.CapitalDemand_Delta.cast<double>().matrix());
// //    writeToCSVfile("p_K_Delta_" + filename + ".csv", CFDeltaVec_Status.p_K_Delta.cast<double>().matrix());
// //    writeToCSVfile("AvgProd_Delta_" + filename + ".csv", CFDeltaVec_Status.AvgProd_Delta.cast<double>().matrix());
// //    writeToCSVfile("WAvgProd_Delta_" + filename + ".csv", CFDeltaVec_Status.WAvgProd_Delta.cast<double>().matrix());
// //    writeToCSVfile("EntrantProd_Delta_" + filename + ".csv", CFDeltaVec_Status.EntrantProd_Delta.cast<double>().matrix());
// //    writeToCSVfile("ExitProd_Delta_" + filename + ".csv", CFDeltaVec_Status.ExitProd_Delta.cast<double>().matrix());
// //    writeToCSVfile("TotEmploy_ur_Delta_" + filename + ".csv", CFDeltaVec_Status.TotEmploy_ur_Delta.cast<double>().matrix());
// //    writeToCSVfile("TotEmploy_uc_Delta_" + filename + ".csv", CFDeltaVec_Status.TotEmploy_uc_Delta.cast<double>().matrix());
// //    writeToCSVfile("TotEmploy_Delta_" + filename + ".csv", CFDeltaVec_Status.TotEmploy_Delta.cast<double>().matrix());
// //    writeToCSVfile("TotOutput_Delta_" + filename + ".csv", CFDeltaVec_Status.TotOutput_Delta.cast<double>().matrix());
// //    writeToCSVfile("Output_per_Capital_Delta_" + filename + ".csv", CFDeltaVec_Status.Output_per_Capital_Delta.cast<double>().matrix());
// //    writeToCSVfile("Output_per_Employ_Delta_" + filename + ".csv", CFDeltaVec_Status.Output_per_Employ_Delta.cast<double>().matrix());
// //    writeToCSVfile("DormLength_Delta_" + filename + ".csv", CFDeltaVec_Status.DormLength.cast<double>().matrix());
// //    writeToCSVfile("ProdLength_Delta_" + filename + ".csv", CFDeltaVec_Status.ProdLength.cast<double>().matrix());
// //    writeToCSVfile("SurvivalLength_Delta_" + filename + ".csv", CFDeltaVec_Status.SurvivalLength.cast<double>().matrix());
// //    writeToCSVfile("PriceIndex_Delta_" + filename + ".csv", CFDeltaVec_Status.PriceIndex_Delta.cast<double>().matrix());
// //    writeToCSVfile("AvgExitRate_Delta_" + filename + ".csv", CFDeltaVec_Status.AvgExitRate.cast<double>().matrix());
// //    writeToCSVfile("K_L_ratio_Delta_" + filename + ".csv", CFDeltaVec_Status.K_L_ratio_Delta.cast<double>().matrix());
// //    writeToCSVfile("firing_share_" + filename + ".csv", CFDeltaVec_Status.firing_share.cast<double>().matrix());
// //    writeToCSVfile("ResidualValue_" + filename + ".csv", CFDeltaVec_Status.ResidualValue.cast<double>().matrix());
// //    writeToCSVfile("EntryCost_" + filename + ".csv", CFDeltaVec_Status.EntryCost.cast<double>().matrix());
// //    writeToCSVfile("ExitSubsidy_" + filename + ".csv", CFDeltaVec_Status.ExitSubsidy.cast<double>().matrix());
// //    writeToCSVfile("AggValueAdded_" + filename + ".csv", CFDeltaVec_Status.AggValueAdded.cast<double>().matrix());
// //    writeToCSVfile("Welfare_Manu_Delta_" + filename + ".csv", CFDeltaVec_Status.Welfare_Manu_Delta.cast<double>().matrix());
// //
// //    return 0;
// //}
// //
// ////
// //
// ///**************************************************************
// //* **************************************************************
// //* Solve the Partial Counterfactural equilibrium
// //* **************************************************************
// ///**************************************************************
// //*** Counterfactual Exercise: Partial Equilibrium; Counterfactual 1.1 (Table): Varying the exit barrier ***
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_1_v1(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    const string & FKL, const ArrayXd & TargetExitRate_vec, MultiThreads::Threads_Management & threadsManagement) {
// //
// //    int N_C = TargetExitRate_vec.size();
// //    ArrayXd ResidualValue_GoodLaborInt(N_C);
// //    ArrayXd ResidualValue_BadLaborInt(N_C);
// //    ArrayXd ResidualValue_GoodCapitalInt(N_C);
// //    ArrayXd ResidualValue_BadCapitalInt(N_C);
// //    //
// //    for (size_t k = 0; k < N_C; ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        ResidualValue_GoodLaborInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 1;
// //        ResidualValue_BadLaborInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //        GoodState = 1; LaborIntensive = 0;
// //        ResidualValue_GoodCapitalInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 0;
// //        ResidualValue_BadCapitalInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //    }
// //    // cout << "ResidualValue_GoodLaborInt(k) = " << ResidualValue_GoodLaborInt(4)
// //    //     << "; ResidualValue_BadLaborInt(k) = " << ResidualValue_BadLaborInt(4)
// //    //     << "; ResidualValue_GoodCapitalInt(k) = " << ResidualValue_GoodCapitalInt(4)
// //    //     << "; ResidualValue_BadCapitalInt(k) = " << ResidualValue_BadCapitalInt(4) << endl;
// //    // throw runtime_error("1149");
// //
// //    ArrayXXd ResidualValue(TargetExitRate_vec.rows(),4);
// //    ResidualValue.col(0) = ResidualValue_GoodLaborInt; //ResidualValue(4,0) = 292.03;
// //    ResidualValue.col(1) = ResidualValue_BadLaborInt; //ResidualValue(4,1) = 457.741;
// //    ResidualValue.col(2) = ResidualValue_GoodCapitalInt; //ResidualValue(4,2) = 4.73785;
// //    ResidualValue.col(3) = ResidualValue_BadCapitalInt; //ResidualValue(4,3) = 168.739;
// //
// //    ArrayXXd VecOnes = ArrayXXd::Ones(TargetExitRate_vec.rows(),4);
// //
// //    ArrayXXd EntryCost(TargetExitRate_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //
// //    tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,VecOnes, EntryCost, FKL, threadsManagement);
// //
// //    return 0;
// //}
// //
// ///**************************************************************
// //*** Counterfactual Exercise: Partial Equilibrium; Counterfactual 1.2 (Graph): Varying the exit barrier ***
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_1_v2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
// //    const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas, const string &FKL,
// //    const ArrayXd &ResVal_vec, MultiThreads::Threads_Management &threadsManagement) {
// //
// //    int N_C = ResVal_vec.size();
// //
// //    ArrayXXd ResidualValue(N_C,4);
// //    ResidualValue.col(0) << ResVal_vec;
// //    ResidualValue.col(1) << ResVal_vec;
// //    ResidualValue.col(2) << ResVal_vec;
// //    ResidualValue.col(3) << ResVal_vec;
// //    cout << "ResidualValue = " << ResidualValue << endl;
// //
// //    ArrayXXd VecOnes = ArrayXXd::Ones(ResidualValue.rows(),4);
// //
// //    ArrayXXd EntryCost(ResidualValue.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(ResidualValue.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(ResidualValue.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(ResidualValue.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(ResidualValue.rows());
// //
// //    tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,VecOnes, EntryCost, FKL, threadsManagement);
// //
// //    return 0;
// //}
// //
// ///**************************************************************
// //*** Counterfactual Exercise: Partial Equilibrium; Counterfactual 1.2 (Table): Varying the firing costs ***
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_2_v1(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    const string & FKL, const ArrayXd & TargetExitRate_vec, MultiThreads::Threads_Management & threadsManagement) {
// //
// //    int N_C = TargetExitRate_vec.size();
// //
// //    ArrayXd firingshare_GoodLaborInt(N_C);
// //    ArrayXd firingshare_BadLaborInt(N_C);
// //    ArrayXd firingshare_GoodCapitalInt(N_C);
// //    ArrayXd firingshare_BadCapitalInt(N_C);
// //    for (size_t k = 0; k < N_C; ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        firingshare_GoodLaborInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 1;
// //        firingshare_BadLaborInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //        GoodState = 1; LaborIntensive = 0;
// //        firingshare_GoodCapitalInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 0;
// //        firingshare_BadCapitalInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_vec(k),threadsManagement);
// //    }
// //
// //    ArrayXXd firingshare(N_C,4);
// //    firingshare.col(0) = firingshare_GoodLaborInt;
// //    firingshare.col(1) = firingshare_BadLaborInt;
// //    firingshare.col(2) = firingshare_GoodCapitalInt;
// //    firingshare.col(3) = firingshare_BadCapitalInt;
// //
// //    ArrayXXd EntryCost(TargetExitRate_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_vec.rows());
// //
// //    ArrayXXd ResidualValue = ArrayXXd::Zero(N_C,4);
// //
// //    tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,firingshare, EntryCost, FKL, threadsManagement);
// //
// //    return 0;
// //}
// //
// ///**************************************************************
// //*** Counterfactual Exercise: Partial Equilibrium; Counterfactual 1.2 (Graph): Varying the firing costs ***
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_2_v2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    const string & FKL, const ArrayXd & firingshare_vec, MultiThreads::Threads_Management & threadsManagement) {
// //
// //    int N_C = firingshare_vec.size();
// //
// //    ArrayXXd firingshare(N_C,4);
// //    firingshare.col(0) << firingshare_vec;
// //    firingshare.col(1) << firingshare_vec;
// //    firingshare.col(2) << firingshare_vec;
// //    firingshare.col(3) << firingshare_vec;
// //
// //    ArrayXXd ResidualValue = ArrayXXd::Zero(firingshare.rows(),4);
// //
// //    ArrayXXd EntryCost(firingshare.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(firingshare.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(firingshare.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(firingshare.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(firingshare.rows());
// //
// //    tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,firingshare, EntryCost, FKL, threadsManagement);
// //    return 0;
// //}
// //
// ///**************************************************************
// //* Counterfactual Exercise 3: Change KsupplyElas and Change the exit barrier to match the US exit rate
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_3_Weighted(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN,
// //    const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd ResidualValue_GoodLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd ResidualValue_BadLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd ResidualValue_GoodCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //    ArrayXd ResidualValue_BadCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //
// //    ArrayXXd EntryCost(TargetExitRate_LaborInt_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //
// //    for (size_t k = 0; k < TargetExitRate_LaborInt_vec.rows(); ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        ResidualValue_GoodLaborInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_LaborInt_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 1;
// //        ResidualValue_BadLaborInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_LaborInt_vec(k),threadsManagement);
// //        GoodState = 1; LaborIntensive = 0;
// //        ResidualValue_GoodCapitalInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_CapitalInt_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 0;
// //        ResidualValue_BadCapitalInt(k) = SolveResidualValueMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_CapitalInt_vec(k),threadsManagement);
// //    }
// //   // throw runtime_error("1149");
// //
// //    ArrayXXd ResidualValue(TargetExitRate_LaborInt_vec.rows(),4);
// //    ResidualValue.col(0) = ResidualValue_GoodLaborInt;
// //    ResidualValue.col(1) = ResidualValue_BadLaborInt;
// //    ResidualValue.col(2) = ResidualValue_GoodCapitalInt;
// //    ResidualValue.col(3) = ResidualValue_BadCapitalInt;
// //
// //    ArrayXXd firingshare = ArrayXXd::Ones(TargetExitRate_LaborInt_vec.rows(),4);
// //
// //    int N_KsupplyElas = 5;
// //    ArrayXd KsupplyElas(N_KsupplyElas);
// //    KsupplyElas << 0.0,0.1,0.2,0.5,0.75;
// //
// //    for(size_t k = 0; k < N_KsupplyElas; ++k) {
// //        string FKL = "_Partial_C3_ResVal_Weighted_KsupplyElas" + to_string(k);
// //        tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //            Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //            Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //            sim_var, SimTbar, SimN, KsupplyElas(k),
// //            ResidualValue,firingshare, EntryCost, FKL, threadsManagement);
// //    }
// ////    throw runtime_error("330");
// //    return 0;
// //}
// //
// //int alias::SolveSimulate_CounterfacturalResult_3_Aggregate(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN,
// //    const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd ResidualValue_GoodLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd ResidualValue_BadLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd ResidualValue_GoodCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //    ArrayXd ResidualValue_BadCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //
// //    ArrayXXd EntryCost(TargetExitRate_LaborInt_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //
// //    for (size_t k = 0; k < TargetExitRate_LaborInt_vec.rows(); ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        ResidualValue_GoodLaborInt(k) = 168.75;
// //        GoodState = 0; LaborIntensive = 1;
// //        ResidualValue_BadLaborInt(k) = 168.75;
// //        GoodState = 1; LaborIntensive = 0;
// //        ResidualValue_GoodCapitalInt(k) = 168.75;
// //        GoodState = 0; LaborIntensive = 0;
// //        ResidualValue_BadCapitalInt(k) = 168.75;
// //    }
// //   // throw runtime_error("1149");
// //
// //    ArrayXXd ResidualValue(TargetExitRate_LaborInt_vec.rows(),4);
// //    ResidualValue.col(0) = ResidualValue_GoodLaborInt;
// //    ResidualValue.col(1) = ResidualValue_BadLaborInt;
// //    ResidualValue.col(2) = ResidualValue_GoodCapitalInt;
// //    ResidualValue.col(3) = ResidualValue_BadCapitalInt;
// //
// //    ArrayXXd firingshare = ArrayXXd::Ones(TargetExitRate_LaborInt_vec.rows(),4);
// //
// //    int N_KsupplyElas = 5;
// //    ArrayXd KsupplyElas(N_KsupplyElas);
// //    KsupplyElas << 0.0,0.1,0.2,0.5,0.75;
// //
// //    for(size_t k = 0; k < N_KsupplyElas; ++k) {
// //        string FKL = "_Partial_C3_ResVal_Aggregate_KsupplyElas" + to_string(k);
// //        tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //            Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //            Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //            sim_var, SimTbar, SimN, KsupplyElas(k),
// //            ResidualValue,firingshare, EntryCost, FKL, threadsManagement);
// //    }
// ////    throw runtime_error("330");
// //    return 0;
// //}
// //
// ///**************************************************************
// //* Counterfactual Exercise 4: Change KsupplyElas and Change firing costs to match the US exit rate
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_4_Weighted(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN,
// //    const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd firingshare_GoodLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd firingshare_BadLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd firingshare_GoodCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //    ArrayXd firingshare_BadCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //
// //    ArrayXXd EntryCost(TargetExitRate_LaborInt_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //
// //    for (size_t k = 0; k < TargetExitRate_LaborInt_vec.rows(); ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        firingshare_GoodLaborInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_LaborInt_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 1;
// //        firingshare_BadLaborInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadLaborInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_LaborInt_vec(k),threadsManagement);
// //        GoodState = 1; LaborIntensive = 0;
// //        firingshare_GoodCapitalInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_GoodCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_CapitalInt_vec(k),threadsManagement);
// //        GoodState = 0; LaborIntensive = 0;
// //        firingshare_BadCapitalInt(k) = SolvefiringshareMatchingTargetExitRate(theta_Est_S,
// //            Value_StatusQuo_BadCapitalInt, sim_var, GoodState, LaborIntensive,
// //            TargetExitRate_CapitalInt_vec(k),threadsManagement);
// //    }
// //
// //    ArrayXXd firingshare(TargetExitRate_LaborInt_vec.rows(),4);
// //    firingshare.col(0) = firingshare_GoodLaborInt;
// //    firingshare.col(1) = firingshare_BadLaborInt;
// //    firingshare.col(2) = firingshare_GoodCapitalInt;
// //    firingshare.col(3) = firingshare_BadCapitalInt;
// //
// //    // ArrayXXd firingshare = 0.5*ArrayXXd::Ones(TargetExitRate_LaborInt_vec.rows(),4);
// //
// //    ArrayXXd ResidualValue = ArrayXXd::Zero(TargetExitRate_LaborInt_vec.rows(),4);
// //
// //    int N_KsupplyElas = 5;
// //    ArrayXd KsupplyElas(N_KsupplyElas);
// //    KsupplyElas << 0.0,0.1,0.2,0.5,0.75;
// //
// //    for(size_t k = 0; k < N_KsupplyElas; ++k) {
// //        string FKL = "_Partial_C4_firingcost_Weighted_KsupplyElas" + to_string(k);
// //        tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //            Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //            Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //            sim_var, SimTbar, SimN, KsupplyElas(k),
// //            ResidualValue,firingshare,EntryCost,FKL,threadsManagement);
// //    }
// //
// //    return 0;
// //}
// //
// //int alias::SolveSimulate_CounterfacturalResult_4_Aggregate(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN,
// //    const ArrayXd & TargetExitRate_LaborInt_vec, const ArrayXd & TargetExitRate_CapitalInt_vec,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd firingshare_GoodLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd firingshare_BadLaborInt(TargetExitRate_LaborInt_vec.rows());
// //    ArrayXd firingshare_GoodCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //    ArrayXd firingshare_BadCapitalInt(TargetExitRate_CapitalInt_vec.rows());
// //
// //    ArrayXXd EntryCost(TargetExitRate_LaborInt_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(TargetExitRate_LaborInt_vec.rows());
// //
// //    for (size_t k = 0; k < TargetExitRate_LaborInt_vec.rows(); ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        firingshare_GoodLaborInt(k) = 0.38125;
// //        GoodState = 0; LaborIntensive = 1;
// //        firingshare_BadLaborInt(k) = 0.38125;
// //        GoodState = 1; LaborIntensive = 0;
// //        firingshare_GoodCapitalInt(k) = 0.38125;
// //        GoodState = 0; LaborIntensive = 0;
// //        firingshare_BadCapitalInt(k) = 0.38125;
// //    }
// //
// //    ArrayXXd firingshare(TargetExitRate_LaborInt_vec.rows(),4);
// //    firingshare.col(0) = firingshare_GoodLaborInt;
// //    firingshare.col(1) = firingshare_BadLaborInt;
// //    firingshare.col(2) = firingshare_GoodCapitalInt;
// //    firingshare.col(3) = firingshare_BadCapitalInt;
// //
// //    // ArrayXXd firingshare = 0.5*ArrayXXd::Ones(TargetExitRate_LaborInt_vec.rows(),4);
// //
// //    ArrayXXd ResidualValue = ArrayXXd::Zero(TargetExitRate_LaborInt_vec.rows(),4);
// //
// //    int N_KsupplyElas = 5;
// //    ArrayXd KsupplyElas(N_KsupplyElas);
// //    KsupplyElas << 0.0,0.1,0.2,0.5,0.75;
// //
// //    for(size_t k = 0; k < N_KsupplyElas; ++k) {
// //        string FKL = "_Partial_C4_firingcost_Aggregate_KsupplyElas" + to_string(k);
// //        tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //            Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //            Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //            sim_var, SimTbar, SimN, KsupplyElas(k),
// //            ResidualValue,firingshare,EntryCost,FKL,threadsManagement);
// //    }
// //
// //    return 0;
// //}
// ///**************************************************************
// //* Counterfactual Exercise 5: Varying both the exit barrier and firing costs at the same time
// //**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_5(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    const double & Residual_up, const double & firingshare_low, const int & N_C,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd ResidualValue_temp = ArrayXd::LinSpaced(N_C, 0.0,Residual_up);
// //    ArrayXXd ResidualValue(N_C,4);
// //    ResidualValue.col(0) = ResidualValue_temp;
// //    ResidualValue.col(1) = ResidualValue_temp;
// //    ResidualValue.col(2) = ResidualValue_temp;
// //    ResidualValue.col(3) = ResidualValue_temp;
// //    cout << "ResidualValue = " << ResidualValue.transpose() << endl;
// //
// //    ArrayXd firing_share_temp = ArrayXd::LinSpaced(N_C, 1,firingshare_low);
// //    ArrayXXd firing_share(N_C,4);
// //    firing_share.col(0) = firing_share_temp;
// //    firing_share.col(1) = firing_share_temp;
// //    firing_share.col(2) = firing_share_temp;
// //    firing_share.col(3) = firing_share_temp;
// //    cout << "firing_share = " << firing_share.transpose() << endl;
// //
// //    ArrayXXd VecZeros = ArrayXXd::Zero(N_C,4);
// //    ArrayXXd VecOnes = ArrayXXd::Ones(N_C,4);
// //
// //    ArrayXXd EntryCost(N_C,4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(N_C);
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(N_C);
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(N_C);
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(N_C);
// //
// ////    //// Only change Residual value
// //    string FKL = "_Partial_C5_OnlyResVal_FiringCost";
// //    tuple<ArrayXXd,ArrayXXd> t_Cal_ResV = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue, VecOnes, EntryCost, FKL, threadsManagement);
// //
// //    //// Only change Firing Cost
// //    FKL = "_Partial_C5_ResVal_OnlyFiringCost";
// //    tuple<ArrayXXd,ArrayXXd> t_Cal_firingcost = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        VecZeros, firing_share, EntryCost, FKL, threadsManagement);
// //
// //    //// Change Residual value and Firing Cost
// //    FKL = "_Partial_C5_Both_ResVal_FiringCost";
// //    tuple<ArrayXXd,ArrayXXd> t_Cal_ResV_firingcost = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue, firing_share, EntryCost, FKL, threadsManagement);
// //
// //    return 0;
// //}
// //
// /////////////**************************************************************
// ///////////* Counterfactual 6: Entry Cost vs Exit cost
// ///////////**************************************************************/
// //int alias::SolveSimulate_CounterfacturalResult_6(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    int N_C = 6;
// //    ArrayXd ResidualValue_temp = ArrayXd::LinSpaced(N_C, 1,500);
// //    ArrayXXd ResidualValue(N_C,4);
// //    ResidualValue.col(0) = ResidualValue_temp;
// //    ResidualValue.col(1) = ResidualValue_temp;
// //    ResidualValue.col(2) = ResidualValue_temp;
// //    ResidualValue.col(3) = ResidualValue_temp;
// //
// //    ArrayXXd firing_share = ArrayXXd::Ones(N_C,4);
// //
// //    ArrayXXd EntryCost(N_C,4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(N_C);
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(N_C);
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(N_C);
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(N_C);
// //
// //    cout << "ResidualValue = " << ResidualValue.transpose() << endl;
// //    cout << "firing_share = " << firing_share.transpose() << endl;
// //    cout << "EntryCost = " << EntryCost.transpose() << endl;
// //    ////
// //    string FKL = "_C6_AdjResVal";
// //    tuple<ArrayXXd,ArrayXXd> t = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue, firing_share, EntryCost, FKL, threadsManagement);
// //    ArrayXXd ExitSubsidy = get<0>(t);
// //    ArrayXXd NumEntrants = get<1>(t);
// //
// //    cout << "ExitSubsidy = " << ExitSubsidy.transpose() << endl;
// //    cout << "NumEntrants = " << NumEntrants.transpose() << endl;
// //    //
// //    // ArrayXXd ExitSubsidy(N_C,4);
// //    // ExitSubsidy.col(0) << 1.21411e+06,7.49741e+06,1.89316e+07;
// //    // ExitSubsidy.col(1) << 52273.2, 1.75283e+07, 5.42987e+07;
// //    // ExitSubsidy.col(2) << 70293.5, 3.20121e+07, 2.47767e+08;
// //    // ExitSubsidy.col(3) << 61284.2, 2.69775e+07, 1.77749e+08;
// //    // //
// //    // ArrayXXd NumEntrants(N_C,4);
// //    // NumEntrants.col(0) << 104336,136933,169558;
// //    // NumEntrants.col(1) << 100058, 125158, 169260;
// //    // NumEntrants.col(2) << 100285, 182953, 667897;
// //    // NumEntrants.col(3) << 100134, 169477, 511835;
// //
// //    ArrayXXd EntrySubsidy_per_Entrants = ArrayXXd::Zero(N_C,4);
// //    for (int k = 0; k < N_C; ++k) {
// //        int GoodState = 1; int LaborIntensive = 1;
// //        cout << "k = " << k << "; GoodState = " << GoodState << "; LaborIntensive = " << LaborIntensive << endl;
// //        EntrySubsidy_per_Entrants(k,0) = CalEntrySubsidy_per_Entrants(theta_Est_S,
// //            Value_StatusQuo_GoodLaborInt, sim_var, SimTbar, SimN, GoodState, LaborIntensive, KsupplyElas,
// //            ExitSubsidy(k,0), NumEntrants(k,0),threadsManagement);
// //        GoodState = 0; LaborIntensive = 1;
// //        EntrySubsidy_per_Entrants(k,1) = CalEntrySubsidy_per_Entrants(theta_Est_S,
// //            Value_StatusQuo_BadLaborInt, sim_var, SimTbar, SimN, GoodState, LaborIntensive, KsupplyElas,
// //            ExitSubsidy(k,1), NumEntrants(k,1),threadsManagement);
// //        GoodState = 1; LaborIntensive = 0;
// //        EntrySubsidy_per_Entrants(k,2) = CalEntrySubsidy_per_Entrants(theta_Est_S,
// //            Value_StatusQuo_GoodCapitalInt, sim_var, SimTbar, SimN, GoodState, LaborIntensive, KsupplyElas,
// //            ExitSubsidy(k,2), NumEntrants(k,2),threadsManagement);
// //        GoodState = 0; LaborIntensive = 0;
// //        EntrySubsidy_per_Entrants(k,3) = CalEntrySubsidy_per_Entrants(theta_Est_S,
// //            Value_StatusQuo_BadCapitalInt, sim_var, SimTbar, SimN, GoodState, LaborIntensive, KsupplyElas,
// //            ExitSubsidy(k,3), NumEntrants(k,3),threadsManagement);
// //    }
// //
// //    cout << "EntrySubsidy_per_Entrants = " << EntrySubsidy_per_Entrants.transpose() << endl;
// //
// //    // EntrySubsidy_per_Entrants.col(0) << 11.8811,59.1564,125.16;
// //    // // EntrySubsidy_per_Entrants.col(1) << 52.1009, 204.395, 328.365;
// //    // // EntrySubsidy_per_Entrants.col(2) << 52.9803, 184.547, 264.177;
// //    // // EntrySubsidy_per_Entrants.col(3) << 52.7393, 186.683, 254.101;
// //    // // EntrySubsidy_per_Entrants = 31.4946 125.756 222.959
// //
// //    ArrayXXd EntryCost_new(N_C,4);
// //    EntryCost_new.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry - EntrySubsidy_per_Entrants.col(0);
// //    EntryCost_new.col(1) = Value_StatusQuo_BadLaborInt.F_Entry - EntrySubsidy_per_Entrants.col(1);
// //    EntryCost_new.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry - EntrySubsidy_per_Entrants.col(2);
// //    EntryCost_new.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry - EntrySubsidy_per_Entrants.col(3);
// //    cout << "EntryCost_new = " << EntryCost_new.transpose() << endl;
// //
// ////     EntryCost_new = 1409.27 1315.01 1217.81
// //
// //    ResidualValue = ArrayXXd::Zero(N_C,4);
// //    FKL = "_C6_AdjEntryCost";
// //    tuple<ArrayXXd,ArrayXXd> t_Cal_EntryCost = CalSimulation_CounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt,
// //        Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue, firing_share, EntryCost_new, FKL, threadsManagement);
// //
// //    return 0;
// //}
// ////
// //// ///////***** Counterfactual 8_1: Calculate total exit subsidy
// //// double alias::CalExitSubsidy(const ArrayXd & theta_Est_S, const CounterFactVariable & Value_StatusQuo_Segment,
// ////     const SimVar & sim_var, const int SimTbar, const int SimN, const int GoodState, const int LaborIntensive,
// ////     const double & KsupplyElas, const ArrayXd & ResidualValue, const ArrayXd & firing_share,
// ////     MultiThreads::Threads_Management & threadsManagement) {
// ////
// ////     ArrayXd theta_C = Assign_thetaC(theta_Est_S,ResidualValue,firing_share,GoodState,LaborIntensive);
// ////     cout << "theta_C = " << theta_C.transpose() << endl;
// ////
// ////     tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// ////         t_C = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C,
// ////         Value_StatusQuo_Segment.F_Entry, Value_StatusQuo_Segment.CapitalDemand,
// ////         Value_StatusQuo_Segment.PriceIndex, Value_StatusQuo_Segment.FirmMass,
// ////         sim_var, GoodState, LaborIntensive, 0.0, threadsManagement);
// ////     EquState0 Equ0_BadCapitalInt_C = get<6>(t_C);
// ////     double FirmMass_C = Equ0_BadCapitalInt_C.FirmMass;
// ////     SimData sim_data_C = get<9>(t_C);
// ////
// ////     int k_ResVFiringShr;
// ////     if (GoodState == 1 & LaborIntensive == 1) {k_ResVFiringShr = 0;}
// ////     else if (GoodState == 0 & LaborIntensive == 1) {k_ResVFiringShr = 1;}
// ////     else if (GoodState == 1 & LaborIntensive == 0) {k_ResVFiringShr = 2;}
// ////     else {k_ResVFiringShr = 3;}
// ////     double ExitSubsidy = ResidualValue(k_ResVFiringShr) * (sim_data_C.Capital * sim_data_C.State_S.cast<double>()).sum() * FirmMass_C;
// ////
// ////     return ExitSubsidy;
// //// }
// ////
// //
// ///**************************************************************
// //* **************************************************************
// //* Partial Equilibrium: functions.
// //* **************************************************************
// //**************************************************************/
// //double alias::SolveResidualValueMatchingTargetExitRate(const ArrayXd & theta_Est_S, const CounterFactVariable & Value_StatusQuo,
// //    const SimVar & sim_var, const int & GoodState, const int & LaborIntensive, const double & TargetExitRate,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd theta_C = theta_Est_S;
// //    double ResidualValue_up = 1000; double ResidualValue_low = -1000;
// //    double ResidualValue_mid;
// //
// //    int n_Est;
// //    if (GoodState == 1) {n_Est = 14;}
// //    else {n_Est = 15;}
// //
// //    for (size_t n = 0; n < 10000; ++n) {
// //        ResidualValue_mid = 0.5*ResidualValue_up + 0.5*ResidualValue_low;
// //        theta_C(n_Est) = theta_Est_S(n_Est) + ResidualValue_mid;
// //        cout << "theta_C = " << theta_C.transpose() << endl;
// //
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData> t_Equ_Counterfactual
// //            = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C, Value_StatusQuo.F_Entry,
// //            Value_StatusQuo.CapitalDemand,Value_StatusQuo.PriceIndex,
// //            Value_StatusQuo.FirmMass, sim_var, GoodState, LaborIntensive, 0.0,
// //            threadsManagement);
// //        EquState0 Equ0_C = get<6>(t_Equ_Counterfactual);
// //        double FirmMass_C = Equ0_C.FirmMass;
// //        SimData sim_data_C = get<9>(t_Equ_Counterfactual);
// //
// //        /*** Calculate Average Exit Rate ***/
// //        double AvgExitRate_C = calAvgExitRate(sim_data_C, FirmMass_C);
// //        if (AvgExitRate_C > TargetExitRate) {
// //            ResidualValue_up = ResidualValue_mid;
// //        }
// //        else {
// //            ResidualValue_low = ResidualValue_mid;
// //        }
// //        cout << "** TargetExitRate = " << TargetExitRate << "; AvgExitRate_C = " << AvgExitRate_C
// //             << "; ResidualValue_up = " << ResidualValue_up << "; ResidualValue_low = " << ResidualValue_low << endl;
// //        if (abs(ResidualValue_up - ResidualValue_low) < 1e-2) {
// //            break;
// //        }
// //    }
// //
// //    return ResidualValue_mid;
// //}
// //
// //double alias::SolvefiringshareMatchingTargetExitRate(const ArrayXd & theta_Est_S, const CounterFactVariable & Value_StatusQuo,
// //    const SimVar & sim_var, const int & GoodState, const int & LaborIntensive, const double & TargetExitRate,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ArrayXd theta_C = theta_Est_S;
// //    double firing_share_up = 10; double firing_share_low = 0.01;
// //    double firing_share_mid;
// //
// //    int n_Est_high_ur; int n_Est_low_ur; int n_Est_uc;
// //    if (GoodState == 1) {
// //        n_Est_high_ur = 4; n_Est_low_ur = 5; n_Est_uc = 9;
// //    }
// //    else {
// //        n_Est_high_ur = 6; n_Est_low_ur = 7; n_Est_uc = 10;
// //    }
// //
// //    for (size_t n = 0; n < 10000; ++n) {
// //
// //        firing_share_mid = 0.5*firing_share_up + 0.5*firing_share_low;
// //
// //        theta_C(n_Est_high_ur) = firing_share_mid * theta_Est_S(n_Est_high_ur);
// //        theta_C(n_Est_low_ur) = firing_share_mid * theta_Est_S(n_Est_low_ur);
// //        theta_C(n_Est_uc) = firing_share_mid * theta_Est_S(n_Est_uc);
// //
// //        cout << "theta_Est_S = " << theta_Est_S.transpose() << endl;
// //        cout << "theta_C = " << theta_C.transpose() << endl;
// //
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData> t_Equ_Counterfactual
// //            = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C, Value_StatusQuo.F_Entry,
// //            Value_StatusQuo.CapitalDemand,Value_StatusQuo.PriceIndex,
// //            Value_StatusQuo.FirmMass, sim_var, GoodState, LaborIntensive, 0.0,
// //            threadsManagement);
// //        EquState0 Equ0_C = get<6>(t_Equ_Counterfactual);
// //        double FirmMass_C = Equ0_C.FirmMass;
// //        SimData sim_data_C = get<9>(t_Equ_Counterfactual);
// //
// //        /*** Calculate Average Exit Rate ***/
// //        double AvgExitRate_C = calAvgExitRate(sim_data_C, FirmMass_C);
// //        if (AvgExitRate_C > TargetExitRate) {
// //            firing_share_low = firing_share_mid;
// //        }
// //        else {
// //            firing_share_up = firing_share_mid;
// //        }
// //        cout << "** TargetExitRate = " << TargetExitRate << "; AvgExitRate_C = " << AvgExitRate_C
// //             << "; firing_share_up = " << firing_share_up << "; firing_share_low = " << firing_share_low << endl;
// //        if (abs(firing_share_up - firing_share_low) < 1e-3) {
// //            break;
// //        }
// //    }
// //
// //    return firing_share_mid;
// //}
// //
// //////
// //tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// //    alias::SolveSimulate_Counterfactural_FEntry_V0_EquPrice( const ArrayXd & theta_C, const double & F_Entry,
// //    const double & CapitalStock, const double & PriceIndex_S, const double & FirmMass_S, const SimVar & sim_var_C,
// //    const int & GoodState, const int & LaborIntensive, const double & KsupplyElas,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    EquState0 Equ0_C;
// //    Equ0_C.F_Entry = F_Entry;
// //
// //    ParaEst para_est_C = constructParaEst(theta_C, GoodState, LaborIntensive);
// //    para_est_C.PI = 1.0;
// //    double p_K_C = Solve_Counterfactural_FEntry_V0_Solve_EquPrice(para_est_C, F_Entry, threadsManagement);
// //    ParaVec para_vec_C = constructParaFull(para_est_C, p_K_C);
// //
// //    //// Simulate the data to get the stationary distribution of capital employment
// //    /** simulate a series of random number **/
// //    int SimN_C = sim_var_C.lnphi_sim.rows();
// //    int SimTbar_C = sim_var_C.lnphi_sim.cols() - 1;
// //
// //    //// solve equilibrium under the counterfactual parameters and a guessed price p_K_mid
// //    //    Rev RevOpt = RevOpt_Prod(para_est, para_vec);
// //    //    ArrayXd VEnd = calVEnd(para_est, para_vec, RevOpt);
// //    ArrayXd VEnd = ArrayXd::Zero(para.N*2);
// //    //// solve value function for t > 1
// //    tuple<EquStateV,EquStateVmat> t_Equ = solveV(para_est_C,para_vec_C,VEnd,threadsManagement);
// //    EquStateV EquV_C = get<0>(t_Equ);
// //    EquStateVmat EValmat_C = get<1>(t_Equ);
// //
// //    //// With free entry condition, solve value/policy function at the initial t = 1
// //    //// not done yet.
// //    tuple<EquStateV0,EquStateV0mat> t_Equ0 = solveV0(para_est_C,para_vec_C,EquV_C.EVal_PD,
// //        threadsManagement);
// //    EquStateV0 EquV0_C = get<0>(t_Equ0);
// //    EquStateV0mat EVal0mat_C = get<1>(t_Equ0);
// //
// //    /** simulate the data under the counterfactual data **/
// //    SimData sim_data_C = Simulation_StatusQuo_Counterfactual(para_est_C, para_vec_C, EquV_C,
// //        EValmat_C, EquV0_C, EVal0mat_C, sim_var_C, SimN_C,
// //        SimTbar_C, GoodState, LaborIntensive, threadsManagement);
// //
// //    double CapitalSupply = CapitalStock * pow(p_K_C,KsupplyElas);
// //    double Capital_Demand = CapitalSupply;
// //    double FirmMass_C = CapitalSupply / (sim_data_C.Capital * sim_data_C.State_S.cast<double>()).sum();
// //    Equ0_C.FirmMass = FirmMass_C;
// //
// //    return tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>(para_est_C,
// //        para_vec_C,EquV_C,EValmat_C,EquV0_C,EVal0mat_C,Equ0_C,p_K_C,Capital_Demand,sim_data_C);
// //}
// //
// //double alias::Solve_Counterfactural_FEntry_V0_Solve_EquPrice(const ParaEst & para_est_C, const double & F_Entry,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    double p_K_low = 0;
// //    double p_K_up = 100;
// //    double p_K_mid = 1.0;
// //
// //    for (size_t n = 0; n < 1000; ++n) {
// //        p_K_mid = 0.5*p_K_low + 0.5*p_K_up;
// //        cout << "p_K_mid = " << p_K_mid << endl;
// //        //// construct the grid of state space
// //        ParaVec para_vec_C = constructParaFull(para_est_C, p_K_mid);
// //
// //        //// solve equilibrium under the counterfactual parameters and a guessed price p_K_mid
// //        //    Rev RevOpt = RevOpt_Prod(para_est, para_vec);
// //        //    ArrayXd VEnd = calVEnd(para_est, para_vec, RevOpt);
// //        ArrayXd VEnd = ArrayXd::Zero(para.N*2);
// //        //// solve value function for t > 1
// //        tuple<EquStateV,EquStateVmat> t_Equ = solveV(para_est_C,para_vec_C,VEnd,threadsManagement);
// //        EquStateV EquV_C = get<0>(t_Equ);
// //        EquStateVmat EValmat_C = get<1>(t_Equ);
// //
// //        //// With free entry condition, solve value/policy function at the initial t = 1
// //        //// not done yet.
// //        tuple<EquStateV0,EquStateV0mat> t_Equ0 = solveV0(para_est_C,para_vec_C,EquV_C.EVal_PD,
// //            threadsManagement);
// //        EquStateV0 EquV0_C = get<0>(t_Equ0);
// //        EquStateV0mat EVal0mat_C = get<1>(t_Equ0);
// //
// //        //// calculate the expected value of entry under the counterfactual parameters and a guessed price p_K_mid
// //        double V_Entry_C = 0.0;
// //        for (size_t k = 0; k < para.N_phi; ++k) {
// //            V_Entry_C = V_Entry_C + para_vec_C.dist0(k) * EquV0_C.EVal_PD0(k * para.N_KL);
// //        }
// //
// //        //// update p_K
// //        if (V_Entry_C >= F_Entry) {
// //            p_K_low = p_K_mid;
// //        } else {
// //            p_K_up = p_K_mid;
// //        }
// //        cout << "p_K_low = " << p_K_low << "; p_K_up = " << p_K_up << "; p_K_mid = " << p_K_mid
// //            << "; diff of p_K_up and p_K_low = " << abs(p_K_up - p_K_low)
// //            << "; V_Entry_C = " << V_Entry_C << "; F_Entry = " << F_Entry << endl;
// //        if (abs(p_K_up - p_K_low) < 1e-3) {
// //            break;
// //        }
// ////        throw runtime_error("273");
// //    }
// //
// //    double p_K = p_K_mid;
// //    return p_K;
// //}
// //
// //tuple<ArrayXXd,ArrayXXd> alias::CalSimulation_CounterfactualResult(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    const ArrayXXd & ResidualValue, const ArrayXXd & firing_share, const ArrayXXd & EntryCost, const string & FKL,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    int N_C = ResidualValue.rows();
// //
// //    CFDeltaVec CFDeltaVec_GoodLaborInt = InitializeCFDeltaVec(N_C);
// //    CFDeltaVec CFDeltaVec_BadLaborInt = InitializeCFDeltaVec(N_C);
// //    CFDeltaVec CFDeltaVec_GoodCapitalInt = InitializeCFDeltaVec(N_C);
// //    CFDeltaVec CFDeltaVec_BadCapitalInt = InitializeCFDeltaVec(N_C);
// //
// //    CFDeltaVec_GoodLaborInt.firing_share = firing_share.col(0);
// //    CFDeltaVec_BadLaborInt.firing_share = firing_share.col(1);
// //    CFDeltaVec_GoodCapitalInt.firing_share = firing_share.col(2);
// //    CFDeltaVec_BadCapitalInt.firing_share = firing_share.col(3);
// //
// //    CFDeltaVec_GoodLaborInt.ResidualValue = ResidualValue.col(0);
// //    CFDeltaVec_BadLaborInt.ResidualValue = ResidualValue.col(1);
// //    CFDeltaVec_GoodCapitalInt.ResidualValue = ResidualValue.col(2);
// //    CFDeltaVec_BadCapitalInt.ResidualValue = ResidualValue.col(3);
// //
// //    CFDeltaVec_GoodLaborInt.EntryCost = EntryCost.col(0);
// //    CFDeltaVec_BadLaborInt.EntryCost = EntryCost.col(1);
// //    CFDeltaVec_GoodCapitalInt.EntryCost = EntryCost.col(2);
// //    CFDeltaVec_BadCapitalInt.EntryCost = EntryCost.col(3);
// //
// //    ArrayXXd ExitSubsidy = ArrayXXd::Zero(N_C,4);
// //    ArrayXXd NumEntrants = ArrayXXd::Zero(N_C,4);
// //
// //    cout << "Start the loop: choose different residual value" << endl;
// //    for (size_t k = 0; k < N_C; ++k) {
// //    // for (size_t k = 4; k < 5; ++k) {
// //
// //        cout << "The counterfactual is simulated for k = " << k << "; ResidualValue = " << ResidualValue.row(k)
// //             << "; firing_share = " << firing_share.row(k) << "; EntryCost = " << EntryCost.row(k) << endl;
// //
// //        ArrayXd theta_C_ResFC_GoodLaborInd = Assign_thetaC(theta_Est_S,ResidualValue.row(k),
// //            firing_share.row(k),1,1);
// //        ArrayXd theta_C_ResFC_BadLaborInd = Assign_thetaC(theta_Est_S,ResidualValue.row(k),
// //            firing_share.row(k),0,1);
// //        ArrayXd theta_C_ResFC_GoodCapitalInd = Assign_thetaC(theta_Est_S,ResidualValue.row(k),
// //            firing_share.row(k),1,0);
// //        ArrayXd theta_C_ResFC_BadCapitalInd = Assign_thetaC(theta_Est_S,ResidualValue.row(k),
// //            firing_share.row(k),0,0);
// //        cout << "theta_Est_S = " << theta_Est_S.transpose() << endl;
// //        cout << "theta_C_ResFC_GoodLaborInd = " << theta_C_ResFC_GoodLaborInd.transpose() << endl;
// //        cout << "theta_C_ResFC_BadLaborInd = " << theta_C_ResFC_BadLaborInd.transpose() << endl;
// //        cout << "theta_C_ResFC_GoodCapitalInd = " << theta_C_ResFC_GoodCapitalInd.transpose() << endl;
// //        cout << "theta_C_ResFC_BadCapitalInd = " << theta_C_ResFC_BadCapitalInd.transpose() << endl;
// //
// //        int LaborIntensive = 1; int GoodState = 1;
// //        cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// //            t_C_GoodLaborInt = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C_ResFC_GoodLaborInd,
// //            EntryCost(k,0), Value_StatusQuo_GoodLaborInt.CapitalDemand,
// //            Value_StatusQuo_GoodLaborInt.PriceIndex, Value_StatusQuo_GoodLaborInt.FirmMass,
// //            sim_var, GoodState, LaborIntensive, KsupplyElas, threadsManagement);
// //        EquState0 Equ0_GoodLaborInt_C = get<6>(t_C_GoodLaborInt);
// //        double p_K_GoodLaborInt_C = get<7>(t_C_GoodLaborInt);
// //        double CapitalDemand_GoodLaborInt_C = get<8>(t_C_GoodLaborInt);
// //        SimData sim_data_GoodLaborInt_C = get<9>(t_C_GoodLaborInt);
// //
// //        ExitSubsidy(k,0) = 0.0;
// //        for (size_t t = 1; t<SimTbar; ++t) {
// //            ExitSubsidy(k,0) = ExitSubsidy(k,0) + pow(para.deltaV,t)*ResidualValue(k,0)
// //                * Equ0_GoodLaborInt_C.FirmMass * sim_data_GoodLaborInt_C.State_E.col(t).cast<double>().sum();
// //        }
// //        NumEntrants(k,0) = Equ0_GoodLaborInt_C.FirmMass * double(SimN);
// //        // AggValueAdded(k,0) = (sim_data_GoodLaborInt_C.State_S.cast<double>() * sim_data_GoodLaborInt_C.Revenue).sum()
// //
// //        LaborIntensive = 1; GoodState = 0;
// //        cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// //            t_C_BadLaborInt = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C_ResFC_BadLaborInd,
// //            EntryCost(k,1), Value_StatusQuo_BadLaborInt.CapitalDemand,
// //            Value_StatusQuo_BadLaborInt.PriceIndex, Value_StatusQuo_BadLaborInt.FirmMass,
// //            sim_var, GoodState, LaborIntensive, KsupplyElas, threadsManagement);
// //        EquState0 Equ0_BadLaborInt_C = get<6>(t_C_BadLaborInt);
// //        double p_K_BadLaborInt_C = get<7>(t_C_BadLaborInt);
// //        double CapitalDemand_BadLaborInt_C = get<8>(t_C_BadLaborInt);
// //        SimData sim_data_BadLaborInt_C = get<9>(t_C_BadLaborInt);
// //
// //        for (size_t t = 1; t<SimTbar; ++t) {
// //            ExitSubsidy(k,1) = ExitSubsidy(k,1) + pow(para.deltaV,t)*ResidualValue(k,1)
// //                * Equ0_BadLaborInt_C.FirmMass * sim_data_BadLaborInt_C.State_E.col(t).cast<double>().sum();
// //        }
// //        NumEntrants(k,1) = Equ0_BadLaborInt_C.FirmMass * double(SimN);
// //
// //        LaborIntensive = 0; GoodState = 1;
// //        cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// //            t_C_GoodCapitalInt = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C_ResFC_GoodCapitalInd,
// //            EntryCost(k,2), Value_StatusQuo_GoodCapitalInt.CapitalDemand,
// //            Value_StatusQuo_GoodCapitalInt.PriceIndex, Value_StatusQuo_GoodCapitalInt.FirmMass,
// //            sim_var, GoodState, LaborIntensive, KsupplyElas, threadsManagement);
// //        EquState0 Equ0_GoodCapitalInt_C = get<6>(t_C_GoodCapitalInt);
// //        double p_K_GoodCapitalInt_C = get<7>(t_C_GoodCapitalInt);
// //        double CapitalDemand_GoodCapitalInt_C = get<8>(t_C_GoodCapitalInt);
// //        SimData sim_data_GoodCapitalInt_C = get<9>(t_C_GoodCapitalInt);
// //
// //        for (size_t t = 1; t<SimTbar; ++t) {
// //            ExitSubsidy(k,2) = ExitSubsidy(k,2) + pow(para.deltaV,t)*ResidualValue(k,2)
// //                * Equ0_GoodCapitalInt_C.FirmMass * sim_data_GoodCapitalInt_C.State_E.col(t).cast<double>().sum();
// //        }
// //        NumEntrants(k,2) = Equ0_GoodCapitalInt_C.FirmMass * double(SimN);
// //
// //        LaborIntensive = 0; GoodState = 0;
// //        cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// //            t_C_BadCapitalInt = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C_ResFC_BadCapitalInd,
// //            EntryCost(k,3), Value_StatusQuo_BadCapitalInt.CapitalDemand,
// //            Value_StatusQuo_BadCapitalInt.PriceIndex, Value_StatusQuo_BadCapitalInt.FirmMass,
// //            sim_var, GoodState, LaborIntensive, KsupplyElas, threadsManagement);
// //        EquState0 Equ0_BadCapitalInt_C = get<6>(t_C_BadCapitalInt);
// //        double p_K_BadCapitalInt_C = get<7>(t_C_BadCapitalInt);
// //        double CapitalDemand_BadCapitalInt_C = get<8>(t_C_BadCapitalInt);
// //        SimData sim_data_BadCapitalInt_C = get<9>(t_C_BadCapitalInt);
// //
// //        for (size_t t = 1; t<SimTbar; ++t) {
// //            ExitSubsidy(k,3) = ExitSubsidy(k,3) + pow(para.deltaV,t)*ResidualValue(k,3)
// //                * Equ0_BadCapitalInt_C.FirmMass * sim_data_BadCapitalInt_C.State_E.col(t).cast<double>().sum();
// //        }
// //        NumEntrants(k,3) = Equ0_BadCapitalInt_C.FirmMass * double(SimN);
// //
// //        CounterFactVariable Value_C_GoodLaborInt = CalEconomicVariable(p_K_GoodLaborInt_C,
// //            Equ0_GoodLaborInt_C.FirmMass,Equ0_GoodLaborInt_C.F_Entry,CapitalDemand_GoodLaborInt_C,
// //            sim_data_GoodLaborInt_C);
// //        CounterFactVariable Value_C_BadLaborInt = CalEconomicVariable(p_K_BadLaborInt_C,
// //            Equ0_BadLaborInt_C.FirmMass,Equ0_BadLaborInt_C.F_Entry,CapitalDemand_BadLaborInt_C,
// //            sim_data_BadLaborInt_C);
// //        CounterFactVariable Value_C_GoodCapitalInt = CalEconomicVariable(p_K_GoodCapitalInt_C,
// //            Equ0_GoodCapitalInt_C.FirmMass,Equ0_GoodCapitalInt_C.F_Entry,CapitalDemand_GoodCapitalInt_C,
// //            sim_data_GoodCapitalInt_C);
// //        CounterFactVariable Value_C_BadCapitalInt = CalEconomicVariable(p_K_BadCapitalInt_C,
// //            Equ0_BadCapitalInt_C.FirmMass,Equ0_BadCapitalInt_C.F_Entry,CapitalDemand_BadCapitalInt_C,
// //            sim_data_BadCapitalInt_C);
// //
// //        cout << "Counterfacutal:: EconVal.K_L_ratio: GoodLabor = " << Value_C_GoodLaborInt.K_L_ratio <<
// //            "; BadLabor " << Value_C_BadLaborInt.K_L_ratio <<
// //            "; GoodCapital " << Value_C_GoodCapitalInt.K_L_ratio <<
// //            "; BadCapital " << Value_C_BadCapitalInt.K_L_ratio << endl;
// //
// //        CounterFactVariable tDelta_GoodLaborInt = CalEconomicVariableDelta(Value_C_GoodLaborInt,
// //            Value_StatusQuo_GoodLaborInt);
// //        CFDeltaVec_GoodLaborInt.p_K_Delta(k) = tDelta_GoodLaborInt.p_K_Delta;
// //        CFDeltaVec_GoodLaborInt.FirmMass_Delta(k) = tDelta_GoodLaborInt.FirmMass_Delta;
// //        CFDeltaVec_GoodLaborInt.CapitalDemand_Delta(k) = tDelta_GoodLaborInt.CapitalDemand_Delta;
// //        CFDeltaVec_GoodLaborInt.NumFirm_Delta(k) = tDelta_GoodLaborInt.NumFirm_Delta;
// //        CFDeltaVec_GoodLaborInt.AvgProd_Delta(k) = tDelta_GoodLaborInt.AvgProd_Delta;
// //        CFDeltaVec_GoodLaborInt.WAvgProd_Delta(k) = tDelta_GoodLaborInt.WAvgProd_Delta;
// //        CFDeltaVec_GoodLaborInt.EntrantProd_Delta(k) = tDelta_GoodLaborInt.EntrantProd_Delta;
// //        CFDeltaVec_GoodLaborInt.ExitProd_Delta(k) = tDelta_GoodLaborInt.ExitProd_Delta;
// //        CFDeltaVec_GoodLaborInt.TotEmploy_ur_Delta(k) = tDelta_GoodLaborInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_GoodLaborInt.TotEmploy_uc_Delta(k) = tDelta_GoodLaborInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_GoodLaborInt.TotEmploy_Delta(k) = tDelta_GoodLaborInt.TotEmploy_Delta;
// //        CFDeltaVec_GoodLaborInt.TotOutput_Delta(k) = tDelta_GoodLaborInt.TotOutput_Delta;
// //        CFDeltaVec_GoodLaborInt.Output_per_Capital_Delta(k) = tDelta_GoodLaborInt.Output_per_Capital_Delta;
// //        CFDeltaVec_GoodLaborInt.Output_per_Employ_Delta(k) = tDelta_GoodLaborInt.Output_per_Employ_Delta;
// //        CFDeltaVec_GoodLaborInt.DormLength(k) = tDelta_GoodLaborInt.DormLength;
// //        CFDeltaVec_GoodLaborInt.ProdLength(k) = tDelta_GoodLaborInt.ProdLength;
// //        CFDeltaVec_GoodLaborInt.SurvivalLength(k) = tDelta_GoodLaborInt.SurvivalLength;
// //        CFDeltaVec_GoodLaborInt.PriceIndex_Delta(k) = tDelta_GoodLaborInt.PriceIndex_Delta;
// //        CFDeltaVec_GoodLaborInt.AvgExitRate(k) = tDelta_GoodLaborInt.AvgExitRate;
// //        CFDeltaVec_GoodLaborInt.K_L_ratio_Delta(k) = tDelta_GoodLaborInt.K_L_ratio_Delta;
// //        CFDeltaVec_GoodLaborInt.AggValueAdded(k) = tDelta_GoodLaborInt.AggValueAdded;
// //        cout << "p_K = " << tDelta_GoodLaborInt.p_K_Delta
// //            << "; tDelta_GoodLaborInt.TotOutput_Delta = " << tDelta_GoodLaborInt.TotOutput_Delta
// //            << "; tDelta_GoodLaborInt.K_L_ratio_Delta = " << tDelta_GoodLaborInt.K_L_ratio_Delta
// //            << "; tDelta_GoodLaborInt.AvgExitRate = " << tDelta_GoodLaborInt.AvgExitRate
// //            << "; Equ0_GoodLaborInt_C.FirmMass = " << Equ0_GoodLaborInt_C.FirmMass
// //            << "; tDelta_GoodLaborInt.TotEmploy_ur_Delta = " << tDelta_GoodLaborInt.TotEmploy_ur_Delta
// //            << "; tDelta_GoodLaborInt.TotEmploy_uc_Delta = " << tDelta_GoodLaborInt.TotEmploy_uc_Delta
// //            << "; tDelta_GoodLaborInt.TotEmploy_Delta = " << tDelta_GoodLaborInt.TotEmploy_Delta
// //            << endl;
// //
// //        CounterFactVariable tDelta_BadLaborInt = CalEconomicVariableDelta(Value_C_BadLaborInt,
// //            Value_StatusQuo_BadLaborInt);
// //        CFDeltaVec_BadLaborInt.p_K_Delta(k) = tDelta_BadLaborInt.p_K_Delta;
// //        CFDeltaVec_BadLaborInt.FirmMass_Delta(k) = tDelta_BadLaborInt.FirmMass_Delta;
// //        CFDeltaVec_BadLaborInt.CapitalDemand_Delta(k) = tDelta_BadLaborInt.CapitalDemand_Delta;
// //        CFDeltaVec_BadLaborInt.NumFirm_Delta(k) = tDelta_BadLaborInt.NumFirm_Delta;
// //        CFDeltaVec_BadLaborInt.AvgProd_Delta(k) = tDelta_BadLaborInt.AvgProd_Delta;
// //        CFDeltaVec_BadLaborInt.WAvgProd_Delta(k) = tDelta_BadLaborInt.WAvgProd_Delta;
// //        CFDeltaVec_BadLaborInt.EntrantProd_Delta(k) = tDelta_BadLaborInt.EntrantProd_Delta;
// //        CFDeltaVec_BadLaborInt.ExitProd_Delta(k) = tDelta_BadLaborInt.ExitProd_Delta;
// //        CFDeltaVec_BadLaborInt.TotEmploy_ur_Delta(k) = tDelta_BadLaborInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_BadLaborInt.TotEmploy_uc_Delta(k) = tDelta_BadLaborInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_BadLaborInt.TotEmploy_Delta(k) = tDelta_BadLaborInt.TotEmploy_Delta;
// //        CFDeltaVec_BadLaborInt.TotOutput_Delta(k) = tDelta_BadLaborInt.TotOutput_Delta;
// //        CFDeltaVec_BadLaborInt.Output_per_Capital_Delta(k) = tDelta_BadLaborInt.Output_per_Capital_Delta;
// //        CFDeltaVec_BadLaborInt.Output_per_Employ_Delta(k) = tDelta_BadLaborInt.Output_per_Employ_Delta;
// //        CFDeltaVec_BadLaborInt.DormLength(k) = tDelta_BadLaborInt.DormLength;
// //        CFDeltaVec_BadLaborInt.ProdLength(k) = tDelta_BadLaborInt.ProdLength;
// //        CFDeltaVec_BadLaborInt.SurvivalLength(k) = tDelta_BadLaborInt.SurvivalLength;
// //        CFDeltaVec_BadLaborInt.PriceIndex_Delta(k) = tDelta_BadLaborInt.PriceIndex_Delta;
// //        CFDeltaVec_BadLaborInt.AvgExitRate(k) = tDelta_BadLaborInt.AvgExitRate;
// //        CFDeltaVec_BadLaborInt.K_L_ratio_Delta(k) = tDelta_BadLaborInt.K_L_ratio_Delta;
// //        CFDeltaVec_BadLaborInt.AggValueAdded(k) = tDelta_BadLaborInt.AggValueAdded;
// //        cout << "p_K = " << tDelta_BadLaborInt.p_K_Delta
// //            << "; tDelta_BadLaborInt.TotOutput_Delta = " << tDelta_BadLaborInt.TotOutput_Delta
// //            << "; tDelta_BadLaborInt.K_L_ratio_Delta = " << tDelta_BadLaborInt.K_L_ratio_Delta
// //            << "; tDelta_BadLaborInt.AvgExitRate = " << tDelta_BadLaborInt.AvgExitRate
// //            << "; Equ0_BadLaborInt_C.FirmMass = " << Equ0_BadLaborInt_C.FirmMass
// //            << "; tDelta_BadLaborInt.TotEmploy_ur_Delta = " << tDelta_BadLaborInt.TotEmploy_ur_Delta
// //            << "; tDelta_BadLaborInt.TotEmploy_uc_Delta = " << tDelta_BadLaborInt.TotEmploy_uc_Delta
// //            << "; tDelta_BadLaborInt.TotEmploy_Delta = " << tDelta_BadLaborInt.TotEmploy_Delta
// //            << endl;
// //
// //        CounterFactVariable tDelta_GoodCapitalInt = CalEconomicVariableDelta(Value_C_GoodCapitalInt,
// //            Value_StatusQuo_GoodCapitalInt);
// //        CFDeltaVec_GoodCapitalInt.p_K_Delta(k) = tDelta_GoodCapitalInt.p_K_Delta;
// //        CFDeltaVec_GoodCapitalInt.FirmMass_Delta(k) = tDelta_GoodCapitalInt.FirmMass_Delta;
// //        CFDeltaVec_GoodCapitalInt.CapitalDemand_Delta(k) = tDelta_GoodCapitalInt.CapitalDemand_Delta;
// //        CFDeltaVec_GoodCapitalInt.NumFirm_Delta(k) = tDelta_GoodCapitalInt.NumFirm_Delta;
// //        CFDeltaVec_GoodCapitalInt.AvgProd_Delta(k) = tDelta_GoodCapitalInt.AvgProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.WAvgProd_Delta(k) = tDelta_GoodCapitalInt.WAvgProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.EntrantProd_Delta(k) = tDelta_GoodCapitalInt.EntrantProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.ExitProd_Delta(k) = tDelta_GoodCapitalInt.ExitProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotEmploy_ur_Delta(k) = tDelta_GoodCapitalInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotEmploy_uc_Delta(k) = tDelta_GoodCapitalInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotEmploy_Delta(k) = tDelta_GoodCapitalInt.TotEmploy_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotOutput_Delta(k) = tDelta_GoodCapitalInt.TotOutput_Delta;
// //        CFDeltaVec_GoodCapitalInt.Output_per_Capital_Delta(k) = tDelta_GoodCapitalInt.Output_per_Capital_Delta;
// //        CFDeltaVec_GoodCapitalInt.Output_per_Employ_Delta(k) = tDelta_GoodCapitalInt.Output_per_Employ_Delta;
// //        CFDeltaVec_GoodCapitalInt.DormLength(k) = tDelta_GoodCapitalInt.DormLength;
// //        CFDeltaVec_GoodCapitalInt.ProdLength(k) = tDelta_GoodCapitalInt.ProdLength;
// //        CFDeltaVec_GoodCapitalInt.SurvivalLength(k) = tDelta_GoodCapitalInt.SurvivalLength;
// //        CFDeltaVec_GoodCapitalInt.PriceIndex_Delta(k) = tDelta_GoodCapitalInt.PriceIndex_Delta;
// //        CFDeltaVec_GoodCapitalInt.AvgExitRate(k) = tDelta_GoodCapitalInt.AvgExitRate;
// //        CFDeltaVec_GoodCapitalInt.K_L_ratio_Delta(k) = tDelta_GoodCapitalInt.K_L_ratio_Delta;
// //        CFDeltaVec_GoodCapitalInt.AggValueAdded(k) = tDelta_GoodCapitalInt.AggValueAdded;
// //        cout << "p_K = " << tDelta_GoodCapitalInt.p_K_Delta
// //            << "; tDelta_GoodCapitalInt.TotOutput_Delta = " << tDelta_GoodCapitalInt.TotOutput_Delta
// //            << "; tDelta_GoodCapitalInt.K_L_ratio_Delta = " << tDelta_GoodCapitalInt.K_L_ratio_Delta
// //            << "; tDelta_GoodCapitalInt.AvgExitRate = " << tDelta_GoodCapitalInt.AvgExitRate
// //            << "; Equ0_GoodCapitalInt_C.FirmMass = " << Equ0_GoodCapitalInt_C.FirmMass
// //            << "; tDelta_GoodCapitalInt.TotEmploy_ur_Delta = " << tDelta_GoodCapitalInt.TotEmploy_ur_Delta
// //            << "; tDelta_GoodCapitalInt.TotEmploy_uc_Delta = " << tDelta_GoodCapitalInt.TotEmploy_uc_Delta
// //            << "; tDelta_GoodCapitalInt.TotEmploy_Delta = " << tDelta_GoodCapitalInt.TotEmploy_Delta
// //            << endl;
// //
// //        CounterFactVariable tDelta_BadCapitalInt = CalEconomicVariableDelta(Value_C_BadCapitalInt,
// //            Value_StatusQuo_BadCapitalInt);
// //        CFDeltaVec_BadCapitalInt.p_K_Delta(k) = tDelta_BadCapitalInt.p_K_Delta;
// //        CFDeltaVec_BadCapitalInt.FirmMass_Delta(k) = tDelta_BadCapitalInt.FirmMass_Delta;
// //        CFDeltaVec_BadCapitalInt.CapitalDemand_Delta(k) = tDelta_BadCapitalInt.CapitalDemand_Delta;
// //        CFDeltaVec_BadCapitalInt.NumFirm_Delta(k) = tDelta_BadCapitalInt.NumFirm_Delta;
// //        CFDeltaVec_BadCapitalInt.AvgProd_Delta(k) = tDelta_BadCapitalInt.AvgProd_Delta;
// //        CFDeltaVec_BadCapitalInt.WAvgProd_Delta(k) = tDelta_BadCapitalInt.WAvgProd_Delta;
// //        CFDeltaVec_BadCapitalInt.EntrantProd_Delta(k) = tDelta_BadCapitalInt.EntrantProd_Delta;
// //        CFDeltaVec_BadCapitalInt.ExitProd_Delta(k) = tDelta_BadCapitalInt.ExitProd_Delta;
// //        CFDeltaVec_BadCapitalInt.TotEmploy_ur_Delta(k) = tDelta_BadCapitalInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_BadCapitalInt.TotEmploy_uc_Delta(k) = tDelta_BadCapitalInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_BadCapitalInt.TotEmploy_Delta(k) = tDelta_BadCapitalInt.TotEmploy_Delta;
// //        CFDeltaVec_BadCapitalInt.TotOutput_Delta(k) = tDelta_BadCapitalInt.TotOutput_Delta;
// //        CFDeltaVec_BadCapitalInt.Output_per_Capital_Delta(k) = tDelta_BadCapitalInt.Output_per_Capital_Delta;
// //        CFDeltaVec_BadCapitalInt.Output_per_Employ_Delta(k) = tDelta_BadCapitalInt.Output_per_Employ_Delta;
// //        CFDeltaVec_BadCapitalInt.DormLength(k) = tDelta_BadCapitalInt.DormLength;
// //        CFDeltaVec_BadCapitalInt.ProdLength(k) = tDelta_BadCapitalInt.ProdLength;
// //        CFDeltaVec_BadCapitalInt.SurvivalLength(k) = tDelta_BadCapitalInt.SurvivalLength;
// //        CFDeltaVec_BadCapitalInt.PriceIndex_Delta(k) = tDelta_BadCapitalInt.PriceIndex_Delta;
// //        CFDeltaVec_BadCapitalInt.AvgExitRate(k) = tDelta_BadCapitalInt.AvgExitRate;
// //        CFDeltaVec_BadCapitalInt.K_L_ratio_Delta(k) = tDelta_BadCapitalInt.K_L_ratio_Delta;
// //        CFDeltaVec_BadCapitalInt.AggValueAdded(k) = tDelta_BadCapitalInt.AggValueAdded;
// //        cout << "p_K = " << tDelta_BadCapitalInt.p_K_Delta
// //            << "; tDelta_BadCapitalInt.TotOutput_Delta = " << tDelta_BadCapitalInt.TotOutput_Delta
// //            << "; tDelta_BadCapitalInt.K_L_ratio_Delta = " << tDelta_BadCapitalInt.K_L_ratio_Delta
// //            << "; tDelta_BadCapitalInt.AvgExitRate = " << tDelta_BadCapitalInt.AvgExitRate
// //            << "; Equ0_BadCapitalInt_C.FirmMass = " << Equ0_BadCapitalInt_C.FirmMass
// //            << "; tDelta_BadCapitalInt.TotEmploy_ur_Delta = " << tDelta_BadCapitalInt.TotEmploy_ur_Delta
// //            << "; tDelta_BadCapitalInt.TotEmploy_uc_Delta = " << tDelta_BadCapitalInt.TotEmploy_uc_Delta
// //            << "; tDelta_BadCapitalInt.TotEmploy_Delta = " << tDelta_BadCapitalInt.TotEmploy_Delta
// //            << endl;
// //
// //        // throw runtime_error("1543");
// //    }
// //
// //    CFDeltaVec_GoodLaborInt.ExitSubsidy = ExitSubsidy.col(0);
// //    CFDeltaVec_BadLaborInt.ExitSubsidy = ExitSubsidy.col(1);
// //    CFDeltaVec_GoodCapitalInt.ExitSubsidy = ExitSubsidy.col(2);
// //    CFDeltaVec_BadCapitalInt.ExitSubsidy = ExitSubsidy.col(3);
// //
// //    string filename = "State" + to_string(1) + "_Sec" + to_string(1) + FKL;
// //    int t_write_GoodLaborInt = writeToCSVfileValueDelta(filename,CFDeltaVec_GoodLaborInt);
// //
// //    filename = "State" + to_string(0) + "_Sec" + to_string(1) + FKL;
// //    int t_write_BadLaborInt = writeToCSVfileValueDelta(filename,CFDeltaVec_BadLaborInt);
// //
// //    filename = "State" + to_string(1) + "_Sec" + to_string(0) + FKL;
// //    int t_write_GoodCapitalInt = writeToCSVfileValueDelta(filename,CFDeltaVec_GoodCapitalInt);
// //
// //    filename = "State" + to_string(0) + "_Sec" + to_string(0) + FKL;
// //    int t_write_BadCapitalInt = writeToCSVfileValueDelta(filename,CFDeltaVec_BadCapitalInt);
// //
// //    return tuple<ArrayXXd,ArrayXXd>(ExitSubsidy,NumEntrants);
// //}
// //
// //
// //double alias::CalEntrySubsidy_per_Entrants(const ArrayXd & theta_Est_S, const CounterFactVariable & Value_StatusQuo_Segment,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const int & GoodState, const int & LaborIntensive,
// //    const double & KsupplyElas, const double & ExitSubsidy, const double & NumEntrants,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    double EntrySubsidy_ini = ExitSubsidy / NumEntrants;
// //    cout << "ExitSubsidy = " << ExitSubsidy << "; NumEntrants = " << NumEntrants
// //        << "; EntrySubsidy_ini = " << EntrySubsidy_ini << endl;
// //    double EntrySubsidy_up = 0;
// //
// //    for (size_t n = 0; n<30; ++n) {
// //
// //        double F_Entry_new = Value_StatusQuo_Segment.F_Entry - EntrySubsidy_ini;
// //
// //        tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat, EquState0, double, double, SimData>
// //            t_C = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_Est_S, F_Entry_new,
// //            Value_StatusQuo_Segment.CapitalDemand,Value_StatusQuo_Segment.PriceIndex,
// //            Value_StatusQuo_Segment.FirmMass, sim_var, GoodState, LaborIntensive, KsupplyElas,
// //            threadsManagement);
// //        EquState0 Equ0_C = get<6>(t_C);
// //        SimData sim_data_C = get<9>(t_C);
// //        double NumEntrants_C = Equ0_C.FirmMass * double(SimN);
// //
// //        EntrySubsidy_up = ExitSubsidy / NumEntrants_C;
// //
// //        double eps = abs(EntrySubsidy_up - EntrySubsidy_ini) / EntrySubsidy_ini;
// //        if (eps < 1e-3) {
// //            break;
// //        }
// //        EntrySubsidy_ini = 0.5*EntrySubsidy_up + 0.5*EntrySubsidy_ini;
// //        cout << "n = " << n << "; EntrySubsidy_ini = " << EntrySubsidy_ini << "; EntrySubsidy_up = " << EntrySubsidy_up
// //            << "; NumEntrants_C = " << NumEntrants_C << "; Equ0_C.FirmMass = " << Equ0_C.FirmMass
// //            << "; sim_data_C.State_S.col(0).cast<double>().sum() = " << sim_data_C.State_S.col(0).cast<double>().sum()
// //            << "; eps = " << eps << endl;
// //    }
// //
// //    double EntrySubsidy_per_Entrants = EntrySubsidy_ini;
// //    return EntrySubsidy_per_Entrants;
// //}
// //
// ///**************************************************************
// //* **************************************************************
// //* Solve the General Counterfactural equilibrium
// //* **************************************************************
// //**************************************************************/
// ///**************************************************************
// //*** Counterfactual Exercise 1: General Equilibrium: Varying the exit barrier ***
// //**************************************************************/
// //int alias::SolveSimulate_GECounterfacturalResult_1(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
// //    const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas, const ArrayXd & ResVal_vec,
// //    MultiThreads::Threads_Management &threadsManagement) {
// //
// //    int N_C = ResVal_vec.size();
// //
// //    ArrayXXd firingshare = ArrayXXd::Ones(ResVal_vec.rows(),4);
// //
// //    ArrayXXd EntryCost(ResVal_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(ResVal_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(ResVal_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(ResVal_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(ResVal_vec.rows());
// //
// //    double p_K_C_low = 1.0; double p_K_C_high = 2.0;
// //    ArrayXXd ResidualValue = ArrayXXd::Zero(ResVal_vec.rows(),4);
// //
// //    /*** change residual value in All Cases ***/
// //    ResidualValue = ArrayXXd::Zero(ResVal_vec.rows(),4);
// //    ResidualValue.col(0) << ResVal_vec;
// //    ResidualValue.col(1) << ResVal_vec;
// //    ResidualValue.col(2) << ResVal_vec;
// //    ResidualValue.col(3) << ResVal_vec;
// //    string FKL = "GE_C1_ResVal_PolicyStateAll";
// //    tuple<ArrayXXd,ArrayXXd> t_GE_C1_ResVal_PolicyStateAll = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,firingshare, EntryCost, FKL, p_K_C_low, p_K_C_high, threadsManagement);
// //    throw runtime_error("1806");
// //    // //
// //    /*** change residual value in Good States ***/
// //    ResidualValue = ArrayXXd::Zero(ResVal_vec.rows(),4);
// //    ResidualValue.col(0) << ResVal_vec;
// //    ResidualValue.col(2) << ResVal_vec;
// //    FKL = "GE_C1_ResVal_PolicyState1";
// //    tuple<ArrayXXd,ArrayXXd> t_GE_C1_ResVal_PolicyState1 = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,firingshare, EntryCost, FKL, p_K_C_low, p_K_C_high, threadsManagement);
// //
// //    /*** change residual value in Bad States ***/
// //    ResidualValue = ArrayXXd::Zero(ResVal_vec.rows(),4);
// //    ResidualValue.col(1) << ResVal_vec;
// //    ResidualValue.col(3) << ResVal_vec;
// //    FKL = "GE_C1_ResVal_PolicyState0";
// //    tuple<ArrayXXd,ArrayXXd> t_GE_C1_ResVal_PolicyState0 = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,firingshare, EntryCost, FKL, p_K_C_low, p_K_C_high, threadsManagement);
// //
// //    return 0;
// //}
// //
// //int alias::SolveSimulate_GECounterfacturalResult_2(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
// //    const SimVar &sim_var, const int SimTbar, const int SimN, const double &KsupplyElas, const ArrayXd &firingshare_vec,
// //    MultiThreads::Threads_Management &threadsManagement) {
// //
// //    int N_C = firingshare_vec.size();
// //
// //    ArrayXXd firingshare = ArrayXXd::Ones(firingshare_vec.rows(),4);
// //
// //    ArrayXXd EntryCost(firingshare_vec.rows(),4);
// //    EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(firingshare_vec.rows());
// //    EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(firingshare_vec.rows());
// //    EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(firingshare_vec.rows());
// //    EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(firingshare_vec.rows());
// //
// //    ArrayXXd ResidualValue = ArrayXXd::Zero(firingshare_vec.rows(),4);
// //
// //    double p_K_C_low = 1.0; double p_K_C_high = 2.0;
// //
// //    // /*** change residual value in All Cases ***/
// //    firingshare = ArrayXXd::Ones(firingshare_vec.rows(),4);
// //    firingshare.col(0) << firingshare_vec;
// //    firingshare.col(1) << firingshare_vec;
// //    firingshare.col(2) << firingshare_vec;
// //    firingshare.col(3) << firingshare_vec;
// //    string FKL = "GE_C2_FiringCost_PolicyStateAll";
// //    tuple<ArrayXXd,ArrayXXd> t_GE_C1_ResVal_PolicyStateAll = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// //        Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //        sim_var, SimTbar, SimN, KsupplyElas,
// //        ResidualValue,firingshare, EntryCost, FKL, p_K_C_low, p_K_C_high, threadsManagement);
// //    // //throw runtime_error("1657");
// //    // //
// //    // /*** change residual value in Good States ***/
// //    // firingshare = ArrayXXd::Ones(firingshare_vec.rows(),4);
// //    // firingshare.col(0) << firingshare_vec;
// //    // firingshare.col(2) << firingshare_vec;
// //    // FKL = "GE_C2_FiringCost_PolicyState1";
// //    // tuple<ArrayXXd,ArrayXXd> t_GE_C1_ResVal_PolicyState1 = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// //    //     Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //    //     sim_var, SimTbar, SimN, KsupplyElas,
// //    //     ResidualValue,firingshare, EntryCost, FKL, p_K_C_low, p_K_C_high, threadsManagement);
// //    // //
// //    // // /*** change residual value in Bad States ***/
// //    // firingshare = ArrayXXd::Ones(firingshare_vec.rows(),4);
// //    // firingshare.col(1) << firingshare_vec;
// //    // firingshare.col(3) << firingshare_vec;
// //    // FKL = "GE_C2_FiringCost_PolicyState0";
// //    // tuple<ArrayXXd,ArrayXXd> t_GE_C1_ResVal_PolicyState0 = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// //    //     Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// //    //     sim_var, SimTbar, SimN, KsupplyElas,
// //    //     ResidualValue,firingshare, EntryCost, FKL, p_K_C_low, p_K_C_high, threadsManagement);
// //
// //    return 0;
// //}
// //
// //tuple<ArrayXXd,ArrayXXd> alias::CalSimulation_GECounterfactualResult(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// //    const CounterFactVariable & Value_StatusQuo_GoodLaborInt, const CounterFactVariable & Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable & Value_StatusQuo_GoodCapitalInt, const CounterFactVariable & Value_StatusQuo_BadCapitalInt,
// //    const SimVar & sim_var, const int SimTbar, const int SimN, const double & KsupplyElas,
// //    const ArrayXXd & ResidualValue, const ArrayXXd & firing_share, const ArrayXXd & EntryCost, const string & FKL,
// //    const double & p_K_C_low, const double & p_K_C_high, MultiThreads::Threads_Management & threadsManagement) {
// //
// //    int N_C = ResidualValue.rows();
// //
// //    CFDeltaVec CFDeltaVec_GoodLaborInt = InitializeCFDeltaVec(N_C);
// //    CFDeltaVec CFDeltaVec_BadLaborInt = InitializeCFDeltaVec(N_C);
// //    CFDeltaVec CFDeltaVec_GoodCapitalInt = InitializeCFDeltaVec(N_C);
// //    CFDeltaVec CFDeltaVec_BadCapitalInt = InitializeCFDeltaVec(N_C);
// //
// //    CFDeltaVec_GoodLaborInt.firing_share = firing_share.col(0);
// //    CFDeltaVec_BadLaborInt.firing_share = firing_share.col(1);
// //    CFDeltaVec_GoodCapitalInt.firing_share = firing_share.col(2);
// //    CFDeltaVec_BadCapitalInt.firing_share = firing_share.col(3);
// //
// //    CFDeltaVec_GoodLaborInt.ResidualValue = ResidualValue.col(0);
// //    CFDeltaVec_BadLaborInt.ResidualValue = ResidualValue.col(1);
// //    CFDeltaVec_GoodCapitalInt.ResidualValue = ResidualValue.col(2);
// //    CFDeltaVec_BadCapitalInt.ResidualValue = ResidualValue.col(3);
// //
// //    CFDeltaVec_GoodLaborInt.EntryCost = EntryCost.col(0);
// //    CFDeltaVec_BadLaborInt.EntryCost = EntryCost.col(1);
// //    CFDeltaVec_GoodCapitalInt.EntryCost = EntryCost.col(2);
// //    CFDeltaVec_BadCapitalInt.EntryCost = EntryCost.col(3);
// //
// //    ArrayXXd ExitSubsidy = ArrayXXd::Zero(N_C,5);
// //    ArrayXXd NumEntrants = ArrayXXd::Zero(N_C,5);
// //
// //    cout << "Start the loop: choose different residual value" << endl;
// //    for (size_t k = 0; k < N_C; ++k) {
// //        cout << "The counterfactual is simulated for k = " << k << "; ResidualValue = " << ResidualValue.row(k)
// //             << "; firing_share = " << firing_share.row(k) << "; EntryCost = " << EntryCost.row(k) << endl;
// //
// //        ArrayXd theta_C_ResFC_GoodLaborInt = Assign_thetaC(theta_Est_S, ResidualValue.row(k),
// //            firing_share.row(k), 1, 1);
// //        ArrayXd theta_C_ResFC_BadLaborInt = Assign_thetaC(theta_Est_S, ResidualValue.row(k),
// //            firing_share.row(k), 0, 1);
// //        ArrayXd theta_C_ResFC_GoodCapitalInt = Assign_thetaC(theta_Est_S, ResidualValue.row(k),
// //            firing_share.row(k), 1, 0);
// //        ArrayXd theta_C_ResFC_BadCapitalInt = Assign_thetaC(theta_Est_S, ResidualValue.row(k),
// //            firing_share.row(k), 0, 0);
// //        cout << "theta_Est_S = " << theta_Est_S.transpose() << endl;
// //        cout << "theta_C_ResFC_GoodLaborInd = " << theta_C_ResFC_GoodLaborInt.transpose() << endl;
// //        cout << "theta_C_ResFC_BadLaborInd = " << theta_C_ResFC_BadLaborInt.transpose() << endl;
// //        cout << "theta_C_ResFC_GoodCapitalInd = " << theta_C_ResFC_GoodCapitalInt.transpose() << endl;
// //        cout << "theta_C_ResFC_BadCapitalInd = " << theta_C_ResFC_BadCapitalInt.transpose() << endl;
// //
// //        EquState0 Equ0_GoodLaborInt_C;
// //        Equ0_GoodLaborInt_C.F_Entry = CFDeltaVec_GoodLaborInt.EntryCost(k);
// //        EquState0 Equ0_BadLaborInt_C;
// //        Equ0_BadLaborInt_C.F_Entry = CFDeltaVec_BadLaborInt.EntryCost(k);
// //        EquState0 Equ0_GoodCapitalInt_C;
// //        Equ0_GoodLaborInt_C.F_Entry = CFDeltaVec_GoodCapitalInt.EntryCost(k);
// //        EquState0 Equ0_BadCapitalInt_C;
// //        Equ0_GoodLaborInt_C.F_Entry = CFDeltaVec_BadCapitalInt.EntryCost(k);
// //
// //        double CapitalStock = shr_RealData.shr_GoodState_LaborInt * Value_StatusQuo_GoodLaborInt.CapitalDemand
// //            + shr_RealData.shr_BadState_LaborInt * Value_StatusQuo_BadLaborInt.CapitalDemand
// //            + shr_RealData.shr_GoodState_CapitalInt * Value_StatusQuo_GoodCapitalInt.CapitalDemand
// //            + shr_RealData.shr_BadState_CapitalInt * Value_StatusQuo_BadCapitalInt.CapitalDemand;
// //        cout << "CapitalStock = " << CapitalStock << endl;
// //
// //        double CapitalDemand_GoodLaborInt_C;
// //        double CapitalDemand_BadLaborInt_C;
// //        double CapitalDemand_GoodCapitalInt_C;
// //        double CapitalDemand_BadCapitalInt_C;
// //
// //        SimData sim_data_C_GoodLaborInt;
// //        SimData sim_data_C_BadLaborInt;
// //        SimData sim_data_C_GoodCapitalInt;
// //        SimData sim_data_C_BadCapitalInt;
// //
// //        double p_K_C_ini = p_K_C_low;
// //        double p_K_C_up = p_K_C_high;
// //        for (size_t n = 0; n<50; ++n) {
// //
// //            int LaborIntensive = 1; int GoodState = 1;
// //            cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //            tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat> t_GoodLaborInt =
// //                Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(theta_C_ResFC_GoodLaborInt,
// //                CFDeltaVec_GoodLaborInt.EntryCost(k), p_K_C_ini, GoodState, LaborIntensive, threadsManagement);
// //            ParaEst para_est_C_GoodLaborInt = get<0>(t_GoodLaborInt);
// //            ParaVec para_vec_C_GoodLaborInt = get<1>(t_GoodLaborInt);
// //            EquStateV EquV_C_GoodLaborInt = get<2>(t_GoodLaborInt);
// //            EquStateVmat EValmat_C_GoodLaborInt = get<3>(t_GoodLaborInt);
// //            EquStateV0 EquV0_C_GoodLaborInt = get<4>(t_GoodLaborInt);
// //            EquStateV0mat EVal0mat_C_GoodLaborInt = get<5>(t_GoodLaborInt);
// //
// //            /** simulate the data under the counterfactual data **/
// //            sim_data_C_GoodLaborInt = Simulation_StatusQuo_Counterfactual(para_est_C_GoodLaborInt,
// //                para_vec_C_GoodLaborInt, EquV_C_GoodLaborInt, EValmat_C_GoodLaborInt,
// //                EquV0_C_GoodLaborInt, EVal0mat_C_GoodLaborInt, sim_var, SimN, SimTbar,
// //                GoodState, LaborIntensive, threadsManagement);
// //
// //            LaborIntensive = 1; GoodState = 0;
// //            cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //            tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat> t_BadLaborInt =
// //                Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(theta_C_ResFC_BadLaborInt,
// //                CFDeltaVec_BadLaborInt.EntryCost(k), p_K_C_ini, GoodState, LaborIntensive, threadsManagement);
// //            ParaEst para_est_C_BadLaborInt = get<0>(t_BadLaborInt);
// //            ParaVec para_vec_C_BadLaborInt = get<1>(t_BadLaborInt);
// //            EquStateV EquV_C_BadLaborInt = get<2>(t_BadLaborInt);
// //            EquStateVmat EValmat_C_BadLaborInt = get<3>(t_BadLaborInt);
// //            EquStateV0 EquV0_C_BadLaborInt = get<4>(t_BadLaborInt);
// //            EquStateV0mat EVal0mat_C_BadLaborInt = get<5>(t_BadLaborInt);
// //
// //            /** simulate the data under the counterfactual data **/
// //            sim_data_C_BadLaborInt = Simulation_StatusQuo_Counterfactual(para_est_C_BadLaborInt,
// //                para_vec_C_BadLaborInt, EquV_C_BadLaborInt, EValmat_C_BadLaborInt,
// //                EquV0_C_BadLaborInt, EVal0mat_C_BadLaborInt, sim_var, SimN, SimTbar,
// //                GoodState, LaborIntensive, threadsManagement);
// //
// //            LaborIntensive = 0; GoodState = 1;
// //            cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //            tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat> t_GoodCapitalInt =
// //                Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(theta_C_ResFC_GoodCapitalInt,
// //                CFDeltaVec_GoodCapitalInt.EntryCost(k), p_K_C_ini, GoodState, LaborIntensive, threadsManagement);
// //            ParaEst para_est_C_GoodCapitalInt = get<0>(t_GoodCapitalInt);
// //            ParaVec para_vec_C_GoodCapitalInt = get<1>(t_GoodCapitalInt);
// //            EquStateV EquV_C_GoodCapitalInt = get<2>(t_GoodCapitalInt);
// //            EquStateVmat EValmat_C_GoodCapitalInt = get<3>(t_GoodCapitalInt);
// //            EquStateV0 EquV0_C_GoodCapitalInt = get<4>(t_GoodCapitalInt);
// //            EquStateV0mat EVal0mat_C_GoodCapitalInt = get<5>(t_GoodCapitalInt);
// //
// //            /** simulate the data under the counterfactual data **/
// //            sim_data_C_GoodCapitalInt = Simulation_StatusQuo_Counterfactual(para_est_C_GoodCapitalInt,
// //                para_vec_C_GoodCapitalInt, EquV_C_GoodCapitalInt, EValmat_C_GoodCapitalInt,
// //                EquV0_C_GoodCapitalInt, EVal0mat_C_GoodCapitalInt, sim_var, SimN, SimTbar,
// //                GoodState, LaborIntensive, threadsManagement);
// //
// //            LaborIntensive = 0; GoodState = 0;
// //            cout << "LaborIntensive = " << LaborIntensive << "; GoodState = " << GoodState << endl;
// //            tuple<ParaEst, ParaVec, EquStateV, EquStateVmat, EquStateV0, EquStateV0mat> t_BadCapitalInt =
// //                Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(theta_C_ResFC_BadCapitalInt,
// //                CFDeltaVec_BadCapitalInt.EntryCost(k), p_K_C_ini, GoodState, LaborIntensive, threadsManagement);
// //            ParaEst para_est_C_BadCapitalInt = get<0>(t_BadCapitalInt);
// //            ParaVec para_vec_C_BadCapitalInt = get<1>(t_BadCapitalInt);
// //            EquStateV EquV_C_BadCapitalInt = get<2>(t_BadCapitalInt);
// //            EquStateVmat EValmat_C_BadCapitalInt = get<3>(t_BadCapitalInt);
// //            EquStateV0 EquV0_C_BadCapitalInt = get<4>(t_BadCapitalInt);
// //            EquStateV0mat EVal0mat_C_BadCapitalInt = get<5>(t_BadCapitalInt);
// //
// //            /** simulate the data under the counterfactual data **/
// //            sim_data_C_BadCapitalInt = Simulation_StatusQuo_Counterfactual(para_est_C_BadCapitalInt,
// //                para_vec_C_BadCapitalInt, EquV_C_BadCapitalInt, EValmat_C_BadCapitalInt, EquV0_C_BadCapitalInt,
// //                EVal0mat_C_BadCapitalInt, sim_var, SimN, SimTbar, GoodState, LaborIntensive,
// //                threadsManagement);
// //
// //            tuple<double, double, double, double> t_FirmMass = CalFirmMass_Given_sim_data(para_est_C_GoodLaborInt,
// //                para_est_C_BadLaborInt, para_est_C_GoodCapitalInt, para_est_C_BadCapitalInt,
// //                Value_StatusQuo_GoodLaborInt, Value_StatusQuo_BadLaborInt, Value_StatusQuo_GoodCapitalInt, Value_StatusQuo_BadCapitalInt,
// //                sim_data_C_GoodLaborInt, sim_data_C_BadLaborInt, sim_data_C_GoodCapitalInt, sim_data_C_BadCapitalInt);
// //            Equ0_GoodLaborInt_C.FirmMass = get<0>(t_FirmMass);
// //            Equ0_BadLaborInt_C.FirmMass = get<1>(t_FirmMass);
// //            Equ0_GoodCapitalInt_C.FirmMass = get<2>(t_FirmMass);
// //            Equ0_BadCapitalInt_C.FirmMass = get<3>(t_FirmMass);
// //
// //            cout << "Equ0_GoodLaborInt_C.FirmMass = " << Equ0_GoodLaborInt_C.FirmMass
// //                 << "; Equ0_BadLaborInt_C.FirmMass = " << Equ0_BadLaborInt_C.FirmMass
// //                 << "; Equ0_GoodCapitalInt_C.FirmMass = " << Equ0_GoodCapitalInt_C.FirmMass
// //                 << "; Equ0_BadCapitalInt_C.FirmMass = " << Equ0_BadCapitalInt_C.FirmMass << endl;
// //            // throw runtime_error("788");
// //
// //            if (abs(p_K_C_high - p_K_C_low) < 1e-10) {
// //                break;
// //            }
// ////
// //            double CapitalSupply = CapitalStock * pow(p_K_C_ini, KsupplyElas);
// //            CapitalDemand_GoodLaborInt_C = Equ0_GoodLaborInt_C.FirmMass * (sim_data_C_GoodLaborInt.Capital *
// //                sim_data_C_GoodLaborInt.State_S.cast<double>()).sum();
// //            CapitalDemand_BadLaborInt_C = Equ0_BadLaborInt_C.FirmMass * (sim_data_C_BadLaborInt.Capital *
// //                sim_data_C_BadLaborInt.State_S.cast<double>()).sum();
// //            CapitalDemand_GoodCapitalInt_C = Equ0_GoodCapitalInt_C.FirmMass * (sim_data_C_GoodCapitalInt.Capital *
// //                sim_data_C_GoodCapitalInt.State_S.cast<double>()).sum();
// //            CapitalDemand_BadCapitalInt_C = Equ0_BadCapitalInt_C.FirmMass * (sim_data_C_BadCapitalInt.Capital *
// //                sim_data_C_BadCapitalInt.State_S.cast<double>()).sum();
// //
// //            double CapitalDemand = shr_RealData.shr_GoodState_LaborInt * CapitalDemand_GoodLaborInt_C
// //                + shr_RealData.shr_BadState_LaborInt * CapitalDemand_BadLaborInt_C
// //                + shr_RealData.shr_GoodState_CapitalInt * CapitalDemand_GoodCapitalInt_C
// //                + shr_RealData.shr_BadState_CapitalInt * CapitalDemand_BadCapitalInt_C;
// //
// //            p_K_C_up = pow(CapitalDemand / CapitalSupply, 0.5) * p_K_C_ini;
// //
// //            double eps = abs(p_K_C_ini - p_K_C_up);
// //            cout << "eps = " << eps << "; CapitalDemand = " << CapitalDemand << "; CapitalSupply = " << CapitalSupply
// //                 << "; p_K_C_ini = " << p_K_C_ini << "; p_K_C_up = " << p_K_C_up << endl;
// //
// //            CounterFactVariable Value_C_GoodLaborInt = CalEconomicVariable(p_K_C_ini,
// //                Equ0_GoodLaborInt_C.FirmMass,Equ0_GoodLaborInt_C.F_Entry,CapitalDemand_GoodLaborInt_C,
// //                sim_data_C_GoodLaborInt);
// //            CounterFactVariable tDelta_GoodLaborInt = CalEconomicVariableDelta(Value_C_GoodLaborInt,
// //                Value_StatusQuo_GoodLaborInt);
// //            cout << "tDelta_GoodLaborInt.TotOutput_Delta = " << tDelta_GoodLaborInt.TotOutput_Delta << endl;
// //            cout << "tDelta_GoodLaborInt.TotEmploy_Delta = " << tDelta_GoodLaborInt.TotEmploy_Delta << endl;
// //            cout << "tDelta_GoodLaborInt.K_L_ratio_Delta = " << tDelta_GoodLaborInt.K_L_ratio_Delta << endl;
// //
// //            CounterFactVariable Value_C_BadLaborInt = CalEconomicVariable(p_K_C_ini,
// //                Equ0_BadLaborInt_C.FirmMass,Equ0_BadLaborInt_C.F_Entry,CapitalDemand_BadLaborInt_C,
// //                sim_data_C_BadLaborInt);
// //            CounterFactVariable tDelta_BadLaborInt = CalEconomicVariableDelta(Value_C_BadLaborInt,
// //                Value_StatusQuo_BadLaborInt);
// //            cout << "tDelta_BadLaborInt.TotOutput_Delta = " << tDelta_BadLaborInt.TotOutput_Delta << endl;
// //            cout << "tDelta_BadLaborInt.TotEmploy_Delta = " << tDelta_BadLaborInt.TotEmploy_Delta << endl;
// //            cout << "tDelta_BadLaborInt.K_L_ratio_Delta = " << tDelta_BadLaborInt.K_L_ratio_Delta << endl;
// //
// //            CounterFactVariable Value_C_GoodCapitalInt = CalEconomicVariable(p_K_C_ini,
// //                Equ0_GoodCapitalInt_C.FirmMass,Equ0_GoodCapitalInt_C.F_Entry,CapitalDemand_GoodCapitalInt_C,
// //                sim_data_C_GoodCapitalInt);
// //            CounterFactVariable tDelta_GoodCapitalInt = CalEconomicVariableDelta(Value_C_GoodCapitalInt,
// //                Value_StatusQuo_GoodCapitalInt);
// //            cout << "tDelta_GoodCapitalInt.TotOutput_Delta = " << tDelta_GoodCapitalInt.TotOutput_Delta << endl;
// //            cout << "tDelta_GoodCapitalInt.TotEmploy_Delta = " << tDelta_GoodCapitalInt.TotEmploy_Delta << endl;
// //            cout << "tDelta_GoodCapitalInt.K_L_ratio_Delta = " << tDelta_GoodCapitalInt.K_L_ratio_Delta << endl;
// //
// //            CounterFactVariable Value_C_BadCapitalInt = CalEconomicVariable(p_K_C_ini,
// //                Equ0_BadCapitalInt_C.FirmMass,Equ0_BadCapitalInt_C.F_Entry,CapitalDemand_BadCapitalInt_C,
// //                sim_data_C_BadCapitalInt);
// //            CounterFactVariable tDelta_BadCapitalInt = CalEconomicVariableDelta(Value_C_BadCapitalInt,
// //                Value_StatusQuo_BadCapitalInt);
// //            cout << "tDelta_BadCapitalInt.TotOutput_Delta = " << tDelta_BadCapitalInt.TotOutput_Delta << endl;
// //            cout << "tDelta_BadCapitalInt.TotEmploy_Delta = " << tDelta_BadCapitalInt.TotEmploy_Delta << endl;
// //            cout << "tDelta_BadCapitalInt.K_L_ratio_Delta = " << tDelta_BadCapitalInt.K_L_ratio_Delta << endl;
// //
// //            if (eps < 1e-3) {
// //                break;
// //            }
// //            p_K_C_ini = 0.5 * p_K_C_ini + 0.5 * p_K_C_up;
// //        }
// //
// //        CounterFactVariable Value_C_GoodLaborInt = CalEconomicVariable(p_K_C_ini,
// //            Equ0_GoodLaborInt_C.FirmMass,Equ0_GoodLaborInt_C.F_Entry,CapitalDemand_GoodLaborInt_C,
// //            sim_data_C_GoodLaborInt);
// //        CounterFactVariable Value_C_BadLaborInt = CalEconomicVariable(p_K_C_ini,
// //            Equ0_BadLaborInt_C.FirmMass,Equ0_BadLaborInt_C.F_Entry,CapitalDemand_BadLaborInt_C,
// //            sim_data_C_BadLaborInt);
// //        CounterFactVariable Value_C_GoodCapitalInt = CalEconomicVariable(p_K_C_ini,
// //            Equ0_GoodCapitalInt_C.FirmMass,Equ0_GoodCapitalInt_C.F_Entry,CapitalDemand_GoodCapitalInt_C,
// //            sim_data_C_GoodCapitalInt);
// //        CounterFactVariable Value_C_BadCapitalInt = CalEconomicVariable(p_K_C_ini,
// //            Equ0_BadCapitalInt_C.FirmMass,Equ0_BadCapitalInt_C.F_Entry,CapitalDemand_BadCapitalInt_C,
// //            sim_data_C_BadCapitalInt);
// //
// //        CounterFactVariable tDelta_GoodLaborInt = CalEconomicVariableDelta(Value_C_GoodLaborInt,
// //            Value_StatusQuo_GoodLaborInt);
// //        CFDeltaVec_GoodLaborInt.p_K_Delta(k) = tDelta_GoodLaborInt.p_K_Delta;
// //        CFDeltaVec_GoodLaborInt.FirmMass_Delta(k) = tDelta_GoodLaborInt.FirmMass_Delta;
// //        CFDeltaVec_GoodLaborInt.CapitalDemand_Delta(k) = tDelta_GoodLaborInt.CapitalDemand_Delta;
// //        CFDeltaVec_GoodLaborInt.NumFirm_Delta(k) = tDelta_GoodLaborInt.NumFirm_Delta;
// //        CFDeltaVec_GoodLaborInt.AvgProd_Delta(k) = tDelta_GoodLaborInt.AvgProd_Delta;
// //        CFDeltaVec_GoodLaborInt.WAvgProd_Delta(k) = tDelta_GoodLaborInt.WAvgProd_Delta;
// //        CFDeltaVec_GoodLaborInt.EntrantProd_Delta(k) = tDelta_GoodLaborInt.EntrantProd_Delta;
// //        CFDeltaVec_GoodLaborInt.ExitProd_Delta(k) = tDelta_GoodLaborInt.ExitProd_Delta;
// //        CFDeltaVec_GoodLaborInt.TotEmploy_ur_Delta(k) = tDelta_GoodLaborInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_GoodLaborInt.TotEmploy_uc_Delta(k) = tDelta_GoodLaborInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_GoodLaborInt.TotEmploy_Delta(k) = tDelta_GoodLaborInt.TotEmploy_Delta;
// //        CFDeltaVec_GoodLaborInt.TotOutput_Delta(k) = tDelta_GoodLaborInt.TotOutput_Delta;
// //        CFDeltaVec_GoodLaborInt.Output_per_Capital_Delta(k) = tDelta_GoodLaborInt.Output_per_Capital_Delta;
// //        CFDeltaVec_GoodLaborInt.Output_per_Employ_Delta(k) = tDelta_GoodLaborInt.Output_per_Employ_Delta;
// //        CFDeltaVec_GoodLaborInt.DormLength(k) = tDelta_GoodLaborInt.DormLength;
// //        CFDeltaVec_GoodLaborInt.ProdLength(k) = tDelta_GoodLaborInt.ProdLength;
// //        CFDeltaVec_GoodLaborInt.SurvivalLength(k) = tDelta_GoodLaborInt.SurvivalLength;
// //        CFDeltaVec_GoodLaborInt.PriceIndex_Delta(k) = tDelta_GoodLaborInt.PriceIndex_Delta;
// //        CFDeltaVec_GoodLaborInt.AvgExitRate(k) = tDelta_GoodLaborInt.AvgExitRate;
// //        CFDeltaVec_GoodLaborInt.K_L_ratio_Delta(k) = tDelta_GoodLaborInt.K_L_ratio_Delta;
// //        CFDeltaVec_GoodLaborInt.AggValueAdded(k) = tDelta_GoodLaborInt.AggValueAdded;
// //        CFDeltaVec_GoodLaborInt.Welfare_Manu_Delta(k) = tDelta_GoodLaborInt.Welfare_Manu_Delta;
// //        cout << "tDelta_GoodLaborInt.TotOutput_Delta = " << tDelta_GoodLaborInt.TotOutput_Delta << endl;
// //        cout << "tDelta_GoodLaborInt.TotEmploy_Delta = " << tDelta_GoodLaborInt.TotEmploy_Delta << endl;
// //        cout << "tDelta_GoodLaborInt.K_L_ratio_Delta = " << tDelta_GoodLaborInt.K_L_ratio_Delta << endl;
// //
// //        CounterFactVariable tDelta_BadLaborInt = CalEconomicVariableDelta(Value_C_BadLaborInt,
// //            Value_StatusQuo_BadLaborInt);
// //        CFDeltaVec_BadLaborInt.p_K_Delta(k) = tDelta_BadLaborInt.p_K_Delta;
// //        CFDeltaVec_BadLaborInt.FirmMass_Delta(k) = tDelta_BadLaborInt.FirmMass_Delta;
// //        CFDeltaVec_BadLaborInt.CapitalDemand_Delta(k) = tDelta_BadLaborInt.CapitalDemand_Delta;
// //        CFDeltaVec_BadLaborInt.NumFirm_Delta(k) = tDelta_BadLaborInt.NumFirm_Delta;
// //        CFDeltaVec_BadLaborInt.AvgProd_Delta(k) = tDelta_BadLaborInt.AvgProd_Delta;
// //        CFDeltaVec_BadLaborInt.WAvgProd_Delta(k) = tDelta_BadLaborInt.WAvgProd_Delta;
// //        CFDeltaVec_BadLaborInt.EntrantProd_Delta(k) = tDelta_BadLaborInt.EntrantProd_Delta;
// //        CFDeltaVec_BadLaborInt.ExitProd_Delta(k) = tDelta_BadLaborInt.ExitProd_Delta;
// //        CFDeltaVec_BadLaborInt.TotEmploy_ur_Delta(k) = tDelta_BadLaborInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_BadLaborInt.TotEmploy_uc_Delta(k) = tDelta_BadLaborInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_BadLaborInt.TotEmploy_Delta(k) = tDelta_BadLaborInt.TotEmploy_Delta;
// //        CFDeltaVec_BadLaborInt.TotOutput_Delta(k) = tDelta_BadLaborInt.TotOutput_Delta;
// //        CFDeltaVec_BadLaborInt.Output_per_Capital_Delta(k) = tDelta_BadLaborInt.Output_per_Capital_Delta;
// //        CFDeltaVec_BadLaborInt.Output_per_Employ_Delta(k) = tDelta_BadLaborInt.Output_per_Employ_Delta;
// //        CFDeltaVec_BadLaborInt.DormLength(k) = tDelta_BadLaborInt.DormLength;
// //        CFDeltaVec_BadLaborInt.ProdLength(k) = tDelta_BadLaborInt.ProdLength;
// //        CFDeltaVec_BadLaborInt.SurvivalLength(k) = tDelta_BadLaborInt.SurvivalLength;
// //        CFDeltaVec_BadLaborInt.PriceIndex_Delta(k) = tDelta_BadLaborInt.PriceIndex_Delta;
// //        CFDeltaVec_BadLaborInt.AvgExitRate(k) = tDelta_BadLaborInt.AvgExitRate;
// //        CFDeltaVec_BadLaborInt.K_L_ratio_Delta(k) = tDelta_BadLaborInt.K_L_ratio_Delta;
// //        CFDeltaVec_BadLaborInt.AggValueAdded(k) = tDelta_BadLaborInt.AggValueAdded;
// //        CFDeltaVec_BadLaborInt.Welfare_Manu_Delta(k) = tDelta_BadLaborInt.Welfare_Manu_Delta;
// //        cout << "tDelta_BadLaborInt.TotOutput_Delta = " << tDelta_BadLaborInt.TotOutput_Delta << endl;
// //        cout << "tDelta_BadLaborInt.TotEmploy_Delta = " << tDelta_BadLaborInt.TotEmploy_Delta << endl;
// //        cout << "tDelta_BadLaborInt.K_L_ratio_Delta = " << tDelta_BadLaborInt.K_L_ratio_Delta << endl;
// //
// //        CounterFactVariable tDelta_GoodCapitalInt = CalEconomicVariableDelta(Value_C_GoodCapitalInt,
// //            Value_StatusQuo_GoodCapitalInt);
// //        CFDeltaVec_GoodCapitalInt.p_K_Delta(k) = tDelta_GoodCapitalInt.p_K_Delta;
// //        CFDeltaVec_GoodCapitalInt.FirmMass_Delta(k) = tDelta_GoodCapitalInt.FirmMass_Delta;
// //        CFDeltaVec_GoodCapitalInt.CapitalDemand_Delta(k) = tDelta_GoodCapitalInt.CapitalDemand_Delta;
// //        CFDeltaVec_GoodCapitalInt.NumFirm_Delta(k) = tDelta_GoodCapitalInt.NumFirm_Delta;
// //        CFDeltaVec_GoodCapitalInt.AvgProd_Delta(k) = tDelta_GoodCapitalInt.AvgProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.WAvgProd_Delta(k) = tDelta_GoodCapitalInt.WAvgProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.EntrantProd_Delta(k) = tDelta_GoodCapitalInt.EntrantProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.ExitProd_Delta(k) = tDelta_GoodCapitalInt.ExitProd_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotEmploy_ur_Delta(k) = tDelta_GoodCapitalInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotEmploy_uc_Delta(k) = tDelta_GoodCapitalInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotEmploy_Delta(k) = tDelta_GoodCapitalInt.TotEmploy_Delta;
// //        CFDeltaVec_GoodCapitalInt.TotOutput_Delta(k) = tDelta_GoodCapitalInt.TotOutput_Delta;
// //        CFDeltaVec_GoodCapitalInt.Output_per_Capital_Delta(k) = tDelta_GoodCapitalInt.Output_per_Capital_Delta;
// //        CFDeltaVec_GoodCapitalInt.Output_per_Employ_Delta(k) = tDelta_GoodCapitalInt.Output_per_Employ_Delta;
// //        CFDeltaVec_GoodCapitalInt.DormLength(k) = tDelta_GoodCapitalInt.DormLength;
// //        CFDeltaVec_GoodCapitalInt.ProdLength(k) = tDelta_GoodCapitalInt.ProdLength;
// //        CFDeltaVec_GoodCapitalInt.SurvivalLength(k) = tDelta_GoodCapitalInt.SurvivalLength;
// //        CFDeltaVec_GoodCapitalInt.PriceIndex_Delta(k) = tDelta_GoodCapitalInt.PriceIndex_Delta;
// //        CFDeltaVec_GoodCapitalInt.AvgExitRate(k) = tDelta_GoodCapitalInt.AvgExitRate;
// //        CFDeltaVec_GoodCapitalInt.K_L_ratio_Delta(k) = tDelta_GoodCapitalInt.K_L_ratio_Delta;
// //        CFDeltaVec_GoodCapitalInt.AggValueAdded(k) = tDelta_GoodCapitalInt.AggValueAdded;
// //        CFDeltaVec_GoodCapitalInt.Welfare_Manu_Delta(k) = tDelta_GoodCapitalInt.Welfare_Manu_Delta;
// //        cout << "tDelta_GoodCapitalInt.TotOutput_Delta = " << tDelta_GoodCapitalInt.TotOutput_Delta << endl;
// //        cout << "tDelta_GoodCapitalInt.TotEmploy_Delta = " << tDelta_GoodCapitalInt.TotEmploy_Delta << endl;
// //        cout << "tDelta_GoodCapitalInt.K_L_ratio_Delta = " << tDelta_GoodCapitalInt.K_L_ratio_Delta << endl;
// //
// //        CounterFactVariable tDelta_BadCapitalInt = CalEconomicVariableDelta(Value_C_BadCapitalInt,
// //            Value_StatusQuo_BadCapitalInt);
// //        CFDeltaVec_BadCapitalInt.p_K_Delta(k) = tDelta_BadCapitalInt.p_K_Delta;
// //        CFDeltaVec_BadCapitalInt.FirmMass_Delta(k) = tDelta_BadCapitalInt.FirmMass_Delta;
// //        CFDeltaVec_BadCapitalInt.CapitalDemand_Delta(k) = tDelta_BadCapitalInt.CapitalDemand_Delta;
// //        CFDeltaVec_BadCapitalInt.NumFirm_Delta(k) = tDelta_BadCapitalInt.NumFirm_Delta;
// //        CFDeltaVec_BadCapitalInt.AvgProd_Delta(k) = tDelta_BadCapitalInt.AvgProd_Delta;
// //        CFDeltaVec_BadCapitalInt.WAvgProd_Delta(k) = tDelta_BadCapitalInt.WAvgProd_Delta;
// //        CFDeltaVec_BadCapitalInt.EntrantProd_Delta(k) = tDelta_BadCapitalInt.EntrantProd_Delta;
// //        CFDeltaVec_BadCapitalInt.ExitProd_Delta(k) = tDelta_BadCapitalInt.ExitProd_Delta;
// //        CFDeltaVec_BadCapitalInt.TotEmploy_ur_Delta(k) = tDelta_BadCapitalInt.TotEmploy_ur_Delta;
// //        CFDeltaVec_BadCapitalInt.TotEmploy_uc_Delta(k) = tDelta_BadCapitalInt.TotEmploy_uc_Delta;
// //        CFDeltaVec_BadCapitalInt.TotEmploy_Delta(k) = tDelta_BadCapitalInt.TotEmploy_Delta;
// //        CFDeltaVec_BadCapitalInt.TotOutput_Delta(k) = tDelta_BadCapitalInt.TotOutput_Delta;
// //        CFDeltaVec_BadCapitalInt.Output_per_Capital_Delta(k) = tDelta_BadCapitalInt.Output_per_Capital_Delta;
// //        CFDeltaVec_BadCapitalInt.Output_per_Employ_Delta(k) = tDelta_BadCapitalInt.Output_per_Employ_Delta;
// //        CFDeltaVec_BadCapitalInt.DormLength(k) = tDelta_BadCapitalInt.DormLength;
// //        CFDeltaVec_BadCapitalInt.ProdLength(k) = tDelta_BadCapitalInt.ProdLength;
// //        CFDeltaVec_BadCapitalInt.SurvivalLength(k) = tDelta_BadCapitalInt.SurvivalLength;
// //        CFDeltaVec_BadCapitalInt.PriceIndex_Delta(k) = tDelta_BadCapitalInt.PriceIndex_Delta;
// //        CFDeltaVec_BadCapitalInt.AvgExitRate(k) = tDelta_BadCapitalInt.AvgExitRate;
// //        CFDeltaVec_BadCapitalInt.K_L_ratio_Delta(k) = tDelta_BadCapitalInt.K_L_ratio_Delta;
// //        CFDeltaVec_BadCapitalInt.AggValueAdded(k) = tDelta_BadCapitalInt.AggValueAdded;
// //        CFDeltaVec_BadCapitalInt.Welfare_Manu_Delta(k) = tDelta_BadCapitalInt.Welfare_Manu_Delta;
// //        cout << "tDelta_BadCapitalInt.TotOutput_Delta = " << tDelta_BadCapitalInt.TotOutput_Delta << endl;
// //        cout << "tDelta_BadCapitalInt.TotEmploy_Delta = " << tDelta_BadCapitalInt.TotEmploy_Delta << endl;
// //        cout << "tDelta_BadCapitalInt.K_L_ratio_Delta = " << tDelta_BadCapitalInt.K_L_ratio_Delta << endl;
// //
// //        // throw runtime_error("722");
// //    }
// ////
// //    CFDeltaVec_GoodLaborInt.ExitSubsidy = ExitSubsidy.col(0);
// //    CFDeltaVec_BadLaborInt.ExitSubsidy = ExitSubsidy.col(1);
// //    CFDeltaVec_GoodCapitalInt.ExitSubsidy = ExitSubsidy.col(2);
// //    CFDeltaVec_BadCapitalInt.ExitSubsidy = ExitSubsidy.col(3);
// //
// //    ExitSubsidy.col(4) = shr_RealData.shr_GoodState_LaborInt * ExitSubsidy.col(0)
// //        + shr_RealData.shr_BadState_LaborInt * ExitSubsidy.col(1)
// //        + shr_RealData.shr_GoodState_CapitalInt * ExitSubsidy.col(2)
// //        + shr_RealData.shr_BadState_CapitalInt * ExitSubsidy.col(3);
// //    NumEntrants.col(4) = shr_RealData.shr_GoodState_LaborInt * NumEntrants.col(0)
// //        + shr_RealData.shr_BadState_LaborInt * NumEntrants.col(1)
// //        + shr_RealData.shr_GoodState_CapitalInt * NumEntrants.col(2)
// //        + shr_RealData.shr_BadState_CapitalInt * NumEntrants.col(3);
// //
// //    string filename = FKL + "_NumState" + to_string(1) + "Sec" + to_string(1);
// //    int t_write_GoodLaborInt = writeToCSVfileValueDelta(filename,CFDeltaVec_GoodLaborInt);
// //
// //    filename = FKL + "_NumState" + to_string(0) + "Sec" + to_string(1);
// //    int t_write_BadLaborInt = writeToCSVfileValueDelta(filename,CFDeltaVec_BadLaborInt);
// //
// //    filename = FKL + "_NumState" + to_string(1) + "Sec" + to_string(0);
// //    int t_write_GoodCapitalInt = writeToCSVfileValueDelta(filename,CFDeltaVec_GoodCapitalInt);
// //
// //    filename = FKL + "_NumState" + to_string(0) + "Sec" + to_string(0);
// //    int t_write_BadCapitalInt = writeToCSVfileValueDelta(filename,CFDeltaVec_BadCapitalInt);
// //    // throw runtime_error("722");
// //    return tuple<ArrayXXd,ArrayXXd>(ExitSubsidy,NumEntrants);
// //}
// //
// //tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat> alias::Solve_Counterfactural_FEntry_V0_Solve_PriceIndex_SingleSegment(
// //    const ArrayXd &theta_C, const double & F_Entry, const double & p_K, const int & GoodState, const int & LaborIntensive,
// //    MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ParaEst para_est_C;
// //    ParaVec para_vec_C;
// //    EquStateV EquV_C;
// //    EquStateVmat EValmat_C;
// //    EquStateV0 EquV0_C;
// //    EquStateV0mat EVal0mat_C;
// //
// //    double PI_ini = 1.0;
// //    double PI_up = 1.0;
// //    for (size_t n = 0; n < 1000; ++n) {
// //
// //        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,double> t = Solve_GECounterfactural_FEntry_V0_Cal_VEntry(
// //            theta_C, PI_ini, F_Entry, p_K, GoodState, LaborIntensive, threadsManagement);
// //        para_est_C = get<0>(t);
// //        para_vec_C = get<1>(t);
// //        EquV_C = get<2>(t);
// //        EValmat_C = get<3>(t);
// //        EquV0_C = get<4>(t);
// //        EVal0mat_C = get<5>(t);
// //        double V_Entry_C = get<6>(t);
// //
// //        PI_up = pow( (F_Entry / V_Entry_C),0.5 ) * PI_ini;
// ////        cout << "PI_ini = " << PI_ini << "; PI_up = " << PI_up << "; abs(PI_up - PI_ini) = " << abs(PI_up - PI_ini)
// ////             << "; V_Entry_C = " << V_Entry_C << "; F_Entry = " << F_Entry << endl;
// //
// //        PI_ini = PI_up*0.5+PI_ini*0.5;
// //        //// update p_K
// //        if (abs(PI_ini - PI_up) < 1e-5) {
// //            break;
// //        }
// //    }
// //
// //    para_est_C.PI = PI_ini;
// //
// //    return tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat>(para_est_C,para_vec_C,
// //        EquV_C,EValmat_C,EquV0_C,EVal0mat_C);
// //}
// //
// //tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,double> alias::Solve_GECounterfactural_FEntry_V0_Cal_VEntry(
// //        const ArrayXd &theta_C, const double & PI, const double & F_Entry, const double & p_K,
// //        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {
// //
// //    ParaEst para_est_C = constructParaEst(theta_C, GoodState, LaborIntensive);
// //    para_est_C.PI = PI;
// //    ParaVec para_vec_C = constructParaFull(para_est_C, p_K);
// //
// //    //// solve equilibrium under the counterfactual parameters and a guessed price p_K_mid
// //    //    Rev RevOpt = RevOpt_Prod(para_est, para_vec);
// //    //    ArrayXd VEnd = calVEnd(para_est, para_vec, RevOpt);
// //    ArrayXd VEnd = ArrayXd::Zero(para.N*2);
// //    //// solve value function for t > 1
// //    tuple<EquStateV,EquStateVmat> t_Equ = solveV(para_est_C,para_vec_C,VEnd,threadsManagement);
// //    EquStateV EquV_C = get<0>(t_Equ);
// //    EquStateVmat EValmat_C = get<1>(t_Equ);
// //
// //    //// With free entry condition, solve value/policy function at the initial t = 1
// //    //// not done yet.
// //    tuple<EquStateV0,EquStateV0mat> t_Equ0 = solveV0(para_est_C,para_vec_C,EquV_C.EVal_PD,
// //        threadsManagement);
// //    EquStateV0 EquV0_C = get<0>(t_Equ0);
// //    EquStateV0mat EVal0mat_C = get<1>(t_Equ0);
// //
// //    //// calculate the expected value of entry under the counterfactual parameters and a guessed price p_K_mid
// //    double V_Entry_C = 0.0;
// //    for (size_t k = 0; k < para.N_phi; ++k) {
// //        V_Entry_C = V_Entry_C + para_vec_C.dist0(k) * EquV0_C.EVal_PD0(k * para.N_KL);
// //    }
// //
// //    return tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,double>(
// //            para_est_C,para_vec_C,EquV_C,EValmat_C,EquV0_C,EVal0mat_C,V_Entry_C);
// //}
// //
// //tuple<double,double,double,double> alias::CalFirmMass_Given_sim_data(const ParaEst &para_est_C_GoodLaborInt,
// //    const ParaEst &para_est_C_BadLaborInt, const ParaEst &para_est_C_GoodCapitalInt, const ParaEst &para_est_C_BadCapitalInt,
// //    const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
// //    const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
// //    SimData & sim_data_C_GoodLaborInt,SimData & sim_data_C_BadLaborInt,
// //    SimData & sim_data_C_GoodCapitalInt,SimData & sim_data_C_BadCapitalInt) {
// //
// //
// //    double powratio = (1.0 - para.alpha_M_G_L * (para.sigma-1)/para.sigma) / (1.0/para.sigma);
// //    double PIsigma_GoodLaborInt = pow(para_est_C_GoodLaborInt.PI,powratio);
// //    powratio = (1.0 - para.alpha_M_B_L * (para.sigma-1)/para.sigma) / (1.0/para.sigma);
// //    double PIsigma_BadLaborInt = pow(para_est_C_BadLaborInt.PI,powratio);
// //    powratio = (1.0 - para.alpha_M_G_K * (para.sigma-1)/para.sigma) / (1.0/para.sigma);
// //    double PIsigma_GoodCapitalInt = pow(para_est_C_GoodCapitalInt.PI,powratio);
// //    powratio = (1.0 - para.alpha_M_B_K * (para.sigma-1)/para.sigma) / (1.0/para.sigma);
// //    double PIsigma_BadCapitalInt = pow(para_est_C_BadCapitalInt.PI,powratio);
// //    cout << "PI = " << para_est_C_GoodLaborInt.PI << "; " << para_est_C_BadLaborInt.PI << "; "
// //        << para_est_C_GoodCapitalInt.PI << "; " << para_est_C_BadCapitalInt.PI << endl;
// //    cout << "PIsigma = " << PIsigma_GoodLaborInt << "; " << PIsigma_BadLaborInt << "; "
// //        << PIsigma_GoodCapitalInt << "; " << PIsigma_BadCapitalInt << endl;
// //
// //    double EP_S_GoodLaborInt = Value_StatusQuo_GoodLaborInt.TotOutput / pow(Value_StatusQuo_GoodLaborInt.PriceIndex,1.0-para.sigma);
// //    double EP_S_BadLaborInt = Value_StatusQuo_BadLaborInt.TotOutput / pow(Value_StatusQuo_BadLaborInt.PriceIndex,1.0-para.sigma);
// //
// //    double lambdaL_GG = EP_S_GoodLaborInt / (EP_S_GoodLaborInt + pow(para.tau_GB,1.0-para.sigma) * EP_S_BadLaborInt);
// //    double lambdaL_GB = pow(para.tau_GB,1.0-para.sigma) * EP_S_BadLaborInt / (EP_S_GoodLaborInt + pow(para.tau_GB,1.0-para.sigma) * EP_S_BadLaborInt);
// //    double lambdaL_BG = pow(para.tau_BG,1.0-para.sigma) * EP_S_GoodLaborInt / (pow(para.tau_BG,1.0-para.sigma) * EP_S_GoodLaborInt + EP_S_BadLaborInt);
// //    double lambdaL_BB = EP_S_BadLaborInt / (pow(para.tau_BG,1.0-para.sigma) * EP_S_GoodLaborInt + EP_S_BadLaborInt);
// //
// //    double tempL = lambdaL_BB*lambdaL_GG - lambdaL_BG*lambdaL_GB;
// //
// //    double EP_GoodLaborInt = (PIsigma_GoodLaborInt*lambdaL_BB - PIsigma_BadLaborInt*lambdaL_GB)/tempL;
// //    double EP_BadLaborInt = (PIsigma_BadLaborInt*lambdaL_GG - PIsigma_GoodLaborInt*lambdaL_BG)/tempL;
// //
// //    double EP_S_GoodCapitalInt = Value_StatusQuo_GoodCapitalInt.TotOutput / pow(Value_StatusQuo_GoodCapitalInt.PriceIndex,1.0-para.sigma);
// //    double EP_S_BadCapitalInt = Value_StatusQuo_BadCapitalInt.TotOutput / pow(Value_StatusQuo_BadCapitalInt.PriceIndex,1.0-para.sigma);
// //
// //    double lambdaK_GG = EP_S_GoodCapitalInt / (EP_S_GoodCapitalInt + pow(para.tau_GB,1.0-para.sigma) * EP_S_BadCapitalInt);
// //    double lambdaK_GB = pow(para.tau_GB,1.0-para.sigma) * EP_S_BadCapitalInt / (EP_S_GoodCapitalInt + pow(para.tau_GB,1.0-para.sigma) * EP_S_BadCapitalInt);
// //    double lambdaK_BG = pow(para.tau_BG,1.0-para.sigma) * EP_S_GoodCapitalInt / (pow(para.tau_BG,1.0-para.sigma) * EP_S_GoodCapitalInt + EP_S_BadCapitalInt);
// //    double lambdaK_BB = EP_S_BadCapitalInt / (pow(para.tau_BG,1.0-para.sigma) * EP_S_GoodCapitalInt + EP_S_BadCapitalInt);
// //
// //    double tempK = lambdaK_BB*lambdaK_GG - lambdaK_BG*lambdaK_GB;
// //    double EP_GoodCapitalInt = (PIsigma_GoodCapitalInt*lambdaK_BB - PIsigma_BadCapitalInt*lambdaK_GB)/tempK;
// //    double EP_BadCapitalInt = (PIsigma_BadCapitalInt*lambdaK_GG - PIsigma_GoodCapitalInt*lambdaK_BG)/tempK;
// //
// //     EP_GoodLaborInt = PIsigma_GoodLaborInt;
// //     EP_BadLaborInt = PIsigma_BadLaborInt;
// //     EP_GoodCapitalInt = PIsigma_GoodCapitalInt;
// //     EP_BadCapitalInt = PIsigma_BadCapitalInt;
// //
// //    cout << "EP_GoodLaborInt = " << EP_GoodLaborInt << "; EP_BadLaborInt = " << EP_BadLaborInt
// //        << "; EP_GoodCapitalInt = " << EP_GoodCapitalInt << "; EP_BadCapitalInt = " << EP_BadCapitalInt << endl;
// //    double a_G = 1.0 - para.IncShrManu;
// //    double a_B = 1.0 - para.IncShrManu;
// //
// //    /**** b_lK ****/
// //    double IncManu_GoodLaborInt_C = calTotalOutput(sim_data_C_GoodLaborInt, 1.0);
// //    double IncManu_GoodCapitalInt_C = calTotalOutput(sim_data_C_GoodCapitalInt, 1.0);
// //    double IncManu_BadLaborInt_C = calTotalOutput(sim_data_C_BadLaborInt, 1.0);
// //    double IncManu_BadCapitalInt_C = calTotalOutput(sim_data_C_BadCapitalInt, 1.0);
// //
// //    double IncManu_GoodLaborInt_hat = IncManu_GoodLaborInt_C / Value_StatusQuo_GoodLaborInt.TotOutput;
// //    double IncManu_GoodCapitalInt_hat = IncManu_GoodCapitalInt_C / Value_StatusQuo_GoodCapitalInt.TotOutput;
// //    double IncManu_BadLaborInt_hat = IncManu_BadLaborInt_C / Value_StatusQuo_BadLaborInt.TotOutput;
// //    double IncManu_BadCapitalInt_hat = IncManu_BadCapitalInt_C / Value_StatusQuo_BadCapitalInt.TotOutput;
// //
// //    double IncManu_GoodLaborInt_StatusQuo_Shr = Value_StatusQuo_GoodLaborInt.TotOutput
// //        / (Value_StatusQuo_GoodLaborInt.TotOutput + Value_StatusQuo_GoodCapitalInt.TotOutput) * para.IncShrManu;
// //    double IncManu_GoodCapitalInt_StatusQuo_Shr = Value_StatusQuo_GoodCapitalInt.TotOutput
// //        / (Value_StatusQuo_GoodLaborInt.TotOutput + Value_StatusQuo_GoodCapitalInt.TotOutput) * para.IncShrManu;
// //    double IncManu_BadLaborInt_StatusQuo_Shr = Value_StatusQuo_BadLaborInt.TotOutput
// //        / (Value_StatusQuo_BadLaborInt.TotOutput + Value_StatusQuo_BadCapitalInt.TotOutput) * para.IncShrManu;
// //    double IncManu_BadCapitalInt_StatusQuo_Shr = Value_StatusQuo_BadCapitalInt.TotOutput
// //        / (Value_StatusQuo_BadLaborInt.TotOutput + Value_StatusQuo_BadCapitalInt.TotOutput) * para.IncShrManu;
// //
// //    double b_LG = IncManu_GoodLaborInt_StatusQuo_Shr * IncManu_GoodLaborInt_hat;
// //    double b_KG = IncManu_GoodCapitalInt_StatusQuo_Shr * IncManu_GoodCapitalInt_hat;
// //    double b_LB = IncManu_BadLaborInt_StatusQuo_Shr * IncManu_BadLaborInt_hat;
// //    double b_KB = IncManu_BadCapitalInt_StatusQuo_Shr * IncManu_BadCapitalInt_hat;
// //
// //    cout << "b_LG = " << b_LG << "; b_LB = " << b_LB
// //        << "; b_KG = " << b_KG << "; b_KB = " << b_KB << endl;
// //
// //    /**** c_lK ****/
// //    double PriceIndex_GoodLaborInt_StatusQuo = Value_StatusQuo_GoodLaborInt.PriceIndex;
// //    double PriceIndex_GoodCapitalInt_StatusQuo = Value_StatusQuo_GoodCapitalInt.PriceIndex;
// //    double PriceIndex_BadLaborInt_StatusQuo = Value_StatusQuo_BadLaborInt.PriceIndex;
// //    double PriceIndex_BadCapitalInt_StatusQuo = Value_StatusQuo_BadCapitalInt.PriceIndex;
// //
// //    double ImpShr_GGLaborInt = pow(PriceIndex_GoodLaborInt_StatusQuo,1.0-para.sigma)
// //        / ( pow(PriceIndex_GoodLaborInt_StatusQuo,1.0-para.sigma) + pow(para.tau_BG * PriceIndex_BadLaborInt_StatusQuo,1.0-para.sigma) );
// //    double ImpShr_BGLaborInt = pow(para.tau_BG * PriceIndex_BadLaborInt_StatusQuo,1.0-para.sigma)
// //        / ( pow(PriceIndex_GoodLaborInt_StatusQuo,1.0-para.sigma) + pow(para.tau_BG * PriceIndex_BadLaborInt_StatusQuo,1.0-para.sigma) );
// //
// //    double ImpShr_GBLaborInt = pow(para.tau_GB * PriceIndex_GoodLaborInt_StatusQuo,1.0-para.sigma)
// //        / ( pow(para.tau_GB * PriceIndex_GoodLaborInt_StatusQuo,1.0-para.sigma) + pow(PriceIndex_BadLaborInt_StatusQuo,1.0-para.sigma) );
// //    double ImpShr_BBLaborInt = pow(PriceIndex_BadLaborInt_StatusQuo,1.0-para.sigma)
// //        / ( pow(para.tau_GB * PriceIndex_GoodLaborInt_StatusQuo,1.0-para.sigma) + pow(PriceIndex_BadLaborInt_StatusQuo,1.0-para.sigma) );
// //
// //    double ImpShr_GGCapitalInt = pow(PriceIndex_GoodCapitalInt_StatusQuo,1-para.sigma)
// //        / ( pow(PriceIndex_GoodCapitalInt_StatusQuo,1-para.sigma) + pow(para.tau_BG * PriceIndex_BadCapitalInt_StatusQuo,1-para.sigma) );
// //    double ImpShr_BGCapitalInt = pow(para.tau_BG * PriceIndex_BadCapitalInt_StatusQuo,1-para.sigma)
// //        / ( pow(PriceIndex_GoodCapitalInt_StatusQuo,1-para.sigma) + pow(para.tau_BG * PriceIndex_BadCapitalInt_StatusQuo,1-para.sigma) );
// //
// //    double ImpShr_GBCapitalInt = pow(para.tau_GB * PriceIndex_GoodCapitalInt_StatusQuo,1-para.sigma)
// //        / ( pow(para.tau_GB * PriceIndex_GoodCapitalInt_StatusQuo,1-para.sigma) + pow(PriceIndex_BadCapitalInt_StatusQuo,1-para.sigma) );
// //    double ImpShr_BBCapitalInt = pow(PriceIndex_BadCapitalInt_StatusQuo,1-para.sigma)
// //        / ( pow(para.tau_GB * PriceIndex_GoodCapitalInt_StatusQuo,1-para.sigma) + pow(PriceIndex_BadCapitalInt_StatusQuo,1-para.sigma) );
// //
// //    double PriceIndex_GoodLaborInt_C = calPriceIndex(sim_data_C_GoodLaborInt, 1.0);
// //    double PriceIndex_GoodCapitalInt_C = calPriceIndex(sim_data_C_GoodCapitalInt, 1.0);
// //    double PriceIndex_BadLaborInt_C = calPriceIndex(sim_data_C_BadLaborInt, 1.0);
// //    double PriceIndex_BadCapitalInt_C = calPriceIndex(sim_data_C_BadCapitalInt, 1.0);
// //
// //    double PriceIndex_GoodLaborInt_hat = PriceIndex_GoodLaborInt_C / PriceIndex_GoodLaborInt_StatusQuo;
// //    double PriceIndex_GoodCapitalInt_hat = PriceIndex_GoodCapitalInt_C / PriceIndex_GoodCapitalInt_StatusQuo;
// //    double PriceIndex_BadLaborInt_hat = PriceIndex_BadLaborInt_C / PriceIndex_BadLaborInt_StatusQuo;
// //    double PriceIndex_BadCapitalInt_hat = PriceIndex_BadCapitalInt_C / PriceIndex_BadCapitalInt_StatusQuo;
// //
// //    double c_GGL = ImpShr_GGLaborInt * pow(PriceIndex_GoodLaborInt_hat,1.0-para.sigma) * EP_GoodLaborInt;
// //    double c_BGL = ImpShr_BGLaborInt * pow(PriceIndex_BadLaborInt_hat,1.0-para.sigma) * EP_GoodLaborInt;
// //    double c_GBL = ImpShr_GBLaborInt * pow(PriceIndex_GoodLaborInt_hat,1.0-para.sigma) * EP_BadLaborInt;
// //    double c_BBL = ImpShr_BBLaborInt * pow(PriceIndex_BadLaborInt_hat,1.0-para.sigma) * EP_BadLaborInt;
// //    cout << "c_GGL = " << c_GGL << "; c_BGL = " << c_BGL << "; c_GBL = " << c_GBL << "; c_BBL = " << c_BBL << endl;
// //
// //    double c_GGK = ImpShr_GGCapitalInt * pow(PriceIndex_GoodCapitalInt_hat,1.0-para.sigma) * EP_GoodCapitalInt;
// //    double c_BGK = ImpShr_BGCapitalInt * pow(PriceIndex_BadCapitalInt_hat,1.0-para.sigma) * EP_GoodCapitalInt;
// //    double c_GBK = ImpShr_GBCapitalInt * pow(PriceIndex_GoodCapitalInt_hat,1.0-para.sigma) * EP_BadCapitalInt;
// //    double c_BBK = ImpShr_BBCapitalInt * pow(PriceIndex_BadCapitalInt_hat,1.0-para.sigma) * EP_BadCapitalInt;
// //    cout << "c_GGK = " << c_GGK << "; c_BGK = " << c_BGK << "; c_GBK = " << c_GBK << "; c_BBK = " << c_BBK << endl;
// //
// //    double denominator = b_KB*b_KG*c_BBL*c_GGL - b_KB*b_KG*c_BGL*c_GBL + b_KB*b_LG*c_BBL*c_GGK - b_KB*b_LG*c_BGL*c_GBK + b_KG*b_LB*c_BBK*c_GGL
// //        - b_KG*b_LB*c_BGK*c_GBL + b_LB*b_LG*c_BBK*c_GGK - b_LB*b_LG*c_BGK*c_GBK - b_KG*c_BBK*c_BBL*c_GGL + b_KG*c_BBK*c_BGL*c_GBL
// //        - b_LG*c_BBK*c_BBL*c_GGK + b_LG*c_BBL*c_BGK*c_GBK - b_KB*c_BBL*c_GGK*c_GGL + b_KB*c_BGL*c_GBL*c_GGK - b_LB*c_BBK*c_GGK*c_GGL
// //        + b_LB*c_BGK*c_GBK*c_GGL + c_BBK*c_BBL*c_GGK*c_GGL - c_BBK*c_BGL*c_GBL*c_GGK - c_BBL*c_BGK*c_GBK*c_GGL + c_BGK*c_BGL*c_GBK*c_GBL;
// //    cout << "denominator = " << denominator << endl;
// //
// //    double FirmMass_GL = ( a_B*(b_KB*b_KG*c_BGL + b_KG*b_LB*c_BGK - b_KG*c_BBL*c_BGK - b_KB*c_BGL*c_GGK) ) / denominator
// //        - ( a_G*(b_KB*b_KG*c_BBL + b_KG*b_LB*c_BBK - b_KG*c_BBK*c_BBL - b_KB*c_BGL*c_GBK) ) / denominator
// //        - ( a_B*(b_KB*b_KG*c_BGL + b_KG*b_LB*c_BGK - b_KG*c_BBK*c_BGL - b_KB*c_BGL*c_GGK + c_BBK*c_BGL*c_GGK - c_BGK*c_BGL*c_GBK) ) / denominator
// //        + ( a_G*(b_KB*b_KG*c_BBL + b_KG*b_LB*c_BBK - b_KG*c_BBK*c_BBL - b_KB*c_BBL*c_GGK - b_LB*c_BBK*c_GGK + b_LB*c_BGK*c_GBK + c_BBK*c_BBL*c_GGK - c_BBL*c_BGK*c_GBK) ) / denominator;
// //
// //    double FirmMass_BL = ( a_G*(b_KB*b_KG*c_GBL + b_KB*b_LG*c_GBK - b_KG*c_BBK*c_GBL - b_KB*c_GBK*c_GGL) ) / denominator
// //        - ( a_B*(b_KB*b_KG*c_GGL + b_KB*b_LG*c_GGK - b_KG*c_BGK*c_GBL - b_KB*c_GGK*c_GGL) ) / denominator
// //        - ( a_G*(b_KB*b_KG*c_GBL + b_KB*b_LG*c_GBK - b_KG*c_BBK*c_GBL - b_KB*c_GBL*c_GGK + c_BBK*c_GBL*c_GGK - c_BGK*c_GBK*c_GBL) ) / denominator
// //        + ( a_B*(b_KB*b_KG*c_GGL + b_KB*b_LG*c_GGK - b_KG*c_BBK*c_GGL - b_LG*c_BBK*c_GGK + b_LG*c_BGK*c_GBK - b_KB*c_GGK*c_GGL + c_BBK*c_GGK*c_GGL - c_BGK*c_GBK*c_GGL) ) / denominator;
// //
// //    double FirmMass_GK = ( a_B*(b_KB*b_LG*c_BGL + b_LB*b_LG*c_BGK - b_LG*c_BBK*c_BGL - b_LB*c_BGK*c_GGL) ) / denominator
// //        - ( a_G*(b_KB*b_LG*c_BBL + b_LB*b_LG*c_BBK - b_LG*c_BBK*c_BBL - b_LB*c_BGK*c_GBL) ) / denominator
// //        - ( a_B*(b_KB*b_LG*c_BGL + b_LB*b_LG*c_BGK - b_LG*c_BBL*c_BGK - b_LB*c_BGK*c_GGL + c_BBL*c_BGK*c_GGL - c_BGK*c_BGL*c_GBL) ) / denominator
// //        + ( a_G*(b_KB*b_LG*c_BBL + b_LB*b_LG*c_BBK - b_LG*c_BBK*c_BBL - b_KB*c_BBL*c_GGL + b_KB*c_BGL*c_GBL - b_LB*c_BBK*c_GGL + c_BBK*c_BBL*c_GGL - c_BBK*c_BGL*c_GBL) ) / denominator;
// //
// //    double FirmMass_BK = ( a_G*(b_KG*b_LB*c_GBL + b_LB*b_LG*c_GBK - b_LG*c_BBL*c_GBK - b_LB*c_GBL*c_GGK) ) / denominator
// //        - ( a_B*(b_KG*b_LB*c_GGL + b_LB*b_LG*c_GGK - b_LG*c_BGL*c_GBK - b_LB*c_GGK*c_GGL) ) / denominator
// //        - ( a_G*(b_KG*b_LB*c_GBL + b_LB*b_LG*c_GBK - b_LG*c_BBL*c_GBK - b_LB*c_GBK*c_GGL + c_BBL*c_GBK*c_GGL - c_BGL*c_GBK*c_GBL) ) / denominator
// //        + ( a_B*(b_KG*b_LB*c_GGL + b_LB*b_LG*c_GGK - b_KG*c_BBL*c_GGL + b_KG*c_BGL*c_GBL - b_LG*c_BBL*c_GGK - b_LB*c_GGK*c_GGL + c_BBL*c_GGK*c_GGL - c_BGL*c_GBL*c_GGK) ) / denominator;
// //
// //    cout << "FirmMass_GL = " << FirmMass_GL << "; FirmMass_BL = " << FirmMass_BL << "; FirmMass_GK = " << FirmMass_GK
// //        << "; FirmMass_BK = " << FirmMass_BK << endl;
// //    // throw runtime_error("2345");
// //    return tuple<double,double,double,double>(FirmMass_GL,FirmMass_BL,FirmMass_GK,FirmMass_BK);
// //}
// //
// //
// //
// //
// ////////
// ////////
// //
// ////
// //
// //
// //
// /////***********************************************************************************************************
// ////*** Counterfactual Exercise 3: General Equilibrium: Varying the exit barrier to match bad / good states; ***
// ////*** Change elasticity of Ksupply elasticity
// ////************************************************************************************************************/
// ////int alias::SolveSimulate_GECounterfacturalResult_3(const ArrayXd & theta_Est_S, const SegmentShr &shr_RealData,
// ////    const CounterFactVariable &Value_StatusQuo_GoodLaborInt, const CounterFactVariable &Value_StatusQuo_BadLaborInt,
// ////    const CounterFactVariable &Value_StatusQuo_GoodCapitalInt, const CounterFactVariable &Value_StatusQuo_BadCapitalInt,
// ////    const SimVar &sim_var, const int SimTbar, const int SimN,
// ////    const int &N_C, MultiThreads::Threads_Management &threadsManagement) {
// ////
// ////    int N_KsupplyElas = 4;
// ////    ArrayXd KsupplyElas(N_KsupplyElas);
// ////    KsupplyElas << 0.0,0.1,0.2,0.5;
// ////
// ////    for(size_t k = 0; k < N_KsupplyElas; ++k) {
// ////        ArrayXXd ResidualValue = ArrayXXd::Zero(2,4);
// ////        ResidualValue.row(0) << theta_Est_S(14),theta_Est_S(14),theta_Est_S(14),theta_Est_S(14);
// ////        ResidualValue.row(1) << theta_Est_S(15),theta_Est_S(15),theta_Est_S(15),theta_Est_S(15);
// ////
// ////        ArrayXXd firingshare = ArrayXXd::Ones(2,4);
// ////
// ////        ArrayXXd EntryCost(N_C*2,4);
// ////        EntryCost.col(0) = Value_StatusQuo_GoodLaborInt.F_Entry * ArrayXd::Ones(2);
// ////        EntryCost.col(1) = Value_StatusQuo_BadLaborInt.F_Entry * ArrayXd::Ones(2);
// ////        EntryCost.col(2) = Value_StatusQuo_GoodCapitalInt.F_Entry * ArrayXd::Ones(2);
// ////        EntryCost.col(3) = Value_StatusQuo_BadCapitalInt.F_Entry * ArrayXd::Ones(2);
// ////
// ////        string FKL = "_GE_C3_ResVal_KsupplyElas" + to_string(k);
// ////        tuple<ArrayXXd,ArrayXXd> t = CalSimulation_GECounterfactualResult(theta_Est_S, shr_RealData,
// ////            Value_StatusQuo_GoodLaborInt,Value_StatusQuo_BadLaborInt,Value_StatusQuo_GoodCapitalInt,Value_StatusQuo_BadCapitalInt,
// ////            sim_var, SimTbar, SimN, KsupplyElas(k),
// ////            ResidualValue,firingshare, EntryCost, FKL, 1.0,1000.0, threadsManagement);
// ////    }
// ////
// ////    return 0;
// ////}
// ////
// ////
// ////////////
// /////////////**************************************************************
// ////////////* Counterfactual 6: With a target exit rate, Draw the PPF of exit barrier and firing costs, with employment loss as the cost curve
// ////////////**************************************************************/
// ////////////int alias::SolveSimulate_CounterfacturalResult_6(const ArrayXd & theta_Est_S, const CounterFactVariable & Val_StatusQuo,
// ////////////    const SimVar & sim_var, const int SimTbar, const int SimN, const int & GoodState, const int & LaborIntensive,
// ////////////    const IndustryPara & para_ind, const double & KsupplyElas, const string & State, const double & TargetExitRate,
// ////////////    const int & N_C, MultiThreads::Threads_Management & threadsManagement) {
// ////////////
// ////////////    ArrayXd temp1 = ArrayXd::LinSpaced(N_C, theta_Est_S(para.dim1+para.dim2+para.dim3-2)-500,
// ////////////        theta_Est_S(para.dim1+para.dim2+para.dim3-2));
// ////////////    ArrayXd temp2 = ArrayXd::LinSpaced(N_C, theta_Est_S(para.dim1+para.dim2+para.dim3-2),
// ////////////        theta_Est_S(para.dim1+para.dim2+para.dim3-2)+500);
// ////////////    ArrayXd ResidualValue(2*N_C - 1);
// ////////////    ResidualValue << temp1, temp2(seqN(1,N_C-1));
// ////////////    cout << "ResidualValue = " << ResidualValue.transpose() << endl;
// ////////////
// ////////////    cout << "Start the loop: choose different residual value" << endl;
// ////////////    string FKL = "_PPF_ExitRate";
// ////////////    ArrayXd TotEmploy_Delta(2*N_C - 1);
// ////////////    ArrayXd TotOutput_Delta(2*N_C - 1);
// ////////////    ArrayXd firing_share(2*N_C - 1);
// ////////////    ArrayXd AvgExitRate(2*N_C - 1);
// ////////////    for (size_t k = 0; k < 2*N_C - 1; ++k) {
// ////////////        ArrayXd theta_C = theta_Est_S;
// ////////////        theta_C(para.dim1+para.dim2+para.dim3-2) = ResidualValue(k);
// ////////////
// ////////////        double firing_share_up = 1; double firing_share_low = 0.01;
// ////////////        double firing_share_mid;
// ////////////
// ////////////        for (size_t n = 0; n < 10000; ++n) {
// ////////////            firing_share_mid = 0.5*firing_share_up + 0.5*firing_share_low;
// ////////////            // Capital Firing Cost
// //////////////        theta_C(para.dim1+para.dim2-1+3) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+3);
// //////////////        theta_C(para.dim1+para.dim2-1+4) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+4);
// ////////////            // Regular Workers Firing Cost
// //////////////            theta_C(para.dim1+para.dim2-1+7) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+7);
// //////////////            theta_C(para.dim1+para.dim2-1+8) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+8);
// ////////////            theta_C(para.dim1+para.dim2-1+9) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+9);
// ////////////            theta_C(para.dim1+para.dim2-1+10) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+10);
// ////////////            // Contract Workers Firing Cost
// ////////////            theta_C(para.dim1+para.dim2-1+13) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+13);
// ////////////            theta_C(para.dim1+para.dim2-1+14) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+14);
// ////////////
// ////////////            tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData> t_Equ_Counterfactual
// ////////////                = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C, Val_StatusQuo.Equ0_S.F_Entry,
// ////////////                Val_StatusQuo.CapitalStock_S,Val_StatusQuo.PriceIndex_S,
// ////////////                Val_StatusQuo.FirmMass_S, sim_var, GoodState, LaborIntensive, 0.0,
// ////////////                threadsManagement);
// ////////////            EquState0 Equ0_C = get<6>(t_Equ_Counterfactual);
// ////////////            double FirmMass_C = Equ0_C.FirmMass;
// ////////////            SimData sim_data_C = get<9>(t_Equ_Counterfactual);
// ////////////
// ////////////            /*** Calculate Average Exit Rate ***/
// ////////////            double AvgExitRate_C = calAvgExitRate(sim_data_C, FirmMass_C);
// ////////////            if (AvgExitRate_C > TargetExitRate) {
// ////////////                firing_share_low = firing_share_mid;
// ////////////            }
// ////////////            else {
// ////////////                firing_share_up = firing_share_mid;
// ////////////            }
// ////////////
// ////////////            cout << "** TargetExitRate = " << TargetExitRate << "; AvgExitRate_C = " << AvgExitRate_C
// ////////////                 << "; firing_share_up = " << firing_share_up << "; firing_share_low = " << firing_share_low << endl;
// ////////////            if (abs(firing_share_up - firing_share_low) < 1e-3) {
// ////////////                break;
// ////////////            }
// ////////////        }
// ////////////
// ////////////        firing_share(k) = firing_share_mid;
// ////////////        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData> t_Equ_Counterfactual
// ////////////            = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C, Val_StatusQuo.Equ0_S.F_Entry,
// ////////////            Val_StatusQuo.CapitalStock_S,Val_StatusQuo.PriceIndex_S,
// ////////////            Val_StatusQuo.FirmMass_S, sim_var, GoodState, LaborIntensive, KsupplyElas,
// ////////////            threadsManagement);
// ////////////        SimData sim_data_C = get<9>(t_Equ_Counterfactual);
// ////////////        EquState0 Equ0_C = get<6>(t_Equ_Counterfactual);
// ////////////        double FirmMass_C = Equ0_C.FirmMass;
// ////////////
// ////////////        /*** Calculate Average Exit Rate ***/
// ////////////        AvgExitRate(k) = calAvgExitRate(sim_data_C, FirmMass_C);
// ////////////
// ////////////        /*** Calculate Total Labor Employment ***/
// ////////////        tuple<double,double> t_TotEmp_C = calTotalEmployment(sim_data_C, FirmMass_C);
// ////////////        double TotEmploy_ur_C = get<0>(t_TotEmp_C);
// ////////////        double TotEmploy_uc_C = get<1>(t_TotEmp_C);
// ////////////        TotEmploy_Delta(k) = (TotEmploy_ur_C+TotEmploy_uc_C) / Val_StatusQuo.TotEmploy_S;
// ////////////
// ////////////        /*** Calculate Total Output ***/
// ////////////        double TotOutput_C = calTotalOutput(sim_data_C, FirmMass_C);
// ////////////        TotOutput_Delta(k) = TotOutput_C / Val_StatusQuo.TotOutput_S;
// ////////////    }
// ////////////
// ////////////    cout << "TotOutput_Delta = " << TotOutput_Delta.transpose() << endl;
// ////////////    cout << "TotEmploy_Delta = " << TotEmploy_Delta.transpose() << endl;
// ////////////    cout << "ResidualValue = " << ResidualValue.transpose() << endl;
// ////////////    cout << "firing_share = " << firing_share.transpose() << endl;
// ////////////    cout << "AvgExitRate = " << AvgExitRate.transpose() << endl;
// ////////////
// ////////////    writeToCSVfile("TotOutput_Delta_" + State + FKL + ".csv", TotOutput_Delta.cast<double>().matrix());
// ////////////    writeToCSVfile("TotEmploy_Delta_" + State + FKL + ".csv", TotEmploy_Delta.cast<double>().matrix());
// ////////////    writeToCSVfile("ResidualValue_" + State + FKL + ".csv", ResidualValue.cast<double>().matrix());
// ////////////    writeToCSVfile("firing_share_" + State + FKL + ".csv", firing_share.cast<double>().matrix());
// ////////////    writeToCSVfile("AvgExitRate_" + State + FKL + ".csv", AvgExitRate.cast<double>().matrix());
// ////////////
// ////////////    return 0;
// ////////////
// ////////////}
// ////////////
// ////////////int alias::SolveSimulate_CounterfacturalResult_7(const ArrayXd & theta_Est_S, const CounterFactVariable & Val_StatusQuo,
// ////////////    const SimVar & sim_var, const int SimTbar, const int SimN, const int & GoodState, const int & LaborIntensive,
// ////////////    const IndustryPara & para_ind, const double & KsupplyElas, const string & State, const double & TargetOutput_Delta,
// ////////////    const int & N_C, MultiThreads::Threads_Management & threadsManagement) {
// ////////////
// ////////////    double TargetOutput = TargetOutput_Delta * Val_StatusQuo.TotOutput_S;
// ////////////
// ////////////    ArrayXd temp1 = ArrayXd::LinSpaced(N_C, theta_Est_S(para.dim1+para.dim2+para.dim3-2)-200,
// ////////////        theta_Est_S(para.dim1+para.dim2+para.dim3-2));
// ////////////    ArrayXd temp2 = ArrayXd::LinSpaced(N_C, theta_Est_S(para.dim1+para.dim2+para.dim3-2),
// ////////////        theta_Est_S(para.dim1+para.dim2+para.dim3-2)+200);
// ////////////    ArrayXd ResidualValue(2*N_C - 1);
// ////////////    ResidualValue << temp1, temp2(seqN(1,N_C-1));
// ////////////    cout << "ResidualValue = " << ResidualValue.transpose() << endl;
// ////////////
// ////////////    cout << "Start the loop: choose different residual value" << endl;
// ////////////    string FKL = "_PPF_Output";
// ////////////    ArrayXd TotEmploy_Delta(2*N_C - 1);
// ////////////    ArrayXd TotOutput_Delta(2*N_C - 1);
// ////////////    ArrayXd firing_share(2*N_C - 1);
// ////////////    ArrayXd AvgExitRate(2*N_C - 1);
// ////////////    for (size_t k = 0; k < 2*N_C - 1; ++k) {
// ////////////        cout << "******* k = " << k << endl;
// ////////////        ArrayXd theta_C = theta_Est_S;
// ////////////        theta_C(para.dim1+para.dim2+para.dim3-2) = ResidualValue(k);
// ////////////        cout << "theta_C = " << theta_C.transpose() << endl;
// ////////////        double firing_share_up = 10; double firing_share_low = 0.01;
// ////////////        double firing_share_mid;
// ////////////
// ////////////        for (size_t n = 0; n < 10000; ++n) {
// ////////////            firing_share_mid = 0.5*firing_share_up + 0.5*firing_share_low;
// ////////////            // Capital Firing Cost
// //////////////        theta_C(para.dim1+para.dim2-1+3) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+3);
// //////////////        theta_C(para.dim1+para.dim2-1+4) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+4);
// ////////////            // Regular Workers Firing Cost
// ////////////            theta_C(para.dim1+para.dim2-1+7) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+7);
// ////////////            theta_C(para.dim1+para.dim2-1+8) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+8);
// ////////////            theta_C(para.dim1+para.dim2-1+9) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+9);
// ////////////            theta_C(para.dim1+para.dim2-1+10) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+10);
// ////////////            // Contract Workers Firing Cost
// ////////////            theta_C(para.dim1+para.dim2-1+13) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+13);
// ////////////            theta_C(para.dim1+para.dim2-1+14) = firing_share_mid*theta_Est_S(para.dim1+para.dim2-1+14);
// ////////////
// ////////////            tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData> t_Equ_Counterfactual
// ////////////                = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C, Val_StatusQuo.Equ0_S.F_Entry,
// ////////////                Val_StatusQuo.CapitalStock_S,Val_StatusQuo.PriceIndex_S,
// ////////////                Val_StatusQuo.FirmMass_S, sim_var, GoodState, LaborIntensive, 0.0,
// ////////////                threadsManagement);
// ////////////            EquState0 Equ0_C = get<6>(t_Equ_Counterfactual);
// ////////////            double FirmMass_C = Equ0_C.FirmMass;
// ////////////            SimData sim_data_C = get<9>(t_Equ_Counterfactual);
// ////////////
// ////////////            /*** Calculate Total Output ***/
// ////////////            double TotOutput_C = calTotalOutput(sim_data_C, FirmMass_C);
// ////////////
// ////////////            if (TotOutput_C > TargetOutput) {
// ////////////                firing_share_low = firing_share_mid;
// ////////////            }
// ////////////            else {
// ////////////                firing_share_up = firing_share_mid;
// ////////////            }
// ////////////
// ////////////            cout << "** TargetOutput = " << TargetOutput << "; TotOutput_C = " << TotOutput_C
// ////////////                 << "; firing_share_up = " << firing_share_up << "; firing_share_low = " << firing_share_low << endl;
// ////////////            if (abs(firing_share_up - firing_share_low) < 1e-3) {
// ////////////                break;
// ////////////            }
// ////////////        }
// ////////////        firing_share(k) = firing_share_mid;
// ////////////        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData> t_Equ_Counterfactual
// ////////////            = SolveSimulate_Counterfactural_FEntry_V0_EquPrice(theta_C, Val_StatusQuo.Equ0_S.F_Entry,
// ////////////            Val_StatusQuo.CapitalStock_S,Val_StatusQuo.PriceIndex_S,
// ////////////            Val_StatusQuo.FirmMass_S, sim_var, GoodState, LaborIntensive, KsupplyElas,
// ////////////            threadsManagement);
// ////////////        SimData sim_data_C = get<9>(t_Equ_Counterfactual);
// ////////////        EquState0 Equ0_C = get<6>(t_Equ_Counterfactual);
// ////////////        double FirmMass_C = Equ0_C.FirmMass;
// ////////////
// ////////////        /*** Calculate Average Exit Rate ***/
// ////////////        AvgExitRate(k) = calAvgExitRate(sim_data_C, FirmMass_C);
// ////////////
// ////////////        /*** Calculate Total Labor Employment ***/
// ////////////        tuple<double,double> t_TotEmp_C = calTotalEmployment(sim_data_C, FirmMass_C);
// ////////////        double TotEmploy_ur_C = get<0>(t_TotEmp_C);
// ////////////        double TotEmploy_uc_C = get<1>(t_TotEmp_C);
// ////////////        TotEmploy_Delta(k) = (TotEmploy_ur_C+TotEmploy_uc_C) / Val_StatusQuo.TotEmploy_S;
// ////////////
// ////////////        /*** Calculate Total Output ***/
// ////////////        double TotOutput_C = calTotalOutput(sim_data_C, FirmMass_C);
// ////////////        TotOutput_Delta(k) = TotOutput_C / Val_StatusQuo.TotOutput_S;
// ////////////
// ////////////        cout << "TotOutput_Delta = " << TotOutput_Delta(k) << endl;
// ////////////        cout << "TotEmploy_Delta = " << TotEmploy_Delta(k) << endl;
// ////////////        cout << "ResidualValue = " << ResidualValue(k) << endl;
// ////////////        cout << "firing_share = " << firing_share(k) << endl;
// ////////////        cout << "AvgExitRate = " << AvgExitRate(k) << endl;
// ////////////    }
// ////////////
// ////////////    cout << "TotOutput_Delta = " << TotOutput_Delta.transpose() << endl;
// ////////////    cout << "TotEmploy_Delta = " << TotEmploy_Delta.transpose() << endl;
// ////////////    cout << "ResidualValue = " << ResidualValue.transpose() << endl;
// ////////////    cout << "firing_share = " << firing_share.transpose() << endl;
// ////////////    cout << "AvgExitRate = " << AvgExitRate.transpose() << endl;
// ////////////
// ////////////    writeToCSVfile("TotOutput_Delta_" + State + FKL + ".csv", TotOutput_Delta.cast<double>().matrix());
// ////////////    writeToCSVfile("TotEmploy_Delta_" + State + FKL + ".csv", TotEmploy_Delta.cast<double>().matrix());
// ////////////    writeToCSVfile("ResidualValue_" + State + FKL + ".csv", ResidualValue.cast<double>().matrix());
// ////////////    writeToCSVfile("firing_share_" + State + FKL + ".csv", firing_share.cast<double>().matrix());
// ////////////    writeToCSVfile("AvgExitRate_" + State + FKL + ".csv", AvgExitRate.cast<double>().matrix());
// ////////////
// ////////////    return 0;
// ////////////
// ////////////}
// //////////
// //////////
// ////////////
// /////////////**************************************************************
// ////////////* Solve the counterfactural equilibrium and simulate firm distribution of counterfactual
// ////////////* (1) the value function at the initial period
// ////////////* (2) the value function for t >= 1
// ////////////**************************************************************/
// ////////////tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>
// ////////////    alias::SolveSimulate_Counterfactural_FEntry_V0( const ArrayXd & theta_C, const double & F_Entry,
// ////////////    const double & CapitalStock, const double & PriceIndex_S, const double & FirmMass_S, const SimVar & sim_var_C,
// ////////////    const int & GoodState, const int & LaborIntensive, const IndustryPara & para_ind, const double & KsupplyElas,
// ////////////    MultiThreads::Threads_Management & threadsManagement) {
// ////////////
// ////////////    EquState0 Equ0_C;
// ////////////    Equ0_C.F_Entry = F_Entry;
// ////////////
// ////////////    double Capital_Demand;
// ////////////
// ////////////    double p_K_C = 1;
// ////////////    double p_K_low = 0;
// ////////////    double p_K_up = 10;
// ////////////
// ////////////    double PI_C = 1;
// ////////////    double PI_low = 0;
// ////////////    double PI_up = 10;
// ////////////
// ////////////    ParaEst para_est_C;
// ////////////    ParaVec para_vec_C;
// ////////////    EquStateV0 EquV0_C;
// ////////////    EquStateV0mat EVal0mat_C;
// ////////////    EquStateV EquV_C;
// ////////////    EquStateVmat EValmat_C;
// ////////////    SimData sim_data_C;
// ////////////
// ////////////    for (size_t k = 0; k < 1000; ++k) {
// ////////////        cout << "p_K_C = " << p_K_C << endl;
// ////////////        double PI_C = 1;
// ////////////        double PI_low = 0;
// ////////////        double PI_up = 5;
// ////////////        for (size_t n = 0; n < 1; ++n) {
// ////////////            //// construct the set of parameters
// ////////////            para_est_C = constructParaEst(theta_C, GoodState, LaborIntensive);
// ////////////            para_est_C.PI = PI_C;
// ////////////            //// construct the grid of state space
// ////////////            para_vec_C = constructParaFull(para_est_C, p_K_C);
// ////////////            //// solve equilibrium under the counterfactual parameters and a guessed price p_K_mid
// ////////////            //    Rev RevOpt_C = RevOpt_Prod(para_est_C, para_vec_C);
// ////////////            //    ArrayXd VEnd_C = calVEnd(para_est_C, para_vec_C, RevOpt_C);
// ////////////            ArrayXd VEnd_C = ArrayXd::Zero(para.N);
// ////////////            //// solve value function for t > 1
// ////////////            tuple<EquStateV, EquStateVmat> t_Equ = solveV(para_est_C, para_vec_C, VEnd_C, threadsManagement);
// ////////////            EquV_C = get<0>(t_Equ);
// ////////////            EValmat_C = get<1>(t_Equ);
// ////////////            //// With free entry condition, solve value/policy function at the initial t = 1
// ////////////            tuple<EquStateV0, EquStateV0mat> t_Equ0 = solveV0(para_est_C,para_vec_C,EquV_C.EVal,threadsManagement);
// ////////////            EquV0_C = get<0>(t_Equ0);
// ////////////            EVal0mat_C = get<1>(t_Equ0);
// ////////////
// ////////////            //// calculate the expected value of entry under the counterfactual parameters and a guessed price p_K_mid
// ////////////            double V_Entry_C = 0.0;
// ////////////            for (size_t k = 0; k < para.N_phi; ++k) {
// ////////////                V_Entry_C = V_Entry_C + para_vec_C.dist0(k) * EquV0_C.EVal0(k * para.N_KL);
// ////////////            }
// ////////////            //// update PI
// ////////////            if (V_Entry_C >= F_Entry) {
// ////////////                PI_up = PI_C;
// ////////////            } else {
// ////////////                PI_low = PI_C;
// ////////////            }
// ////////////            cout << "V_Entry_C = " << V_Entry_C << "; F_Entry = " << F_Entry << "; PI_up = " << PI_up << "; PI_low = "
// ////////////                 << PI_low
// ////////////                 << "; diff of PI_up and PI_low = " << abs(PI_up - PI_low) << endl;
// ////////////            if (abs(PI_up - PI_low) < 1e-4) {
// ////////////                break;
// ////////////            }
// //////////////            PI_C = 0.5 * PI_low + 0.5 * PI_up;
// ////////////        }
// ////////////        //// Simulate the data to get the stationary distribution of capital employment
// ////////////        /** simulate a series of random number **/
// ////////////        int SimN_C = sim_var_C.lnphi_sim.rows();
// ////////////        int SimTbar_C = sim_var_C.lnphi_sim.cols() - 1;
// ////////////        /** simulate the data under the counterfactual data **/
// ////////////        sim_data_C = Simulation_StatusQuo_Counterfactual(para_est_C, para_vec_C, EquV_C,
// ////////////            EValmat_C, EquV0_C, EVal0mat_C, sim_var_C, SimN_C,
// ////////////            SimTbar_C, GoodState, LaborIntensive, threadsManagement);
// ////////////        double PI_temp = calPriceIndex(sim_data_C, 1.0);
// ////////////        double FirmMass_C = pow(PI_C,1-para.sigma) * pow(PriceIndex_S, 1.0 - para.sigma) / pow(PI_temp, 1.0 - para.sigma);
// ////////////        Capital_Demand = FirmMass_C * (sim_data_C.Capital * sim_data_C.State_S.cast<double>()).sum();
// ////////////        Equ0_C.FirmMass = FirmMass_C;
// ////////////        double tot_output = calTotalOutput(sim_data_C, FirmMass_C);
// ////////////
// ////////////        double AvgExitRate_C = calAvgExitRate(sim_data_C, FirmMass_C);
// //////////////        cout << "AvgExitRate_C = " << AvgExitRate_C << endl;
// ////////////
// ////////////        //// update p_K
// ////////////        if (Capital_Demand >= CapitalStock) {
// ////////////            p_K_low = p_K_C;
// ////////////        } else {
// ////////////            p_K_up = p_K_C;
// ////////////        }
// ////////////        cout << "PI_up = " << PI_up << "; para_est_C.PI = " << para_est_C.PI << "; para.sigma-1.0 = " << para.sigma-1.0
// ////////////            << "; PI_sigma = " << pow(para_est_C.PI,para.sigma-1.0)
// ////////////            << "; FirmMass_C = " << FirmMass_C << "; tot_output = " << tot_output
// ////////////            << "; Capital Demand = " << Capital_Demand << "; Capital Supply = " << CapitalStock
// ////////////            << "; p_K_up = " << p_K_up << "; p_K_low = " << p_K_low
// ////////////            << "; diff of p_K_up and p_K_low = " << abs(p_K_up - p_K_low) << endl;
// ////////////        if (abs(p_K_up - p_K_low) < 1e-3) {
// ////////////            break;
// ////////////        }
// ////////////        p_K_C = 0.5 * p_K_low + 0.5 * p_K_up;
// ////////////    }
// ////////////        throw runtime_error("163");
// ////////////
// ////////////    return tuple<ParaEst,ParaVec,EquStateV,EquStateVmat,EquStateV0,EquStateV0mat,EquState0,double,double,SimData>(para_est_C,
// ////////////        para_vec_C,EquV_C,EValmat_C,EquV0_C,EVal0mat_C,Equ0_C,p_K_C,Capital_Demand,sim_data_C);
// ////////////}