#include "SolveModel.h"

using namespace Eigen;
using namespace std;
using namespace alias;
using namespace Ipopt_Wrapper;

/**************************************************************
* solve value function for period >=1
**************************************************************/
tuple<EquStateV,EquStateVmat> alias::solveV(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & VEnd_PD,
    MultiThreads::Threads_Management & threadsManagement) {

    EquStateV EquV_PD;
    // exit value
    ArrayXd ResVal_E = calResidual_Value_Exit(para_est,para_vec);
    ArrayXd Vprime_PD_ini = ArrayXd::Zero(VEnd_PD.size());

    for (size_t n = 0; n < 50; ++n) {
        EquV_PD = solveV_OneLoop(para_est,para_vec,Vprime_PD_ini,ResVal_E,threadsManagement);
        Vprime_PD_ini = EquV_PD.EVal_PD;
    }
//    cout << "value function iteration: eps = " << eps << "; n = " << n << endl;
    // throw runtime_error("45");

    EquStateVmat Evalmat;
    Evalmat.EVal_PD_Lur_mat = ArrayXXd::Zero(para.N_Lur*2,para.N_phi*para.N_K);
    Evalmat.EVal_PD_SE_mat = ArrayXXd::Zero(para.N_Lur*2,para.N_phi*para.N_K);
    for (size_t i_PD = 0; i_PD < 2; ++i_PD) {
        ArrayXd EVal_Lur = EquV_PD.EVal_PD_Lur.segment(i_PD*para.N,para.N);
        Eigen::Map<const Eigen::ArrayXXd> temp_Lur_mat(EVal_Lur.data(),para.N_Lur,para.N_phi*para.N_K);
        Evalmat.EVal_PD_Lur_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) =  temp_Lur_mat;

        ArrayXd EVal_SE = EquV_PD.EVal_PD_SE.segment(i_PD*para.N,para.N);
        Eigen::Map<const Eigen::ArrayXXd> temp_mat(EVal_SE.data(),para.N_Lur,para.N_phi*para.N_K);
        Evalmat.EVal_PD_SE_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) =  temp_mat;
    }

    return tuple<EquStateV,EquStateVmat>(EquV_PD,Evalmat);
}

/**************************************************************
* Calculate the flow value added for every grids in the space
**************************************************************/
Rev alias::RevOpt_Prod_vec(const ParaEst & para_est, const ParaVec & para_vec) {

    /*** CES production function: Cobb Douglas ***/
    Rev RevOpt;

    ArrayXd vec_LK_Full = para_vec.vec_Lur_Full.pow(para_est.alpha_Lr) * para_vec.vec_K_Full.pow(para_est.alpha_K);

    // given the productivity phi, the expected phi in this period
    const double pi_term = pow(para_est.PI, para_est.sigma_tilde);
    RevOpt.RevOpt_Prod = para_vec.vec_phi_Full * vec_LK_Full.pow(para_est.sigma_tilde) * pi_term;

    Eigen::Map<const Eigen::ArrayXXd> RevOpt_Prod_mat_temp(RevOpt.RevOpt_Prod.data(),para.N_Lur,
        para.N_phi*para.N_K);
    RevOpt.RevOpt_Prod_mat = RevOpt_Prod_mat_temp;

    return RevOpt;
}


/**************************************************************
* Calculate the residual value of firms if exiting
**************************************************************/
//// calculation for every grids in the space
ArrayXd alias::calResidual_Value_Exit(const ParaEst & para_est, const ParaVec & para_vec) {

    ArrayXd log_Kbase = (1.0 + para_vec.vec_K_Full).log();
    ArrayXd log_Lbase = (1.0 + para_vec.vec_Lur_Full).log();
    ArrayXd F_Kbase = para_est.F_E_c_FK * log_Kbase;
    ArrayXd F_Lbase_ur = (para_est.c_low_F_ur * log_Lbase)
            * (para_vec.vec_Lur_Full < para.Lur_cutoff1 ).cast<double>()
            + (para_est.c_high_F_ur * log_Lbase)
            * (para_vec.vec_Lur_Full >= para.Lur_cutoff1 ).cast<double>();

    ArrayXd ResVal_E = - F_Lbase_ur - F_Kbase + para_vec.vec_K_Full;
    return ResVal_E;
}


/**************************************************************
* Define the form of hiring and firing costs
**************************************************************/
double alias::CalHFcost_Lbase(const double & c, const double & Lbase) {
    const double inv_Lbase2 = 1.0 / (Lbase * Lbase);
    double HFcost_Lbase =  c * log(Lbase+1.0) * inv_Lbase2;
    return HFcost_Lbase;
}

/**************************************************************
* Per loop, update the value function
**************************************************************/
EquStateV alias::solveV_OneLoop(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & Vprime_PD_ini,
    const ArrayXd & ResVal_E, MultiThreads::Threads_Management & threadsManagement) {

    EquStateV EquV;
    //// the expected value function with the productivity transition
    int PD = 2;
    ArrayXd EVprime_PD = para.deltaV * CalEVprime_mat(PD,Vprime_PD_ini,para_vec.tran_phi);

    /// Step 1. solve L_ur
    ArrayXd EVprime_PD_ELur = CalEVal_Lur_error(PD, para_vec, EVprime_PD, para_est.sigma_Lerror_ur,threadsManagement); // With a chosen target, firms' expected value depends on the shocks on employment
    int InitialPeriod = 0; // whether this is the initial period.
    tuple<ArrayXd,ArrayXd,ArrayXd,ArrayXd> t_Lur = solveOptLur(PD, para_est, para_vec, EVprime_PD_ELur, InitialPeriod,threadsManagement);
    EquV.OptLur_P = get<0>(t_Lur);
    EquV.OptLur_D = get<1>(t_Lur);
    EquV.EVal_PD_OptLurP = get<2>(t_Lur);
    EquV.EVal_PD_OptLurD = get<3>(t_Lur);
    // Eigen::Map<const Eigen::ArrayXXd> EVal_PD_SE_mat1(EquV.EVal_PD_SE.data(),para.N_Lur,para.N_phi*para.N_K);
    // cout << "EVal_PD_SE_mat1 = " << EVal_PD_SE_mat1 << endl;
    // Eigen::Map<const Eigen::ArrayXXd> EVal_PD_SE_mat2(EquV.EVal_PD_SE.data() + para.N,para.N_Lur,para.N_phi*para.N_K);
    // cout << "EVal_PD_SE_mat2 = " << EVal_PD_SE_mat2 << endl;
    // throw runtime_error("Not implemented");

    //// Step 3. Firms decide if produce or be dormant
    tuple<ArrayXd,ArrayXd,ArrayXd> t_PD = calEV_ProbStatus_PD(para_est,para_vec,EquV.EVal_PD_OptLurP,EquV.EVal_PD_OptLurD,
        threadsManagement);
    ArrayXd EVal_PD = get<0>(t_PD);
    EquV.ProbPD_P = get<1>(t_PD);
    EquV.ProbPD_D = get<2>(t_PD);

    //// Step 1. Entry and exit
    tuple<ArrayXd,ArrayXd,ArrayXd> t_SE = calEV_ProbStatus_SE_logit(para_est, para_vec, ResVal_E, EquV.EVal_PD_SE,
        threadsManagement);
    EquV.EVal_PD = get<0>(t_SE);
    EquV.ProbPD_S = get<1>(t_SE);
    EquV.ProbPD_E = get<2>(t_SE);
    // Eigen::Map<const Eigen::ArrayXXd> EVal_PD_mat1(EquV.EVal_PD.data(),para.N_Lur,para.N_phi*para.N_K);
    // cout << "EVal_PD_mat1 = " << EVal_PD_mat1 << endl;
    // Eigen::Map<const Eigen::ArrayXXd> EVal_PD_mat2(EquV.EVal_PD.data() + para.N,para.N_Lur,para.N_phi*para.N_K);
    // cout << "EVal_PD_mat2 = " << EVal_PD_mat2 << endl;
    // throw runtime_error("Not implemented");
    return EquV;
}
//
/**************************************************************
* Within each loop: we solve it by backward induction
**************************************************************/
// the expected value function with the productivity transition
 ArrayXd alias::CalEVprime_mat(const int & PD, const ArrayXd & Vprime_PD, const ArrayXXd & tran_phi) {

     ArrayXd EVprime_PD(para.N * PD);
     const MatrixXd tran_phi_T = tran_phi.matrix().transpose();

     for (int i_PD = 0; i_PD < PD; ++i_PD) {
         Eigen::Map<const MatrixXd> Vprime_mat(Vprime_PD.data() + i_PD * para.N, para.N_KL, para.N_phi);
         MatrixXd EVprime_mat(para.N_KL, para.N_phi);
         EVprime_mat.noalias() = Vprime_mat * tran_phi_T;
         Eigen::Map<ArrayXd>(EVprime_PD.data() + i_PD * para.N, para.N) =
             Eigen::Map<const ArrayXd>(EVprime_mat.data(), para.N);
     }

     return EVprime_PD;
 }

////// Step 3. solve Lur
///* With a chosen target, firms' expected value depends on the shocks on employment */
ArrayXd alias::CalEVal_Lur_error(const int & PD, const ParaVec & para_vec, const ArrayXd & EVal_PD,
    const double & sigmaL_ur, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd EVprime_PD_ELur(para.N*PD);
    const int N_Lur = para.N_Lur;
    const int N_phi_K = para.N_phi * para.N_K;
    const double tail_upper = 1e10;

    for (int i_PD = 0; i_PD < PD; ++i_PD) {
        Eigen::Map<const Eigen::ArrayXXd> EVal_mat(EVal_PD.data() + i_PD * para.N, N_Lur, N_phi_K);
        ArrayXd EVprime_ELur(para.N);

        auto worker = [&](size_t i, unsigned thread_id) {
            ArrayXd tempY = EVal_mat.col(i);
            ArrayXXd coef_vec = LinearInterpolation1D_Coeff(para_vec.vec_Lur, tempY);
            int i_base = static_cast<int>(i) * N_Lur;

            for (int i_Lur = 0; i_Lur < N_Lur; ++i_Lur) {
                double mu = log(para_vec.vec_Lur(i_Lur));
                double value = coef_vec(0, 0) * lognormCDF(0.0, para_vec.vec_Lur(0), mu, sigmaL_ur)
                               + coef_vec(0, 1) * lognormConditionalEpectation(0.0, para_vec.vec_Lur(0), mu, sigmaL_ur);

                for (int kk = 1; kk < N_Lur; ++kk) {
                    value += coef_vec(kk, 0)
                             * lognormCDF(para_vec.vec_Lur(kk - 1), para_vec.vec_Lur(kk), mu, sigmaL_ur);
                    value += coef_vec(kk, 1)
                             * lognormConditionalEpectation(para_vec.vec_Lur(kk - 1), para_vec.vec_Lur(kk), mu, sigmaL_ur);
                }

                value += tempY(N_Lur - 1) * lognormCDF(para_vec.vec_Lur(N_Lur - 1), tail_upper, mu, sigmaL_ur);
                EVprime_ELur(i_base + i_Lur) = value;
            }
        };
        MultiThreads::simple_parallel_for(worker, N_phi_K, threadsManagement);

        EVprime_PD_ELur.segment(i_PD*para.N,para.N) = EVprime_ELur;

    }

    return EVprime_PD_ELur;
}

///* Choose the targeted Luc to maximize the expected value */
tuple<ArrayXd,ArrayXd,ArrayXd,ArrayXd> alias::solveOptLur(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
    const ArrayXd & EVal_PD_Lur, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd TotVOpt_PD = ArrayXd::Zero(para.N*PD);
    ArrayXd OptLur_PD = ArrayXd::Zero(para.N*PD);
    const double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));

    ArrayXXd EVal_Lur_PD_mat( para.N_Lur*PD, para.N_phi * para.N_K);
    for (int i_PD = 0; i_PD < PD; ++i_PD) {
        Eigen::Map<const Eigen::ArrayXXd> EVal_Lur_mat(EVal_PD_Lur.data() + i_PD * para.N, para.N_Lur, para.N_phi*para.N_K);
        EVal_Lur_PD_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) = EVal_Lur_mat;
    }

    auto worker = [&](size_t ii, unsigned thread_id) {
    // size_t thread_id = 0;
    // for (size_t ii = 0; ii < para.N*PD; ++ii) {
        size_t i_PD = ii / para.N;
        size_t i = ii - i_PD * para.N;

        const ArrayXXd EVal_Lur_mat = EVal_Lur_PD_mat(seqN(i_PD * para.N_Lur,para.N_Lur), all);

        Pos pos;
        pos.i_phi = i / (para.N_K * para.N_Lur);
        pos.i_K = (i - pos.i_phi * para.N_K * para.N_Lur ) / para.N_Lur;
        pos.i_Lur = i - pos.i_phi * para.N_K * para.N_Lur - pos.i_K * para.N_Lur;
        int i_state = static_cast<int>(pos.i_phi * para.N_K + pos.i_K);

        // cout << "i = " << i << "; pos.i_phi = " << pos.i_phi << "; pos.i_K = " << pos.i_K
        //     << "; pos.i_Lur = " << pos.i_Lur << "; int i_state = " << i_state << endl;

        double Lbase = para_vec.vec_Lur(pos.i_Lur);

        double H_Lbase_ur = 0.0; double F_Lbase_ur = 0.0;
        if (InitialPeriod == 0) {
            H_Lbase_ur = CalHFcost_Lbase(para_est.c_H_ur, Lbase);
            if (Lbase < para.Lur_cutoff1) {
                F_Lbase_ur = CalHFcost_Lbase(para_est.c_low_F_ur, Lbase);
            }
            else {
                F_Lbase_ur = CalHFcost_Lbase(para_est.c_high_F_ur, Lbase);
            }
        }

        tuple<double,double,int> t_V = solve1D_L_LinearSpline(H_Lbase_ur, F_Lbase_ur, wstar, Lbase,
            para_vec.vec_Lur, EVal_Lur_mat.col(i_state));
        OptLur_PD(ii) = get<0>(t_V);
        TotVOpt_PD(ii) = get<1>(t_V);
    // }
    };
    MultiThreads::simple_parallel_for(worker, para.N*PD, threadsManagement);
//    cout << "699 = " << TotVOpt_PD.transpose() << endl;
    // throw runtime_error("Error in SolveModel function");
    return tuple<ArrayXd,ArrayXd>(OptLur_PD,TotVOpt_PD);
}


//// Step 2. Firms decide if produce or be dormant
tuple<ArrayXd,ArrayXd,ArrayXd> alias::calEV_ProbStatus_PD(const ParaEst & para_est, const ParaVec & para_vec,
    const ArrayXd & Vprime_LurP, const ArrayXd & Vprime_LurD, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd log_diffV_PD = (Vprime_LurP - Vprime_LurD).log();
    ArrayXd dEVal_PD_sigma(para.N * 2);
    dEVal_PD_sigma.segment(0, para.N) = (log_diffV_PD - para_est.F_P_PP) / para_est.sigma_PD;
    dEVal_PD_sigma.segment(para.N, para.N) = (log_diffV_PD - para_est.F_P_DP) / para_est.sigma_PD;
    ArrayXd ProbPD_P = normCDF_vec(dEVal_PD_sigma, 0, 1);
    ArrayXd ProbPD_D = 1.0 - ProbPD_P;

    ArrayXd ExpProdF_PD(para.N*2);
    auto worker = [&](size_t i, unsigned thread_id){
        const double rev = RevOpt.RevOpt_Prod(i);
        const double vP = Vprime_P(i);
        const double vD = Vprime_D(i);
        const double cutoff = rev + vP - vD;

        const double prob_pp = ProbPD_P(i);
        const double prob_dp = ProbPD_P(i + para.N);

        const double expF_pp = lognormConditionalEpectation(0.0, cutoff, para_est.F_P_PP, para_est.sigma_PD);
        const double expF_dp = lognormConditionalEpectation(0.0, cutoff, para_est.F_P_DP, para_est.sigma_PD);

        const double expRev_pp = (rev + vP) * prob_pp + vD * (1.0 - prob_pp);
        const double expRev_dp = (rev + vP) * prob_dp + vD * (1.0 - prob_dp);

        ExpProdF_PD(i) = expRev_pp - expF_pp;
        ExpProdF_PD(i + para.N) = expRev_dp - expF_dp;
    };
    MultiThreads::simple_parallel_for(worker, para.N, threadsManagement);
    ArrayXd Vprime_PD_up = ExpProdF_PD;

    return tuple<ArrayXd,ArrayXd,ArrayXd>(Vprime_PD_up,ProbPD_P,ProbPD_D);
}


//
//// Step 1. Entry and exit ()
tuple<ArrayXd,ArrayXd,ArrayXd> alias::calEV_ProbStatus_SE_logit(const ParaEst & para_est,
    const ArrayXd & ResVal_E, const ArrayXd & EVal_PD, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd ResVal_PD_E(para.N*2);
    ResVal_PD_E << ResVal_E,ResVal_E;
    /***********************************************************************/
    ArrayXd EVal_PD_S_sigma = EVal_PD / para_est.sigma_SE;
    ArrayXd EVal_PD_E_sigma = (ResVal_PD_E + para_est.F_E) / para_est.sigma_SE;
    ArrayXd dEVal_PD_SE_sigma = EVal_PD_S_sigma - EVal_PD_E_sigma;

    ArrayXd Vprime_PD_up(para.N*2); ArrayXd ProbPD_S(para.N*2); ArrayXd ProbPD_E(para.N*2);
    auto worker = [&](size_t i, unsigned thread_id) {
        //    size_t thread_id = 0;
        //    for (size_t i = 0; i < para.N*2; ++i) {
        double x = dEVal_PD_SE_sigma(i);
        if (x <= 700) {
            double exp_x = exp(x);
            ProbPD_S(i) = exp(x - log1p(exp_x));
            ProbPD_E(i) = 1.0 - ProbPD_S(i);
            Vprime_PD_up(i) = para_est.sigma_SE * log1p(exp_x) + para_est.sigma_SE * EVal_PD_E_sigma(i);
        }
        else {
            ProbPD_S(i) = 1.0;
            ProbPD_E(i) = 0.0;
            Vprime_PD_up(i) = para_est.sigma_SE * x + para_est.sigma_SE * EVal_PD_E_sigma(i);
        }
        //    }
    };
    MultiThreads::simple_parallel_for(worker, para.N*2, threadsManagement);
    return tuple<ArrayXd,ArrayXd,ArrayXd>(Vprime_PD_up,ProbPD_S,ProbPD_E);
}

//

//
// /* Choose the targeted Luc to maximize the expected value by iteration */
// tuple<ArrayXd,ArrayXd> alias::solveOptLur_iteration(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
//     const ArrayXd & EVal_PD_Lur, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement) {
//
//     ArrayXd TotVOpt_PD = ArrayXd::Zero(para.N*PD);
//     ArrayXd OptLur_PD = ArrayXd::Zero(para.N*PD);
//     const double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
//
//     ArrayXXd EVal_Lur_PD_mat( para.N_Lur*PD, para.N_phi * para.N_K);
//     for (int i_PD = 0; i_PD < PD; ++i_PD) {
//         Eigen::Map<const Eigen::ArrayXXd> EVal_Lur_mat(EVal_PD_Lur.data() + i_PD * para.N, para.N_Lur, para.N_phi*para.N_K);
//         EVal_Lur_PD_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) = EVal_Lur_mat;
//     }
//
//     auto worker = [&](size_t ii, unsigned thread_id) {
//         //    size_t thread_id = 0;
//         //    for (size_t ii = 0; ii < para.N*PD; ++ii) {
//         size_t i_PD = ii / para.N;
//         size_t i = ii - i_PD * para.N;
//
//         const auto EVal_Lur_mat = EVal_Lur_PD_mat(seqN(i_PD * para.N_Lur,para.N_Lur), all);
//
//         Pos pos;
//         pos.i_phi = i / (para.N_K * para.N_Lur);
//         pos.i_K = (i - pos.i_phi * para.N_K * para.N_Lur ) / para.N_Lur;
//         pos.i_Lur = (i - pos.i_phi * para.N_K * para.N_Lur - pos.i_K * para.N_Lur ) / para.N_Lur;
//         int i_state = static_cast<int>(pos.i_phi * para.N_K + pos.i_K);
//
//         //        cout << "i = " << i << "; pos.i_phi = " << pos.i_phi << "; pos.i_K = " << pos.i_K
//         // << "; pos.i_Lur = " << pos.i_Lur << "; int i_N = " << i_N << endl;
//
//         double Lbase = para_vec.vec_Lur(pos.i_Lur);
//
//         double H_Lbase_ur = 0.0; double F_Lbase_ur = 0.0;
//         if (InitialPeriod == 0) {
//             H_Lbase_ur = CalHFcost_Lbase(para_est.c_H_ur, Lbase);
//             if (Lbase < para.Lur_cutoff1) {
//                 F_Lbase_ur = CalHFcost_Lbase(para_est.c_low_F_ur, Lbase);
//             }
//             else {
//                 F_Lbase_ur = CalHFcost_Lbase(para_est.c_high_F_ur, Lbase);
//             }
//         }
//
//         ArrayXd diff_vec_L = para_vec.vec_Lur - Lbase;
//         ArrayXd TotV = EVal_Lur_mat.col(i_state)
//             - H_Lbase_ur * diff_vec_L.pow(2) * (diff_vec_L > 0).cast<double>()
//             - F_Lbase_ur * diff_vec_L.pow(2) * (diff_vec_L <= 0).cast<double>()
//             - wstar * para_vec.vec_Lur;
//
//         ArrayXd::Index max_Index;
//         double Vdefault = TotV.maxCoeff(&max_Index);
//         double xdefault =  para_vec.vec_Lur(max_Index);
//
//         OptLur_PD(ii) = xdefault;
//         TotVOpt_PD(ii) = Vdefault;
// //    }
//     };
//     MultiThreads::simple_parallel_for(worker, para.N*PD, threadsManagement);
//
//     return tuple<ArrayXd,ArrayXd>(OptLur_PD,TotVOpt_PD);
// }

//// Find the optimal solution (Lur) in Step 2
tuple<double,double,int> alias::solve1D_L_LinearSpline(const double & H, const double & F, const double & wstar,
    const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec) {

    int N_L = vec_L.size();

    ArrayXd diff_vec_L = vec_L - Lbase;
    ArrayXd TotV = EVal_L_vec - H * diff_vec_L.pow(2) * (diff_vec_L > 0).cast<double>()
                   - F * diff_vec_L.pow(2) * (diff_vec_L <= 0).cast<double>()
                   - wstar * vec_L;

    ArrayXd::Index max_Index;
    double Vdefault = TotV.maxCoeff(&max_Index);
    double xdefault = vec_L(max_Index);
    int i_L = max_Index;

    double L_opt = vec_L(i_L); double v_max = Vdefault;
    int L_index_max = i_L;

    for (int offset = -1; offset <= 0; ++offset) {
        int L_index = i_L + offset;
        if (L_index >= 0 && L_index <= N_L-2) {

            ArrayXd Lvec = vec_L.segment(L_index,2);
            ArrayXd coef_EVal = CalLinearspline1D_Coeff(vec_L(L_index),
                                                        vec_L(L_index+1),EVal_L_vec(L_index),
                                                        EVal_L_vec(L_index+1));

            tuple<double, double> t_max = SolveMax_Linear1D_L(coef_EVal, vec_L(L_index),
                                                              vec_L(L_index+1), TotV(L_index),
                                                              TotV(L_index+1), H, F, wstar, Lbase);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
                L_index_max = L_index;
            }
        }
    }

    return tuple<double,double,int>(L_opt,v_max,L_index_max);
}

double alias::solve1D_L_LinearSpline_Diff(const double & H, const double & F, const double & wstar,
    const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const int & L_index_max) {

    int N_L = vec_L.size();

    ArrayXd diff_vec_L = vec_L - Lbase;
    ArrayXd TotV = EVal_L_vec - H * diff_vec_L.pow(2) * (diff_vec_L > 0).cast<double>()
        - F * diff_vec_L.pow(2) * (diff_vec_L <= 0).cast<double>()
        - wstar * vec_L;

    double Vdefault = TotV(L_index_max);
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

    return L_opt;
}

//
// /**************************************************************
// * Solve the value/policy function in the first period
// **************************************************************/
// tuple<EquStateV0,EquStateV0mat> alias::solveV0(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & Vprime_PD,
//     MultiThreads::Threads_Management & threadsManagement) {
//
//     Rev RevOpt = RevOpt_Prod(para_est, para_vec);
//
//     EquStateV0 EquV0;
//     //// the expected value function with the productivity transition
//     int PD = 1;
//     ArrayXd EVprime_PD = para.deltaV * CalEVprime_mat(PD,Vprime_PD.segment(0,para.N), para_vec.tran_phi);
//
//     //// Firms always produce at the first period: EVal_PD_temp only have a dimension of para.N
//     EquV0.ProbPD_P0 = ArrayXd::Ones(para.N);
//     EquV0.ProbPD_D0 = ArrayXd::Zero(para.N);
//
//     EquV0.EVal_PD_Lur0 = CalEVal_Lur_error(PD, para_vec, EVprime_PD, para_est.sigma_Lerror_ur,threadsManagement); // With a chosen target, firms' expected value depends on the shocks on employment
//     int InitialPeriod = 1; // whether this is the initial period.
//     //// solve L_ur
//     tuple<ArrayXd,ArrayXd> t_Lur = solveOptLur(PD, para_est, para_vec,EquV0.EVal_PD_Lur0,InitialPeriod,threadsManagement);
//     EquV0.OptLur_PD0 = get<0>(t_Lur);
//     EquV0.EVal_PD_K0 = get<1>(t_Lur);
//
//     //// solve K
//     tuple<ArrayXd,ArrayXd> t_K = solveOptK(PD, para_est,para_vec,EquV0.EVal_PD_K0,InitialPeriod,threadsManagement);
//     EquV0.OptK_PD0 = get<0>(t_K);
//     EquV0.EVal_PD0 = get<1>(t_K);
//
//     EquStateV0mat EVal0mat;
//     Eigen::Map<const Eigen::ArrayXXd> temp_Lur0_mat(EquV0.EVal_PD_Lur0.data(),para.N_Lur,para.N_phi*para.N_K);
//     EVal0mat.EVal_PD_Lur0_mat =  temp_Lur0_mat;
//
//     Eigen::Map<const Eigen::ArrayXXd> temp_K0_mat(EquV0.EVal_PD_K0.data(),para.N_Lur,para.N_phi*para.N_K);
//     EVal0mat.EVal_PD_K0_mat = temp_K0_mat;
//
//     Eigen::Map<const Eigen::ArrayXXd> temp_0_mat(EquV0.EVal_PD0.data(),para.N_Lur,para.N_phi*para.N_K);
//     EVal0mat.EVal_PD0_mat = temp_0_mat;
//
//     return tuple<EquStateV0,EquStateV0mat>(EquV0,EVal0mat);
// }
// //
// /* Choose the targeted K to maximize the expected value */
// tuple<ArrayXd,ArrayXd> alias::solveOptK(const int & PD, const ParaEst & para_est,const ParaVec & para_vec,
//     const ArrayXd & EVal_PD_K, MultiThreads::Threads_Management & threadsManagement) {
//
//     ArrayXd TotVOpt_PD = ArrayXd::Zero(para.N*PD);
//     ArrayXd OptK_PD = ArrayXd::Zero(para.N*PD);
//
//     ArrayXXd EVal_K_PD_mat(para.N_Lur*PD, para.N_phi * para.N_K);
//     for (size_t i_PD = 0; i_PD < PD; ++i_PD) {
//         Eigen::Map<const Eigen::ArrayXXd> EVal_K_mat(EVal_PD_K.data() + i_PD * para.N, para.N_Lur, para.N_phi * para.N_K);
//         EVal_K_PD_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) = EVal_K_mat;
//     }
//
//     auto worker = [&](size_t ii, unsigned thread_id) {
// //    size_t thread_id = 0;
// //    for (size_t ii = 0; ii < para.N*PD; ++ii) {
//         size_t i_PD = ii / para.N;
//         size_t i = ii - i_PD * para.N;
//
//         Pos pos;
//         pos.i_phi = i / (para.N_K * para.N_Lur);
//         pos.i_K = (i - pos.i_phi * para.N_K * para.N_Lur) / para.N_Lur;
//         pos.i_Lur = i - pos.i_phi * para.N_K * para.N_Lur - pos.i_K * para.N_Lur;
//
//         ArrayXXd EVal_K_mat = EVal_K_PD_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all);
//         ArrayXd EVal_K_vec = EVal_K_mat(pos.i_Lur, seqN(pos.i_phi * para.N_K, para.N_K));
//
//         double p_K = 1.0;
//
//         ArrayXd diff_vec_K = vec_K - Kbase;
//         ArrayXd TotV = EVal_K_vec - H * diff_vec_K.pow(2) * (diff_vec_K > 0).cast<double>()
//             - F * diff_vec_K.pow(2) * (diff_vec_K <= 0).cast<double>()
//             - rstar * vec_K
//             - p_K * (vec_K - Kbase) * (diff_vec_K > 0).cast<double>()
//             - c_K * p_K * (vec_K - Kbase) * (diff_vec_K <= 0).cast<double>();
//
//         ArrayXd::Index max_Index;
//         double Vdefault = TotV.maxCoeff(&max_Index);
//         double xdefault = vec_K(max_Index);
//         int i_K = max_Index;
//
//
//         tuple<double,double,int> t_K = solve1D_K_LinearSpline(H_Kbase,F_Kbase,
//             c_K,rstar,p_K,Kbase,para_vec.vec_K,EVal_K_vec);
//         OptK_PD(ii) = get<0>(t_K);
//         TotVOpt_PD(ii) = get<1>(t_K);
// //        cout << "i = " << i << "; Kbase = " << Kbase <<"; OptK(i) = " <<  OptK(i) << "; TotVOpt(i) = " << TotVOpt(i) << endl;
// //    }
//     };
//     MultiThreads::simple_parallel_for(worker, para.N*PD, threadsManagement);
// //    throw runtime_error("311");
//
//     return tuple<ArrayXd,ArrayXd>(OptK_PD,TotVOpt_PD);
// }

//
//
/////**************************************************************
//// * Calculate the value function of the end period:
//// * Approximation: Assume that firms stop change capital or employment after 40 period. This is for computation only
////**************************************************************/
////ArrayXd alias::calVEnd(const ParaEst & para_est, const ParaVec & para_vec, const Rev & RevOpt) {
////
////    ArrayXd VEnd = RevOpt.RevOpt_Prod - para_est.r_K*para_vec.p_K*para_vec.vec_K_Full
////        - para_est.w_ur*para_vec.vec_Lur_Full - para_est.w_uc*para_vec.vec_Luc_Full;
////
////    return VEnd;
////}
////
//
//
////
/////**************************************************************
////* Testing the concacvity of value functions
////**************************************************************/
////void alias::testConcavity_L(const ArrayXd & Y, const ArrayXd & vec_Luc, const ArrayXd & vec_Lur) {
////
////    Eigen::Map<const Eigen::ArrayXXd> Y_mat(Y.data(),para.N_Luc,para.N_Lur);
////
////    ArrayXXd dY_dLuc(para.N_Luc,para.N_Lur);
////    ArrayXXd dY_dLur(para.N_Luc,para.N_Lur);
////
////    for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
////        dY_dLuc.col(i_Lur) = construct_dv_hermite_Schumaker(vec_Luc,Y_mat.col(i_Lur));
////    }
//////    cout << "Y_mat = \n" << Y_mat << endl;
//////    cout << "dY_dLuc = \n" << dY_dLuc << endl;
////    ArrayXXi dY_dLuc_mat = ( dY_dLuc(seqN(1,para.N_Luc-1),all)
////                           > dY_dLuc(seqN(0,para.N_Luc-1),all)+1e-5 ).cast<int>();
////    if (dY_dLuc_mat.sum() > 0) {
////        cout << "!!!! non concave dY_dLuc_mat" << endl;
////        cout << "Y_mat = \n" << Y_mat << endl;
////        cout << "dY_dLuc = \n" << dY_dLuc << endl;
////        cout << "diff dY_dLuc = \n" << dY_dLuc(seqN(1,para.N_Luc-1),all)
////            - dY_dLuc(seqN(0,para.N_Luc-1),all) << endl;
////        cout << "dY_dLuc_mat = \n" << dY_dLuc_mat << endl;
////        cout << "max dY_dLuc(seqN(1,para.N_Luc-1),all) - dY_dLuc(seqN(0,para.N_Luc-1),all) = "
////             << (dY_dLuc(seqN(1,para.N_Luc-1),all)
////                 - dY_dLuc(seqN(0,para.N_Luc-1),all)).maxCoeff() << endl;
////        throw runtime_error("1034");
////    }
////
////
////    for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
////        dY_dLur.row(i_Luc) = construct_dv_hermite_Schumaker(vec_Lur,Y_mat.row(i_Luc)).transpose();
////    }
////
////    ArrayXXi dY_dLur_mat = ( dY_dLur(all,seqN(1,para.N_Lur-1))
////                             > dY_dLur(all,seqN(0,para.N_Lur-1))+1e-5 ).cast<int>();
////    if (dY_dLur_mat.sum() > 0) {
////        cout << "!!!! non concave dY_dLur_mat" << endl;
////        cout << "Y_mat = \n" << Y_mat << endl;
////        cout << "dY_dLur = \n" << dY_dLur << endl;
////        cout << "diff dY_dLur = \n" << dY_dLur(all,seqN(1,para.N_Lur-1))
////            - dY_dLur(all,seqN(0,para.N_Lur-1)) << endl;
////        cout << "dY_dLur_mat = \n" << dY_dLur_mat << endl;
////        cout << "max dY_dLur(all,seqN(1,para.N_Lur-1)) - dY_dLur(all,seqN(0,para.N_Lur-1)) = "
////             << (dY_dLur(all,seqN(1,para.N_Lur-1))
////                 - dY_dLur(all,seqN(0,para.N_Lur-1))).maxCoeff() << endl;
////        throw runtime_error("1041");
////    }
////}
////
////void alias::testConcavity_K(const ArrayXXd & Y, const ArrayXd & vec_K) {
////
////    int N = Y.rows();
////    ArrayXXd dY_dK(N,para.N_K);
////    for (size_t i = 0; i < N; ++i) {
////        dY_dK.row(i) = construct_dv_hermite_Schumaker(vec_K,Y.row(i)).transpose();
////    }
////
////    ArrayXXi dY_dK_mat = ( dY_dK(all,seqN(1,para.N_K-1))
////            > dY_dK(all,seqN(0,para.N_K-1))+1e-5 ).cast<int>();
////
//////    cout << "Y = \n" << Y << endl;
//////    cout << "dY_dK = \n" << dY_dK << endl;
////    if (dY_dK_mat.sum() > 0) {
////        cout << "!!!! non concave dY_dK_mat" << endl;
////        cout << "Y = \n" << Y << endl;
////        cout << "dY_dK = \n" << dY_dK << endl;
////        cout << "diff dY_dK = \n" << dY_dK(all,seqN(1,para.N_K-1))
////            - dY_dK(all,seqN(0,para.N_K-1)) << endl;
////        cout << "dY_dK_mat = \n" << dY_dK_mat << endl;
////        cout << "max dY_dK(all,seqN(1,para.N_K-1)) - dY_dK(all,seqN(0,para.N_K-1)) = "
////            << (dY_dK(all,seqN(1,para.N_K-1))
////            - dY_dK(all,seqN(0,para.N_K-1))).maxCoeff() << endl;
////        throw runtime_error("1059");
////    }
////}
////
////void alias::testConcavity_L_loop(const ArrayXd & Eval, const ArrayXd & vec_Luc, const ArrayXd & vec_Lur) {
////
////    cout << "******************* testConcavity_L_loop ************************" << endl;
////    Eigen::Map<const Eigen::ArrayXXd> EVal_mat(Eval.data(),para.N_Luc*para.N_Lur,para.N_phi*para.N_K);
////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
////        for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
//////            cout << "i_phi = " << i_phi << "; i_K = " << i_K << endl;
////            ArrayXd temp1 = EVal_mat.col(i_phi*para.N_K+i_K);
////            testConcavity_L(EVal_mat.col(i_phi*para.N_K+i_K), vec_Luc, vec_Lur);
////        }
////    }
////    cout << "******************* testConcavity_L_loop::Success!! ************************" << endl;
////}
////
////void alias::testConcavity_K_loop(const ArrayXd & Eval, const ArrayXd & vec_K) {
////
////    cout << "******************* testConcavity_K_loop ************************" << endl;
////    Eigen::Map<const Eigen::ArrayXXd> Eval_mat(Eval.data(),para.N_Luc*para.N_Lur,para.N_phi*para.N_K);
////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
////        for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
//////            cout << "i_phi = " << i_phi  << "i_Lur = " << i_Lur << endl;
////            testConcavity_K(Eval_mat(seqN(i_Lur * para.N_Luc, para.N_Luc),
////                seqN(i_phi * para.N_K, para.N_K)), vec_K);
////        }
////    }
////    cout << "******************* testConcavity_K_loop::Success!! ************************" << endl;
////}
////
////void alias::testMonotonicity_by_K_loop(const ArrayXd & OptKL) {
////
////    cout << "******************* testMonotonicity_by_K_loop ************************" << endl;
////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
////        for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
////            for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
////                for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
////                    if (i_K > 0) {
////                        int i_plus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
////                                     i_Lur * para.N_Luc + i_Luc;
////                        int i_minus = i_phi * para.N_K * para.N_Lur * para.N_Luc + (i_K - 1) * para.N_Lur * para.N_Luc +
////                                      i_Lur * para.N_Luc + i_Luc;
////
////                        if (OptKL(i_plus) < OptKL(i_minus) - 1e-9) {
////                            cout << "K::: i_phi = " << i_phi << "; i_K = " << i_K << "; i_Luc = " << i_Luc << "; i_Lur = " << i_Lur
////                                 << "; OptKL(i_plus) = " << OptKL(i_plus) << "; OptKL(i_minus) = " << OptKL(i_minus)
////                                 << "; OptKL(i_plus) - OptKL(i_minus)= " << OptKL(i_plus) - OptKL(i_minus) << endl;
//////                            throw runtime_error("340");
////                        }
////                    }
////                }
////            }
////        }
////    }
////    cout << "******************* testMonotonicity_by_K_loop::Success!! ************************" << endl;
////}
////
////void alias::testMonotonicity_by_Lur_loop(const ArrayXd & OptKL) {
////
////    cout << "******************* testMonotonicity_by_Lur_loop ************************" << endl;
////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
////        for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
////            for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
////                for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
////                    if (i_Lur > 0) {
////                        int i_plus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
////                                     i_Lur * para.N_Luc + i_Luc;
////                        int i_minus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
////                                      (i_Lur - 1) * para.N_Luc + i_Luc;
////
////                        if (OptKL(i_plus) < OptKL(i_minus) - 1e-9) {
////                            cout << "Lur:::i_K = " << i_K << "; i_Luc = " << i_Luc << "; i_Lur = " << i_Lur
////                                 << "; OptKL(i_plus) = " << OptKL(i_plus) << "; OptKL(i_minus) = " << OptKL(i_minus)
////                                 << "; OptKL(i_plus) - OptK(i_minus)= " << OptKL(i_plus) - OptKL(i_minus) << endl;
//////                            throw runtime_error("340");
////                        }
////                    }
////                }
////            }
////        }
////    }
////    cout << "******************* testMonotonicity_by_Lur_loop::Success!! ************************" << endl;
////
////}
////
////void alias::testMonotonicity_by_Luc_loop(const ArrayXd & OptKL) {
////
////    cout << "******************* testMonotonicity_by_Luc_loop ************************" << endl;
////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
////        for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
////            for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
////                for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
////                    if (i_Luc > 0) {
////                        int i_plus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
////                                     i_Lur * para.N_Luc + i_Luc;
////                        int i_minus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
////                                      i_Lur * para.N_Luc + i_Luc - 1;
////
////                        if (OptKL(i_plus) < OptKL(i_minus) - 1e-9) {
////                            cout << "Luc:::i_K = " << i_K << "; i_Luc = " << i_Luc << "; i_Lur = " << i_Lur
////                                 << "; OptKL(i_plus) = " << OptKL(i_plus) << "; OptKL(i_minus) = " << OptKL(i_minus)
////                                 << "; OptKL(i_plus) - OptKL(i_minus)= " << OptKL(i_plus) - OptKL(i_minus) << endl;
//////                            throw runtime_error("340");
////                        }
////                    }
////                }
////            }
////        }
////    }
////    cout << "******************* testMonotonicity_by_Luc_loop::Success!! ************************" << endl;
////
////}
////
////
////
//////
//////tuple<ArrayXd,ArrayXd,ArrayXd> alias::calEV_ProbStatus_SE_lognormal(const ParaEst & para_est, const ParaVec & para_vec,
//////    const ArrayXd & ResVal_E, const ArrayXd & EVal, MultiThreads::Threads_Management & threadsManagement) {
//////
//////    ArrayXd Ones_d = ArrayXd::Ones(para.N);
//////    /***********************************************************************/
//////    ArrayXd dEVal_SE_sigma = EVal - ResVal_E;
//////
//////    ArrayXd Vprime_up(para.N);
//////    ArrayXd Prob_S(para.N);
//////    ArrayXd Prob_E(para.N);
////////    auto worker = [&](size_t i, unsigned thread_id) {
//////    for (size_t i = 0; i < para.N; ++i) {
//////
//////        if (dEVal_SE_sigma(i) <= 0) {
//////            Prob_S(i) = 0.0;
//////            Prob_E(i) = 1.0;
//////            Vprime_up(i) = ResVal_E(i) + exp(para_est.F_E);
//////        } else {
//////            Prob_S(i) = normCDF((log(dEVal_SE_sigma(i)) - para_est.F_E) / para_est.sigma_SE, 0.0, 1.0);
//////            Prob_E(i) = normCDF(-(log(dEVal_SE_sigma(i)) - para_est.F_E) / para_est.sigma_SE, 0.0, 1.0);;
//////            Vprime_up(i) = EVal(i) * Prob_S(i) + (ResVal_E(i) + exp(para_est.F_E)) * Prob_E(i);
//////        }
////////        cout << "dEVal_SE_sigma = " << dEVal_SE_sigma(i) << "; Prob_S(i) = " << Prob_S(i) << endl;
//////    }
////////    };
////////throw runtime_error("205");
//////    return tuple<ArrayXd,ArrayXd,ArrayXd>(Vprime_up,Prob_S,Prob_E);
//////}
////
////
////
