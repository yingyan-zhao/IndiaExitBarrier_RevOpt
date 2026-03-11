#include "SolveModel.h"
#include <limits>

using namespace Eigen;
using namespace std;
using namespace alias;
// using namespace Ipopt_Wrapper;

namespace {
MatrixXd BuildLurTransitionWeights(const ArrayXd& vec_Lur, const double sigmaL_ur) {
    const int N_Lur = static_cast<int>(vec_Lur.size());
    const ArrayXd log_Lur = vec_Lur.log();
    MatrixXd p_lur_transpose(N_Lur, N_Lur);
    for (int i_Lur = 0; i_Lur < N_Lur; ++i_Lur) {
        p_lur_transpose.col(i_Lur) =
            NormalGridProbabilities(log_Lur, std::log(vec_Lur(i_Lur)), sigmaL_ur).matrix();
    }
    return p_lur_transpose;
}

ArrayXd ApplyLurTransitionWeights(const int PD, const ArrayXd& EVal_PD, const MatrixXd& p_lur_transpose) {
    const int N_Lur = static_cast<int>(p_lur_transpose.rows());
    const int N_phi_K = para.N_phi * para.N_K;

    ArrayXd EVprime_PD_ELur(para.N * PD);
    for (int i_PD = 0; i_PD < PD; ++i_PD) {
        Eigen::Map<const MatrixXd> EVal_mat(EVal_PD.data() + i_PD * para.N, N_Lur, N_phi_K);
        MatrixXd EVprime_mat(N_Lur, N_phi_K);
        EVprime_mat.noalias() = p_lur_transpose.transpose() * EVal_mat;
        Eigen::Map<ArrayXd>(EVprime_PD_ELur.data() + i_PD * para.N, para.N) =
            Eigen::Map<const ArrayXd>(EVprime_mat.data(), para.N);
    }
    return EVprime_PD_ELur;
}
}

/**************************************************************
* solve value function for period >=1
**************************************************************/
EquStateV alias::solveV(const ParaEst & para_est, const ParaVec & para_vec, MultiThreads::Threads_Management & threadsManagement) {

    EquStateV EquV_PD;
    // exit value
    ArrayXd ResVal_E = calResidual_Value_Exit(para_est,para_vec);
    // Revenue per period
    ArrayXXd RevOpt_Prod_mat_val = RevOpt_Prod_mat(para_est, para_vec);
    const MatrixXd p_lur_transpose = BuildLurTransitionWeights(para_vec.vec_Lur, para_est.sigma_Lerror_ur);
    // value function of the ending period
    ArrayXd Vprime_PD_ini = ArrayXd::Zero(para.N+para.N);
    const int max_iter = 50;
    const double vf_tol = 1e-6;
    for (int n = 0; n < max_iter; ++n) {
        //// the expected value function with the productivity transition
        int PD = 2;
        ArrayXd EVprime_PD = para.deltaV * CalEVprime_mat(PD,Vprime_PD_ini,para_vec.tran_phi);

        /// Step 3. solve L_ur
        EquV_PD.EVal_PD_Lur_error = ApplyLurTransitionWeights(PD, EVprime_PD, p_lur_transpose); // With a chosen target, firms' expected value depends on the shocks on employment

        tuple<ArrayXd,ArrayXd,ArrayXd,ArrayXd> t_Lur = solveOptLur(PD, para_est, para_vec, EquV_PD.EVal_PD_Lur_error,
            RevOpt_Prod_mat_val, threadsManagement);
        EquV_PD.OptLur_P = get<0>(t_Lur);
        EquV_PD.OptLur_D = get<1>(t_Lur);
        EquV_PD.EVal_PD_OptLurP = get<2>(t_Lur);
        EquV_PD.EVal_PD_OptLurD = get<3>(t_Lur);

        // //// Step 2. Firms decide if produce or be dormant
        tuple<ArrayXd,ArrayXd,ArrayXd> t_PD = calEV_ProbStatus_PD(para_est,EquV_PD.EVal_PD_OptLurP,
        EquV_PD.EVal_PD_OptLurD, threadsManagement);
        EquV_PD.EVal_PD_ProdDorm = get<0>(t_PD);
        EquV_PD.ProbPD_P = get<1>(t_PD);
        EquV_PD.ProbPD_D = get<2>(t_PD);

        //// Step 1. Entry and exit
        tuple<ArrayXd,ArrayXd,ArrayXd> t_SE = calEV_ProbStatus_SE_logit(para_est, ResVal_E, EquV_PD.EVal_PD_ProdDorm,
        threadsManagement);
        EquV_PD.EVal_PD = get<0>(t_SE);
        EquV_PD.ProbPD_S = get<1>(t_SE);
        EquV_PD.ProbPD_E = get<2>(t_SE);

        const double err = (EquV_PD.EVal_PD - Vprime_PD_ini).abs().maxCoeff();
        Vprime_PD_ini = EquV_PD.EVal_PD;
        if (err < vf_tol) { break; }
    }

    //
    // Eigen::Map<const Eigen::ArrayXXd> TotVOpt_P_mat3(EquV_PD.EVal_PD.data(),para.N_Lur,para.N_phi*para.N_K);
    // cout << "TotVOpt_P_mat3 = " << endl;
    // cout << TotVOpt_P_mat3 << endl;
    //
    // Eigen::Map<const Eigen::ArrayXXd> ProbPD_S_mat3(EquV_PD.ProbPD_S.data(),para.N_Lur,para.N_phi*para.N_K);
    // cout << "ProbPD_S_mat3 = " << endl;
    // cout << ProbPD_S_mat3 << endl;
    // // throw runtime_error("54");

//
//     EquStateVmat Evalmat;
//     Evalmat.EVal_PD_Lur_mat = ArrayXXd::Zero(para.N_Lur*2,para.N_phi*para.N_K);
//     Evalmat.EVal_PD_SE_mat = ArrayXXd::Zero(para.N_Lur*2,para.N_phi*para.N_K);
//     for (size_t i_PD = 0; i_PD < 2; ++i_PD) {
//         ArrayXd EVal_Lur = EquV_PD.EVal_PD_Lur.segment(i_PD*para.N,para.N);
//         Eigen::Map<const Eigen::ArrayXXd> temp_Lur_mat(EVal_Lur.data(),para.N_Lur,para.N_phi*para.N_K);
//         Evalmat.EVal_PD_Lur_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) =  temp_Lur_mat;
//
//         ArrayXd EVal_SE = EquV_PD.EVal_PD_SE.segment(i_PD*para.N,para.N);
//         Eigen::Map<const Eigen::ArrayXXd> temp_mat(EVal_SE.data(),para.N_Lur,para.N_phi*para.N_K);
//         Evalmat.EVal_PD_SE_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) =  temp_mat;
//     }

    return EquV_PD;
}

/**************************************************************
* Some basic functions
**************************************************************/
//// Calculate the residual value of firms if exiting
ArrayXd alias::calResidual_Value_Exit(const ParaEst & para_est, const ParaVec & para_vec) {

    ArrayXd log_Kbase = (0.001 + para_vec.vec_K_Full);
    ArrayXd log_Lbase = (1.0 + para_vec.vec_Lur_Full);
    ArrayXd F_Kbase = para_est.F_E_c_FK * log_Kbase;
    ArrayXd F_Lbase_ur = (para_est.c_low_F_ur * log_Lbase)
            * (para_vec.vec_Lur_Full < para.Lur_cutoff1 ).cast<double>()
            + (para_est.c_high_F_ur * log_Lbase)
            * (para_vec.vec_Lur_Full >= para.Lur_cutoff1 ).cast<double>();

    ArrayXd ResVal_E = - F_Lbase_ur - F_Kbase + para_vec.vec_K_Full;

    return ResVal_E;
}

//// Define the form of hiring and firing costs
tuple<double,double> alias::CalHFcost_Lbase(const double & c_H_ur, const double & c_high_F_ur, const double & c_low_F_ur,
    const double & Lbase) {

    const double inv_Lbase2 = (Lbase+1.0) / (Lbase * Lbase);
    double H_Lbase_ur = c_H_ur * inv_Lbase2;
    const double F_adj = (Lbase > para.Lur_cutoff1) ? c_high_F_ur : c_low_F_ur;
    double F_Lbase_ur = F_adj * inv_Lbase2;
    return tuple<double,double>(H_Lbase_ur,F_Lbase_ur);
}

//// Calculate the flow value added for every grids in the space
ArrayXXd alias::RevOpt_Prod_mat(const ParaEst & para_est, const ParaVec & para_vec) {

    /*** CES production function: Cobb Douglas ***/
    ArrayXd vec_LK_Full = para_vec.vec_Lur_Full.pow(para_est.alpha_tilde_Lr) * para_vec.vec_K_Full.pow(para_est.alpha_tilde_K);

    // given the productivity phi, the expected phi in this period
    const double pi_term = pow(para_est.PI, para_est.sigma_tilde);
    ArrayXd RevOpt_Prod = para_vec.vec_phi_Full * vec_LK_Full * pi_term;

    Eigen::Map<const Eigen::ArrayXXd> RevOpt_Prod_mat_temp(RevOpt_Prod.data(),para.N_Lur,para.N_phi*para.N_K);
    ArrayXXd RevOpt_Prod_mat = RevOpt_Prod_mat_temp;

    return RevOpt_Prod_mat;
}

ArrayXXd alias::RevOpt_Prod_vec(const ParaEst & para_est, const ParaVec & para_vec) {

    /*** CES production function: Cobb Douglas ***/
    ArrayXd vec_LK_Full = para_vec.vec_Lur_Full_PD.pow(para_est.alpha_tilde_Lr) * para_vec.vec_K_Full_PD.pow(para_est.alpha_tilde_K);

    // given the productivity phi, the expected phi in this period
    const double pi_term = pow(para_est.PI, para_est.sigma_tilde);
    ArrayXd RevOpt_Prod = para_vec.vec_phi_Full_PD * vec_LK_Full * pi_term;

    return RevOpt_Prod;
}

//// Calculate the flow value added for every grids of Lur in the space
tuple<double,double,double> alias::RevOpt_Prod_ELur(const ParaEst & para_est, const double & phi, const double & K, const double & Lur) {

    const double pi_term = pow(para_est.PI, para_est.sigma_tilde);
    double A = pow(K,para_est.alpha_tilde_K) * pi_term * phi;

    double logL_alpha = para_est.alpha_tilde_Lr * log(Lur)
        + 0.5 * para_est.alpha_tilde_Lr * para_est.alpha_tilde_Lr * para_est.sigma_Lerror_ur * para_est.sigma_Lerror_ur;
    double RevOpt = A * exp(logL_alpha);

    double dRevOpt_dL = para_est.alpha_tilde_Lr * RevOpt / Lur;
    double dRevOpt_dL_coef = A * para_est.alpha_tilde_Lr
        * exp(0.5 * para_est.alpha_tilde_Lr * para_est.alpha_tilde_Lr * para_est.sigma_Lerror_ur * para_est.sigma_Lerror_ur);
    return tuple<double,double,double>(RevOpt,dRevOpt_dL,dRevOpt_dL_coef);
}

tuple<double,double,double> alias::RevOpt_Prod_K(const ParaEst & para_est, const double & phi, const double & K, const double & Lur) {

    const double pi_term = pow(para_est.PI, para_est.sigma_tilde);
    double logL_alpha = para_est.alpha_tilde_Lr * log(Lur)
        + 0.5 * para_est.alpha_tilde_Lr * para_est.alpha_tilde_Lr * para_est.sigma_Lerror_ur * para_est.sigma_Lerror_ur;
    double RevOpt = pi_term * phi * exp(logL_alpha) * pow(K,para_est.alpha_tilde_K);

    double dRevOpt_dK = para_est.alpha_tilde_K * RevOpt / K;
    double dRevOpt_dK_coef = para_est.alpha_tilde_K * pi_term * phi * exp(logL_alpha);

    return tuple<double,double,double>(RevOpt,dRevOpt_dK,dRevOpt_dK_coef);
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
        Eigen::Map<ArrayXd>(EVprime_PD.data() + i_PD * para.N, para.N) = Eigen::Map<const ArrayXd>(EVprime_mat.data(), para.N);
    }

    return EVprime_PD;
}

////// Step 3. solve Lur
///* With a chosen target, firms' expected value depends on the shocks on employment */
ArrayXd alias::CalEVal_Lur_error(const int & PD, const ParaVec & para_vec, const ArrayXd & EVal_PD, const double & sigmaL_ur) {
    const MatrixXd p_lur_transpose = BuildLurTransitionWeights(para_vec.vec_Lur, sigmaL_ur);
    return ApplyLurTransitionWeights(PD, EVal_PD, p_lur_transpose);
}

///* Choose the targeted Lur to maximize the expected value */
tuple<ArrayXd,ArrayXd,ArrayXd,ArrayXd> alias::solveOptLur(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
    const ArrayXd & EVal_PD, const ArrayXXd & RevOpt_Prod_mat, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd OptLur_P = ArrayXd::Zero(para.N*PD);
    ArrayXd OptLur_D = ArrayXd::Zero(para.N*PD);
    ArrayXd TotVOpt_P = ArrayXd::Zero(para.N*PD);
    ArrayXd TotVOpt_D = ArrayXd::Zero(para.N*PD);

    const double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
    const int N_Lur = para.N_Lur;
    const int N_K = para.N_K;
    const int N_state = para.N_phi * para.N_K;

    // Hiring/firing adjustments depend only on Lbase grid; precompute once.
    ArrayXd H_Lbase_ur_vec = ArrayXd::Zero(N_Lur);
    ArrayXd F_Lbase_ur_vec = ArrayXd::Zero(N_Lur);
    for (int i_Lur = 0; i_Lur < N_Lur; ++i_Lur) {
        const double Lbase = para_vec.vec_Lur(i_Lur);
        const tuple<double,double> t_HF = CalHFcost_Lbase(para_est.c_H_ur, para_est.c_high_F_ur, para_est.c_low_F_ur, Lbase);
        H_Lbase_ur_vec(i_Lur) = get<0>(t_HF);
        F_Lbase_ur_vec(i_Lur) = get<1>(t_HF);
    }

    for (int i_PD = 0; i_PD < PD; ++i_PD) {
        const int pd_offset = i_PD * para.N;
        Eigen::Map<const Eigen::ArrayXXd> EVal_Lur_mat(EVal_PD.data() + pd_offset, N_Lur, N_state);

        auto worker = [&](size_t i_raw, unsigned thread_id) {
        // for (int i_raw = 0; i_raw < para.N; ++i_raw) {
            const int i = static_cast<int>(i_raw);
            const int i_state = i / N_Lur;
            const int i_Lur = i - i_state * N_Lur;
            const int i_phi = i_state / N_K;
            const int i_K = i_state - i_phi * N_K;
            const int out_idx = pd_offset + i;

            const double Lbase = para_vec.vec_Lur(i_Lur);
            const double H_Lbase_ur = H_Lbase_ur_vec(i_Lur);
            const double F_Lbase_ur = F_Lbase_ur_vec(i_Lur);

            // Lur choices and value if firms choose dormancy.
            tuple<double,double,int> t_V_D = solve1D_L_LinearSpline(
                H_Lbase_ur, F_Lbase_ur, wstar, Lbase, para_vec.vec_Lur, EVal_Lur_mat.col(i_state));
            OptLur_D(out_idx) = get<0>(t_V_D);
            TotVOpt_D(out_idx) = get<1>(t_V_D);

            // Lur choices and value if firms choose production.
            tuple<double,double,int> t_V_P = solve1D_L_LinearSpline_RevOpt(
                H_Lbase_ur, F_Lbase_ur, wstar, Lbase,
                para_vec.vec_phi(i_phi), para_vec.vec_K(i_K), para_est, para_vec.vec_Lur,
                EVal_Lur_mat.col(i_state), RevOpt_Prod_mat.col(i_state));
            OptLur_P(out_idx) = get<0>(t_V_P);
            TotVOpt_P(out_idx) = get<1>(t_V_P);
        };
        // }
        MultiThreads::simple_parallel_for(worker, para.N, threadsManagement);
    }
    //
    // Eigen::Map<const Eigen::ArrayXXd> TotVOpt_P_mat1(TotVOpt_P.data(),para.N_Lur,para.N_phi*para.N_K);
    // cout << "TotVOpt_P_mat1 = " << endl;
    // cout << TotVOpt_P_mat1 << endl;
    // Eigen::Map<const Eigen::ArrayXXd> TotVOpt_P_mat2(TotVOpt_P.data()+para.N,para.N_Lur,para.N_phi*para.N_K);
    // cout << "TotVOpt_P_mat2 = " << endl;
    // cout << TotVOpt_P_mat2 << endl;
    //
    // Eigen::Map<const Eigen::ArrayXXd> TotVOpt_D_mat1(TotVOpt_D.data(),para.N_Lur,para.N_phi*para.N_K);
    // cout << "TotVOpt_D_mat1 = " << endl;
    // cout << TotVOpt_D_mat1 << endl;
    // Eigen::Map<const Eigen::ArrayXXd> TotVOpt_D_mat2(TotVOpt_D.data()+para.N,para.N_Lur,para.N_phi*para.N_K);
    // cout << "TotVOpt_D_mat2 = " << endl;
    // cout << TotVOpt_D_mat2 << endl;
    // //    cout << "699 = " << TotVOpt_PD.transpose() << endl;
    // throw runtime_error("Error in SolveModel function");
    return tuple<ArrayXd,ArrayXd,ArrayXd,ArrayXd>(OptLur_P,OptLur_D,TotVOpt_P,TotVOpt_D);
}

/// solve optimal L given firms decides to be dormant
tuple<double,double,int> alias::solve1D_L_LinearSpline(const double & H, const double & F, const double & wstar,
    const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec) {

    const int N_L = static_cast<int>(vec_L.size());

    auto total_value = [&](int idx) {
        const double L = vec_L(idx);
        const double d = L - Lbase;
        const double adj = (d > 0.0) ? H : F;
        return EVal_L_vec(idx) - adj * d * d - wstar * L;
    };

    int i_L = 0;
    double v_max = total_value(0);
    for (int i = 1; i < N_L; ++i) {
        const double v = total_value(i);
        if (v > v_max) {
            v_max = v;
            i_L = i;
        }
    }
    double L_opt = vec_L(i_L);
    int L_index_max = i_L;

    // derivative at i_L: right and left derivatives
    if (i_L >= 1 && i_L <= N_L-2) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_L(i_L-1), vec_L(i_L), EVal_L_vec(i_L-1), EVal_L_vec(i_L) );
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_L(i_L), vec_L(i_L+1), EVal_L_vec(i_L), EVal_L_vec(i_L+1) );
        double adj_L = (L_opt > Lbase) ? H : F;
        double slope_Left = coef_Left(1) - 2*adj_L*(L_opt - Lbase) - wstar;
        double slope_Right = coef_Right(1) - 2*adj_L*(L_opt - Lbase) - wstar;

        if (slope_Left>0 and slope_Right>0) { // maximum point lies in the right interval
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double> t_max = SolveMax_Linear1D_L(coef_Right, vec_L(i_L), vec_L(i_L+1), TotV0, TotV1,
                H, F, wstar, Lbase);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right<0){ // maximum point lies in the left interval
            double TotV0 = total_value(i_L-1);
            double TotV1 = total_value(i_L);
            tuple<double, double> t_max = SolveMax_Linear1D_L(coef_Left, vec_L(i_L-1), vec_L(i_L),
                TotV0, TotV1, H, F, wstar, Lbase);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right>0) { // this case should not happen
            // calculate the right interval first.
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double> t_maxRight = SolveMax_Linear1D_L(coef_Right, vec_L(i_L), vec_L(i_L+1),
                TotV0, TotV1, H, F, wstar, Lbase);
            double L_opt_temp_right = get<0>(t_maxRight);
            double v_max_temp_right = get<1>(t_maxRight);

            if (v_max_temp_right > v_max) {
                L_opt = L_opt_temp_right; v_max = v_max_temp_right;
            }

            // calculate the left interval then, and compare.
            TotV0 = total_value(i_L-1);
            TotV1 = total_value(i_L);
            tuple<double, double> t_maxLeft = SolveMax_Linear1D_L(coef_Left, vec_L(i_L-1), vec_L(i_L),
                TotV0, TotV1, H, F, wstar, Lbase);
            const double L_opt_temp_left = get<0>(t_maxLeft);
            const double v_max_temp_left = get<1>(t_maxLeft);
            if (v_max_temp_left > v_max) {
                L_opt = L_opt_temp_left;
                v_max = v_max_temp_left;
            }
        }
        else {L_opt = L_opt; v_max = v_max;}  // maximum point is right at the point.
    }
    else if (i_L == 0) {
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_L(i_L), vec_L(i_L+1), EVal_L_vec(i_L), EVal_L_vec(i_L+1) );
        double adj_L = (L_opt > Lbase) ? H : F;
        double slope_Right = coef_Right(1) - 2*adj_L*(L_opt - Lbase) - wstar;

        if (slope_Right>0) { // maximum point lies in the right interval
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double> t_max = SolveMax_Linear1D_L(coef_Right, vec_L(i_L), vec_L(i_L+1),
                TotV0, TotV1, H, F, wstar, Lbase);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else {L_opt = L_opt; v_max = v_max;}
    }
    else if (i_L == N_L-1) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_L(i_L-1), vec_L(i_L), EVal_L_vec(i_L-1), EVal_L_vec(i_L) );
        double adj_L = (L_opt > Lbase) ? H : F;
        double slope_Left = coef_Left(1) - 2*adj_L*(L_opt - Lbase) - wstar;

        if (slope_Left<0){ // maximum point lies in the left interval
            double TotV0 = total_value(i_L-1);
            double TotV1 = total_value(i_L);
            tuple<double, double> t_max = SolveMax_Linear1D_L(coef_Left, vec_L(i_L-1), vec_L(i_L),
                TotV0, TotV1, H, F, wstar, Lbase);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else {L_opt = L_opt; v_max = v_max;}
    }
    else {L_opt = L_opt; v_max = v_max;}
    return tuple<double,double,int>(L_opt,v_max,L_index_max);
}

tuple<double, double> alias::SolveMax_Linear1D_L_part(const ArrayXd & coef,const double & Adj, const double & wstar,
    const double & Lbase) {

    const double x_opt = (coef(1) - wstar) / (2.0 * Adj) + Lbase;
    const double dx = x_opt - Lbase;
    const double v_max = coef(0) + coef(1) * x_opt - Adj * dx * dx - wstar * x_opt;

    return tuple<double,double>(x_opt,v_max);
}

tuple<double, double> alias::SolveMax_Linear1D_L(const ArrayXd & coef, const double & x0, const double & x1,
    const double & TotV0,const double & TotV1, const double & H, const double & F, const double & wstar, const double & Lbase) {

    double x_opt = x0;
    double v_max = TotV0;

    if (Lbase <= x0 || Lbase >= x1) {
        double Adj = (Lbase <= x0) ? H : F;

        tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef,Adj,wstar,Lbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);
        if (x_opt_temp < x0) { x_opt = x0; v_max = TotV0; }
        else if (x_opt_temp > x1) { x_opt = x1; v_max = TotV1; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else {
        double VLbase = coef(0) + coef(1) * Lbase - wstar * Lbase;
        double dVLbase = coef(1) - wstar;
        if (dVLbase >= 0) {
            double Adj = H;

            tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef, Adj, wstar, Lbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < Lbase) { x_opt = Lbase; v_max = VLbase; }
            else if (x_opt_temp > x1) { x_opt = x1; v_max = TotV1; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        } else {
            double Adj = F;

            tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef, Adj, wstar, Lbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < x0) { x_opt = x0; v_max = TotV0; }
            else if (x_opt_temp > Lbase) { x_opt = Lbase; v_max = VLbase; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        }
    }

    return tuple<double, double>(x_opt,v_max);
}

/// solve optimal L given firms decides to produce
tuple<double,double,int> alias::solve1D_L_LinearSpline_RevOpt(const double & H, const double & F, const double & wstar,
    const double & Lbase, const double & phi, const double & K, const ParaEst & para_est,
    const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const ArrayXd & RevOpt_Prod_vec) {

    const int N_L = static_cast<int>(vec_L.size());

    auto total_value = [&](int idx) {
        const double L = vec_L(idx);
        const double d = L - Lbase;
        const double adj = (d > 0.0) ? H : F;
        return RevOpt_Prod_vec(idx) + EVal_L_vec(idx) - adj * d * d - wstar * L;
    };

    int i_L = 0;
    double v_max = total_value(0);
    for (int i = 1; i < N_L; ++i) {
        const double v = total_value(i);
        if (v > v_max) {
            v_max = v;
            i_L = i;
        }
    }

    double L_opt = vec_L(i_L);
    int L_index_max = i_L;

    tuple<double,double,double> t_RevOpt_Prod = RevOpt_Prod_ELur(para_est, phi, K, L_opt);
    double dRevOpt_dL = get<1>(t_RevOpt_Prod);
    double dRevOpt_dL_coef = get<2>(t_RevOpt_Prod);

    if (L_index_max >= 1 && L_index_max <= N_L-2) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_L(i_L-1), vec_L(i_L), EVal_L_vec(i_L-1), EVal_L_vec(i_L) );
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_L(i_L), vec_L(i_L+1), EVal_L_vec(i_L), EVal_L_vec(i_L+1) );
        double adj_L = (L_opt > Lbase) ? H : F;
        double slope_Left = coef_Left(1) - 2*adj_L*(L_opt - Lbase) - wstar + dRevOpt_dL;
        double slope_Right = coef_Right(1) - 2*adj_L*(L_opt - Lbase) - wstar + dRevOpt_dL;

        if (slope_Left>0 and slope_Right>0) { // maximum point lies in the right interval
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double> t_max = SolveMax_Linear1D_L_RevOpt(coef_Right, dRevOpt_dL_coef,
                vec_L(i_L), vec_L(i_L+1), TotV0, TotV1, H, F, wstar, Lbase, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right<0){ // maximum point lies in the left interval
            double TotV0 = total_value(i_L-1);
            double TotV1 = total_value(i_L);
            tuple<double, double> t_max = SolveMax_Linear1D_L_RevOpt(coef_Left, dRevOpt_dL_coef,
                vec_L(i_L-1), vec_L(i_L), TotV0, TotV1, H, F, wstar, Lbase, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right>0) { // this case should not happen
            // calculate the right interval first.
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double> t_maxRight = SolveMax_Linear1D_L_RevOpt(coef_Right, dRevOpt_dL_coef,
                vec_L(i_L), vec_L(i_L+1), TotV0, TotV1, H, F, wstar, Lbase, phi, K, para_est);
            double L_opt_temp_right = get<0>(t_maxRight);
            double v_max_temp_right = get<1>(t_maxRight);

            if (v_max_temp_right > v_max) {
                L_opt = L_opt_temp_right; v_max = v_max_temp_right;
            }

            // calculate the left interval then, and compare.
            TotV0 = total_value(i_L-1);
            TotV1 = total_value(i_L);
            tuple<double, double> t_maxLeft = SolveMax_Linear1D_L_RevOpt(coef_Left, dRevOpt_dL_coef,
                vec_L(i_L-1), vec_L(i_L), TotV0, TotV1, H, F, wstar, Lbase, phi, K, para_est);
            const double L_opt_temp_left = get<0>(t_maxLeft);
            const double v_max_temp_left = get<1>(t_maxLeft);
            if (v_max_temp_left > v_max) {
                L_opt = L_opt_temp_left;
                v_max = v_max_temp_left;
            }
        }
        else {L_opt = L_opt; v_max = v_max;}  // maximum point is right at the point.
    }
    else if (i_L == 0) {
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_L(i_L), vec_L(i_L+1), EVal_L_vec(i_L), EVal_L_vec(i_L+1) );
        double adj_L = (L_opt > Lbase) ? H : F;
        double slope_Right = coef_Right(1) - 2*adj_L*(L_opt - Lbase) - wstar + dRevOpt_dL;
        if (slope_Right>0) {
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double> t_max = SolveMax_Linear1D_L_RevOpt(coef_Right, dRevOpt_dL_coef,
                vec_L(i_L), vec_L(i_L+1), TotV0, TotV1, H, F, wstar, Lbase, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else {L_opt = L_opt; v_max = v_max;}
    }
    else if (i_L == N_L-1) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_L(i_L-1), vec_L(i_L), EVal_L_vec(i_L-1), EVal_L_vec(i_L) );
        double adj_L = (L_opt > Lbase) ? H : F;
        double slope_Left = coef_Left(1) - 2*adj_L*(L_opt - Lbase) - wstar + dRevOpt_dL;
        if (slope_Left < 0) {
            double TotV0 = total_value(i_L-1);
            double TotV1 = total_value(i_L);
            tuple<double, double> t_max = SolveMax_Linear1D_L_RevOpt(coef_Left, dRevOpt_dL_coef,
                vec_L(i_L-1), vec_L(i_L), TotV0, TotV1, H, F, wstar, Lbase, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp;
            }
        }
        else {L_opt = L_opt; v_max = v_max;}
    }
    else {L_opt = L_opt; v_max = v_max;}
    return tuple<double,double,int>(L_opt,v_max,L_index_max);
}

tuple<double, double> alias::SolveMax_Linear1D_L_RevOpt(const ArrayXd & coef, const double & dRevOpt_dL_coef,
    const double & x0,const double & x1, const double & TotV0,const double & TotV1, const double & H, const double & F,
    const double & wstar, const double & Lbase, const double & phi, const double & K, const ParaEst & para_est) {

    double L_opt = (TotV0 > TotV1) ? x0 : x1;
    double v_max = (TotV0 > TotV1) ? TotV0 : TotV1;

    auto eval_value = [&](const double L, const double Adj) {
        tuple<double,double,double> t_rev = RevOpt_Prod_ELur(para_est, phi, K, L);
        const double rev = get<0>(t_rev);
        return coef(0) + coef(1) * L - Adj * (L - Lbase) * (L - Lbase) - wstar * L + rev;
    };

    if (Lbase <= x0 || Lbase >= x1) {
        double Adj = H;
        if (Lbase >= x1) {Adj = F;}

        double a = dRevOpt_dL_coef;
        double b = -2*Adj;
        double c = coef(1) + 2*Adj*Lbase - wstar;
        double alpha = para_est.alpha_tilde_Lr-1;

        NewtonResult x_result = SolvePowerEqPositiveNewton(a, alpha, b, c, Lbase);
        double L_opt_temp = x_result.x;
        if (x_result.converged && L_opt_temp < x1 && L_opt_temp > x0) {
            const double v_max_temp = eval_value(L_opt_temp, Adj);
            if (v_max_temp > v_max) {
                L_opt = L_opt_temp;
                v_max = v_max_temp;
            }
        }
    }
    else {
        tuple<double,double,double> t_revLbase = RevOpt_Prod_ELur(para_est, phi, K, Lbase);
        double revLbase = get<0>(t_revLbase);
        double dRevOpt_dLur_Lbase = get<1>(t_revLbase);

        double VLbase = coef(0) + coef(1) * Lbase - wstar * Lbase + revLbase;
        if (VLbase > v_max) {
            L_opt = Lbase;
            v_max = VLbase;
        }
        double dVLbase = coef(1) - wstar + dRevOpt_dLur_Lbase;
        if (dVLbase >= 0) {
            double Adj = H;

            double a = dRevOpt_dL_coef;
            double b = -2*Adj;
            double c = coef(1) + 2*Adj*Lbase - wstar;
            double alpha = para_est.alpha_tilde_Lr-1;

            NewtonResult x_result = SolvePowerEqPositiveNewton(a, alpha, b, c, Lbase);
            double L_opt_temp = x_result.x;
            if (x_result.converged && L_opt_temp <= x1 && L_opt_temp >= Lbase) {
                const double v_max_temp = eval_value(L_opt_temp, Adj);
                if (v_max_temp > v_max) {
                    L_opt = L_opt_temp;
                    v_max = v_max_temp;
                }
            }
        } else {
            double Adj = F;

            double a = dRevOpt_dL_coef;
            double b = -2*Adj;
            double c = coef(1) + 2*Adj*Lbase - wstar;
            double alpha = para_est.alpha_tilde_Lr-1;

            NewtonResult x_result = SolvePowerEqPositiveNewton(a, alpha, b, c, Lbase);
            double L_opt_temp = x_result.x;
            if (x_result.converged && L_opt_temp <= Lbase && L_opt_temp >= x0) {
                const double v_max_temp = eval_value(L_opt_temp, Adj);
                if (v_max_temp > v_max) {
                    L_opt = L_opt_temp;
                    v_max = v_max_temp;
                }
            }
        }
    }

    return tuple<double, double>(L_opt,v_max);
}


//// Step 2. Firms decide if produce or be dormant
tuple<ArrayXd,ArrayXd,ArrayXd> alias::calEV_ProbStatus_PD(const ParaEst & para_est,
    const ArrayXd & Vprime_LurP, const ArrayXd & Vprime_LurD, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd diffV_PD = Vprime_LurP - Vprime_LurD;
    ArrayXd ProbPD_P = ArrayXd::Zero(para.N * 2);
    ArrayXd ProbPD_D = ArrayXd::Ones(para.N * 2);
    ArrayXd expF = ArrayXd::Zero(para.N * 2);

    const double sigma_PD = para_est.sigma_PD;
    const double sigma2_PD = sigma_PD * sigma_PD;
    static const boost::math::normal stdn(0.0, 1.0);

    auto worker = [&](size_t i_raw, unsigned thread_id) {
        const int i = static_cast<int>(i_raw);
        const double diff = diffV_PD(i);
        if (diff <= 0.0) { return; }

        const double mu = (i < para.N) ? para_est.F_P_PP : para_est.F_P_DP;
        const double log_diff = std::log(diff);
        const double z0 = (log_diff - mu) / sigma_PD;
        const double p = boost::math::cdf(stdn, z0);
        ProbPD_P(i) = p;
        ProbPD_D(i) = 1.0 - p;

        if (p > 0.0) {
            const double z1 = (log_diff - mu - sigma2_PD) / sigma_PD;
            const double scale = std::exp(mu + 0.5 * sigma2_PD);
            expF(i) = scale * boost::math::cdf(stdn, z1) / p;
        }
    };
    MultiThreads::simple_parallel_for(worker, para.N * 2, threadsManagement);

    ArrayXd EVal_PD_ProdDorm = (Vprime_LurP - expF) * ProbPD_P + Vprime_LurD * ProbPD_D;

    return tuple<ArrayXd,ArrayXd,ArrayXd>(EVal_PD_ProdDorm,ProbPD_P,ProbPD_D);
}


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
// //

//
// //
//
// //
// // /* Choose the targeted Luc to maximize the expected value by iteration */
// // tuple<ArrayXd,ArrayXd> alias::solveOptLur_iteration(const int & PD, const ParaEst & para_est, const ParaVec & para_vec,
// //     const ArrayXd & EVal_PD_Lur, const int & InitialPeriod, MultiThreads::Threads_Management & threadsManagement) {
// //
// //     ArrayXd TotVOpt_PD = ArrayXd::Zero(para.N*PD);
// //     ArrayXd OptLur_PD = ArrayXd::Zero(para.N*PD);
// //     const double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
// //
// //     ArrayXXd EVal_Lur_PD_mat( para.N_Lur*PD, para.N_phi * para.N_K);
// //     for (int i_PD = 0; i_PD < PD; ++i_PD) {
// //         Eigen::Map<const Eigen::ArrayXXd> EVal_Lur_mat(EVal_PD_Lur.data() + i_PD * para.N, para.N_Lur, para.N_phi*para.N_K);
// //         EVal_Lur_PD_mat(seqN(i_PD*para.N_Lur,para.N_Lur),all) = EVal_Lur_mat;
// //     }
// //
// //     auto worker = [&](size_t ii, unsigned thread_id) {
// //         //    size_t thread_id = 0;
// //         //    for (size_t ii = 0; ii < para.N*PD; ++ii) {
// //         size_t i_PD = ii / para.N;
// //         size_t i = ii - i_PD * para.N;
// //
// //         const auto EVal_Lur_mat = EVal_Lur_PD_mat(seqN(i_PD * para.N_Lur,para.N_Lur), all);
// //
// //         Pos pos;
// //         pos.i_phi = i / (para.N_K * para.N_Lur);
// //         pos.i_K = (i - pos.i_phi * para.N_K * para.N_Lur ) / para.N_Lur;
// //         pos.i_Lur = (i - pos.i_phi * para.N_K * para.N_Lur - pos.i_K * para.N_Lur ) / para.N_Lur;
// //         int i_state = static_cast<int>(pos.i_phi * para.N_K + pos.i_K);
// //
// //         //        cout << "i = " << i << "; pos.i_phi = " << pos.i_phi << "; pos.i_K = " << pos.i_K
// //         // << "; pos.i_Lur = " << pos.i_Lur << "; int i_N = " << i_N << endl;
// //
// //         double Lbase = para_vec.vec_Lur(pos.i_Lur);
// //
// //         double H_Lbase_ur = 0.0; double F_Lbase_ur = 0.0;
// //         if (InitialPeriod == 0) {
// //             H_Lbase_ur = CalHFcost_Lbase(para_est.c_H_ur, Lbase);
// //             if (Lbase < para.Lur_cutoff1) {
// //                 F_Lbase_ur = CalHFcost_Lbase(para_est.c_low_F_ur, Lbase);
// //             }
// //             else {
// //                 F_Lbase_ur = CalHFcost_Lbase(para_est.c_high_F_ur, Lbase);
// //             }
// //         }
// //
// //         ArrayXd diff_vec_L = para_vec.vec_Lur - Lbase;
// //         ArrayXd TotV = EVal_Lur_mat.col(i_state)
// //             - H_Lbase_ur * diff_vec_L.pow(2) * (diff_vec_L > 0).cast<double>()
// //             - F_Lbase_ur * diff_vec_L.pow(2) * (diff_vec_L <= 0).cast<double>()
// //             - wstar * para_vec.vec_Lur;
// //
// //         ArrayXd::Index max_Index;
// //         double Vdefault = TotV.maxCoeff(&max_Index);
// //         double xdefault =  para_vec.vec_Lur(max_Index);
// //
// //         OptLur_PD(ii) = xdefault;
// //         TotVOpt_PD(ii) = Vdefault;
// // //    }
// //     };
// //     MultiThreads::simple_parallel_for(worker, para.N*PD, threadsManagement);
// //
// //     return tuple<ArrayXd,ArrayXd>(OptLur_PD,TotVOpt_PD);
// // }
//
// //// Find the optimal solution (Lur) in Step 2

//
// double alias::solve1D_L_LinearSpline_Diff(const double & H, const double & F, const double & wstar,
//     const double & Lbase, const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const int & L_index_max) {
//
//     int N_L = vec_L.size();
//
//     ArrayXd diff_vec_L = vec_L - Lbase;
//     ArrayXd TotV = EVal_L_vec - H * diff_vec_L.pow(2) * (diff_vec_L > 0).cast<double>()
//         - F * diff_vec_L.pow(2) * (diff_vec_L <= 0).cast<double>()
//         - wstar * vec_L;
//
//     double Vdefault = TotV(L_index_max);
//     double L_opt = vec_L(L_index_max);
//     double v_max = Vdefault;
//
//     int L_index = L_index_max;
//     if (L_index_max == N_L-1) { L_index = L_index_max - 1; }
//
//     ArrayXd Lvec = vec_L.segment(L_index, 2);
//     ArrayXd coef_EVal = CalLinearspline1D_Coeff(vec_L(L_index),vec_L(L_index+1),
//         EVal_L_vec(L_index),EVal_L_vec(L_index+1));
//
//     tuple<double, double> t_max = SolveMax_Linear1D_L(coef_EVal, vec_L(L_index),
//         vec_L(L_index+1), TotV(L_index),TotV(L_index+1), H, F,
//         wstar, Lbase);
//     double L_opt_temp = get<0>(t_max);
//     double v_max_temp = get<1>(t_max);
//
//     if (v_max_temp > v_max) {
//         L_opt = L_opt_temp; v_max = v_max_temp;
//     }
//
//     return L_opt;
// }
//
// //


// //
// /////**************************************************************
// //// * Calculate the value function of the end period:
// //// * Approximation: Assume that firms stop change capital or employment after 40 period. This is for computation only
// ////**************************************************************/
// ////ArrayXd alias::calVEnd(const ParaEst & para_est, const ParaVec & para_vec, const Rev & RevOpt) {
// ////
// ////    ArrayXd VEnd = RevOpt.RevOpt_Prod - para_est.r_K*para_vec.p_K*para_vec.vec_K_Full
// ////        - para_est.w_ur*para_vec.vec_Lur_Full - para_est.w_uc*para_vec.vec_Luc_Full;
// ////
// ////    return VEnd;
// ////}
// ////
// //
// //
// ////
// /////**************************************************************
// ////* Testing the concacvity of value functions
// ////**************************************************************/
// ////void alias::testConcavity_L(const ArrayXd & Y, const ArrayXd & vec_Luc, const ArrayXd & vec_Lur) {
// ////
// ////    Eigen::Map<const Eigen::ArrayXXd> Y_mat(Y.data(),para.N_Luc,para.N_Lur);
// ////
// ////    ArrayXXd dY_dLuc(para.N_Luc,para.N_Lur);
// ////    ArrayXXd dY_dLur(para.N_Luc,para.N_Lur);
// ////
// ////    for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
// ////        dY_dLuc.col(i_Lur) = construct_dv_hermite_Schumaker(vec_Luc,Y_mat.col(i_Lur));
// ////    }
// //////    cout << "Y_mat = \n" << Y_mat << endl;
// //////    cout << "dY_dLuc = \n" << dY_dLuc << endl;
// ////    ArrayXXi dY_dLuc_mat = ( dY_dLuc(seqN(1,para.N_Luc-1),all)
// ////                           > dY_dLuc(seqN(0,para.N_Luc-1),all)+1e-5 ).cast<int>();
// ////    if (dY_dLuc_mat.sum() > 0) {
// ////        cout << "!!!! non concave dY_dLuc_mat" << endl;
// ////        cout << "Y_mat = \n" << Y_mat << endl;
// ////        cout << "dY_dLuc = \n" << dY_dLuc << endl;
// ////        cout << "diff dY_dLuc = \n" << dY_dLuc(seqN(1,para.N_Luc-1),all)
// ////            - dY_dLuc(seqN(0,para.N_Luc-1),all) << endl;
// ////        cout << "dY_dLuc_mat = \n" << dY_dLuc_mat << endl;
// ////        cout << "max dY_dLuc(seqN(1,para.N_Luc-1),all) - dY_dLuc(seqN(0,para.N_Luc-1),all) = "
// ////             << (dY_dLuc(seqN(1,para.N_Luc-1),all)
// ////                 - dY_dLuc(seqN(0,para.N_Luc-1),all)).maxCoeff() << endl;
// ////        throw runtime_error("1034");
// ////    }
// ////
// ////
// ////    for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
// ////        dY_dLur.row(i_Luc) = construct_dv_hermite_Schumaker(vec_Lur,Y_mat.row(i_Luc)).transpose();
// ////    }
// ////
// ////    ArrayXXi dY_dLur_mat = ( dY_dLur(all,seqN(1,para.N_Lur-1))
// ////                             > dY_dLur(all,seqN(0,para.N_Lur-1))+1e-5 ).cast<int>();
// ////    if (dY_dLur_mat.sum() > 0) {
// ////        cout << "!!!! non concave dY_dLur_mat" << endl;
// ////        cout << "Y_mat = \n" << Y_mat << endl;
// ////        cout << "dY_dLur = \n" << dY_dLur << endl;
// ////        cout << "diff dY_dLur = \n" << dY_dLur(all,seqN(1,para.N_Lur-1))
// ////            - dY_dLur(all,seqN(0,para.N_Lur-1)) << endl;
// ////        cout << "dY_dLur_mat = \n" << dY_dLur_mat << endl;
// ////        cout << "max dY_dLur(all,seqN(1,para.N_Lur-1)) - dY_dLur(all,seqN(0,para.N_Lur-1)) = "
// ////             << (dY_dLur(all,seqN(1,para.N_Lur-1))
// ////                 - dY_dLur(all,seqN(0,para.N_Lur-1))).maxCoeff() << endl;
// ////        throw runtime_error("1041");
// ////    }
// ////}
// ////
// ////void alias::testConcavity_K(const ArrayXXd & Y, const ArrayXd & vec_K) {
// ////
// ////    int N = Y.rows();
// ////    ArrayXXd dY_dK(N,para.N_K);
// ////    for (size_t i = 0; i < N; ++i) {
// ////        dY_dK.row(i) = construct_dv_hermite_Schumaker(vec_K,Y.row(i)).transpose();
// ////    }
// ////
// ////    ArrayXXi dY_dK_mat = ( dY_dK(all,seqN(1,para.N_K-1))
// ////            > dY_dK(all,seqN(0,para.N_K-1))+1e-5 ).cast<int>();
// ////
// //////    cout << "Y = \n" << Y << endl;
// //////    cout << "dY_dK = \n" << dY_dK << endl;
// ////    if (dY_dK_mat.sum() > 0) {
// ////        cout << "!!!! non concave dY_dK_mat" << endl;
// ////        cout << "Y = \n" << Y << endl;
// ////        cout << "dY_dK = \n" << dY_dK << endl;
// ////        cout << "diff dY_dK = \n" << dY_dK(all,seqN(1,para.N_K-1))
// ////            - dY_dK(all,seqN(0,para.N_K-1)) << endl;
// ////        cout << "dY_dK_mat = \n" << dY_dK_mat << endl;
// ////        cout << "max dY_dK(all,seqN(1,para.N_K-1)) - dY_dK(all,seqN(0,para.N_K-1)) = "
// ////            << (dY_dK(all,seqN(1,para.N_K-1))
// ////            - dY_dK(all,seqN(0,para.N_K-1))).maxCoeff() << endl;
// ////        throw runtime_error("1059");
// ////    }
// ////}
// ////
// ////void alias::testConcavity_L_loop(const ArrayXd & Eval, const ArrayXd & vec_Luc, const ArrayXd & vec_Lur) {
// ////
// ////    cout << "******************* testConcavity_L_loop ************************" << endl;
// ////    Eigen::Map<const Eigen::ArrayXXd> EVal_mat(Eval.data(),para.N_Luc*para.N_Lur,para.N_phi*para.N_K);
// ////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
// ////        for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
// //////            cout << "i_phi = " << i_phi << "; i_K = " << i_K << endl;
// ////            ArrayXd temp1 = EVal_mat.col(i_phi*para.N_K+i_K);
// ////            testConcavity_L(EVal_mat.col(i_phi*para.N_K+i_K), vec_Luc, vec_Lur);
// ////        }
// ////    }
// ////    cout << "******************* testConcavity_L_loop::Success!! ************************" << endl;
// ////}
// ////
// ////void alias::testConcavity_K_loop(const ArrayXd & Eval, const ArrayXd & vec_K) {
// ////
// ////    cout << "******************* testConcavity_K_loop ************************" << endl;
// ////    Eigen::Map<const Eigen::ArrayXXd> Eval_mat(Eval.data(),para.N_Luc*para.N_Lur,para.N_phi*para.N_K);
// ////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
// ////        for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
// //////            cout << "i_phi = " << i_phi  << "i_Lur = " << i_Lur << endl;
// ////            testConcavity_K(Eval_mat(seqN(i_Lur * para.N_Luc, para.N_Luc),
// ////                seqN(i_phi * para.N_K, para.N_K)), vec_K);
// ////        }
// ////    }
// ////    cout << "******************* testConcavity_K_loop::Success!! ************************" << endl;
// ////}
// ////
// ////void alias::testMonotonicity_by_K_loop(const ArrayXd & OptKL) {
// ////
// ////    cout << "******************* testMonotonicity_by_K_loop ************************" << endl;
// ////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
// ////        for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
// ////            for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
// ////                for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
// ////                    if (i_K > 0) {
// ////                        int i_plus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
// ////                                     i_Lur * para.N_Luc + i_Luc;
// ////                        int i_minus = i_phi * para.N_K * para.N_Lur * para.N_Luc + (i_K - 1) * para.N_Lur * para.N_Luc +
// ////                                      i_Lur * para.N_Luc + i_Luc;
// ////
// ////                        if (OptKL(i_plus) < OptKL(i_minus) - 1e-9) {
// ////                            cout << "K::: i_phi = " << i_phi << "; i_K = " << i_K << "; i_Luc = " << i_Luc << "; i_Lur = " << i_Lur
// ////                                 << "; OptKL(i_plus) = " << OptKL(i_plus) << "; OptKL(i_minus) = " << OptKL(i_minus)
// ////                                 << "; OptKL(i_plus) - OptKL(i_minus)= " << OptKL(i_plus) - OptKL(i_minus) << endl;
// //////                            throw runtime_error("340");
// ////                        }
// ////                    }
// ////                }
// ////            }
// ////        }
// ////    }
// ////    cout << "******************* testMonotonicity_by_K_loop::Success!! ************************" << endl;
// ////}
// ////
// ////void alias::testMonotonicity_by_Lur_loop(const ArrayXd & OptKL) {
// ////
// ////    cout << "******************* testMonotonicity_by_Lur_loop ************************" << endl;
// ////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
// ////        for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
// ////            for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
// ////                for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
// ////                    if (i_Lur > 0) {
// ////                        int i_plus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
// ////                                     i_Lur * para.N_Luc + i_Luc;
// ////                        int i_minus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
// ////                                      (i_Lur - 1) * para.N_Luc + i_Luc;
// ////
// ////                        if (OptKL(i_plus) < OptKL(i_minus) - 1e-9) {
// ////                            cout << "Lur:::i_K = " << i_K << "; i_Luc = " << i_Luc << "; i_Lur = " << i_Lur
// ////                                 << "; OptKL(i_plus) = " << OptKL(i_plus) << "; OptKL(i_minus) = " << OptKL(i_minus)
// ////                                 << "; OptKL(i_plus) - OptK(i_minus)= " << OptKL(i_plus) - OptKL(i_minus) << endl;
// //////                            throw runtime_error("340");
// ////                        }
// ////                    }
// ////                }
// ////            }
// ////        }
// ////    }
// ////    cout << "******************* testMonotonicity_by_Lur_loop::Success!! ************************" << endl;
// ////
// ////}
// ////
// ////void alias::testMonotonicity_by_Luc_loop(const ArrayXd & OptKL) {
// ////
// ////    cout << "******************* testMonotonicity_by_Luc_loop ************************" << endl;
// ////    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) {
// ////        for (size_t i_K = 0; i_K < para.N_K; ++i_K) {
// ////            for (size_t i_Luc = 0; i_Luc < para.N_Luc; ++i_Luc) {
// ////                for (size_t i_Lur = 0; i_Lur < para.N_Lur; ++i_Lur) {
// ////                    if (i_Luc > 0) {
// ////                        int i_plus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
// ////                                     i_Lur * para.N_Luc + i_Luc;
// ////                        int i_minus = i_phi * para.N_K * para.N_Lur * para.N_Luc + i_K * para.N_Lur * para.N_Luc +
// ////                                      i_Lur * para.N_Luc + i_Luc - 1;
// ////
// ////                        if (OptKL(i_plus) < OptKL(i_minus) - 1e-9) {
// ////                            cout << "Luc:::i_K = " << i_K << "; i_Luc = " << i_Luc << "; i_Lur = " << i_Lur
// ////                                 << "; OptKL(i_plus) = " << OptKL(i_plus) << "; OptKL(i_minus) = " << OptKL(i_minus)
// ////                                 << "; OptKL(i_plus) - OptKL(i_minus)= " << OptKL(i_plus) - OptKL(i_minus) << endl;
// //////                            throw runtime_error("340");
// ////                        }
// ////                    }
// ////                }
// ////            }
// ////        }
// ////    }
// ////    cout << "******************* testMonotonicity_by_Luc_loop::Success!! ************************" << endl;
// ////
// ////}
// ////
// ////
// ////
// //////
// //////tuple<ArrayXd,ArrayXd,ArrayXd> alias::calEV_ProbStatus_SE_lognormal(const ParaEst & para_est, const ParaVec & para_vec,
// //////    const ArrayXd & ResVal_E, const ArrayXd & EVal, MultiThreads::Threads_Management & threadsManagement) {
// //////
// //////    ArrayXd Ones_d = ArrayXd::Ones(para.N);
// //////    /***********************************************************************/
// //////    ArrayXd dEVal_SE_sigma = EVal - ResVal_E;
// //////
// //////    ArrayXd Vprime_up(para.N);
// //////    ArrayXd Prob_S(para.N);
// //////    ArrayXd Prob_E(para.N);
// ////////    auto worker = [&](size_t i, unsigned thread_id) {
// //////    for (size_t i = 0; i < para.N; ++i) {
// //////
// //////        if (dEVal_SE_sigma(i) <= 0) {
// //////            Prob_S(i) = 0.0;
// //////            Prob_E(i) = 1.0;
// //////            Vprime_up(i) = ResVal_E(i) + exp(para_est.F_E);
// //////        } else {
// //////            Prob_S(i) = normCDF((log(dEVal_SE_sigma(i)) - para_est.F_E) / para_est.sigma_SE, 0.0, 1.0);
// //////            Prob_E(i) = normCDF(-(log(dEVal_SE_sigma(i)) - para_est.F_E) / para_est.sigma_SE, 0.0, 1.0);;
// //////            Vprime_up(i) = EVal(i) * Prob_S(i) + (ResVal_E(i) + exp(para_est.F_E)) * Prob_E(i);
// //////        }
// ////////        cout << "dEVal_SE_sigma = " << dEVal_SE_sigma(i) << "; Prob_S(i) = " << Prob_S(i) << endl;
// //////    }
// ////////    };
// ////////throw runtime_error("205");
// //////    return tuple<ArrayXd,ArrayXd,ArrayXd>(Vprime_up,Prob_S,Prob_E);
// //////}
// ////
// ////
// ////
