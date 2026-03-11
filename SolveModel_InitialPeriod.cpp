#include "SolveModel_InitialPeriod.h"
#include <limits>

using namespace Eigen;
using namespace std;
using namespace alias;
// using namespace Ipopt_Wrapper;

/**************************************************************
* Solve the value/policy function in the first period
**************************************************************/
EquStateV0 alias::solveV0(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & Vprime_PD,
    MultiThreads::Threads_Management & threadsManagement) {

    EquStateV0 EquV0;
    //// the expected value function with the productivity transition
    int PD = 1;
    ArrayXd EVprime_PD = para.deltaV * CalEVprime_mat(PD,Vprime_PD.segment(0,para.N), para_vec.tran_phi);

    /// Step 3. solve L_ur
    EquV0.EVal_Lur0_error = CalEVal_Lur_error(PD, para_vec, EVprime_PD, para_est.sigma_Lerror_ur); // With a chosen target, firms' expected value depends on the shocks on employment
    ArrayXXd RevOpt_Prod_mat_val = RevOpt_Prod_mat(para_est, para_vec);

    tuple<ArrayXd,ArrayXd,ArrayXd> t_Lur = solveOptLur0(para_est, para_vec, EquV0.EVal_Lur0_error, RevOpt_Prod_mat_val, threadsManagement);
    EquV0.OptLur0 = get<0>(t_Lur); // dim = N_phi * N_K
    EquV0.EVal_Lur0 = get<1>(t_Lur); // dim = N_phi * N_K
    EquV0.EVpart_Lur0 = get<2>(t_Lur); // dim = N_phi * N_K

    /// Step 3. solve K
    tuple<ArrayXd,ArrayXd> t_K = solveOptK0(para_est, para_vec, EquV0.EVpart_Lur0, EquV0.OptLur0, threadsManagement);
    EquV0.OptK0 = get<0>(t_K); // dim = N_phi
    EquV0.EVal0 = get<1>(t_K); // dim = N_phi

    return EquV0;
}


/*****************************************************************
///* Choose the targeted Lur to maximize the expected value
*****************************************************************/
tuple<ArrayXd,ArrayXd,ArrayXd> alias::solveOptLur0(const ParaEst & para_est, const ParaVec & para_vec,
    const ArrayXd & EVal_PD, const ArrayXXd & RevOpt_Prod_mat, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd OptLur_P = ArrayXd::Zero(para.N_phi*para.N_K);
    ArrayXd TotVOpt_P = ArrayXd::Zero(para.N_phi*para.N_K);
    ArrayXd EVLur_P = ArrayXd::Zero(para.N_phi*para.N_K);

    const double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
    const int N_K = para.N_K;

    Eigen::Map<const Eigen::ArrayXXd> EVal_Lur_mat(EVal_PD.data(), para.N_Lur, para.N_phi*para.N_K);

    auto worker = [&](size_t i_raw, unsigned thread_id) {
        const int i = static_cast<int>(i_raw);
        const int i_phi = i / N_K;
        const int i_K = i - i_phi * N_K;
        const int i_state = i_phi * para.N_K + i_K;

        tuple<double,double,double,int> t_V_P = solve1D_L0_LinearSpline_RevOpt(
            wstar, para_vec.vec_phi(i_phi), para_vec.vec_K(i_K), para_est, para_vec.vec_Lur,
            EVal_Lur_mat.col(i_state), RevOpt_Prod_mat.col(i_state));
        OptLur_P(i) = get<0>(t_V_P);
        TotVOpt_P(i) = get<1>(t_V_P);
        EVLur_P(i) = get<2>(t_V_P);
    };
    MultiThreads::simple_parallel_for(worker, para.N_phi * para.N_K, threadsManagement);

    return tuple<ArrayXd,ArrayXd,ArrayXd>(OptLur_P,TotVOpt_P,EVLur_P);
}

tuple<double,double,double,int> alias::solve1D_L0_LinearSpline_RevOpt(const double & wstar,
    const double & phi, const double & K, const ParaEst & para_est,
    const ArrayXd & vec_L, const ArrayXd & EVal_L_vec, const ArrayXd & RevOpt_Prod_vec) {

    const int N_L = static_cast<int>(vec_L.size());

    auto total_value = [&](int idx) {
        const double L = vec_L(idx);
        return RevOpt_Prod_vec(idx) + EVal_L_vec(idx) - wstar * L;
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
    double EV_max = EVal_L_vec(i_L);

    tuple<double,double,double> t_RevOpt_Prod = RevOpt_Prod_ELur(para_est, phi, K, L_opt);
    double dRevOpt_dL = get<1>(t_RevOpt_Prod);
    double dRevOpt_dL_coef = get<2>(t_RevOpt_Prod);

    if (L_index_max >= 1 && L_index_max <= N_L-2) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_L(i_L-1), vec_L(i_L), EVal_L_vec(i_L-1), EVal_L_vec(i_L) );
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_L(i_L), vec_L(i_L+1), EVal_L_vec(i_L), EVal_L_vec(i_L+1) );

        double slope_Left = coef_Left(1) - wstar + dRevOpt_dL;
        double slope_Right = coef_Right(1) - wstar + dRevOpt_dL;

        if (slope_Left>0 and slope_Right>0) { // maximum point lies in the right interval
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double, double> t_max = SolveMax_Linear1D_L0_RevOpt(coef_Right, dRevOpt_dL_coef,
                vec_L(i_L), vec_L(i_L+1), TotV0, TotV1,
                wstar, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            double EV_max_temp = get<2>(t_max);
            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp; EV_max = EV_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right<0){ // maximum point lies in the left interval
            double TotV0 = total_value(i_L-1);
            double TotV1 = total_value(i_L);
            tuple<double, double, double> t_max = SolveMax_Linear1D_L0_RevOpt(coef_Left, dRevOpt_dL_coef,
                vec_L(i_L-1), vec_L(i_L), TotV0, TotV1,
                wstar, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            double EV_max_temp = get<2>(t_max);
            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp; EV_max = EV_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right>0) { // this case should not happen
            // calculate the right interval first.
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double, double> t_maxRight = SolveMax_Linear1D_L0_RevOpt(coef_Right, dRevOpt_dL_coef,
                vec_L(i_L), vec_L(i_L+1), TotV0, TotV1,
                wstar, phi, K, para_est);
            double L_opt_temp_right = get<0>(t_maxRight);
            double v_max_temp_right = get<1>(t_maxRight);
            double EV_max_temp_right = get<2>(t_maxRight);
            if (v_max_temp_right > v_max) {
                L_opt = L_opt_temp_right; v_max = v_max_temp_right; EV_max = EV_max_temp_right;
            }

            // calculate the left interval then, and compare.
            TotV0 = total_value(i_L-1);
            TotV1 = total_value(i_L);
            tuple<double, double, double> t_maxLeft = SolveMax_Linear1D_L0_RevOpt(coef_Left, dRevOpt_dL_coef,
                vec_L(i_L-1), vec_L(i_L), TotV0, TotV1,
                wstar, phi, K, para_est);
            double L_opt_temp_left = get<0>(t_maxLeft);
            double v_max_temp_left = get<1>(t_maxLeft);
            double EV_max_temp_left = get<2>(t_maxLeft);
            if (v_max_temp_left > v_max) {
                L_opt = L_opt_temp_left; v_max = v_max_temp_left; EV_max = EV_max_temp_left;
            }
        }
    }
    else if (i_L == 0) {
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_L(i_L), vec_L(i_L+1), EVal_L_vec(i_L), EVal_L_vec(i_L+1) );
        double slope_Right = coef_Right(1) - wstar + dRevOpt_dL;
        if (slope_Right>0) {
            double TotV0 = total_value(i_L);
            double TotV1 = total_value(i_L+1);
            tuple<double, double, double> t_max = SolveMax_Linear1D_L0_RevOpt(coef_Right, dRevOpt_dL_coef,
                vec_L(i_L), vec_L(i_L+1), TotV0, TotV1,
                wstar, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            double EV_max_temp = get<2>(t_max);
            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp; EV_max = EV_max_temp;
            }
        }
    }
    else if (i_L == N_L-1) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_L(i_L-1), vec_L(i_L), EVal_L_vec(i_L-1), EVal_L_vec(i_L) );
        double slope_Left = coef_Left(1) - wstar + dRevOpt_dL;
        if (slope_Left < 0) {
            double TotV0 = total_value(i_L-1);
            double TotV1 = total_value(i_L);
            tuple<double, double, double> t_max = SolveMax_Linear1D_L0_RevOpt(coef_Left, dRevOpt_dL_coef,
                vec_L(i_L-1), vec_L(i_L), TotV0, TotV1,
                wstar, phi, K, para_est);
            double L_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            double EV_max_temp = get<2>(t_max);

            if (v_max_temp > v_max) {
                L_opt = L_opt_temp; v_max = v_max_temp; EV_max = EV_max_temp;
            }
        }
    }
    return tuple<double,double,double,int>(L_opt,v_max,EV_max,L_index_max);
}

tuple<double, double, double> alias::SolveMax_Linear1D_L0_RevOpt(const ArrayXd & coef, const double & dRevOpt_dL_coef,
    const double & x0,const double & x1, const double & TotV0,const double & TotV1,
    const double & wstar, const double & phi, const double & K, const ParaEst & para_est) {

    double L_opt = (TotV0 > TotV1) ? x0 : x1;;
    double v_max = (TotV0 > TotV1) ? TotV0 : TotV1;

    double a = dRevOpt_dL_coef;
    double c = coef(1) - wstar;
    double alpha = para_est.alpha_tilde_Lr-1;
    double L_opt_temp = solve_axalpha_plus_c_eq0(a, alpha, c);

    if (L_opt_temp > x0 && L_opt_temp < x1) {
        tuple<double,double,double> t_rev = RevOpt_Prod_ELur(para_est, phi, K, L_opt_temp);
        double rev = get<0>(t_rev);
        double v_max_temp = coef(0) + coef(1) * L_opt_temp - wstar * L_opt_temp + rev;

        if (v_max_temp > v_max) {
            L_opt = L_opt_temp; v_max = v_max_temp;
        }
    }
    double EV_max = coef(0) + coef(1)*L_opt;

    return tuple<double,double,double>(L_opt,v_max,EV_max);
}

/**************************************************************
* At period 0, choose the optimal K
**************************************************************/
/* Choose the targeted K to maximize the expected value */
tuple<ArrayXd,ArrayXd> alias::solveOptK0(const ParaEst & para_est, const ParaVec & para_vec, const ArrayXd & EVpart_Lur0,
    const ArrayXd & OptLur0, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd TotVOpt0 = ArrayXd::Zero(para.N_phi);
    ArrayXd OptK0 = ArrayXd::Zero(para.N_phi);

    Eigen::Map<const Eigen::ArrayXXd> EVpart_Lur0_mat(EVpart_Lur0.data(), para.N_K, para.N_phi);
    Eigen::Map<const Eigen::ArrayXXd> OptLur0_mat(OptLur0.data(), para.N_K, para.N_phi);

    auto worker = [&](size_t ii, unsigned thread_id) {
        const int i_phi = static_cast<int>(ii);

        ArrayXd EVpart_Lur0_vec = EVpart_Lur0_mat.col(i_phi);
        ArrayXd OptLur0_vec = OptLur0_mat.col(i_phi);
        tuple<double,double,int> t_K = solve1D_K_LinearSpline(
            para_vec.vec_K, EVpart_Lur0_vec, OptLur0_vec, para_est, para_vec.vec_phi(i_phi));
        OptK0(i_phi) = get<0>(t_K);
        TotVOpt0(i_phi) = get<1>(t_K);
    };
    MultiThreads::simple_parallel_for(worker, para.N_phi, threadsManagement);

    return tuple<ArrayXd,ArrayXd>(OptK0,TotVOpt0);
}


//// Find the optimal solution (K) in Step 2
tuple<double,double,int> alias::solve1D_K_LinearSpline(const ArrayXd & vec_K, const ArrayXd & EVpart_Lur0_vec,
    const ArrayXd & OptLur0_vec, const ParaEst & para_est, const double & phi) {

    double p_K = 1;
    const int N_K = static_cast<int>(vec_K.size());
    ArrayXd total_value(N_K);
    ArrayXd dRevOpt_dK(N_K);
    for (int ii = 0; ii < N_K; ++ii) {

        tuple<double,double,double> t_rev = RevOpt_Prod_K(para_est, phi, vec_K(ii), OptLur0_vec(ii));
        double RevOpt_K = get<0>(t_rev);
        dRevOpt_dK(ii) = get<1>(t_rev);
        total_value(ii) = RevOpt_K + EVpart_Lur0_vec(ii) - p_K * vec_K(ii);
    }

    int i_K = 0;
    double v_max = total_value(0);
    for (int i = 1; i < N_K; ++i) {
        const double v = total_value(i);
        if (v > v_max) {
            v_max = v;
            i_K = i;
        }
    }

    double K_opt = vec_K(i_K);
    int K_index_max = i_K;
    double dRevOpt_dK_max = dRevOpt_dK(i_K);

    if (K_index_max >= 1 && K_index_max <= N_K-2) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_K(i_K-1), vec_K(i_K), EVpart_Lur0_vec(i_K-1), EVpart_Lur0_vec(i_K) );
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_K(i_K), vec_K(i_K+1), EVpart_Lur0_vec(i_K), EVpart_Lur0_vec(i_K+1) );
        double slope_Left = coef_Left(1) + dRevOpt_dK_max - p_K;
        double slope_Right = coef_Right(1) + dRevOpt_dK_max - p_K;

        if (slope_Left>0 and slope_Right>0) { // maximum point lies in the right interval
            double TotV0 = total_value(i_K);
            double TotV1 = total_value(i_K+1);
            tuple<double, double> t_max = SolveMax_Linear1D_K0_RevOpt(coef_Right, vec_K(i_K), vec_K(i_K+1),
                TotV0, TotV1, EVpart_Lur0_vec(i_K), EVpart_Lur0_vec(i_K+1),
                OptLur0_vec(i_K), OptLur0_vec(i_K+1), phi, para_est);
            double K_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                K_opt = K_opt_temp; v_max = v_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right<0){ // maximum point lies in the left interval
            double TotV0 = total_value(i_K-1);
            double TotV1 = total_value(i_K);
            tuple<double, double> t_max = SolveMax_Linear1D_K0_RevOpt(coef_Left, vec_K(i_K-1), vec_K(i_K),
                TotV0, TotV1, EVpart_Lur0_vec(i_K-1), EVpart_Lur0_vec(i_K),
                OptLur0_vec(i_K-1), OptLur0_vec(i_K), phi, para_est);
            double K_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);

            if (v_max_temp > v_max) {
                K_opt = K_opt_temp; v_max = v_max_temp;
            }
        }
        else if (slope_Left<0 and slope_Right>0) { // this case should not happen
            // calculate the right interval first.
            double TotV0 = total_value(i_K);
            double TotV1 = total_value(i_K+1);
            tuple<double, double> t_maxRight = SolveMax_Linear1D_K0_RevOpt(coef_Right, vec_K(i_K), vec_K(i_K+1),
                TotV0, TotV1, EVpart_Lur0_vec(i_K), EVpart_Lur0_vec(i_K+1),
                OptLur0_vec(i_K), OptLur0_vec(i_K+1), phi, para_est);
            double K_opt_temp_right = get<0>(t_maxRight);
            double v_max_temp_right = get<1>(t_maxRight);

            if (v_max_temp_right > v_max) {
                K_opt = K_opt_temp_right; v_max = v_max_temp_right;
            }

            // calculate the left interval then, and compare.
            TotV0 = total_value(i_K-1);
            TotV1 = total_value(i_K);
            tuple<double, double> t_max_Left = SolveMax_Linear1D_K0_RevOpt(coef_Left, vec_K(i_K-1), vec_K(i_K),
                TotV0, TotV1, EVpart_Lur0_vec(i_K-1), EVpart_Lur0_vec(i_K),
                OptLur0_vec(i_K-1), OptLur0_vec(i_K), phi, para_est);
            double K_opt_temp_left = get<0>(t_max_Left);
            double v_max_temp_left = get<1>(t_max_Left);
            if (v_max_temp_left > v_max) {
                K_opt = K_opt_temp_left;
                v_max = v_max_temp_left;
            }
        }
    }
    else if (i_K == 0) {
        ArrayXd coef_Right = LinearInterpolation1D_Coeff( vec_K(i_K), vec_K(i_K+1), EVpart_Lur0_vec(i_K), EVpart_Lur0_vec(i_K+1) );
        double TotV0 = total_value(i_K);
        double TotV1 = total_value(i_K+1);
        tuple<double, double> t_max = SolveMax_Linear1D_K0_RevOpt(coef_Right, vec_K(i_K), vec_K(i_K+1),
            TotV0, TotV1, EVpart_Lur0_vec(i_K), EVpart_Lur0_vec(i_K+1),
            OptLur0_vec(i_K), OptLur0_vec(i_K+1), phi, para_est);
        double K_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);

        if (v_max_temp > v_max) {
            K_opt = K_opt_temp; v_max = v_max_temp;
        }
    }
    else if (i_K == N_K-1) {
        ArrayXd coef_Left = LinearInterpolation1D_Coeff( vec_K(i_K-1), vec_K(i_K), EVpart_Lur0_vec(i_K-1), EVpart_Lur0_vec(i_K) );
        double TotV0 = total_value(i_K-1);
        double TotV1 = total_value(i_K);
        tuple<double, double> t_max = SolveMax_Linear1D_K0_RevOpt(coef_Left, vec_K(i_K-1), vec_K(i_K),
                TotV0, TotV1, EVpart_Lur0_vec(i_K-1), EVpart_Lur0_vec(i_K),
                OptLur0_vec(i_K-1), OptLur0_vec(i_K), phi, para_est);
        double K_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);

        if (v_max_temp > v_max) {
            K_opt = K_opt_temp; v_max = v_max_temp;
        }
    }

    return tuple<double,double,int>(K_opt,v_max,K_index_max);
}

tuple<double, double> alias::SolveMax_Linear1D_K0_RevOpt(const ArrayXd & coef, const double & x0,const double & x1,
    const double & TotV0, const double & TotV1, const double & EVpart_Lur0_0, const double & EVpart_Lur0_1,
    const double & OptLur0, const double & OptLur1, const double & phi, const ParaEst & para_est) {

    double p_K = 1.0;

    double K_opt = (TotV0 > TotV1) ? x0 : x1;;
    double v_max = (TotV0 > TotV1) ? TotV0 : TotV1;

    double logL_alpha0 = para_est.alpha_tilde_Lr * log(OptLur0)
        + 0.5 * para_est.alpha_tilde_Lr * para_est.alpha_tilde_Lr * para_est.sigma_Lerror_ur * para_est.sigma_Lerror_ur;
    double Rev_K0 = pow(para_est.PI, para_est.sigma_tilde) * phi * exp(logL_alpha0);

    double logL_alpha1 = para_est.alpha_tilde_Lr * log(OptLur1)
        + 0.5 * para_est.alpha_tilde_Lr * para_est.alpha_tilde_Lr * para_est.sigma_Lerror_ur * para_est.sigma_Lerror_ur;
    double Rev_K1 = pow(para_est.PI, para_est.sigma_tilde) * phi * exp(logL_alpha1);

    ArrayXd coef_b = LinearInterpolation1D_Coeff(x0, x1, Rev_K0, Rev_K1);
    ArrayXd coef_c = LinearInterpolation1D_Coeff(x0, x1, EVpart_Lur0_0, EVpart_Lur0_1);

    //Solve: d1*k^a + d0*k^(a-1) + d2 = 0, with k > 0
    double d1 = para_est.alpha_tilde_K * coef_b(1);
    double d0 = para_est.alpha_tilde_K * coef_b(0);
    double d2 = coef_c(1) - p_K;
    double a = para_est.alpha_tilde_K;
    RootResult opt_sol = solve_k_equation_hybrid(d1, d0, d2, a);
    double K_opt_temp = opt_sol.k;

    if (K_opt_temp > x0 && K_opt_temp < x1) {
        double rev = (coef_b(0) + coef_b(1) * K_opt_temp) * pow(K_opt_temp, para_est.alpha_tilde_K);
        double v_max_temp = coef_c(0) + coef_c(1)*K_opt_temp + rev - p_K*K_opt_temp;

        if (v_max_temp > v_max) {
            K_opt = K_opt_temp; v_max = v_max_temp;
        }
    }
    return tuple<double,double>(K_opt,v_max);
}


/**************************************************************
* Solve the entry cost
**************************************************************/
double alias::SolveFEntry_StatusQuoEqu(const ParaVec & para_vec, const EquStateV0 & EquV0) {

    double F_Entry;
    //// calculate entry cost f^e
    F_Entry = 0.0;
    for (size_t k = 0; k < para.N_phi; ++k) {
        F_Entry = F_Entry + para_vec.dist0(k) * EquV0.EVal0(k);
    }

    return F_Entry;
}
