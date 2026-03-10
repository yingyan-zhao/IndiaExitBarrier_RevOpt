#include "Basics.h"
#include <random>
//
using namespace std;
using namespace alias;
using namespace Eigen;

/**************************************************************
* Full Parameter Definition
**************************************************************/
ParaEst alias::constructParaEst(const ArrayXd & theta,const int & GoodState, const int & LaborIntensive) {
    ParaEst para_est;

    // cout << "theta = " << theta.transpose() << endl;
    // cout << "GoodState = " << GoodState << "; LaborIntensive = " << LaborIntensive << endl;

    if (GoodState == 1) { para_est.gamma0 = theta[para.dim3+0]; }
    else { para_est.gamma0 = theta[para.dim3+1]; }
    para_est.gamma1 = theta[para.dim3+2];
    para_est.sigma_phi_eps = theta[para.dim3+3];
    cout << "para_est.gamma0 = " << para_est.gamma0 << "; para_est.gamma1 = " << para_est.gamma1
        << "; para_est.sigma_phi_eps = " << para_est.sigma_phi_eps << endl;

    if (GoodState == 1 & LaborIntensive == 1) {
        para_est.alpha_M = para.Mshr_G_L;
        para_est.alpha_K = para.Kshr_G_L;
        para_est.alpha_Lr = para.Lurshr_G_L;
        para_est.alpha_Lc = para.Lucshr_G_L;

        // cout << "sigma_tilde = " << sigma_tilde << endl;
        // cout << "para_est.alpha_tilde_K + para_est.alpha_tilde_L = " <<  para_est.alpha_tilde_K + para_est.alpha_tilde_L << endl;
    }
    else if (GoodState == 0 & LaborIntensive == 1) {
        para_est.alpha_M = para.Mshr_B_L;
        para_est.alpha_K = para.Kshr_B_L;
        para_est.alpha_Lr = para.Lurshr_B_L;
        para_est.alpha_Lc = para.Lucshr_B_L;
    }
    else if (GoodState == 1 & LaborIntensive == 0) {
        para_est.alpha_M = para.Mshr_G_K;
        para_est.alpha_K = para.Kshr_G_K;
        para_est.alpha_Lr = para.Lurshr_G_K;
        para_est.alpha_Lc = para.Lucshr_G_K;
    }
    else {
        para_est.alpha_M = para.Mshr_B_K;
        para_est.alpha_K = para.Kshr_B_K;
        para_est.alpha_Lr = para.Lurshr_B_K;
        para_est.alpha_Lc = para.Lucshr_B_K;
    }
    para_est.sigma_tilde = (para.sigma-1.0)/para.sigma / (1.0-(para_est.alpha_M+para_est.alpha_Lc)*(para.sigma-1)/para.sigma);
    para_est.alpha_tilde_Lr = para_est.alpha_Lr * para_est.sigma_tilde;
    para_est.alpha_tilde_K = para_est.alpha_K * para_est.sigma_tilde;

    para_est.F_P_PP = theta(para.dim3+para.dim1);
    para_est.F_P_DP = theta(para.dim3+para.dim1+1);
    para_est.sigma_PD = theta(para.dim3+para.dim1+2);
    cout << "para_est.F_P_PP = " << para_est.F_P_PP << "; para_est.F_P_DP = " << para_est.F_P_DP
        << "; para_est.sigma_PD = " << para_est.sigma_PD << endl;

    int nn = 0;
    para_est.c_H_ur = theta(nn+0);

    if (GoodState == 1) {
        para_est.c_high_F_ur = theta(nn+1);
        para_est.c_low_F_ur = theta(nn+2);
    }
    else {
        para_est.c_high_F_ur = theta(nn+3);
        para_est.c_low_F_ur = theta(nn+4);
    }

    para_est.sigma_Lerror_ur = theta(nn+5);

    if (GoodState == 1) {
        para_est.F_E = theta(nn+6);
        para_est.w_ur = para.w_ur_G_L;
        para_est.w_uc = para.w_uc_G_L;
    }
    else if (GoodState == 0) {
        para_est.F_E = theta(nn+7);
        para_est.w_ur = para.w_ur_B_L;
        para_est.w_uc = para.w_uc_B_L;
    }

    para_est.F_E_c_FK = theta(nn+8);
    para_est.sigma_SE = theta(nn+9);

    cout << "para_est.c_H_ur = " << para_est.c_H_ur << "; para_est.c_high_F_ur = " << para_est.c_high_F_ur << "; para_est.c_low_F_ur = " << para_est.c_low_F_ur
        << "; para_est.sigma_Lerror_ur = " << para_est.sigma_Lerror_ur
        << "; para_est.F_E = " << para_est.F_E << "; para_est.F_E_c_FK = " << para_est.F_E_c_FK << "; para_est.sigma_SE = " << para_est.sigma_SE << endl;
    return para_est;
}

/**************************************************************
* Construct the grid for the value function and policy functions
**************************************************************/
//// Construct the grid for the value function and policy functions
ParaVec alias::constructParaFull(const ParaEst & para_est) {
    ParaVec para_vec;

    para_vec.prod_mu0 = para_est.gamma0 / (1.0 - para_est.gamma1);
    para_vec.prod_sigma0 = para_est.sigma_phi_eps / sqrt(1-pow(para_est.gamma1,2));
    cout << "para_vec.prod_mu0 = " << para_vec.prod_mu0 << "; para_vec.prod_sigma0 = " << para_vec.prod_sigma0 << endl;

    //// theta
    para_vec.vec_lnphi = para.vec_norm_phi*para_vec.prod_sigma0 + para_vec.prod_mu0*ArrayXd::Ones(para.N_phi);
    para_vec.vec_phi = exp(para_vec.vec_lnphi);

    para_vec.vec_phi_Full = KroneckerProduct(para_vec.vec_phi,ArrayXd::Ones(para.N_K*para.N_Lur));
    para_vec.vec_phi_Full_PD = ArrayXd::Zero(para.N_phi * para.N_K * para.N_Lur * 2);
    para_vec.vec_phi_Full_PD << para_vec.vec_phi_Full,para_vec.vec_phi_Full;

    // initial distribution
    para_vec.dist0 = NormalGridProbabilities(para_vec.vec_lnphi, para_vec.prod_mu0, para_vec.prod_sigma0);

    //// transition of theta
    para_vec.tran_phi = GenerateTransitionMatrixAR1(para_vec.vec_lnphi, para_est.gamma0, para_est.gamma1, para_est.sigma_phi_eps);
    cout << "para_vec.tran_phi = \n" << para_vec.tran_phi << endl;

    //// grids of Luc
    tuple<ArrayXd,ArrayXd> t = genGridStates();
    para_vec.vec_K =  get<0>(t);
    para_vec.vec_Lur =  get<1>(t);
    cout << "para_vec.vec_Lur = " << para_vec.vec_Lur.transpose() << endl;

    ArrayXd temp_K = KroneckerProduct(para_vec.vec_K,ArrayXd::Ones(para.N_Lur));
    para_vec.vec_K_Full = KroneckerProduct(ArrayXd::Ones(para.N_phi),temp_K);
    para_vec.vec_K_Full_PD = ArrayXd::Zero(para.N_phi * para.N_K * para.N_Lur * 2);
    para_vec.vec_K_Full_PD << para_vec.vec_K_Full,para_vec.vec_K_Full;

    ArrayXd temp_Lur = para_vec.vec_Lur;
    para_vec.vec_Lur_Full = KroneckerProduct(ArrayXd::Ones(para.N_phi*para.N_K),temp_Lur);
    para_vec.vec_Lur_Full_PD = ArrayXd::Zero(para.N_phi * para.N_K * para.N_Lur * 2);
    para_vec.vec_Lur_Full_PD << para_vec.vec_Lur_Full,para_vec.vec_Lur_Full;

    return para_vec;
}

//// get the Grids of K and L_ur
tuple<ArrayXd,ArrayXd> alias::genGridStates() {

    ArrayXd K = exp( ArrayXd::LinSpaced(para.N_K,log(0.001),log(3500)));

    int N_temp = para.N_Lur/2;
    ArrayXd L_ur_temp1 = ArrayXd::LinSpaced(N_temp,log(1),log(para.Lur_cutoff1));
    ArrayXd L_ur_temp2 = ArrayXd::LinSpaced(para.N_Lur - N_temp + 1,log(para.Lur_cutoff1),log(1500));
    ArrayXd L_ur(para.N_Lur);
    L_ur << exp(L_ur_temp1),exp(L_ur_temp2.segment(1,para.N_Lur - N_temp));

    return tuple<ArrayXd,ArrayXd>(K,L_ur);
}
//

// //
// ///**************************************************************
// //* Parameter Partition: Part 1: Production function, productivity evolution
// //**************************************************************/
// //ParaEst1 alias::constructParaEst_Part1(const std::vector<double> & theta1, const int & GoodState,
// //    const int & LaborIntensive) {
// //
// //    ParaEst1 para_est;
// //
// //    if (GoodState == 1) {
// //        para_est.gamma0 = theta1[0];
// //    }
// //    else
// //    {
// //        para_est.gamma0 = theta1[1];
// //    }
// //
// //    para_est.gamma1 = theta1[2];
// //    para_est.sigma_phi_eps = theta1[3];
// //
// //    if (GoodState == 1 & LaborIntensive == 1) {
// //        double alpha_M = para.alpha_M_G_L;
// //        double sigma_tilde = (para.sigma-1.0)/para.sigma / (1.0-alpha_M*(para.sigma-1)/para.sigma);
// //        para_est.alpha_tilde_K = sigma_tilde*para.Kshr_G_L;
// //        para_est.alpha_tilde_L = sigma_tilde*para.Lshr_G_L;
// //
// //        para_est.alpha_K = para.Kshr_G_L;
// //        para_est.alpha_L = para.Lshr_G_L;
// //
// //        para_est.alpha_Lr = para.Lurshr_G_L;
// //        para_est.alpha_Lc = 1 - para_est.alpha_Lr;
// //    }
// //    else if (GoodState == 0 & LaborIntensive == 1) {
// //        double alpha_M = para.alpha_M_B_L;
// //        double sigma_tilde = (para.sigma-1.0)/para.sigma / (1.0-alpha_M*(para.sigma-1)/para.sigma);
// //        para_est.alpha_tilde_K = sigma_tilde*para.Kshr_B_L;
// //        para_est.alpha_tilde_L = sigma_tilde*para.Lshr_B_L;
// //
// //        para_est.alpha_K = para.Kshr_B_L;
// //        para_est.alpha_L = para.Lshr_B_L;
// //
// //        para_est.alpha_Lr = para.Lurshr_B_L;
// //        para_est.alpha_Lc = 1 - para_est.alpha_Lr;
// //    }
// //    else if (GoodState == 1 & LaborIntensive == 0) {
// //        double alpha_M = para.alpha_M_G_K;
// //        double sigma_tilde = (para.sigma-1.0)/para.sigma / (1.0-alpha_M*(para.sigma-1)/para.sigma);
// //        para_est.alpha_tilde_K = sigma_tilde*para.Kshr_G_K;
// //        para_est.alpha_tilde_L = sigma_tilde*para.Lshr_G_K;
// //
// //        para_est.alpha_K = para.Kshr_G_K;
// //        para_est.alpha_L = para.Lshr_G_K;
// //
// //        para_est.alpha_Lr = para.Lurshr_G_K;
// //        para_est.alpha_Lc = 1 - para_est.alpha_Lr;
// //    }
// //    else if (GoodState == 0 & LaborIntensive == 0) {
// //        double alpha_M = para.alpha_M_B_K;
// //        double sigma_tilde = (para.sigma-1.0)/para.sigma / (1.0-alpha_M*(para.sigma-1)/para.sigma);
// //        para_est.alpha_tilde_K = sigma_tilde*para.Kshr_B_K;
// //        para_est.alpha_tilde_L = sigma_tilde*para.Lshr_B_K;
// //
// //        para_est.alpha_K = para.Kshr_B_K;
// //        para_est.alpha_L = para.Lshr_B_K;
// //
// //        para_est.alpha_Lr = para.Lurshr_B_K;
// //        para_est.alpha_Lc = 1 - para_est.alpha_Lr;
// //    }
// //    else
// //    {
// //        throw std::runtime_error("para_est not initiated!");
// //    }
// //    return para_est;
// //}
// //
// ///**************************************************************
// //* Parameter Partition: Part 2: Fixed cost of production
// //**************************************************************/
// //ParaEst2 alias::constructParaEst_Part2(const std::vector<double> & theta2) {
// //    ParaEst2 para_est;
// //
// //    para_est.F_P_PP = theta2[0];
// //    para_est.F_P_DP = theta2[1];
// //    para_est.sigma_PD = theta2[2];
// //
// //    return para_est;
// //}
// //


//

// ///**************************************************************
// //* ProbMissing: no log
// //**************************************************************/
// //double alias::ProbMissing(const double & lagK, const double & lagLur, const double & lagLuc, const double & lagRev,
// //    const int & lagCapital_miss, const int & lagEmploy_ur_miss, const int & lagEmploy_uc_miss,
// //    const int & lagState_D_NonMiss, const int & lagState_DP_Miss, const int & lagState_P_NonMiss, const int & lagyr_DP,
// //    const int & state, const int & industry, const int & year, const int & gap_yr,
// //    const int & GoodState, const int & LaborIntensive) {
// //
// //    double lagL = lagLuc + lagLur;
// //
// //    double KK = log(lagK+1);
// //    double KK_2 = pow(KK,2);
// //    double KK_3 = pow(KK,3);
// //
// //    double LLur = log(lagLur+1);
// //    double LLur_2 = pow(LLur,2);
// //    double LLur_3 = pow(LLur,3);
// //
// //    double LLuc = log(lagLuc+1);
// //    double LLuc_2 = pow(LLuc,2);
// //    double LLuc_3 = pow(LLuc,3);
// //
// ////    double RR = log(lagRev+1);
// ////    double RR_2 = pow(RR,2);
// ////    double RR_3 = pow(RR,3);
// //    double RR = lagRev;
// //    double RR_2 = pow(RR,2);
// //    double RR_3 = pow(RR,3);
// //    double temp;
// //
// //    if (lagState_DP_Miss == 1) {
// //
// //        ArrayXd coef_state(18);
// //        coef_state << 0, .1684732, .0384604, .0815654, .1188896, .1217691, .0015481, .0982617, .0130511, .0235926,
// //                .0367733, -.0196394, .1855818, .1733737, .0518239, .0538207, .0842114, .1517285;
// //
// //        ArrayXd coef_Ind(22);
// //        coef_Ind << 0, -.0435508, -.0371459, -.0894159, -.0187104, -.0315826, -.0361409, -.0614661, -.1221007, -.0524369,
// //                .006573, .0113729, -.0151031, -.0327795, -.0858529, -.3160268, -.1044498, -.2094656, -.2616107, -.0414388,
// //                -.1542498, -.093947;
// //
// //        ArrayXd coef_Yr(20);
// //        coef_Yr << 0,0,0, .1308678, .2018673, -.1544541, .0901728, -.0714419, -.0032234, .2872135,
// //                .0206356, -.0121229, .0234932, .0060042, -.2292563, .0789542, -.055694, -.3597297, -.8246444, 0;
// //
// //        double cons = .5205318;
// //
// ////        ArrayXd coef_yrDP(18);
// ////        coef_yrDP << 0,0, -.0309584, .0090439, .0106635, .0085234, .0007484, -.0038698, .0366868, .0646716,
// ////                .0371894, .0475036, .0248395, .0416352, .042862, .0625848, .0255384, .1390512, 0;
// //
// ////        double coef_K = -.0132032;
// ////        double coef_Ksq = .0012106;
// ////        double coef_K3 = .0000542;
// ////
// ////        double coef_Lur = .3439815;
// ////        double coef_Lursq = -.1065881;
// ////        double coef_Lur3 = .0080506;
// ////
// ////        double coef_Luc = .0718589;
// ////        double coef_Lucsq = -.02036;
// ////        double coef_Luc3 = -.0000104;
// ////
// ////        double coef_R = .0718589;
// ////        double coef_Rsq = -.02036;
// ////        double coef_R3 = -.0000104;
// //
// //        temp =  coef_state(state) + coef_Ind(industry) + coef_Yr(year) + cons;
// //    }
// //    else if (lagState_D_NonMiss == 1) {
// //
// //        ArrayXd coef_state(18);
// //        coef_state << 0, .231846, -.2074491, .1877177, .4068985, .2014489, -.2053168, .1948572, -.1743698, -.0471848,
// //                -.083257, .1409839, .3230155, .279519, .2105926, -.0637654, .3114099, .2180873;
// //
// //        ArrayXd coef_Ind(22);
// //        coef_Ind << 0, -.1596978, -.3363319, -.5464795, -.2720573, .1240994, -.0856572, -.0009448, -.1053783, .0273086,
// //                .1974939, .156131, .0279396, .0354309, .0375032, -.3689668, -.0926224, -.3903284, -.4291488, -.0118739,
// //                -.2264286, -.3288074;
// //
// //        ArrayXd coef_Yr(20);
// //        coef_Yr << 0,0,0, -.1165843, -.0785034, -.2168914, .4194869, .13495, .0533095, -.0042959,
// //                .4502382, .3940486, .6669637, .6007829, 0, -.4843299, -.9706835, -.6605941, -.9692766, 0;
// //
// //        double cons = -.3257867;
// //
// ////        double coef_K = .0217803;
// ////        double coef_Ksq = -.0403614;
// ////        double coef_K3 = .0055843;
// ////
// ////        double coef_Lur = 1.47291;
// ////        double coef_Lursq = -.4362544;
// ////        double coef_Lur3 = .0337552;
// ////
// ////        double coef_Luc = .1736863;
// ////        double coef_Lucsq = -.053351;
// ////        double coef_Luc3 = .0016209;
// //
// ////        double coef_R = -.0070573;
// ////        double coef_Rsq = .0000139;
// ////        double coef_R3 = -6.68e-09;
// //
// ////        temp = coef_K * KK + coef_Ksq * KK_2 + coef_K3 * KK_3
// ////            + coef_Lur * LLur + coef_Lursq * LLur_2 + coef_Lur3 * LLur_3
// ////            + coef_Luc * LLuc + coef_Lucsq * LLuc_2 + coef_Luc3 * LLuc_3
// ////            + coef_state(state) + coef_Ind(industry) + coef_Yr(year) + cons;
// //        temp = coef_state(state) + coef_Ind(industry) + coef_Yr(year) + cons;
// //    }
// //    else if (lagState_P_NonMiss == 1) {
// //
// //        ArrayXd coef_state(18);
// //        coef_state << 0, .3020706, -.0852999, .3431962, .4142807, .3035544, -.2365343, .2541077, -.1361434, -.0797838,
// //                -.0441754, .1008614, .4665753, .4401032, .3393725, .2091128, .1065, .2183392;
// //
// //        ArrayXd coef_Ind(22);
// //        coef_Ind << 0, -.3102829, -.2136383, -.2980829, -.2232223, .0346608, .0947834, -.0095096, .0364475, .0118376,
// //                .2025745, .1549559, .1268542, .0540658, .0495195, -.4026941, .0466926, -.2006017, -.4083589, -.0323983,
// //                -.2762689, -.3949701;
// //
// //        ArrayXd coef_Yr(20);
// //        coef_Yr << 0,0, -.2388566, -.4234949, -.4560446, -.6847614, -.2167606, -.1271036, -.240545, -.5509735,
// //                -.2639105, -.3474888, -.4870865, -.4499028, -.6498845, -1.140154, -1.281872, -1.182892, -1.594506, 0;
// //
// //        double cons = .2323881;
// ////
// ////        double coef_K = -.0009083;
// ////        double coef_Ksq = 1.75e-06;
// ////        double coef_K3 = -7.14e-10;
// ////
// ////        double coef_Lur = -.0149415;
// ////        double coef_Lursq = .0000317;
// ////        double coef_Lur3 = -1.88e-08;
// //////
// ////        double coef_Luc = -.0124904;
// ////        double coef_Lucsq = .0000301;
// ////        double coef_Luc3 = -2.09e-08;
// //////
// //        double coef_R = -.0122483;
// //        double coef_Rsq = .0000201;
// //        double coef_R3 = -8.89e-09;
// //
// ////        temp = coef_K * KK + coef_Ksq * KK_2 + coef_K3 * KK_3
// ////            + coef_Lur * LLur + coef_Lursq * LLur_2 + coef_Lur3 * LLur_3
// ////            + coef_Luc * LLuc + coef_Lucsq * LLuc_2 + coef_Luc3 * LLuc_3
// ////            + coef_state(state) + coef_Ind(industry) + coef_Yr(year) + cons;
// //        temp = coef_R * RR + coef_Rsq * RR_2 + coef_R3 * RR_3
// //               + coef_state(state) + coef_Ind(industry) + coef_Yr(year) + cons;
// //    }
// //
// //    double temp_final = min(temp,7.5);
// //    temp_final = max(temp_final,-7.5);
// //
// //    double prob = normCDF(temp_final, 0, 1);
// //    return prob;
// //}
// ////
// /////**************************************************************
// ////// generate 10 initial values for parameter guess; Not a good programming format. Random numbers are generated differently every time
// ////// I record a result of one try and then uses these numbers as my initial guess.
// ////**************************************************************/
// ////ArrayXXd alias::GenerateInitialGuessParaRandom() {
// ////    std::vector<double> lb(para.dim3); std::vector<double> ub(para.dim3);
// ////    lb[0] = 0; ub[0] = 2000;
// ////    lb[1] = 0; ub[1] = 2000;
// ////
// ////    lb[2] = 0; ub[2] = 2000;
// ////    lb[3] = 0; ub[3] = 2000;
// ////    lb[4] = 0; ub[4] = 2000;
// ////    lb[5] = 0; ub[5] = 2000;
// ////    lb[6] = 0; ub[6] = 2000;
// ////
// ////    lb[7] = 0; ub[7] = 2000;
// ////    lb[8] = 0; ub[8] = 2000;
// ////    lb[9] = 0; ub[9] = 2000;
// ////
// ////    lb[10] = 0.1; ub[10] = 2;
// ////    lb[11] = 0.1; ub[11] = 2;
// ////    lb[12] = 0.1; ub[12] = 2;
// ////
// ////    lb[13] = -1000; ub[13] = 1000;
// ////    lb[14] = -1000; ub[14] = 1000;
// ////    lb[15] = -1; ub[15] = 20;
// ////    lb[16] = 10; ub[16] = 1000;
// ////
// ////    std::default_random_engine generator;
// ////    std::uniform_real_distribution<double> distribution(0.0,1.0);
// ////    ArrayXXd theta_RandomInitial(para.dim3,100);
// ////    for (size_t n = 0; n < 100; ++n) {
// ////        for (size_t j = 0; j < para.dim3; ++j) {
// ////            double number = distribution(generator);
// //////            cout << "number = " << number << endl;
// ////            theta_RandomInitial(j,n) = lb[j] + number * (ub[j]-lb[j]);
// ////        }
// ////    }
// ////
// ////    ArrayXXd theta_RandomInitial_select(para.dim3,10);
// ////    theta_RandomInitial_select = theta_RandomInitial(all,seqN(50,10));
// ////    cout << "theta_RandomInitial_select = " << theta_RandomInitial_select.transpose() << endl;
// //////    theta_RandomInitial_select =
// ////////    165.695,1734.02,95.4231,436,214.284,433.189,234.701,4.6839,73.0999,188.557,0.392284,0.562564,0.741573,-4.94866,-506.188,2.34148,468.115
// //////    f = 2.30221, error = 2
// //////            --------------------------------------
// //////    para_est = 181.086, 1734.02, 145.423, 368.875, 335.284, 349.189, 257.201, 6.39484, 68.428, 117.526, 0.520092, 0.609673, 0.839814, -75.5112, -524.313, -0.777799, 563.428, ;
// ////
// //////    141.204,1256.53,172.944,341.123,38.5633,308.938,459.815,96.1575,217.896,45.1865,0.742439,0.839686,0.882827,-463.223,-522.822,7.51294,,67.363
// //////    f = 2.30303, error = 0.968527
// //////            --------------------------------------
// //////    para_est = 182.227, 1256.54, 188.546, 186.631, 200.641, 252.934, 279.608, 50.8424, 67.9429, 50.6623, 0.519495, 0.614937, 0.838952, -307.246, -443.45, -0.467955, 475.683, ;
// ////
// ////
// //////    471.791,1843.69,382.766,,314.23,261.379,309.654,225.081,26.4898,231.206,203.855,0.312579,0.408445,0.401305,-21.7791,-97.6034,5.74385,82.5234
// //////    f = 2.30689, error = 4
// //////            --------------------------------------
// //////    para_est = 401.791, 1843.69, 383.641, 235.23, 222.379, 229.654, 213.081, 6.22417, 69.206, 139.855, 0.521372, 0.609875, 0.841571, -59.7791, -157.603, -0.573953, 416.523, ;
// ////
// ////
// //////    236.391,1592.12,425.625,216.294,425.68,368.668,279.011,29.8349,39.2397,37.2754,0.303744,0.761219,0.118109,-702.122,-644.828,5.2496,12.8638
// //////    f = 2.3003, error = 8
// //////            --------------------------------------
// //////    para_est = 203.891, 1592.12, 426.031, 241.294, 207.68, 331.668, 282.511, 4.8349, 99.2397, 131.275, 0.519687, 0.609012, 0.835764, -645.122, -654.828, -0.728645, 629.364, ;
// //////
// //////    316.86,1310.08,145.557,183.784,,314.01,154.475,43.1124,49.2609,212.907,29.2005,0.722439,0.225717,0.603862,-683.238,-223.118,3.81805,153.158
// //////    f = 2.3013, error = 2
// //////            --------------------------------------
// //////    para_est = 275.86, 1310.08, 167.557, 199.909, 169.51, 218.475, 255.112, 7.5109, 60.8445, 95.2005, 0.520664, 0.610376, 0.835618, -583.301, -293.118, -0.473026, 492.158, ;
// ////
// ////
// //////    214.249,,1222.4,463.171,272.305,260.181,332.672,416.221,39.8707,173.096,111.397,0.516444,0.690335,0.416208,-431.344,-466.489,8.38435,,,59.01
// //////    f = 2.30154, error = 4
// //////            --------------------------------------
// //////    para_est = 141.249, 1222.41, 462.343, 290.305, 266.181, 315.672, 321.221, 4.3707, 141.096, 145.397, 0.5201, 0.608889, 0.836903, -231.344, -442.989, -0.450136, 589.885, ;
// ////
// //////    291.412,1395.24,390.701,365.128,335.027,479.146,491.055,,78.283,66.0733,155.758,0.498284,0.688665,0.351458,-675.272,-521.916,19.0441,166.247
// //////    f = 2.30247, error = 0.91419
// //////            --------------------------------------
// //////    para_est = 296.678, 1683.28, 153.815, 286.614, 281.24, 248.057, 286.548, 4.99791, 115.777, 88.227, 0.520436, 0.60976, 0.840115, -223.571, -731.139, -0.619191, 581.707, ;
// ////
// ////
// //////    367.428,1683.28,127.784,294.739,,401.24,258.604,397.618,,34.994,139.781,208.227,0.44392,0.843794,0.25351,-175.509,-772.108,2.41685,297.848
// //////    f = 2.30247, error = 0.91419
// //////            --------------------------------------
// //////    para_est = 296.912, 1683.28, 153.94, 286.739, 281.24, 248.182, 286.532, 4.99815, 115.773, 88.226, 0.520436, 0.60976, 0.840115, -223.571, -731.163, -0.619191, 581.723, ;
// ////
// //////    10.4981,1747.97,221.667,253.837,313.701,442.798,,495.81,25.1659,11.3174,86.1474,0.155603,0.543941,0.306581,-793.72,-137.103,11.7828,102.549
// //////    f = 2.30372, error = 2
// //////            --------------------------------------
// //////    para_est = 1014.5, 1811.97, 220.917, 198.837, 109.701, 406.798, 381.81, 5.2909, 67.3174, 162.147, 0.524843, 0.609012, 0.835671, -1051.72, -225.103, -0.868758, 622.549, ;
// ////
// //////    449.654,1963.32,285.625,208.442,219.636,237.765,197.065,92.4639,135.026,127.049,0.526747,0.410587,0.375975,-420.641,-675.59,18.9007,421.395
// //////    f = 2.30139, error = 0.0562989
// //////            --------------------------------------
// //////    para_est = 355.654, 1963.33, 290.765, 228.08, 192.636, 202.418, 182.525, 4.4639, 71.026, 129.727, 0.521282, 0.608913, 0.839076, -466.766, -706.682, -0.512928, 525.395, ;
// //////
// ////////    theta_RandomInitial_select =   635.49  1468.05  348.666  1738.77  833.812   1727.3  917.147  74.4222  579.116  1506.48 0.794174  1.19859  1.62374  987.628 -265.469  2.34148  935.579
// //////    -obj = 2.32265
// //////    para_est = 589.49, 1468.06, 350.041, 1646.77, 882.812, 1695.3, 841.147, 10.4222, 467.116, 1390.48, 0.522728, 0.608746, 0.840537, 952.628, -247.469, -1, 1287.58, ;
// ////
// ////////    535.525  513.065  665.079  1351.52  116.585  1220.16  1835.98  1922.37  1742.14  354.911  1.62579  1.85676  1.95921 -158.059 -307.056  7.51294  125.897
// //////    f = 2.30935, error = 4
// //////            --------------------------------------
// //////    para_est = 256.025, 513.659, 665.423, 1055.52, 660.585, 1181.72, 1497.48, 1922.04, 1002.14, 151.817, 0.52084, 0.609873, 0.839952, -236.059, -248.415, -1, 1706.9, ;
// ////
// //////    1884.86  1687.38  1521.49  1241.75  1026.04  1223.08   877.88  514.946  1849.05  1629.36 0.604875 0.832558 0.815599  945.552  755.991  5.74385  156.527
// //////    924.044  1184.23  1696.43  842.015  1696.65  1463.95  1098.01  582.524  307.146  291.369 0.583891   1.6704 0.143009 -755.306 -612.071   5.2496  15.7861
// //////    1252.49  620.165  553.293  709.322  1240.86  589.695  135.152  974.968  1702.07   226.51  1.57829 0.398579  1.29667 -708.095  442.206  3.81805  299.237
// //////    833.671  444.804  1849.68  1070.63  1021.15  1317.03  1658.04  785.266   1382.3  886.722  1.08906  1.50204 0.850994 -78.3589 -166.222  8.38435   109.02
// //////    1148.62  790.471  1553.88   1449.5  1326.64  1914.88  1963.49  1561.27  522.677  1243.04  1.04592  1.49808 0.697214 -688.179 -304.789  19.0441  325.684
// //////    1458.89  1366.57  480.753   1162.2   1596.9  1014.71  1582.11  686.746  1114.71  1664.48 0.916809  1.86651 0.464586  561.228 -930.271  2.41685  591.571
// //////    2.03312  1495.95  863.948  995.252   1239.6  1766.52   1982.9  488.199   82.871  683.915 0.232058  1.15436 0.590629 -984.301  657.244  11.7828  196.987
// //////    1794.51  1926.64     1125  809.967  855.658  929.652  763.532  1847.75  1076.52  1012.44  1.11353 0.837644 0.755442 -51.6031 -688.976  18.9007  841.186
// ////
// //////    throw runtime_error("39");
// ////    return theta_RandomInitial_select;
// ////}
// ////
// //
// //
// ////
// //
