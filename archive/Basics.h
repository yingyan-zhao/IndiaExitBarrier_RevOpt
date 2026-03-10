#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#include "Auxillaries.h"
//#include <unsupported/Eigen/MatrixFunctions>
//#include "Interpolation_spline.h"
////
namespace alias {
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;
//
    /**************************************************************
    * Pre calibrated Parameters
    **************************************************************/
    struct ConstPara{
        double sigma;
        double deltaV; //discount factor

        // Number of grid for theta and alpha
        int N_phi;
        int N_plus;
        int N_Lur;
        int N_K;
        int N;
        int N_KL;

        // small vs big
        double Lur_cutoff1;

        ArrayXd vec_norm_phi;
        ArrayXd vec_norm_phi_cut;

        int SimTbar;
        int EstTbar;

        int dim1; // number of parameters in group 1
        int dim2; // number of parameters in group 2
        int dim3; // number of parameters in group 3

        double Mshr_G_L;
        double Mshr_B_L;
        double Mshr_G_K;
        double Mshr_B_K;

        double Kshr_G_L;
        double Kshr_B_L;
        double Kshr_G_K;
        double Kshr_B_K;

        double Lucshr_G_L;
        double Lucshr_B_L;
        double Lucshr_G_K;
        double Lucshr_B_K;

        double Lurshr_G_L;
        double Lurshr_B_L;
        double Lurshr_G_K;
        double Lurshr_B_K;

        double w_ur_G_L;
        double w_uc_G_L;
        double w_ur_B_L;
        double w_uc_B_L;
//
//        double IncShrManu;
//        double tau_GB;
//        double tau_BG;
//        double tau_ld;
//
//        double sigma_Lrc;
//        double sigma_LK;

        ConstPara(){
            sigma = 3.94; // reference: For the markup we use μ = 1.34, which is the median markup across all Indian sectors reported by (De Loecker et al., 2016). This markup implies an elasticity of substitution across suppliers σ = 3.94.
            // sigma = 3.5;

            deltaV = 0.9; // discount factor

            // Number of grid for theta and alpha
            N_phi = 7;
            N_plus = 5;

            N_Lur = N_phi+N_plus+1;
            N_K = N_phi+N_plus;

            N = N_phi * N_K * N_Lur;
            N_KL = N_K * N_Lur;

            // Lur_cutoff
            Lur_cutoff1 = 100;

            tuple<ArrayXd,ArrayXd> t1 = GenerateNormalGrid(N_phi);
            vec_norm_phi = get<0>(t1);
            vec_norm_phi_cut = get<1>(t1);

            dim1 = 4;  // number of parameters in group 1
            dim2 = 3;  // number of parameters in group 2

            dim3 = 10; // number of parameters in group 3

            // number of periods
            SimTbar = 20; // number of years for simulation 1999-2018
            EstTbar = 11; // number of years for estimation 1999-2018

            // Share of intermediate inputs
            Mshr_G_L = .5684316;
            Mshr_B_L = .5684316;
            Mshr_G_K = .6460927;
            Mshr_B_K = .6460927;

            Kshr_G_L = .2450455;
            Kshr_B_L = .2450455;
            Kshr_G_K = .2055914;
            Kshr_B_K = .2055914;

            Lurshr_G_L = 0.1;
            Lurshr_B_L = 0.1;
            Lurshr_G_K = 0.1;
            Lurshr_B_K = 0.1;

            Lucshr_G_L = 1 - Lurshr_G_L - Kshr_G_L - Mshr_G_L;
            Lucshr_B_L = 1 - Lurshr_B_L - Kshr_B_L - Mshr_B_L;
            Lucshr_G_K = 1 - Lurshr_G_K - Kshr_G_K - Mshr_G_K;
            Lucshr_B_K = 1 - Lurshr_B_K - Kshr_B_K - Mshr_B_K;

            w_ur_G_L = 0.04485519;
            w_ur_B_L = 0.03934232;
            w_uc_G_L = 0.1203235;
            w_uc_B_L = 0.1078408;

//            IncShrManu = 0.15; //https://economictimes.indiatimes.com/manufacturings-share-in-gdp-growth-rises-to-22-in-04-05/articleshow/1019843.cms?from=mdr

            // tau_GB = 2.0;
            // tau_BG = 2.0;
            // tau_ld = 2.0;

            // tau_GB = 1.5;
            // tau_BG = 1.5;
            // tau_ld = 1.5;
            //
//            tau_GB = 2.0;
//            tau_BG = 2.0;
//            tau_ld = 2.0;
//
//            sigma_Lrc = 1.1;
//            sigma_LK = 1.5;
        }
    };
    const ConstPara para;
////
//    struct SegmentShr {
//        double shr_GoodState_LaborInt;
//        double shr_BadState_LaborInt;
//        double shr_GoodState_CapitalInt;
//        double shr_BadState_CapitalInt;
//    };
//
    // /**************************************************************
    // * Parameter Partition: Part 1: Production function, productivity evolution
    // **************************************************************/
    // struct ParaEst1 {
    //     double gamma0;
    //     double gamma1;
    //     double sigma_phi_eps;
    //
    //     double alpha_tilde_K;
    //     double alpha_tilde_L;
    //
    //     double alpha_K;
    //     double alpha_L;
    //
    //     double alpha_Lc;
    //     double alpha_Lr;
    //
    //     double alpha_M;
    // };
    // ParaEst1 constructParaEst_Part1(const std::vector<double> & theta1, const int & GoodState, const int & LaborIntensive);

//    /**************************************************************
//    * Parameter Partition: Part 2: Fixed cost of production
//    **************************************************************/
//    struct ParaEst2 {
//        double F_P_PP;
//        double F_P_DP;
////        double F_P_Kslope;
//        double sigma_PD;
//    };
//    ParaEst2 constructParaEst_Part2(const std::vector<double> & theta2);
//
    /**************************************************************
    * Parameter: the full set of parameters
    **************************************************************/
    struct ParaEst{

        double gamma0;
        double gamma1;
        double sigma_phi_eps;

        double alpha_M;
        double alpha_K;
        double alpha_Lc;

        double alpha_Lr;
        double alpha_tilde_Lr;

        double sigma_tilde;

        double F_P_PP;
        double F_P_DP;
        double sigma_PD;

        double c_H_ur;
        double c_high_F_ur;
        double c_low_F_ur;

        double sigma_Lerror_ur;

        double F_E;
        double F_E_c_FK;

        double sigma_SE;

        double w_ur;
        double w_uc;
        double PI;
    };
    ParaEst constructParaEst(const ArrayXd & theta,const int & GoodState, const int & LaborIntensive);

    /**************************************************************
    * Data Structure to store the value function and policy functions
    **************************************************************/
    //// Equilibrium value/policy functions for period t >= 1
    struct EquStateV {
        ArrayXd OptLur_P;
        ArrayXd OptLur_D;

        ArrayXd EVal_PD_OptLurP;
        ArrayXd EVal_PD_OptLurD;

        ArrayXd EVal_PD_SE;
        ArrayXd EVal_PD;

        ArrayXd ProbPD_P;
        ArrayXd ProbPD_D;

        ArrayXd ProbPD_S;
        ArrayXd ProbPD_E;
    };

    //// Within each period, firms (Step1) exit or stay (step2) choose Employ_ur (step3) choose production or dormancy.
    //// EquStateVmat includes the value function at each interim steps.
    struct EquStateVmat {
        ArrayXXd EVal_PD_Lur_mat;
        ArrayXXd EVal_PD_SE_mat;
    };

    /**************************************************************
    * Construct grids for value function and policy function
    **************************************************************/
    //// The grid for the value function and policy functions
    struct ParaVec{
        double prod_mu0;
        double prod_sigma0;

        ArrayXd vec_lnphi;
        ArrayXd vec_lnphi_cut;
        ArrayXd vec_phi;
        ArrayXXd tran_phi;

        ArrayXd vec_Lur;
        ArrayXd vec_K;

        ArrayXd vec_phi_Full;
        ArrayXd vec_Lur_Full;
        ArrayXd vec_K_Full;

        ArrayXd vec_phi_Full_PD;
        ArrayXd vec_Lur_Full_PD;
        ArrayXd vec_K_Full_PD;

        ArrayXd dist0;
    };

    //// Construct the grid for the value function and policy functions
    ParaVec constructParaFull(const ParaEst & para_est);
    //// get the Grids of L_uc and L_ur
    tuple<ArrayXd,ArrayXd> genGridStates();

//    /**************************************************************
//    * Probability of Missing
//    **************************************************************/
//    double ProbMissing(const double & lagK, const double & lagLur, const double & lagLuc, const double & lagRev,
//        const int & lagCapital_miss, const int & lagEmploy_ur_miss, const int & lagEmploy_uc_miss,
//        const int & lagState_D_NonMiss, const int & lagState_DP_Miss, const int & lagState_P_NonMiss,
//        const int & lagyr_DP,
//        const int & state, const int & industry, const int & year, const int & gap_yr,
//        const int & GoodState, const int & LaborIntensive);
//
    /**************************************************************
    * Value function in the initial period
    **************************************************************/
    //// Equilibrium value/policy functions for period t == 0
    struct EquStateV0 {
        ArrayXd EVal_PD_K0;
        ArrayXd EVal_PD_Lur0;

        ArrayXd OptK_PD0;
        ArrayXd OptLur_PD0;

        ArrayXd EVal_PD0;

        ArrayXd ProbPD_P0;
        ArrayXd ProbPD_D0;
    };
//
    //// Within each period, firms (Step1) entry or not (step2) choose capital (step3) choose Employ_ur (step4) choose Employ_uc.
    //// DiffEVal0 includes the value function at each interim steps.
    struct EquStateV0mat {
        ArrayXXd EVal_PD_K0_mat;
        ArrayXXd EVal_PD_Lur0_mat;
        ArrayXXd EVal_PD0_mat;
    };
//
    //// The equilibrium state of the initial period
    struct EquState0 {
        double F_Entry;
        double FirmMass;
    };
////
//
////    /**************************************************************
////    * Industry and State specific parameters
////    **************************************************************/
////    struct IndustryPara {
////        double alpha_M;
////        double w_ur;
////        double w_uc;
////    };
////
//
////    /**************************************************************
////     * randomly generate initial guess for parameters
////    **************************************************************/
////    ArrayXXd GenerateInitialGuessParaRandom();
////
//
////
//
//
}
