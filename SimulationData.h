#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "RandomDraw.h"

// #include "multithread_loop.h"
// //#include "Auxillaries.h"
// //
// ////#include "Basics.h"
// //
// ////
// ////
// //////#include "SolveModel.h"
// //////#include "SolveModel_V0_FEntry.h"
// ////
// ////
namespace alias {
    // using namespace std;
    // using ArrayXd = Eigen::ArrayXd;
    // using ArrayXXd = Eigen::ArrayXXd;
    // using ArrayXi = Eigen::ArrayXi;
    // using ArrayXXi = Eigen::ArrayXXi;

    /**************************************************************
    * Sequences of random numbers
    **************************************************************/
    struct SimVar {
        ArrayXXd lnphi_sim;
        ArrayXXd lnLur_sim;
        ArrayXXd lnF_PD;
        ArrayXXd lnF_SE;
        ArrayXXd missing_sim;

        ArrayXXd PDsim_Empirical;
        ArrayXXd SEsim_Empirical;
        ArrayXXd Misssim_Empirical;

        ArrayXXd StateCategory_Sim;
    };

    /**************************************************************
    * Simulate sequences of random numbers
    **************************************************************/
    SimVar SimulationBase(const int & SimN, const int & SimTbar);


}

// ////
// //    /**************************************************************
// //    * Simulated/Real Data Structure including all relevant variables
// //    **************************************************************/
// //    struct SimData {
// //        ArrayXXd phi; ArrayXXd Revenue; ArrayXXi Revenue_Miss;
// //        ArrayXXd Capital; ArrayXXi CapitalMiss;
// //        ArrayXXd Employ_ur; ArrayXXi Employ_urMiss;
// //        ArrayXXd Employ_uc; ArrayXXi Employ_ucMiss;
// //
// //        ArrayXXd Revenue_P; ArrayXXd Revenue_D;
// //        //// Capital/Employment if the previous period is production
// //        ArrayXXd Capital_P; ArrayXXd Employ_ur_P; ArrayXXd Employ_uc_P;
// //        //// Capital/Employment if the previous period is dormancy
// //        ArrayXXd Capital_D; ArrayXXd Employ_ur_D; ArrayXXd Employ_uc_D;
// //
// ////        ArrayXXd Employ_ur_cutoff;
// ////        ArrayXXd profit; ArrayXXd cum_profit;
// //        ArrayXXd Price; ArrayXXd Price_P; ArrayXXd Price_D;
// //        ArrayXXd cost; ArrayXXd cost_P; ArrayXXd cost_D;
// //
// //        ArrayXi good; ArrayXi labour_intensive_ind; ArrayXi state_code; ArrayXi industry; ArrayXi first_year_in_Data;
// //        ArrayXi last_year_in_Data;
// ////        ArrayXXi count_DP; ArrayXXi MissCount_DP;
// //
// //        ArrayXXi State_P; ArrayXXi State_D; ArrayXXi State_S; ArrayXXi State_E;
// //        ArrayXXi State_P_nonMiss; ArrayXXi State_D_nonMiss; ArrayXXi State_DP_Miss; ArrayXXi State_E_Miss;
// //        ArrayXXi yr_DP;
// //        ArrayXXi Employ_ur_Extreme; ArrayXXi Employ_uc_Extreme; ArrayXXi Capital_Extreme;
// //
// //        //        ArrayXXd Prob_Missing;
// //        //// The distribution of states at each moment without missing cases
// //        // Prob_P + Prob_D + Prob_E = 1
// //        // Prob_S = Prob_P + Prob_D
// //        ArrayXXd Prob_P; ArrayXXd Prob_D; ArrayXXd Prob_S; ArrayXXd Prob_E;
// //
// //        // These are the conditional transition probability: prod/dorm/exit/stay probability given the state in the previous period
// //        ArrayXXd ProbP_P; ArrayXXd ProbP_D; ArrayXXd ProbD_P; ArrayXXd ProbD_D;
// //        ArrayXXd ProbP_S; ArrayXXd ProbP_E; ArrayXXd ProbD_S; ArrayXXd ProbD_E;
// //
// //        //// the distribution of states at each moment over time taking the missing probability into account
// //        ArrayXXd Prob_P_nonMiss; ArrayXXd Prob_D_nonMiss; ArrayXXd Prob_DP_Miss; ArrayXXd Prob_E_Miss;
// //        //// conditional transition probability given the firm is production in the previous period
// //        ArrayXXd ProbP_P_nonMiss; ArrayXXd ProbP_D_nonMiss; ArrayXXd ProbP_DP_Miss_UnObs; ArrayXXd ProbP_E_Miss_UnObs;
// //        ArrayXXd ProbP_DP_Miss; ArrayXXd ProbP_E_Miss;
// //        //// conditional transition probability given the firm is dormant in the previous period
// //        ArrayXXd ProbD_P_nonMiss; ArrayXXd ProbD_D_nonMiss; ArrayXXd ProbD_DP_Miss_UnObs; ArrayXXd ProbD_E_Miss_UnObs;
// //        ArrayXXd ProbD_DP_Miss; ArrayXXd ProbD_E_Miss;
// //
// //        ArrayXXd Prob_Missing_P1; ArrayXXd Prob_Missing_D1; ArrayXXd Prob_Missing_DP1;
// //        // probability of each states without considering exiting probability;
// //        // This is the actual states of firms but not completely observable
// //        // as researcher may see a firm never reappears in the data but it is actually in state Prob_DP_Miss_UnObs
// //        // Note that sim_data_missing.Prob_P_nonMiss + sim_data_missing.Prob_D_nonMiss + Prob_DP_Miss_UnObs + Prob_E_Miss_UnObs = 1
// //        ArrayXXd Prob_DP_Miss_UnObs; ArrayXXd Prob_E_Miss_UnObs;
// //
// //        ArrayXXd K_sol; ArrayXXd Lur_sol; ArrayXXd Luc_sol;
// //        ArrayXXd K_sol_P; ArrayXXd Lur_sol_P; ArrayXXd Luc_sol_P;
// //        ArrayXXd K_sol_D; ArrayXXd Lur_sol_D; ArrayXXd Luc_sol_D;
// //
// //        ArrayXXd CutoffVal;
// //        ArrayXXd F_PD;
// //
// //    };
// ////
// //    /**************************************************************
// //    * Import Read Data
// //    **************************************************************/
// //    SimData createRealData(const string & RealSim, const int & GoodState, const int & LaborIntensive);
// //

//

// }
