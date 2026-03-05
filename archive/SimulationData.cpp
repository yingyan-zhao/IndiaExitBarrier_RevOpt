#include "SimulationData.h"

using namespace Eigen;
using namespace std;
using namespace alias;
//
/**************************************************************
* Import Read Data
**************************************************************/
SimData alias::createRealData(const string & RealSim, const int & GoodState, const int & LaborIntensive) {

    string source_path(SOURCE_FOLDER);
    string path = source_path + "/StructuralData_InflationAdj_NormalizeAll/";
    string cat = "Good" + to_string(GoodState) + "_Sector" + to_string(LaborIntensive) + RealSim + ".csv";
    cout << "path: " << path << endl;

    SimData RealData;
    /* Production, dormancy etc Status */
    load_csv(RealData.State_P_nonMiss,path + "ASI_Data_Estimation_State_P_NonMiss_" + cat);
    load_csv(RealData.State_D_nonMiss,path + "ASI_Data_Estimation_State_D_NonMiss_" + cat);
    load_csv(RealData.State_DP_Miss,path + "ASI_Data_Estimation_State_DP_Miss_" + cat);
    load_csv(RealData.State_E_Miss,path + "ASI_Data_Estimation_State_E_Miss_" + cat);

    /* industry and state code index */
    ArrayXXd Indicator(RealData.State_P_nonMiss.rows(),6);
    load_csv(Indicator,path + "ASI_Data_Estimation_Indicator_" + cat);
    RealData.good = Indicator.col(0).cast<int>();
    RealData.labour_intensive_ind = Indicator.col(1).cast<int>();
    RealData.state_code = Indicator.col(2).cast<int>();
    RealData.industry = Indicator.col(3).cast<int>();
    RealData.first_year_in_Data = Indicator.col(4).cast<int>();
    RealData.last_year_in_Data = Indicator.col(5).cast<int>();
    cout << "34" << endl;

    /* value of value added; capital and Employment */
    load_csv(RealData.Revenue,path + "ASI_Data_Estimation_va_" + cat);
    load_csv(RealData.Capital,path + "ASI_Data_Estimation_Capital_" + cat);
//    load_csv(RealData.Employ_ur,path + "ASI_Data_Estimation_payroll_regular_" + cat);
//    load_csv(RealData.Employ_uc,path + "ASI_Data_Estimation_payroll_contract_" + cat);
    load_csv(RealData.Employ_ur,path + "ASI_Data_Estimation_Employ_ur_" + cat);
    load_csv(RealData.Employ_uc,path + "ASI_Data_Estimation_Employ_uc_" + cat);

    /* Indicator matrix: = 1 if value added / capital / Employment are missing */
    load_csv(RealData.Revenue_Miss,path + "ASI_Data_Estimation_va_Miss_" + cat);
    load_csv(RealData.CapitalMiss,path + "ASI_Data_Estimation_Capital_Miss_" + cat);
//    load_csv(RealData.Employ_urMiss,path + "ASI_Data_Estimation_payroll_regular_Miss_" + cat);
//    load_csv(RealData.Employ_ucMiss,path + "ASI_Data_Estimation_payroll_contract_Miss_" + cat);
    load_csv(RealData.Employ_urMiss,path + "ASI_Data_Estimation_Employ_ur_Miss_" + cat);
    load_csv(RealData.Employ_ucMiss,path + "ASI_Data_Estimation_Employ_uc_Miss_" + cat);

    load_csv(RealData.yr_DP,path + "ASI_Data_Estimation_yr_DP_" + cat);
    load_csv(RealData.Employ_ur_Extreme,path + "ASI_Data_Estimation_Employ_ur_Extreme_" + cat);
    load_csv(RealData.Employ_uc_Extreme,path + "ASI_Data_Estimation_Employ_uc_Extreme_" + cat);
    load_csv(RealData.Capital_Extreme,path + "ASI_Data_Estimation_Capital_Extreme_" + cat);
    cout << "56" << endl;
    /* Note */
    /* the function "load_csv" is in Auxillaries.h */
    /* The csv file cannot contain missing data, otherwise, the program will have error. */
    /* Replace all missing to 0, and use an indicator matrix to indicate missing status. */

    return RealData;
}

/**************************************************************
* Simulate sequences of random numbers
**************************************************************/
SimVar alias::SimulationBase(const int & SimN, const int & SimTbar) {
    SimVar sim_var;

    ArrayXXd randomdraw = quasi_random_uniform(SimTbar*11+1, SimN);
    // initial period phi
    // phi
    ArrayXXd lnphi_sim_temp = randomdraw(all,seqN(0,SimTbar+1));
    sim_var.lnphi_sim = quasi_random_uniform_to_normal(lnphi_sim_temp, 0, 1);
//    sim_var.lnphi_sim = ArrayXXd::Zero(para.SimN,para.SimTbar);

    // Labor shock
    ArrayXXd lnLuc_sim_temp = randomdraw(all,seqN(SimTbar+1,SimTbar));
    sim_var.lnLuc_sim = quasi_random_uniform_to_normal(lnLuc_sim_temp, 0, 1);
    ArrayXXd lnLur_sim_temp = randomdraw(all,seqN(SimTbar*2+1,SimTbar));
    sim_var.lnLur_sim = quasi_random_uniform_to_normal(lnLur_sim_temp, 0, 1);
    ArrayXXd lnK_sim_temp = randomdraw(all,seqN(SimTbar*3+1,SimTbar));
    sim_var.lnK_sim = quasi_random_uniform_to_normal(lnK_sim_temp, 0, 1);
    // Fixed production cost / exit cost
    sim_var.lnF_PD = randomdraw(all,seqN(SimTbar*4+1,SimTbar));
    sim_var.lnF_SE = randomdraw(all,seqN(SimTbar*5+1,SimTbar));
    // missing probability
    sim_var.missing_sim = randomdraw(all,seqN(SimTbar*6+1,SimTbar));

    // probability of status
    sim_var.PDsim_Empirical = randomdraw(all,seqN(SimTbar*7+1,SimTbar));
    sim_var.SEsim_Empirical = randomdraw(all,seqN(SimTbar*8+1,SimTbar));
    sim_var.Misssim_Empirical = randomdraw(all,seqN(SimTbar*9+1,SimTbar));

    sim_var.StateCategory_Sim = randomdraw(all,seqN(SimTbar*10+1,SimTbar));

    return sim_var;
}
