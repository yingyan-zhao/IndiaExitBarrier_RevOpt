#include "EstimationIndirectInference_StandardDeviation.h"

using namespace Eigen;
using namespace std;
using namespace alias;


/**************************************************************
* Standard deviation
**************************************************************/
ArrayXd alias::EstimationIndirectInference_StandardDeviation(const ArrayXd & theta_Est,const ArrayXd & theta_a,
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    cout << "theta_Est = " << theta_Est.transpose() << endl;
    cout << "theta_a = " << theta_a.transpose() << endl;

    int GoodState; int LaborIntensive;
    GoodState = 1; LaborIntensive = 1;
    SimData simdata_panel_GoodState_LaborInt = createRealData("_ForSim",GoodState,LaborIntensive);
    GoodState = 0; LaborIntensive = 1;
    SimData simdata_panel_BadState_LaborInt = createRealData("_ForSim",GoodState,LaborIntensive);
    GoodState = 1; LaborIntensive = 0;
    SimData simdata_panel_GoodState_CapitalInt = createRealData("_ForSim",GoodState,LaborIntensive);
    GoodState = 0; LaborIntensive = 0;
    SimData simdata_panel_BadState_CapitalInt = createRealData("_ForSim",GoodState,LaborIntensive);
    cout << "Data imported for simulation" << endl;

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


    /******* JMatrix *******/
    // to prepare for the calculation of the derivative at theta_a, we calculate the value/policy function at theta_a and theta_a + epsilon
    cout << "Start calculate JMatrix " << endl;
    ArrayXXd JArray = EstimationIndirectInference_SD_JMatrix_RealData(theta_a,
        ReadData_GoodState_LaborInt, ReadData_BadState_LaborInt,
        ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
        threadsManagement);
    writeToCSVfile("JArray.csv", JArray.cast<double>().matrix());
    // ArrayXXd JArray = ArrayXXd::Zero( theta_Est.rows(), theta_Est.rows() );
    MatrixXd JMatrix = JArray.matrix();
    cout << "JMatrix done" << endl;
    cout << "JArray = " << JArray << endl;
    MatrixXd JMatrix_inv = JMatrix.inverse();
    cout << "JMatrix_inv done" << endl;
    cout << "JMatrix_inv = " << JMatrix_inv << endl;
    writeToCSVfile("JMatrix_inv.csv", JMatrix_inv.cast<double>().matrix());
    // throw runtime_error("30");
    // // // //
    // ArrayXXd JArray;
    // load_csv(JArray,"JArray.csv");
    // MatrixXd JMatrix = JArray.matrix();
    // cout << "JArray = " << JArray << endl;
    //
    cout << "Start calculate IArray " << endl;
    ArrayXXd IArray = EstimationIndirectInference_SD_IMatrix_RealData(theta_a,
        ReadData_GoodState_LaborInt, ReadData_BadState_LaborInt,
        ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
        threadsManagement);
    // // ArrayXXd IArray = ArrayXXd::Zero( theta_Est.rows(), theta_Est.rows() );
    MatrixXd IMatrix = IArray.matrix();
    cout << "IMatrix done" << endl;
    cout << "IArray = " << IArray << endl;
    writeToCSVfile("IArray.csv", IArray.cast<double>().matrix());
    // throw runtime_error("97");
    //
    // ArrayXXd IArray;
    // load_csv(IArray,"IArray.csv");
    // MatrixXd IMatrix = IArray.matrix();
    // cout << "IArray = " << IArray << endl;
    // // throw runtime_error("38");
    // // //
    // // //
    //
    cout << "Start calculate IStarArray " << endl;
    ArrayXXd IStarArray = EstimationIndirectInference_SD_IStarMatrix_SimData(theta_Est,theta_a,
        simdata_panel_GoodState_LaborInt,simdata_panel_BadState_LaborInt,
        simdata_panel_GoodState_CapitalInt,simdata_panel_BadState_CapitalInt,
        sim_var_2step_GoodState_LaborInt,sim_var_2step_BadState_LaborInt,
        sim_var_2step_GoodState_CapitalInt,sim_var_2step_BadState_CapitalInt,
        threadsManagement);
    // ArrayXXd IStarArray = ArrayXXd::Zero( theta_Est.rows(), theta_Est.rows() );
    MatrixXd IStarMatrix = IStarArray.matrix();
    cout << "IStarMatrix done" << endl;
    cout << "IStarArray = " << IStarArray << endl;
    writeToCSVfile("IStarArray.csv", IStarArray.cast<double>().matrix());
    // throw runtime_error("38");
    // // // //
    // ArrayXXd IStarArray;
    // load_csv(IStarArray,"IStarArray.csv");
    // MatrixXd IStarMatrix = IStarArray.matrix();
    // cout << "IStarMatrix = " << IStarMatrix << endl;
    // // throw runtime_error("71");
    // //
    cout << "Start calculate db_dtheta_a " << endl;
    ArrayXXd db_dtheta_a = Cal_db_dtheta_a(theta_Est,theta_a,
        simdata_panel_GoodState_LaborInt, simdata_panel_BadState_LaborInt,
        simdata_panel_GoodState_CapitalInt, simdata_panel_BadState_CapitalInt,
        sim_var_2step_GoodState_LaborInt, sim_var_2step_BadState_LaborInt,
        sim_var_2step_GoodState_CapitalInt, sim_var_2step_BadState_CapitalInt,
        threadsManagement);
    MatrixXd db_dtheta_a_matrix = db_dtheta_a.matrix();
    writeToCSVfile("db_dtheta_a.csv", db_dtheta_a.cast<double>().matrix());
    // throw runtime_error("76");
    //

    // ArrayXXd db_dtheta_a;
    // load_csv(db_dtheta_a,"db_dtheta_a.csv");
    // MatrixXd db_dtheta_a_Matrix = db_dtheta_a.matrix();
    // cout << "db_dtheta_a = " << db_dtheta_a << endl;

    // cout << "Start calculate db_dtheta0_v1 " << endl;
    // ArrayXXd db_dtheta0_v1 = Cal_db_dtheta0_v1(theta_Est,theta_a,
    //     simdata_panel_GoodState_LaborInt, simdata_panel_BadState_LaborInt,
    //     simdata_panel_GoodState_CapitalInt, simdata_panel_BadState_CapitalInt,
    //     sim_var_2step_GoodState_LaborInt, sim_var_2step_BadState_LaborInt,
    //     sim_var_2step_GoodState_CapitalInt, sim_var_2step_BadState_CapitalInt,
    //     threadsManagement);
    //     MatrixXd db_dtheta0_matrix = db_dtheta0.matrix();
    // writeToCSVfile("db_dtheta0.csv", db_dtheta0.cast<double>().matrix());
    // MatrixXd db_dtheta0_Matrix = db_dtheta0.matrix();
    // cout << "db_dtheta0 = " << db_dtheta0 << endl;
    // // throw runtime_error("127");
    //
    // ArrayXXd db_dtheta0_v1;
    // load_csv(db_dtheta0_v1,"db_dtheta0.csv");
    // cout << "db_dtheta0_v1 = " << db_dtheta0_v1 << endl;
    // //
    //
    cout << "Start calculate db_dtheta0_v2 " << endl;
    ArrayXXd db_dtheta0_v2 = Cal_db_dtheta0_v2(theta_Est,theta_a,
        simdata_panel_GoodState_LaborInt, simdata_panel_BadState_LaborInt,
        simdata_panel_GoodState_CapitalInt, simdata_panel_BadState_CapitalInt,
        sim_var_2step_GoodState_LaborInt, sim_var_2step_BadState_LaborInt,
        sim_var_2step_GoodState_CapitalInt, sim_var_2step_BadState_CapitalInt,
        threadsManagement);
    writeToCSVfile("db_dtheta0_v2.csv", db_dtheta0_v2.cast<double>().matrix());
    cout << "db_dtheta0_v2 = " << db_dtheta0_v2 << endl;
    throw runtime_error("127");
    //
    // MatrixXd db_dtheta0_Matrix = db_dtheta0_v2.matrix();
    //
    // double S = 1.0;
    // // MatrixXd theta_Var = dtheta_a_dtheta0_Matrix.inverse() * JMatrix.inverse() * (IMatrix + 1.0/S * IStarMatrix)
    // //     * JMatrix.inverse() * dtheta_a_dtheta0_Matrix.inverse().transpose()
    // //     / double(SimN_2step_GoodState_LaborInt + SimN_2step_BadState_LaborInt
    // //       + SimN_2step_GoodState_CapitalInt + SimN_2step_BadState_CapitalInt);
    // MatrixXd theta_Var = db_dtheta0_Matrix.inverse() * db_dtheta_a_Matrix * JMatrix.inverse() * (IMatrix + 1.0/S * IStarMatrix)
    //     * JMatrix.inverse() * db_dtheta_a_Matrix * db_dtheta0_Matrix.inverse();
    //
    // cout << "theta_Var = " << theta_Var << endl;
    // ArrayXd theta_SD = theta_Var.diagonal();
    // cout << "theta_SD = " << theta_SD.transpose() << endl;
     ArrayXd theta_SD;
     throw runtime_error("143");
     return theta_SD;
}

//
/** ****************************************************/
/*** Calculate the JMatrix ***/
/** ****************************************************/
ArrayXXd alias::EstimationIndirectInference_SD_JMatrix_RealData(const ArrayXd & theta,
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    /*** Part 1 ***/
    ArrayXXd JMatrix_1_temp = ArrayXXd::Zero(para.dim1, para.dim1);
    ArrayXi NumFirms_1_Good = ArrayXi::Zero(para.dim1);
    ArrayXi NumFirms_1_Bad = ArrayXi::Zero(para.dim1);
    for (size_t k = 0; k < para.dim1; k++) {
        cout << "k = " << k << endl;
        double eps = abs(theta(para.dim3 + k)) * 1e-8;

        ArrayXd theta_plus = theta;
        theta_plus(para.dim3 + k) = theta(para.dim3 + k) + eps;
        ArrayXd theta_minus = theta;
        theta_minus(para.dim3 + k) = theta(para.dim3 + k) - eps;
        cout << theta_plus.transpose() << endl;
        cout << theta_minus.transpose() << endl;

        tuple<ArrayXXd,int,int,int,int,int> t_dL_part1_plus = EstimationIndirectInference_SD_dl_dtheta_Part1_RealData(theta_plus,
            ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part1_plus = get<0>(t_dL_part1_plus);
        int count_part1_plus = get<1>(t_dL_part1_plus);
        int count_firms_part1_GoodState_LaborInt_plus = get<2>(t_dL_part1_plus);
        int count_firms_part1_BadState_LaborInt_plus = get<3>(t_dL_part1_plus);
        int count_firms_part1_GoodState_CapitalInt_plus = get<4>(t_dL_part1_plus);
        int count_firms_part1_BadState_CapitalInt_plus = get<5>(t_dL_part1_plus);
        tuple<ArrayXXd,int,int,int,int,int> t_dL_part1_minus = EstimationIndirectInference_SD_dl_dtheta_Part1_RealData(theta_minus,
            ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part1_minus = get<0>(t_dL_part1_minus);
        int count_part1_minus = get<1>(t_dL_part1_minus);
        int count_firms_part1_GoodState_LaborInt_minus = get<2>(t_dL_part1_minus);
        int count_firms_part1_BadState_LaborInt_minus = get<3>(t_dL_part1_minus);
        int count_firms_part1_GoodState_CapitalInt_minus = get<4>(t_dL_part1_minus);
        int count_firms_part1_BadState_CapitalInt_minus = get<5>(t_dL_part1_minus);
        ArrayXXd ddL_part1 = ( -dL_part1_plus - (-dL_part1_minus) ) / 2.0 / eps;
        cout << "count_part1_plus = " << count_part1_plus << "; count_part1_minus = " << count_part1_minus << endl;
        // cout << "dL_part1_plus = " << dL_part1_plus.rows() << "; " << dL_part1_plus.cols() << endl;
        // cout << "dL_part1_minus = " << dL_part1_minus.rows() << "; " << dL_part1_minus.cols() << endl;

        ArrayXXd dL_part1_plus_vec = Vectorize_dl_dtheta_part(dL_part1_plus, para.dim1);
        ArrayXXd dL_part1_minus_vec = Vectorize_dl_dtheta_part(dL_part1_minus, para.dim1);
        // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;
        // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;

        NumFirms_1_Good(k) = count_firms_part1_GoodState_LaborInt_plus + count_firms_part1_GoodState_CapitalInt_plus;
        NumFirms_1_Bad(k) = count_firms_part1_BadState_LaborInt_plus + count_firms_part1_BadState_CapitalInt_plus;

        JMatrix_1_temp.row(k) = ( -dL_part1_plus_vec - (-dL_part1_minus_vec) ).colwise().sum() / 2.0 / eps;
        // cout << "dL_part1_plus.rows() = " << dL_part1_plus.rows() << endl;
        // cout << "Num = " << Num << endl;
        cout << "JMatrix_1.row(k) = " << JMatrix_1_temp.row(k) << endl;
        cout << "count_part1_plus = " << count_part1_plus << endl;
        cout << "dL_part1_plus.rows = " << dL_part1_plus.rows() << endl;
        cout << "abs(dL_part1_plus_vec - dL_part1_minus_vec) = "
            << ((dL_part1_plus_vec - dL_part1_minus_vec).abs() > 0).cast<int>().colwise().sum() << endl;
    }
    cout << "NumFirms_1_Bad = " << NumFirms_1_Bad.transpose() << endl;
    cout << "NumFirms_1_Good = " << NumFirms_1_Good.transpose() << endl;
    for (size_t k = 0; k < para.dim1; k++) {
        for (size_t kk = 0; kk < para.dim1; kk++) {
            if (k == 0 or kk == 0) {
                JMatrix_1_temp(k,kk) = JMatrix_1_temp(k,kk) / double(NumFirms_1_Good(k));
            }
            else if (k == 1 or kk == 1) {
                JMatrix_1_temp(k,kk) = JMatrix_1_temp(k,kk) / double(NumFirms_1_Bad(k));
            }
            else {
                JMatrix_1_temp(k,kk) = JMatrix_1_temp(k,kk) / double(NumFirms_1_Good(k) + NumFirms_1_Bad(k));
            }
        }
    }
    cout << "JMatrix_1_temp = " << JMatrix_1_temp << endl;
    ArrayXXd JMatrix_1_tran = JMatrix_1_temp.transpose();
    ArrayXXd JMatrix_1 = JMatrix_1_temp*0.5+JMatrix_1_tran*0.5;


    /*** Part 2 ***/
    ArrayXXd JMatrix_2_temp = ArrayXXd::Zero(para.dim2, para.dim2);
    ArrayXi NumFirms_2_P = ArrayXi::Zero(para.dim2);
    ArrayXi NumFirms_2_D = ArrayXi::Zero(para.dim2);
    for (size_t k = 0; k < para.dim2; k++) {
        cout << "k = " << k << endl;
        double eps = abs(theta(para.dim3 + para.dim1 + k)) * 1e-8;

        ArrayXd theta_plus = theta;
        theta_plus(para.dim3 + para.dim1 + k) = theta(para.dim3 + para.dim1 + k) + eps;
        ArrayXd theta_minus = theta;
        theta_minus(para.dim3 + para.dim1 + k) = theta(para.dim3 + para.dim1 + k) - eps;
        cout << theta_plus.transpose() << endl;
        cout << theta_minus.transpose() << endl;

        tuple<ArrayXXd,int,int,int,int> t_dL_part2_plus = EstimationIndirectInference_SD_dl_dtheta_Part2_RealData(theta_plus,
            ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part2_plus = get<0>(t_dL_part2_plus);
        int count_part2_plus = get<1>(t_dL_part2_plus);
        int countP_firms_part2_plus = get<2>(t_dL_part2_plus);
        int countD_firms_part2_plus = get<3>(t_dL_part2_plus);
        int count_firms_part2_plus = get<4>(t_dL_part2_plus);

        tuple<ArrayXXd,int,int,int,int> t_dL_part2_minus = EstimationIndirectInference_SD_dl_dtheta_Part2_RealData(theta_minus,
            ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part2_minus = get<0>(t_dL_part2_minus);
        int count_part2_minus = get<1>(t_dL_part2_minus);
        int countP_firms_part2_minus = get<2>(t_dL_part2_minus);
        int countD_firms_part2_minus = get<3>(t_dL_part2_minus);
        int count_firms_part2_minus = get<4>(t_dL_part2_minus);

        ArrayXXd ddL_part2 = ( -dL_part2_plus - (-dL_part2_minus) ) / 2.0 / eps;
        // cout << "dL_part2_plus = " << dL_part2_plus.rows() << "; " << dL_part2_plus.cols() << endl;
        // cout << "dL_part2_minus = " << dL_part2_minus.rows() << "; " << dL_part2_minus.cols() << endl;

        ArrayXXd dL_part2_plus_vec = Vectorize_dl_dtheta_part(dL_part2_plus, para.dim2);
        ArrayXXd dL_part2_minus_vec = Vectorize_dl_dtheta_part(dL_part2_minus, para.dim2);
        // cout << "dL_part1_plus_vec = " << dL_part2_plus_vec.rows() << "; " << dL_part2_plus_vec.cols() << endl;
        // cout << "dL_part2_minus_vec = " << dL_part2_minus_vec.rows() << "; " << dL_part2_minus_vec.cols() << endl;

        NumFirms_2_P(k) = countP_firms_part2_plus;
        NumFirms_2_D(k) = countD_firms_part2_plus;

        JMatrix_2_temp.row(k) = ( -dL_part2_plus_vec - (-dL_part2_minus_vec) ).colwise().sum() / 2.0 / eps;

        cout << "JMatrix_2.row(k) = " << JMatrix_2_temp.row(k) << endl;
        cout << "count_part2_plus = " << count_part2_plus << endl;
        cout << "dL_part2_plus.rows = " << dL_part2_plus.rows() << endl;
        cout << "abs(dL_part2_plus_vec - dL_part2_minus_vec) = "
            << (abs(dL_part2_plus_vec - dL_part2_minus_vec) > 0).cast<int>().colwise().sum() << endl;
    }
    cout << "NumFirms_2_P = " << NumFirms_2_P.transpose() << endl;
    cout << "NumFirms_2_D = " << NumFirms_2_D.transpose() << endl;
    for (size_t k = 0; k < para.dim2; k++) {
        for (size_t kk = 0; kk < para.dim2; kk++) {
            if (k == 0 or kk == 0) {
                JMatrix_2_temp(k,kk) = JMatrix_2_temp(k,kk) / double(NumFirms_2_P(k));
            }
            else if (k == 1 or kk == 1) {
                JMatrix_2_temp(k,kk) = JMatrix_2_temp(k,kk) / double(NumFirms_2_D(k));
            }
            else {
                JMatrix_2_temp(k,kk) = JMatrix_2_temp(k,kk) / double(NumFirms_2_P(k) + NumFirms_2_D(k));
            }
        }
    }
    cout << "JMatrix_2_temp = " << JMatrix_2_temp << endl;
    ArrayXXd JMatrix_2_tran = JMatrix_2_temp.transpose();
    ArrayXXd JMatrix_2 = JMatrix_2_temp*0.5+JMatrix_2_tran*0.5;
    // throw runtime_error("232");

    /*** Part 3 ***/
    // ArrayXXd JMatrix_3_temp = ArrayXXd::Zero(para.dim3, para.dim3);
    // for (size_t k = 0; k < para.dim3; k++) {
    //     cout << "k = " << k << endl;
    //
    //     double eps = abs(theta(k)) * 1e-2;
    //
    //     ArrayXd theta_plus = theta;
    //     theta_plus(k) = theta(k) + eps;
    //     ArrayXd theta_minus = theta;
    //     theta_minus(k) = theta(k) - eps;
    //     cout << theta_plus.transpose() << endl;
    //     cout << theta_minus.transpose() << endl;
    //
    //     int GoodState; int LaborIntensive;
    //     GoodState = 1; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_GoodState_LaborInt_plus = CalEquV_plus_minus(theta_plus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_BadState_LaborInt_plus = CalEquV_plus_minus(theta_plus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 1; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_GoodState_CapitalInt_plus = CalEquV_plus_minus(theta_plus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_BadState_CapitalInt_plus = CalEquV_plus_minus(theta_plus, GoodState, LaborIntensive, threadsManagement);
    //     tuple<ArrayXXd,int,int> t_dL_part3_plus = EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(
    //         ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
    //         EquV_GoodState_LaborInt_plus,EquV_BadState_LaborInt_plus,
    //         EquV_GoodState_CapitalInt_plus,EquV_BadState_CapitalInt_plus,
    //         threadsManagement);
    //     ArrayXXd dL_part3_plus = get<0>(t_dL_part3_plus);
    //     int count_Good_part3_plus = get<1>(t_dL_part3_plus);
    //     int count_Bad_part3_plus = get<2>(t_dL_part3_plus);
    //
    //     GoodState = 1; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_GoodState_LaborInt_minus = CalEquV_plus_minus(theta_minus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_BadState_LaborInt_minus = CalEquV_plus_minus(theta_minus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 1; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_GoodState_CapitalInt_minus = CalEquV_plus_minus(theta_minus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_BadState_CapitalInt_minus = CalEquV_plus_minus(theta_minus, GoodState, LaborIntensive, threadsManagement);
    //     tuple<ArrayXXd,int,int> t_dL_part3_minus = EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(
    //         ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
    //         EquV_GoodState_LaborInt_minus,EquV_BadState_LaborInt_minus,
    //         EquV_GoodState_CapitalInt_minus,EquV_BadState_CapitalInt_minus,
    //         threadsManagement);
    //     ArrayXXd dL_part3_minus = get<0>(t_dL_part3_minus);
    //
    //     ArrayXXd ddL_part3 = ( -dL_part3_plus -  (-dL_part3_minus) ) / 2.0 / eps;
    //     // cout << "ddL_part3 = " << ddL_part3.colwise().sum().cols() << "; " << ddL_part3.rows() << "; EstTbar = " << para.EstTbar << endl;
    //     // cout << "ddL_part3 = " << ddL_part3.colwise().sum() << endl;
    //     // throw runtime_error("ddL_part3 done");
    //
    //     ArrayXi NumObs = ArrayXi::Ones(para.dim3);
    //     for (size_t kk = 0; kk < para.dim3; ++kk) {
    //         ArrayXXi ddL_part3_index_mat = ( ( dL_part3_plus.middleCols(kk*para.EstTbar,para.EstTbar)
    //             - dL_part3_minus.middleCols(kk*para.EstTbar,para.EstTbar) ).abs() > 0 ).cast<int>();
    //         ArrayXi ddL_part3_index_vec = ddL_part3_index_mat.rowwise().maxCoeff();
    //
    //         if (ddL_part3_index_vec.sum() > 1) {NumObs(kk) = ddL_part3_index_vec.sum();}
    //     }
    //     cout << "NumObs = " << NumObs.transpose() << endl;
    //
    //     ArrayXXd dL_part3_plus_vec = Vectorize_dl_dtheta_part(dL_part3_plus,para.dim3);
    //     ArrayXXd dL_part3_minus_vec = Vectorize_dl_dtheta_part(dL_part3_minus,para.dim3);
    //     // cout << "dL_part3_plus_vec size = " << dL_part3_plus_vec.colwise().sum() << endl;
    //     // cout << "dL_part3_minus_vec size = " << dL_part3_minus_vec.colwise().sum() << endl;
    //     //
    //     JMatrix_3_temp.row(k) = ( -dL_part3_plus_vec - (-dL_part3_minus_vec) ).colwise().sum() / 2.0 / eps / NumObs.cast<double>().transpose();
    //     // cout << "double(ddL_part3.rows()) = " << double(ddL_part3.rows()) << endl;
    //     cout << "JMatrix_3_temp.row(k) = " << JMatrix_3_temp.row(k) << endl;
    //     // cout << "count_part3_plus = " << count_part3_plus << endl;
    //     // throw runtime_error("JMatrix done: k = 0");
    //     // cout << "dL_part3_plus.rows = " << dL_part3_plus.rows() << endl;
    // }

    ArrayXXd JMatrix_3_temp = ArrayXXd::Zero(para.dim3, para.dim3);
    for (size_t k = 0; k < para.dim3; k++) {
        cout << "k = " << k << endl;

        double eps_theta1 = abs(theta(k)) * 1e-8;

        ArrayXd theta_a_plus = theta;
        theta_a_plus(k) = theta(k) + eps_theta1;
        ArrayXd theta_a_minus = theta;
        theta_a_minus(k) = theta(k) - eps_theta1;
        cout << theta_a_plus.transpose() << endl;
        cout << theta_a_minus.transpose() << endl;

        for (size_t kk = 0; kk < para.dim3; kk++) {
            cout << "kk = " << kk << endl;
            double eps_theta2 = abs(theta(kk)) * 1e-8;

            ArrayXd theta_a_plus_plus = theta_a_plus;
            theta_a_plus_plus(kk) = theta_a_plus(kk) + eps_theta2;
            ArrayXd theta_a_plus_minus = theta_a_plus;
            theta_a_plus_minus(kk) = theta_a_plus(kk) - eps_theta2;

            ArrayXd theta_a_minus_plus = theta_a_minus;
            theta_a_minus_plus(kk) = theta_a_minus(kk) + eps_theta2;
            ArrayXd theta_a_minus_minus = theta_a_minus;
            theta_a_minus_minus(kk) = theta_a_minus(kk) - eps_theta2;

            ArrayXXd logL_part3_plus_plus = EstimationIndirectInference_SD_l_Part3_RealData(theta_a_plus_plus,
                ReadData_GoodState_LaborInt, ReadData_BadState_LaborInt, ReadData_GoodState_CapitalInt, ReadData_BadState_CapitalInt,
                threadsManagement);

            ArrayXXd logL_part3_plus_minus = EstimationIndirectInference_SD_l_Part3_RealData(theta_a_plus_minus,
                ReadData_GoodState_LaborInt, ReadData_BadState_LaborInt, ReadData_GoodState_CapitalInt, ReadData_BadState_CapitalInt,
                threadsManagement);

            ArrayXXd logL_part3_minus_plus = EstimationIndirectInference_SD_l_Part3_RealData(theta_a_minus_plus,
                ReadData_GoodState_LaborInt, ReadData_BadState_LaborInt, ReadData_GoodState_CapitalInt, ReadData_BadState_CapitalInt,
                threadsManagement);

            ArrayXXd logL_part3_minus_minus = EstimationIndirectInference_SD_l_Part3_RealData(theta_a_minus_minus,
                ReadData_GoodState_LaborInt, ReadData_BadState_LaborInt, ReadData_GoodState_CapitalInt, ReadData_BadState_CapitalInt,
                threadsManagement);

            ArrayXXd ddL_part3 = - ( (logL_part3_plus_plus - logL_part3_plus_minus) / 2.0 / eps_theta2
                - (logL_part3_minus_plus - logL_part3_minus_minus) / 2.0 / eps_theta2 ) / 2.0 / eps_theta1;

            ArrayXd ddL_part3_sum = ddL_part3.abs().rowwise().sum();
            int NumObs = (ddL_part3_sum.abs() > 0).cast<int>().sum();

            JMatrix_3_temp(k,kk) = ddL_part3.sum() / ( double(NumObs)+0.001 );
            cout << "k = " << k << "; kk = " << kk << "; JMatrix_3_temp(k,kk) = " << JMatrix_3_temp(k,kk)
                << "; NumObs = " << NumObs << endl;
        }
    }

    cout << "JMatrix_3_temp = " << JMatrix_3_temp << endl;
    writeToCSVfile("JMatrix_3_temp.csv", JMatrix_3_temp.cast<double>().matrix());
    // throw runtime_error("451");
    ArrayXXd JMatrix_3_tran = JMatrix_3_temp.transpose();
    ArrayXXd JMatrix_3 = JMatrix_3_temp*0.5+JMatrix_3_tran*0.5;

    /*** Full JMatrix ***/
    ArrayXXd JMatrix = ArrayXXd::Zero(para.dim1+para.dim2+para.dim3,para.dim1+para.dim2+para.dim3);
    JMatrix.block(para.dim3,para.dim3,para.dim1,para.dim1) = JMatrix_1;
    JMatrix.block(para.dim3+para.dim1,para.dim3+para.dim1,para.dim2,para.dim2) = JMatrix_2;
    JMatrix.block(0,0,para.dim3,para.dim3) = JMatrix_3;
    cout << "JMatrix done" << endl;
    cout << "JMatrix = " << JMatrix << endl;
    // throw runtime_error("JMatrix done");
    return JMatrix;
}

/**** dl_dtheta with Real Data ****/
tuple<ArrayXXd,int,int,int,int,int> alias::EstimationIndirectInference_SD_dl_dtheta_Part1_RealData(const ArrayXd & theta,
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    // *********
    /** Part 1 of the full model estimation **/
    std::vector<double> theta_part1(para.dim1);
    for (size_t n = 0; n < para.dim1; ++n) {theta_part1[n] = theta(para.dim3+n);}

    int GoodState = 1; int LaborIntensive = 1;
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Good_LaborInt = EstimationAuxiliaryModel_part1_ObjFun(
        theta_part1,ReadData_GoodState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Good_LaborInt = get<0>(t_Part1_Good_LaborInt);
    std::vector<double> grad_a_Part1_sum_Good_LaborInt = get<1>(t_Part1_Good_LaborInt);
    int Ncount_part1_Good_LaborInt = get<2>(t_Part1_Good_LaborInt);
    int Ncount_firms_part1_Good_LaborInt = get<3>(t_Part1_Good_LaborInt);
    ArrayXXd dlogL_d0_part1_Good_LaborInt = get<4>(t_Part1_Good_LaborInt);
    ArrayXXd dlogL_d1_part1_Good_LaborInt = ArrayXXd::Zero(dlogL_d0_part1_Good_LaborInt.rows(),dlogL_d0_part1_Good_LaborInt.cols());
    ArrayXXd dlogL_d2_part1_Good_LaborInt = get<5>(t_Part1_Good_LaborInt);
    ArrayXXd dlogL_d3_part1_Good_LaborInt = get<6>(t_Part1_Good_LaborInt);
    // cout << "part1_Good_LaborInt size = " << dlogL_d0_part1_Good_LaborInt.rows() << "; " << dlogL_d0_part1_Good_LaborInt.cols()
    //     << "; " << dlogL_d2_part1_Good_LaborInt.rows() << "; " << dlogL_d2_part1_Good_LaborInt.cols()
    //     << "; " << dlogL_d3_part1_Good_LaborInt.rows() << "; " << dlogL_d3_part1_Good_LaborInt.cols() << endl;

    GoodState = 0; LaborIntensive = 1;
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Bad_LaborInt = EstimationAuxiliaryModel_part1_ObjFun(
        theta_part1,ReadData_BadState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Bad_LaborInt = get<0>(t_Part1_Bad_LaborInt);
    std::vector<double> grad_a_Part1_sum_Bad_LaborInt = get<1>(t_Part1_Bad_LaborInt);
    int Ncount_part1_Bad_LaborInt = get<2>(t_Part1_Bad_LaborInt);
    int Ncount_firms_part1_Bad_LaborInt = get<3>(t_Part1_Bad_LaborInt);
    ArrayXXd dlogL_d1_part1_Bad_LaborInt = get<4>(t_Part1_Bad_LaborInt);
    ArrayXXd dlogL_d0_part1_Bad_LaborInt = ArrayXXd::Zero(dlogL_d1_part1_Bad_LaborInt.rows(),dlogL_d1_part1_Bad_LaborInt.cols());
    ArrayXXd dlogL_d2_part1_Bad_LaborInt = get<5>(t_Part1_Bad_LaborInt);
    ArrayXXd dlogL_d3_part1_Bad_LaborInt = get<6>(t_Part1_Bad_LaborInt);
    // cout << "part1_Bad_LaborInt size = " << dlogL_d1_part1_Bad_LaborInt.rows() << "; " << dlogL_d1_part1_Bad_LaborInt.cols()
    //     << "; " << dlogL_d2_part1_Bad_LaborInt.rows() << "; " << dlogL_d2_part1_Bad_LaborInt.cols()
    //     << "; " << dlogL_d3_part1_Bad_LaborInt.rows() << "; " << dlogL_d3_part1_Bad_LaborInt.cols() << endl;

     GoodState = 1; LaborIntensive = 0;
     tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Good_CapitalInt = EstimationAuxiliaryModel_part1_ObjFun(
         theta_part1,ReadData_GoodState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
     double y_part1_sum_Good_CapitalInt = get<0>(t_Part1_Good_CapitalInt);
     std::vector<double> grad_a_Part1_sum_Good_CapitalInt = get<1>(t_Part1_Good_CapitalInt);
     int Ncount_part1_Good_CapitalInt = get<2>(t_Part1_Good_CapitalInt);
     int Ncount_firms_part1_Good_CapitalInt = get<3>(t_Part1_Good_CapitalInt);
     ArrayXXd dlogL_d0_part1_Good_CapitalInt = get<4>(t_Part1_Good_CapitalInt);
     ArrayXXd dlogL_d1_part1_Good_CapitalInt = ArrayXXd::Zero(dlogL_d0_part1_Good_CapitalInt.rows(),dlogL_d0_part1_Good_CapitalInt.cols());
     ArrayXXd dlogL_d2_part1_Good_CapitalInt = get<5>(t_Part1_Good_CapitalInt);
     ArrayXXd dlogL_d3_part1_Good_CapitalInt = get<6>(t_Part1_Good_CapitalInt);
     // cout << "part1_Good_CapitalInt size = " << dlogL_d0_part1_Good_CapitalInt.rows() << "; " << dlogL_d0_part1_Good_CapitalInt.cols()
     //     << "; " << dlogL_d2_part1_Good_CapitalInt.rows() << "; " << dlogL_d2_part1_Good_CapitalInt.cols()
     //     << "; " << dlogL_d3_part1_Good_CapitalInt.rows() << "; " << dlogL_d3_part1_Good_CapitalInt.cols() << endl;

     GoodState = 0; LaborIntensive = 0;
     tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Bad_CapitalInt = EstimationAuxiliaryModel_part1_ObjFun(
         theta_part1,ReadData_BadState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
     double y_part1_sum_Bad_CapitalInt = get<0>(t_Part1_Bad_CapitalInt);
     std::vector<double> grad_a_Part1_sum_Bad_CapitalInt = get<1>(t_Part1_Bad_CapitalInt);
     int Ncount_part1_Bad_CapitalInt = get<2>(t_Part1_Bad_CapitalInt);
     int Ncount_firms_part1_Bad_CapitalInt = get<3>(t_Part1_Bad_CapitalInt);
     ArrayXXd dlogL_d1_part1_Bad_CapitalInt = get<4>(t_Part1_Bad_CapitalInt);
     ArrayXXd dlogL_d0_part1_Bad_CapitalInt = ArrayXXd::Zero(dlogL_d1_part1_Bad_CapitalInt.rows(),dlogL_d1_part1_Bad_CapitalInt.cols());
     ArrayXXd dlogL_d2_part1_Bad_CapitalInt = get<5>(t_Part1_Bad_CapitalInt);
     ArrayXXd dlogL_d3_part1_Bad_CapitalInt = get<6>(t_Part1_Bad_CapitalInt);
     // cout << "part1_Bad_CapitalInt size = " << dlogL_d1_part1_Bad_CapitalInt.rows() << "; " << dlogL_d1_part1_Bad_CapitalInt.cols()
     //     << "; " << dlogL_d2_part1_Bad_CapitalInt.rows() << "; " << dlogL_d2_part1_Bad_CapitalInt.cols()
     //     << "; " << dlogL_d3_part1_Bad_CapitalInt.rows() << "; " << dlogL_d3_part1_Bad_CapitalInt.cols() << endl;

    ArrayXXd dlogL_d0_part1 = AppendMatrix(dlogL_d0_part1_Good_LaborInt,dlogL_d0_part1_Bad_LaborInt,
        dlogL_d0_part1_Good_CapitalInt,dlogL_d0_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d1_part1 = AppendMatrix(dlogL_d1_part1_Good_LaborInt,dlogL_d1_part1_Bad_LaborInt,
        dlogL_d1_part1_Good_CapitalInt,dlogL_d1_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d2_part1 = AppendMatrix(dlogL_d2_part1_Good_LaborInt,dlogL_d2_part1_Bad_LaborInt,
        dlogL_d2_part1_Good_CapitalInt,dlogL_d2_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d3_part1 = AppendMatrix(dlogL_d3_part1_Good_LaborInt,dlogL_d3_part1_Bad_LaborInt,
        dlogL_d3_part1_Good_CapitalInt,dlogL_d3_part1_Bad_CapitalInt);
    // cout << "para.EstTbar = " << para.EstTbar << "; dlogL_d0_part1; dlogL_d1_part1; dlogL_d2_part1 done" << endl;
    // cout << "dlogL_d0_part1.rows() = " << dlogL_d0_part1.rows() << "; dlogL_d1_part1.rows() = " << dlogL_d1_part1.rows()
    //     << "; dlogL_d2_part1.rows() = " << dlogL_d2_part1.rows() << "; dlogL_d3_part1.rows() = " << dlogL_d3_part1.rows() << endl;

    ArrayXXd dlogL_part1 = ArrayXXd::Zero(dlogL_d0_part1.rows(),para.EstTbar*para.dim1);
    dlogL_part1.middleCols(0,para.EstTbar) = dlogL_d0_part1;
    dlogL_part1.middleCols(para.EstTbar,para.EstTbar) = dlogL_d1_part1;
    dlogL_part1.middleCols(2*para.EstTbar,para.EstTbar) = dlogL_d2_part1;
    dlogL_part1.middleCols(3*para.EstTbar,para.EstTbar) = dlogL_d3_part1;
    // cout << "dlogL_part1 done; size dlogL_part1 = " << dlogL_part1.rows() << "; " << dlogL_part1.cols() << endl;

    int Ncount_part1 = Ncount_part1_Good_LaborInt + Ncount_part1_Bad_LaborInt + Ncount_part1_Good_CapitalInt
        + Ncount_part1_Bad_CapitalInt;
    return tuple<ArrayXXd,int,int,int,int,int>(dlogL_part1,Ncount_part1,Ncount_firms_part1_Good_LaborInt,
        Ncount_firms_part1_Bad_LaborInt,Ncount_firms_part1_Good_CapitalInt,Ncount_firms_part1_Bad_CapitalInt);
}


tuple<ArrayXXd,int,int,int,int> alias::EstimationIndirectInference_SD_dl_dtheta_Part2_RealData(const ArrayXd & theta,
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;
    // *********
    /** Part 2 of the full model estimation **/
    std::vector<double> theta_part1(para.dim1);
    for (size_t k = 0; k < para.dim1; ++k) {theta_part1[k] = theta(para.dim3+k);}
    std::vector<double> theta_part2(para.dim2);
    for (size_t n = 0; n < para.dim2; ++n) {theta_part2[n] = theta(para.dim3+para.dim1+n);}

    GoodState = 1; LaborIntensive = 1;
    ParaEst1 para_est1_Good_LaborInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_LaborInt = Backout_phi_a_PD(ReadData_GoodState_LaborInt.Revenue_P,
        ReadData_GoodState_LaborInt.Capital_P, ReadData_GoodState_LaborInt.Employ_ur_P,
        ReadData_GoodState_LaborInt.Employ_uc_P,
        para_est1_Good_LaborInt.alpha_tilde_K, para_est1_Good_LaborInt.alpha_tilde_L,
        para_est1_Good_LaborInt.alpha_Lr, para_est1_Good_LaborInt.alpha_Lc);
    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Good_LaborInt
        = EstimationAuxiliaryModel_part2_ObjFun(theta_part2,para_est1_Good_LaborInt,ReadData_GoodState_LaborInt,
        lnphi_a_GoodState_LaborInt);
    double y_part2_sum_Good_LaborInt = get<0>(t_Part2_Good_LaborInt);
    std::vector<double> grad_a_Part2_sum_Good_LaborInt = get<1>(t_Part2_Good_LaborInt);
    int Ncount_part2_Good_LaborInt = get<2>(t_Part2_Good_LaborInt);
    int NcountP_firms_part2_Good_LaborInt = get<3>(t_Part2_Good_LaborInt);
    int NcountD_firms_part2_Good_LaborInt = get<4>(t_Part2_Good_LaborInt);
    int Ncount_firms_part2_Good_LaborInt = get<5>(t_Part2_Good_LaborInt);
    ArrayXXd dlogL_d0_part2_Good_LaborInt = get<6>(t_Part2_Good_LaborInt);
    ArrayXXd dlogL_d1_part2_Good_LaborInt = get<7>(t_Part2_Good_LaborInt);
    ArrayXXd dlogL_d2_part2_Good_LaborInt = get<8>(t_Part2_Good_LaborInt);
     // cout << "part2_Good_LaborInt size = " << dlogL_d0_part2_Good_LaborInt.rows() << "; " << dlogL_d0_part2_Good_LaborInt.cols()
     //     << "; " << dlogL_d1_part2_Good_LaborInt.rows() << "; " << dlogL_d1_part2_Good_LaborInt.cols()
     //     << "; " << dlogL_d2_part2_Good_LaborInt.rows() << "; " << dlogL_d2_part2_Good_LaborInt.cols() << endl;

    GoodState = 0; LaborIntensive = 1;
    ParaEst1 para_est1_Bad_LaborInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_BadState_LaborInt = Backout_phi_a_PD(ReadData_BadState_LaborInt.Revenue_P,
        ReadData_BadState_LaborInt.Capital_P, ReadData_BadState_LaborInt.Employ_ur_P,
        ReadData_BadState_LaborInt.Employ_uc_P,
        para_est1_Bad_LaborInt.alpha_tilde_K, para_est1_Bad_LaborInt.alpha_tilde_L,
        para_est1_Bad_LaborInt.alpha_Lr, para_est1_Bad_LaborInt.alpha_Lc);
    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Bad_LaborInt
        = EstimationAuxiliaryModel_part2_ObjFun(theta_part2,para_est1_Bad_LaborInt,ReadData_BadState_LaborInt,
        lnphi_a_BadState_LaborInt);
    double y_part2_sum_Bad_LaborInt = get<0>(t_Part2_Bad_LaborInt);
    std::vector<double> grad_a_Part2_sum_Bad_LaborInt = get<1>(t_Part2_Bad_LaborInt);
    int Ncount_part2_Bad_LaborInt = get<2>(t_Part2_Bad_LaborInt);
    int NcountP_firms_part2_Bad_LaborInt = get<3>(t_Part2_Bad_LaborInt);
    int NcountD_firms_part2_Bad_LaborInt = get<4>(t_Part2_Bad_LaborInt);
    int Ncount_firms_part2_Bad_LaborInt = get<5>(t_Part2_Bad_LaborInt);
    ArrayXXd dlogL_d0_part2_Bad_LaborInt = get<6>(t_Part2_Bad_LaborInt);
    ArrayXXd dlogL_d1_part2_Bad_LaborInt = get<7>(t_Part2_Bad_LaborInt);
    ArrayXXd dlogL_d2_part2_Bad_LaborInt = get<8>(t_Part2_Bad_LaborInt);
     // cout << "part2_Bad_LaborInt size = " << dlogL_d0_part2_Bad_LaborInt.rows() << "; " << dlogL_d0_part2_Bad_LaborInt.cols()
     //     << "; " << dlogL_d1_part2_Bad_LaborInt.rows() << "; " << dlogL_d1_part2_Bad_LaborInt.cols()
     //     << "; " << dlogL_d2_part2_Bad_LaborInt.rows() << "; " << dlogL_d2_part2_Bad_LaborInt.cols() << endl;

    GoodState = 1; LaborIntensive = 0;
    ParaEst1 para_est1_Good_CapitalInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_CapitalInt = Backout_phi_a_PD(ReadData_GoodState_CapitalInt.Revenue_P,
        ReadData_GoodState_CapitalInt.Capital_P, ReadData_GoodState_CapitalInt.Employ_ur_P,
        ReadData_GoodState_CapitalInt.Employ_uc_P,
        para_est1_Good_CapitalInt.alpha_tilde_K, para_est1_Good_CapitalInt.alpha_tilde_L,
        para_est1_Good_CapitalInt.alpha_Lr, para_est1_Good_CapitalInt.alpha_Lc);
    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Good_CapitalInt
        = EstimationAuxiliaryModel_part2_ObjFun(theta_part2,para_est1_Good_CapitalInt,ReadData_GoodState_CapitalInt,
        lnphi_a_GoodState_CapitalInt);
    double y_part2_sum_Good_CapitalInt = get<0>(t_Part2_Good_CapitalInt);
    std::vector<double> grad_a_Part2_sum_Good_CapitalInt = get<1>(t_Part2_Good_CapitalInt);
    int Ncount_part2_Good_CapitalInt = get<2>(t_Part2_Good_CapitalInt);
    int NcountP_firms_part2_Good_CapitalInt = get<3>(t_Part2_Good_CapitalInt);
    int NcountD_firms_part2_Good_CapitalInt = get<4>(t_Part2_Good_CapitalInt);
    int Ncount_firms_part2_Good_CapitalInt = get<5>(t_Part2_Good_CapitalInt);
    ArrayXXd dlogL_d0_part2_Good_CapitalInt = get<6>(t_Part2_Good_CapitalInt);
    ArrayXXd dlogL_d1_part2_Good_CapitalInt = get<7>(t_Part2_Good_CapitalInt);
    ArrayXXd dlogL_d2_part2_Good_CapitalInt = get<8>(t_Part2_Good_CapitalInt);
     // cout << "part2_Good_CapitalInt size = " << dlogL_d0_part2_Good_CapitalInt.rows() << "; " << dlogL_d0_part2_Good_CapitalInt.cols()
     //     << "; " << dlogL_d1_part2_Good_CapitalInt.rows() << "; " << dlogL_d1_part2_Good_CapitalInt.cols()
     //     << "; " << dlogL_d2_part2_Good_CapitalInt.rows() << "; " << dlogL_d2_part2_Good_CapitalInt.cols() << endl;

    GoodState = 0; LaborIntensive = 0;
    ParaEst1 para_est1_Bad_CapitalInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_BadState_CapitalInt = Backout_phi_a_PD(ReadData_BadState_CapitalInt.Revenue_P,
        ReadData_BadState_CapitalInt.Capital_P, ReadData_BadState_CapitalInt.Employ_ur_P,
        ReadData_BadState_CapitalInt.Employ_uc_P,
        para_est1_Bad_CapitalInt.alpha_tilde_K, para_est1_Bad_CapitalInt.alpha_tilde_L,
        para_est1_Bad_CapitalInt.alpha_Lr, para_est1_Bad_CapitalInt.alpha_Lc);
    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Bad_CapitalInt
        = EstimationAuxiliaryModel_part2_ObjFun(theta_part2,para_est1_Bad_CapitalInt,ReadData_BadState_CapitalInt,
        lnphi_a_BadState_CapitalInt);
    double y_part2_sum_Bad_CapitalInt = get<0>(t_Part2_Bad_CapitalInt);
    std::vector<double> grad_a_Part2_sum_Bad_CapitalInt = get<1>(t_Part2_Bad_CapitalInt);
    int Ncount_part2_Bad_CapitalInt = get<2>(t_Part2_Bad_CapitalInt);
    int NcountP_firms_part2_Bad_CapitalInt = get<3>(t_Part2_Bad_CapitalInt);
    int NcountD_firms_part2_Bad_CapitalInt = get<4>(t_Part2_Bad_CapitalInt);
    int Ncount_firms_part2_Bad_CapitalInt = get<5>(t_Part2_Bad_CapitalInt);
    ArrayXXd dlogL_d0_part2_Bad_CapitalInt = get<6>(t_Part2_Bad_CapitalInt);
    ArrayXXd dlogL_d1_part2_Bad_CapitalInt = get<7>(t_Part2_Bad_CapitalInt);
    ArrayXXd dlogL_d2_part2_Bad_CapitalInt = get<8>(t_Part2_Bad_CapitalInt);
     // cout << "part2_Bad_CapitalInt size = " << dlogL_d0_part2_Bad_CapitalInt.rows() << "; " << dlogL_d0_part2_Bad_CapitalInt.cols()
     //     << "; " << dlogL_d1_part2_Bad_CapitalInt.rows() << "; " << dlogL_d1_part2_Bad_CapitalInt.cols()
     //     << "; " << dlogL_d2_part2_Bad_CapitalInt.rows() << "; " << dlogL_d2_part2_Bad_CapitalInt.cols() << endl;

     ArrayXXd dlogL_d0_part2 = AppendMatrix(dlogL_d0_part2_Good_LaborInt,dlogL_d0_part2_Bad_LaborInt,
         dlogL_d0_part2_Good_CapitalInt,dlogL_d0_part2_Bad_CapitalInt);
     ArrayXXd dlogL_d1_part2 = AppendMatrix(dlogL_d1_part2_Good_LaborInt,dlogL_d1_part2_Bad_LaborInt,
         dlogL_d1_part2_Good_CapitalInt,dlogL_d1_part2_Bad_CapitalInt);
     ArrayXXd dlogL_d2_part2 = AppendMatrix(dlogL_d2_part2_Good_LaborInt,dlogL_d2_part2_Bad_LaborInt,
         dlogL_d2_part2_Good_CapitalInt,dlogL_d2_part2_Bad_CapitalInt);

     // cout << "para.EstTbar = " << para.EstTbar << "; dlogL_d0_part2; dlogL_d1_part2; dlogL_d2_part2 done" << endl;
     // cout << "dlogL_d0_part2.rows() = " << dlogL_d0_part2.rows() << "; dlogL_d1_part2.rows() = " << dlogL_d1_part2.rows()
     //     << "; dlogL_d2_part2.rows() = " << dlogL_d2_part2.rows() << endl;

     ArrayXXd dlogL_part2 = ArrayXXd::Zero(dlogL_d0_part2.rows(),para.EstTbar*para.dim2);
     dlogL_part2 << dlogL_d0_part2,dlogL_d1_part2,dlogL_d2_part2;
     // cout << "dlogL_part2 done; size dlogL_part2 = " << dlogL_part2.rows() << "; " << dlogL_part2.cols() << endl;

    int Ncount_part2 = Ncount_part2_Good_LaborInt + Ncount_part2_Bad_LaborInt + Ncount_part2_Good_CapitalInt
        + Ncount_part2_Bad_CapitalInt;
    int NcountP_firms_part2 = NcountP_firms_part2_Good_LaborInt + NcountP_firms_part2_Bad_LaborInt
        + NcountP_firms_part2_Good_CapitalInt + NcountP_firms_part2_Bad_CapitalInt;
    int NcountD_firms_part2 = NcountD_firms_part2_Good_LaborInt + NcountD_firms_part2_Bad_LaborInt
        + NcountD_firms_part2_Good_CapitalInt + NcountD_firms_part2_Bad_CapitalInt;
    int Ncount_firms_part2 = Ncount_firms_part2_Good_LaborInt + Ncount_firms_part2_Bad_LaborInt
        + Ncount_firms_part2_Good_CapitalInt + Ncount_firms_part2_Bad_CapitalInt;
    return tuple<ArrayXXd,int,int,int,int>(dlogL_part2,Ncount_part2,NcountP_firms_part2,NcountD_firms_part2,Ncount_firms_part2);
}

//
// ArrayXd alias::EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(const ArrayXd & theta_a,
//     const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
//     const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
//     MultiThreads::Threads_Management &threadsManagement) {
//
//     ArrayXd dlogL_Part3_dtheta_a(para.dim3);
//
//     int GoodState; int LaborIntensive;
//     for (size_t kk = 0; kk < 1; kk++) {
//         double eps = abs(theta_a(kk)) * 1e-2;
//
//         ArrayXd theta_a_plus = theta_a;
//         theta_a_plus(kk) = theta_a(kk) + eps;
//
//         GoodState = 1; LaborIntensive = 1;
//         tuple<ArrayXXd,int> t_part3_Good_LaborInt_plus = Cal_logL_Part3_RealData(ReadData_GoodState_LaborInt,
//             theta_a_plus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Good_LaborInd_plus = get<0>(t_part3_Good_LaborInt_plus);
//         GoodState = 0; LaborIntensive = 1;
//         tuple<ArrayXXd,int> t_part3_Bad_LaborInt_plus = Cal_logL_Part3_RealData(ReadData_BadState_LaborInt,
//             theta_a_plus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Bad_LaborInd_plus = get<0>(t_part3_Bad_LaborInt_plus);
//         GoodState = 1; LaborIntensive = 0;
//         tuple<ArrayXXd,int> t_part3_Good_CapitalInt_plus = Cal_logL_Part3_RealData(ReadData_GoodState_CapitalInt,
//             theta_a_plus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Good_CapitalInd_plus = get<0>(t_part3_Good_CapitalInt_plus);
//         GoodState = 0; LaborIntensive = 0;
//         tuple<ArrayXXd,int> t_part3_Bad_CapitalInt_plus = Cal_logL_Part3_RealData(ReadData_BadState_CapitalInt,
//             theta_a_plus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Bad_CapitalInd_plus = get<0>(t_part3_Bad_CapitalInt_plus);
//
//         ArrayXXd logL_part3_plus = AppendMatrix(logL_Part3_Good_LaborInd_plus,logL_Part3_Bad_LaborInd_plus,
//             logL_Part3_Good_CapitalInd_plus,logL_Part3_Bad_CapitalInd_plus);
//
//         ArrayXd theta_a_minus = theta_a;
//         theta_a_minus(kk) = theta_a(kk) - eps;
//
//         GoodState = 1; LaborIntensive = 1;
//         tuple<ArrayXXd,int> t_part3_Good_LaborInt_minus = Cal_logL_Part3_RealData(ReadData_GoodState_LaborInt,
//             theta_a_minus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Good_LaborInd_minus = get<0>(t_part3_Good_LaborInt_minus);
//         GoodState = 0; LaborIntensive = 1;
//         tuple<ArrayXXd,int> t_part3_Bad_LaborInt_minus = Cal_logL_Part3_RealData(ReadData_BadState_LaborInt,
//             theta_a_minus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Bad_LaborInd_minus = get<0>(t_part3_Bad_LaborInt_minus);
//         GoodState = 1; LaborIntensive = 0;
//         tuple<ArrayXXd,int> t_part3_Good_CapitalInt_minus = Cal_logL_Part3_RealData(ReadData_GoodState_CapitalInt,
//             theta_a_minus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Good_CapitalInd_minus = get<0>(t_part3_Good_CapitalInt_minus);
//         GoodState = 0; LaborIntensive = 0;
//         tuple<ArrayXXd,int> t_part3_Bad_CapitalInt_minus = Cal_logL_Part3_RealData(ReadData_BadState_CapitalInt,
//             theta_a_minus, GoodState, LaborIntensive, threadsManagement);
//         ArrayXXd logL_Part3_Bad_CapitalInd_minus = get<0>(t_part3_Bad_CapitalInt_minus);
//
//         ArrayXXd logL_part3_minus = AppendMatrix(logL_Part3_Good_LaborInd_minus,logL_Part3_Bad_LaborInd_minus,
//             logL_Part3_Good_CapitalInd_minus,logL_Part3_Bad_CapitalInd_minus);
//
//         cout << ( (logL_part3_plus - logL_part3_minus).abs() > 0.0 ).cast<int>().sum() << endl;
//         cout << ( (logL_part3_plus - logL_part3_minus).abs() > 0.0 ).sum() << endl;
//         ArrayXi index_temp = ( (logL_part3_plus - logL_part3_minus).abs() > 0.0 ).cast<int>().rowwise().sum();
//         ArrayXi index = (index_temp > 0).cast<int>();
//
//         dlogL_Part3_dtheta_a(kk) = (logL_part3_plus - logL_part3_minus).sum() / 2.0 / eps;
//         cout << "dlogL_Part3_dtheta_a(kk) = " << dlogL_Part3_dtheta_a(kk) << endl;
//         cout << "index sum = " << index.sum() << endl;
//         cout << "logL_part3_minus = " << logL_part3_minus.rows() << endl;
//         // throw runtime_error("638");
//     }
//     return dlogL_Part3_dtheta_a;
// }
//
//
tuple<ArrayXXd,int> alias::Cal_logL_Part3_RealData(const SimData & ReadData, const ArrayXd & theta_a,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement) {

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

    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
    ParaEst para_est = get<0>(t_Equ);
    ParaVec para_vec = get<1>(t_Equ);
    EquStateV EquV_a = get<2>(t_Equ);
    EquStateVmat Evalmat_a = get<3>(t_Equ);

    tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi> t_logL = EstimationAuxiliaryModel_part3_likelihood(
        ReadData, para_est, para_vec, EquV_a, Evalmat_a, GoodState, LaborIntensive, threadsManagement);
    int N_logL_part3 = get<1>(t_logL);
    ArrayXXd logL_part3 = get<2>(t_logL);

    return tuple<ArrayXXd,int> (logL_part3,N_logL_part3);
}

ArrayXXd alias::EstimationIndirectInference_SD_l_Part3_RealData(const ArrayXd & theta_a,
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;
    GoodState = 1; LaborIntensive = 1;
    tuple<ArrayXXd,int> t_part3_Good_LaborInt = Cal_logL_Part3_RealData(ReadData_GoodState_LaborInt, theta_a,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_Part3_Good_LaborInd = get<0>(t_part3_Good_LaborInt);
    GoodState = 0; LaborIntensive = 1;
    tuple<ArrayXXd,int> t_part3_Bad_LaborInt = Cal_logL_Part3_RealData(ReadData_BadState_LaborInt, theta_a,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_Part3_Bad_LaborInd = get<0>(t_part3_Bad_LaborInt);
    GoodState = 1; LaborIntensive = 0;
    tuple<ArrayXXd,int> t_part3_Good_CapitalInt = Cal_logL_Part3_RealData(ReadData_GoodState_CapitalInt, theta_a,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_Part3_Good_CapitalInd = get<0>(t_part3_Good_CapitalInt);
    GoodState = 0; LaborIntensive = 0;
    tuple<ArrayXXd,int> t_part3_Bad_CapitalInt = Cal_logL_Part3_RealData(ReadData_BadState_CapitalInt, theta_a,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_Part3_Bad_CapitalInd = get<0>(t_part3_Bad_CapitalInt);

    ArrayXXd logL_part3 = AppendMatrix(logL_Part3_Good_LaborInd,logL_Part3_Bad_LaborInd,
        logL_Part3_Good_CapitalInd,logL_Part3_Bad_CapitalInd);

    return logL_part3;
}



tuple<ArrayXXd,int,int> alias::EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    const EquV_Auxiliary & EquV_GoodState_LaborInt, const EquV_Auxiliary & EquV_BadState_LaborInt,
    const EquV_Auxiliary & EquV_GoodState_CapitalInt, const EquV_Auxiliary & EquV_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;
    // *********
    /** Part 3 of the full model estimation **/
    GoodState = 1; LaborIntensive = 1;
    tuple<ArrayXXd,int> t_part3_Good_LaborInt = Cal_dlogL_Part3_RealData(ReadData_GoodState_LaborInt,
        EquV_GoodState_LaborInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Good_LaborInt = get<0>(t_part3_Good_LaborInt);
    int Ncount_part3_Good_LaborInt = dlogL_part3_Good_LaborInt.rows();
    // cout << "part3_Good_LaborInt size = " << dlogL_part3_Good_LaborInt.rows() << "; " << dlogL_part3_Good_LaborInt.cols() << endl;
    GoodState = 0; LaborIntensive = 1;
    tuple<ArrayXXd,int> t_part3_Bad_LaborInt = Cal_dlogL_Part3_RealData(ReadData_BadState_LaborInt,
        EquV_BadState_LaborInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Bad_LaborInt = get<0>(t_part3_Bad_LaborInt);
    int Ncount_part3_Bad_LaborInt = dlogL_part3_Bad_LaborInt.rows();
     // cout << "part3_Bad_LaborInt size = " << dlogL_part3_Bad_LaborInt.rows() << "; " << dlogL_part3_Bad_LaborInt.cols() << endl;
    GoodState = 1; LaborIntensive = 0;
    tuple<ArrayXXd,int> t_part3_Good_CapitalInt = Cal_dlogL_Part3_RealData(ReadData_GoodState_CapitalInt,
        EquV_GoodState_CapitalInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Good_CapitalInt = get<0>(t_part3_Good_CapitalInt);
    int Ncount_part3_Good_CapitalInt = dlogL_part3_Good_CapitalInt.rows();
    // cout << "part3_Good_CapitalInt size = " << dlogL_part3_Good_CapitalInt.rows() << "; " << dlogL_part3_Good_CapitalInt.cols() << endl;
    GoodState = 0; LaborIntensive = 0;
    tuple<ArrayXXd,int> t_part3_Bad_CapitalInt = Cal_dlogL_Part3_RealData(ReadData_BadState_CapitalInt,
        EquV_BadState_CapitalInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Bad_CapitalInt = get<0>(t_part3_Bad_CapitalInt);
    int Ncount_part3_Bad_CapitalInt = dlogL_part3_Bad_CapitalInt.rows();
    // cout << "part3_Bad_CapitalInt size = " << dlogL_part3_Bad_CapitalInt.rows() << "; " << dlogL_part3_Bad_CapitalInt.cols() << endl;

    ArrayXXd dlogL_part3 = AppendMatrix(dlogL_part3_Good_LaborInt,dlogL_part3_Bad_LaborInt,
        dlogL_part3_Good_CapitalInt,dlogL_part3_Bad_CapitalInt);
    int Ncount_part3 = Ncount_part3_Good_LaborInt + Ncount_part3_Bad_LaborInt + Ncount_part3_Good_CapitalInt
        + Ncount_part3_Bad_CapitalInt;
    int Ncount_part3_Good = Ncount_part3_Good_LaborInt + Ncount_part3_Good_CapitalInt;
    int Ncount_part3_Bad = Ncount_part3_Bad_LaborInt + Ncount_part3_Bad_CapitalInt;

    return tuple<ArrayXXd,int,int>(dlogL_part3,Ncount_part3_Good,Ncount_part3_Bad);
}

tuple<ArrayXXd,int> alias::Cal_dlogL_Part3_RealData(const SimData & ReadData, const EquV_Auxiliary & EquV,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement) {

    tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi> t_Part3 = EstimationAuxiliaryModel_part3_likelihood(
        ReadData, EquV.para_est_a_point, EquV.para_vec_a_point, EquV.EquV_a_point, EquV.Evalmat_a_point,
        GoodState, LaborIntensive, threadsManagement);
    double obj_point_sum = get<0>(t_Part3);
    int Ncount_part3 = get<1>(t_Part3);
    ArrayXXd logL_part3 = get<2>(t_Part3);
    ArrayXXi K_index_max_mat = get<3>(t_Part3);
    ArrayXXi Lur_index_max_mat = get<4>(t_Part3);
    ArrayXXi Luc_index_max_mat = get<5>(t_Part3);

    ArrayXXd dlogL_part3 = ArrayXXd::Zero(ReadData.good.rows(),para.EstTbar * para.dim3);
    ArrayXd grad_a3_sum = ArrayXd::Zero(para.dim3);
    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
    //        cout << "i_a = " << i_a << endl;
        double obj_plus_sum = 0.0;
        if ( ( (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15) and GoodState == 1 ) or
            ( (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14) and GoodState == 0 ) ) {
            obj_plus_sum = 0.0;
            grad_a3_sum(i_a) = 0.0;
        }
        else {
            tuple<double,int,ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_likelihood_Diff(ReadData,
                EquV.para_est_a_plus[i_a], EquV.para_vec_a_plus[i_a], EquV.EquV_a_plus[i_a], EquV.Evalmat_a_plus[i_a],
                K_index_max_mat, Lur_index_max_mat, Luc_index_max_mat,
                GoodState, LaborIntensive, threadsManagement);
            obj_plus_sum = get<0>(t_plus);
            grad_a3_sum(i_a) = (obj_plus_sum - obj_point_sum) / EquV.eps_diff(i_a);
            ArrayXXd logL_part3_plus = get<2>(t_plus);
            dlogL_part3.middleCols(i_a*para.EstTbar,para.EstTbar) = (-logL_part3_plus - (-logL_part3)) / EquV.eps_diff(i_a);
        }
    }

    return tuple<ArrayXXd,int>(dlogL_part3,Ncount_part3);
}

ArrayXXd alias::EstimationIndirectInference_SD_IMatrix_RealData(const ArrayXd & theta,
    const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    tuple<ArrayXXd,int,int,int,int,int> t_dL_part1 = EstimationIndirectInference_SD_dl_dtheta_Part1_RealData(theta,
        ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd dL_part1 = get<0>(t_dL_part1);
    tuple<ArrayXXd,int,int,int,int> t_dL_part2 = EstimationIndirectInference_SD_dl_dtheta_Part2_RealData(theta,
        ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd dL_part2 = get<0>(t_dL_part2);

    int GoodState; int LaborIntensive;
    GoodState = 1; LaborIntensive = 1;
    EquV_Auxiliary EquV_GoodState_LaborInt = CalEquV_plus_minus(theta, GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 1;
    EquV_Auxiliary EquV_BadState_LaborInt = CalEquV_plus_minus(theta, GoodState, LaborIntensive, threadsManagement);
    GoodState = 1; LaborIntensive = 0;
    EquV_Auxiliary EquV_GoodState_CapitalInt = CalEquV_plus_minus(theta, GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 0;
    EquV_Auxiliary EquV_BadState_CapitalInt = CalEquV_plus_minus(theta, GoodState, LaborIntensive, threadsManagement);
    tuple<ArrayXXd,int,int> t_dL_part3 = EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(
        ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,
        ReadData_GoodState_CapitalInt,ReadData_BadState_CapitalInt,
        EquV_GoodState_LaborInt,EquV_BadState_LaborInt,EquV_GoodState_CapitalInt,EquV_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd dL_part3 = get<0>(t_dL_part3);

    ArrayXXd dL_AllPart = Vectorize_dl_dtheta(dL_part1,dL_part2,dL_part3);
    cout << "dL_AllPart size = " << dL_AllPart.rows() << "; " << dL_AllPart.cols() << endl;
    cout << "dL_AllPart = " << dL_AllPart.colwise().sum() << endl;

    ArrayXXd Var_dL_AllPart = VarianceMatrix_dL(dL_AllPart);
    cout << "Var_dL_AllPart size = " << Var_dL_AllPart.rows() << "; " << Var_dL_AllPart.cols() << endl;
    cout << "Var_dL_AllPart = " << Var_dL_AllPart << endl;
    //
    // for (size_t i_a = 0; i_a < Var_dL_AllPart.rows(); i_a++) {
    //     for (size_t j_a = 0; j_a < Var_dL_AllPart.rows(); j_a++) {
    //         if ( (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14 or i_a == 17) and (j_a == 6 or j_a == 7 or j_a == 10 or j_a == 15 or j_a == 18) ) {
    //             Var_dL_AllPart(i_a,j_a) = 0.0;
    //         }
    //         if ( (j_a == 4 or j_a == 5 or j_a == 9 or j_a == 14 or j_a == 17) and (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15 or i_a == 18) ) {
    //             Var_dL_AllPart(i_a,j_a) = 0.0;
    //         }
    //     }
    // }

    // throw runtime_error("Not implemented");
    return Var_dL_AllPart;
}
// //

ArrayXXd alias::EstimationIndirectInference_SD_IStarMatrix_SimData(const ArrayXd & theta_Est,const ArrayXd & theta_a,
    const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
    const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;
    GoodState = 1; LaborIntensive = 1;
    SimData sim_data_Para_GoodState_LaborInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 1;
    SimData sim_data_Para_BadState_LaborInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);
    GoodState = 1; LaborIntensive = 0;
    SimData sim_data_Para_GoodState_CapitalInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 0;
    SimData sim_data_Para_BadState_CapitalInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    cout << "EstimationIndirectInference_SD_IStarMatrix_SimData: Simulate the full model" << endl;

    tuple<ArrayXXd,ArrayXXd,double> t_part1 = EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(theta_a,
        sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd dL_part1 = get<0>(t_part1);
    double Nprob_part1 = get<2>(t_part1);
    tuple<ArrayXXd,ArrayXXd,double> t_part2 = EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(theta_a,
        sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd dL_part2 = get<0>(t_part2);
    double Nprob_part2 = get<2>(t_part2);

    GoodState = 1; LaborIntensive = 1;
    EquV_Auxiliary EquV_GoodState_LaborInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 1;
    EquV_Auxiliary EquV_BadState_LaborInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    GoodState = 1; LaborIntensive = 0;
    EquV_Auxiliary EquV_GoodState_CapitalInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 0;
    EquV_Auxiliary EquV_BadState_CapitalInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
    tuple<ArrayXXd,ArrayXXd,double> t_part3 = EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
        sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,
        sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
        EquV_GoodState_LaborInt,EquV_BadState_LaborInt,
        EquV_GoodState_CapitalInt,EquV_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd dL_part3 = get<0>(t_part3);

    ArrayXXd dL_AllPart = Vectorize_dl_dtheta(dL_part1,dL_part2,dL_part3);
    cout << "dL_AllPart = " << dL_AllPart.colwise().sum() << endl;

    ArrayXXd Var_dL_AllPart = VarianceMatrix_dL(dL_AllPart);

    // for (size_t i_a = 0; i_a < Var_dL_AllPart.rows(); i_a++) {
    //     for (size_t j_a = 0; j_a < Var_dL_AllPart.rows(); j_a++) {
    //         if ( (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14 or i_a == 17) and (j_a == 6 or j_a == 7 or j_a == 10 or j_a == 15 or j_a == 18) ) {
    //             Var_dL_AllPart(i_a,j_a) = 0.0;
    //         }
    //         if ( (j_a == 4 or j_a == 5 or j_a == 9 or j_a == 14 or j_a == 17) and (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15 or i_a == 18) ) {
    //             Var_dL_AllPart(i_a,j_a) = 0.0;
    //         }
    //     }
    // }
    // throw runtime_error("1103");
    return Var_dL_AllPart;
}


tuple<ArrayXXd,ArrayXXd,double> alias::EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(const ArrayXd & theta,
    const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
    const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    // *********
    /** Part 1 of the full model estimation **/
    std::vector<double> theta_part1(para.dim1);
    for (size_t n = 0; n < para.dim1; ++n) {theta_part1[n] = theta(para.dim3+n);}

    int GoodState = 1; int LaborIntensive = 1;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Good_LaborInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
        theta_part1,sim_data_Para_GoodState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Good_LaborInt = get<0>(t_Part1_Good_LaborInt);
    std::vector<double> grad_a_Part1_sum_Good_LaborInt = get<1>(t_Part1_Good_LaborInt);
    int Ncount_part1_Good_LaborInt = get<2>(t_Part1_Good_LaborInt);
    double Ncount_prob_part1_Good_LaborInt = get<3>(t_Part1_Good_LaborInt);
    ArrayXXd logL_part1_Good_LaborInt = get<4>(t_Part1_Good_LaborInt);
    ArrayXXd dlogL_d0_part1_Good_LaborInt = get<5>(t_Part1_Good_LaborInt);
    ArrayXXd dlogL_d1_part1_Good_LaborInt = ArrayXXd::Zero(dlogL_d0_part1_Good_LaborInt.rows(),dlogL_d0_part1_Good_LaborInt.cols());
    ArrayXXd dlogL_d2_part1_Good_LaborInt = get<6>(t_Part1_Good_LaborInt);
    ArrayXXd dlogL_d3_part1_Good_LaborInt = get<7>(t_Part1_Good_LaborInt);
    // cout << "part1_Good_LaborInt size = " << dlogL_d0_part1_Good_LaborInt.rows() << "; " << dlogL_d0_part1_Good_LaborInt.cols()
    //     << "; " << dlogL_d2_part1_Good_LaborInt.rows() << "; " << dlogL_d2_part1_Good_LaborInt.cols()
    //     << "; " << dlogL_d3_part1_Good_LaborInt.rows() << "; " << dlogL_d3_part1_Good_LaborInt.cols() << endl;
    // cout << "grad_a_Part1_sum_Good_LaborInt = " << grad_a_Part1_sum_Good_LaborInt[0]
    //     << "; " << grad_a_Part1_sum_Good_LaborInt[1] << "; " << grad_a_Part1_sum_Good_LaborInt[2] << endl;

    GoodState = 0; LaborIntensive = 1;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Bad_LaborInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
        theta_part1,sim_data_Para_BadState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Bad_LaborInt = get<0>(t_Part1_Bad_LaborInt);
    std::vector<double> grad_a_Part1_sum_Bad_LaborInt = get<1>(t_Part1_Bad_LaborInt);
    int Ncount_part1_Bad_LaborInt = get<2>(t_Part1_Bad_LaborInt);
    double Ncount_prob_part1_Bad_LaborInt = get<3>(t_Part1_Bad_LaborInt);
    ArrayXXd logL_part1_Bad_LaborInt = get<4>(t_Part1_Bad_LaborInt);
    ArrayXXd dlogL_d1_part1_Bad_LaborInt = get<5>(t_Part1_Bad_LaborInt);
    ArrayXXd dlogL_d0_part1_Bad_LaborInt = ArrayXXd::Zero(dlogL_d1_part1_Bad_LaborInt.rows(),dlogL_d1_part1_Bad_LaborInt.cols());
    ArrayXXd dlogL_d2_part1_Bad_LaborInt = get<6>(t_Part1_Bad_LaborInt);
    ArrayXXd dlogL_d3_part1_Bad_LaborInt = get<7>(t_Part1_Bad_LaborInt);
    // cout << "part1_Bad_LaborInt size = " << dlogL_d1_part1_Bad_LaborInt.rows() << "; " << dlogL_d1_part1_Bad_LaborInt.cols()
    //     << "; " << dlogL_d2_part1_Bad_LaborInt.rows() << "; " << dlogL_d2_part1_Bad_LaborInt.cols()
    //     << "; " << dlogL_d3_part1_Bad_LaborInt.rows() << "; " << dlogL_d3_part1_Bad_LaborInt.cols() << endl;
    // cout << "grad_a_Part1_sum_Bad_LaborInt = " << grad_a_Part1_sum_Bad_LaborInt[0]
    //     << "; " << grad_a_Part1_sum_Bad_LaborInt[1] << "; " << grad_a_Part1_sum_Bad_LaborInt[2] << endl;

    GoodState = 1; LaborIntensive = 0;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Good_CapitalInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
        theta_part1,sim_data_Para_GoodState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Good_CapitalInt = get<0>(t_Part1_Good_CapitalInt);
    std::vector<double> grad_a_Part1_sum_Good_CapitalInt = get<1>(t_Part1_Good_CapitalInt);
    int Ncount_part1_Good_CapitalInt = get<2>(t_Part1_Good_CapitalInt);
    double Ncount_prob_part1_Good_CapitalInt = get<3>(t_Part1_Good_CapitalInt);
    ArrayXXd logL_part1_Good_CapitalInt = get<4>(t_Part1_Good_CapitalInt);
    ArrayXXd dlogL_d0_part1_Good_CapitalInt = get<5>(t_Part1_Good_CapitalInt);
    ArrayXXd dlogL_d1_part1_Good_CapitalInt = ArrayXXd::Zero(dlogL_d0_part1_Good_CapitalInt.rows(),dlogL_d0_part1_Good_CapitalInt.cols());
    ArrayXXd dlogL_d2_part1_Good_CapitalInt = get<6>(t_Part1_Good_CapitalInt);
    ArrayXXd dlogL_d3_part1_Good_CapitalInt = get<7>(t_Part1_Good_CapitalInt);
     // cout << "part1_Good_CapitalInt size = " << dlogL_d0_part1_Good_CapitalInt.rows() << "; " << dlogL_d0_part1_Good_CapitalInt.cols()
     //     << "; " << dlogL_d2_part1_Good_CapitalInt.rows() << "; " << dlogL_d2_part1_Good_CapitalInt.cols()
     //     << "; " << dlogL_d3_part1_Good_CapitalInt.rows() << "; " << dlogL_d3_part1_Good_CapitalInt.cols() << endl;
    // cout << "grad_a_Part1_sum_Good_CapitalInt = " << grad_a_Part1_sum_Good_CapitalInt[0]
    //     << "; " << grad_a_Part1_sum_Good_CapitalInt[1] << "; " << grad_a_Part1_sum_Good_CapitalInt[2] << endl;

    GoodState = 0; LaborIntensive = 0;
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part1_Bad_CapitalInt
        = EstimationAuxiliaryModel_part1_ObjFun_FullVersion_Step2(
        theta_part1,sim_data_Para_BadState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
    double y_part1_sum_Bad_CapitalInt = get<0>(t_Part1_Bad_CapitalInt);
    std::vector<double> grad_a_Part1_sum_Bad_CapitalInt = get<1>(t_Part1_Bad_CapitalInt);
    int Ncount_part1_Bad_CapitalInt = get<2>(t_Part1_Bad_CapitalInt);
    double Ncount_prob_part1_Bad_CapitalInt = get<3>(t_Part1_Bad_CapitalInt);
    ArrayXXd logL_part1_Bad_CapitalInt = get<4>(t_Part1_Bad_CapitalInt);
    ArrayXXd dlogL_d1_part1_Bad_CapitalInt = get<5>(t_Part1_Bad_CapitalInt);
    ArrayXXd dlogL_d0_part1_Bad_CapitalInt = ArrayXXd::Zero(dlogL_d1_part1_Bad_CapitalInt.rows(),dlogL_d1_part1_Bad_CapitalInt.cols());
    ArrayXXd dlogL_d2_part1_Bad_CapitalInt = get<6>(t_Part1_Bad_CapitalInt);
    ArrayXXd dlogL_d3_part1_Bad_CapitalInt = get<7>(t_Part1_Bad_CapitalInt);
    // cout << "part1_Bad_CapitalInt size = " << dlogL_d1_part1_Bad_CapitalInt.rows() << "; " << dlogL_d1_part1_Bad_CapitalInt.cols()
    //     << "; " << dlogL_d2_part1_Bad_CapitalInt.rows() << "; " << dlogL_d2_part1_Bad_CapitalInt.cols()
    //     << "; " << dlogL_d3_part1_Bad_CapitalInt.rows() << "; " << dlogL_d3_part1_Bad_CapitalInt.cols() << endl;
    // cout << "grad_a_Part1_sum_Bad_CapitalInt = " << grad_a_Part1_sum_Bad_CapitalInt[0]
    //     << "; " << grad_a_Part1_sum_Bad_CapitalInt[1] << "; " << grad_a_Part1_sum_Bad_CapitalInt[2] << endl;

    double Nobs = double(Ncount_part1_Good_LaborInt+Ncount_part1_Bad_LaborInt+Ncount_part1_Good_CapitalInt+Ncount_part1_Bad_CapitalInt);
    double Ncount_prob = double(Ncount_prob_part1_Good_LaborInt+Ncount_prob_part1_Bad_LaborInt+Ncount_prob_part1_Good_CapitalInt+Ncount_prob_part1_Bad_CapitalInt);
    ArrayXd grad_a1(para.dim1);
    grad_a1(0) = (grad_a_Part1_sum_Good_LaborInt[0]+grad_a_Part1_sum_Good_CapitalInt[0]) / Nobs;
    grad_a1(1) = (grad_a_Part1_sum_Bad_LaborInt[0]+grad_a_Part1_sum_Bad_CapitalInt[0])  / Nobs;
    grad_a1(2) = (grad_a_Part1_sum_Good_LaborInt[1]+grad_a_Part1_sum_Bad_LaborInt[1]
        +grad_a_Part1_sum_Good_CapitalInt[1]+grad_a_Part1_sum_Bad_CapitalInt[1]) / Nobs;
    grad_a1(3) = (grad_a_Part1_sum_Good_LaborInt[2]+grad_a_Part1_sum_Bad_LaborInt[2]
        +grad_a_Part1_sum_Good_CapitalInt[2]+grad_a_Part1_sum_Bad_CapitalInt[2]) / Nobs;
    // cout << "grad_a1 = " << grad_a1.transpose() << "; sum = " << (grad_a1*grad_a1).sum() << endl;

    ArrayXXd logL_part1 = AppendMatrix(logL_part1_Good_LaborInt,logL_part1_Bad_LaborInt,
        logL_part1_Good_CapitalInt,logL_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d0_part1 = AppendMatrix(dlogL_d0_part1_Good_LaborInt,dlogL_d0_part1_Bad_LaborInt,
        dlogL_d0_part1_Good_CapitalInt,dlogL_d0_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d1_part1 = AppendMatrix(dlogL_d1_part1_Good_LaborInt,dlogL_d1_part1_Bad_LaborInt,
        dlogL_d1_part1_Good_CapitalInt,dlogL_d1_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d2_part1 = AppendMatrix(dlogL_d2_part1_Good_LaborInt,dlogL_d2_part1_Bad_LaborInt,
        dlogL_d2_part1_Good_CapitalInt,dlogL_d2_part1_Bad_CapitalInt);
    ArrayXXd dlogL_d3_part1 = AppendMatrix(dlogL_d3_part1_Good_LaborInt,dlogL_d3_part1_Bad_LaborInt,
        dlogL_d3_part1_Good_CapitalInt,dlogL_d3_part1_Bad_CapitalInt);
    // cout << "para.EstTbar = " << para.EstTbar << "; dlogL_d0_part1; dlogL_d1_part1; dlogL_d2_part1 done" << endl;
    // cout << "dlogL_d0_part1.rows() = " << dlogL_d0_part1.rows() << "; dlogL_d1_part1.rows() = " << dlogL_d1_part1.rows()
    //     << "; dlogL_d2_part1.rows() = " << dlogL_d2_part1.rows() << "; dlogL_d3_part1.rows() = " << dlogL_d3_part1.rows() << endl;

    ArrayXXd dlogL_part1 = ArrayXXd::Zero(dlogL_d0_part1.rows(),para.EstTbar*para.dim1);
    dlogL_part1.middleCols(0,para.EstTbar) = dlogL_d0_part1;
    dlogL_part1.middleCols(para.EstTbar,para.EstTbar) = dlogL_d1_part1;
    dlogL_part1.middleCols(2*para.EstTbar,para.EstTbar) = dlogL_d2_part1;
    dlogL_part1.middleCols(3*para.EstTbar,para.EstTbar) = dlogL_d3_part1;
    // cout << "dlogL_part1 done; size dlogL_part1 = " << dlogL_part1.rows() << "; " << dlogL_part1.cols() << endl;
    // cout << "Ncount_prob = " << Ncount_prob << "; " << endl;

    return tuple<ArrayXXd,ArrayXXd,double>(dlogL_part1,logL_part1,Ncount_prob);
}


tuple<ArrayXXd,ArrayXXd,double> alias::EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(const ArrayXd & theta,
    const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
    const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;
    // *********
    /** Part 2 of the full model estimation **/
    std::vector<double> theta_part1(para.dim1);
    for (size_t k = 0; k < para.dim1; ++k) {theta_part1[k] = theta(para.dim3+k);}
    std::vector<double> theta_part2(para.dim2);
    for (size_t n = 0; n < para.dim2; ++n) {theta_part2[n] = theta(para.dim3+para.dim1+n);}

    GoodState = 1; LaborIntensive = 1;
    ParaEst1 para_est1_Good_LaborInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_LaborInt = Backout_phi_a_PD(sim_data_Para_GoodState_LaborInt.Revenue_P,
        sim_data_Para_GoodState_LaborInt.Capital_P, sim_data_Para_GoodState_LaborInt.Employ_ur_P,
        sim_data_Para_GoodState_LaborInt.Employ_uc_P,
        para_est1_Good_LaborInt.alpha_tilde_K, para_est1_Good_LaborInt.alpha_tilde_L,
        para_est1_Good_LaborInt.alpha_Lr, para_est1_Good_LaborInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Good_LaborInt
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_part2,sim_data_Para_GoodState_LaborInt,
        threadsManagement);
    double y_part2_sum_Good_LaborInt = get<0>(t_Part2_Good_LaborInt);
    std::vector<double> grad_a_Part2_sum_Good_LaborInt = get<1>(t_Part2_Good_LaborInt);
    int Ncount_part2_Good_LaborInt = get<2>(t_Part2_Good_LaborInt);
    double Nprob_part2_Good_LaborInt = get<3>(t_Part2_Good_LaborInt);
    ArrayXXd logL_part2_Good_LaborInt = get<4>(t_Part2_Good_LaborInt);
    ArrayXXd dlogL_d0_part2_Good_LaborInt = get<5>(t_Part2_Good_LaborInt);
    ArrayXXd dlogL_d1_part2_Good_LaborInt = get<6>(t_Part2_Good_LaborInt);
    ArrayXXd dlogL_d2_part2_Good_LaborInt = get<7>(t_Part2_Good_LaborInt);
    // cout << "part2_Good_LaborInt size = " << dlogL_d0_part2_Good_LaborInt.rows() << "; " << dlogL_d0_part2_Good_LaborInt.cols()
    //     << "; " << dlogL_d1_part2_Good_LaborInt.rows() << "; " << dlogL_d1_part2_Good_LaborInt.cols()
    //     << "; " << dlogL_d2_part2_Good_LaborInt.rows() << "; " << dlogL_d2_part2_Good_LaborInt.cols() << endl;

    GoodState = 0; LaborIntensive = 1;
    ParaEst1 para_est1_Bad_LaborInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_BadState_LaborInt = Backout_phi_a_PD(sim_data_Para_BadState_LaborInt.Revenue_P,
        sim_data_Para_BadState_LaborInt.Capital_P, sim_data_Para_BadState_LaborInt.Employ_ur_P,
        sim_data_Para_BadState_LaborInt.Employ_uc_P,
        para_est1_Bad_LaborInt.alpha_tilde_K, para_est1_Bad_LaborInt.alpha_tilde_L,
        para_est1_Bad_LaborInt.alpha_Lr, para_est1_Bad_LaborInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Bad_LaborInt
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_part2,sim_data_Para_BadState_LaborInt,
        threadsManagement);
    double y_part2_sum_Bad_LaborInt = get<0>(t_Part2_Bad_LaborInt);
    std::vector<double> grad_a_Part2_sum_Bad_LaborInt = get<1>(t_Part2_Bad_LaborInt);
    int Ncount_part2_Bad_LaborInt = get<2>(t_Part2_Bad_LaborInt);
    double Nprob_part2_Bad_LaborInt = get<3>(t_Part2_Bad_LaborInt);
    ArrayXXd logL_part2_Bad_LaborInt = get<4>(t_Part2_Bad_LaborInt);
    ArrayXXd dlogL_d0_part2_Bad_LaborInt = get<5>(t_Part2_Bad_LaborInt);
    ArrayXXd dlogL_d1_part2_Bad_LaborInt = get<6>(t_Part2_Bad_LaborInt);
    ArrayXXd dlogL_d2_part2_Bad_LaborInt = get<7>(t_Part2_Bad_LaborInt);
    // cout << "part2_Bad_LaborInt size = " << dlogL_d0_part2_Bad_LaborInt.rows() << "; " << dlogL_d0_part2_Bad_LaborInt.cols()
    //     << "; " << dlogL_d1_part2_Bad_LaborInt.rows() << "; " << dlogL_d1_part2_Bad_LaborInt.cols()
    //     << "; " << dlogL_d2_part2_Bad_LaborInt.rows() << "; " << dlogL_d2_part2_Bad_LaborInt.cols() << endl;

    GoodState = 1; LaborIntensive = 0;
    ParaEst1 para_est1_Good_CapitalInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_CapitalInt = Backout_phi_a_PD(sim_data_Para_GoodState_CapitalInt.Revenue_P,
        sim_data_Para_GoodState_CapitalInt.Capital_P, sim_data_Para_GoodState_CapitalInt.Employ_ur_P,
        sim_data_Para_GoodState_CapitalInt.Employ_uc_P,
        para_est1_Good_CapitalInt.alpha_tilde_K, para_est1_Good_CapitalInt.alpha_tilde_L,
        para_est1_Good_CapitalInt.alpha_Lr, para_est1_Good_CapitalInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Good_CapitalInt
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_part2,sim_data_Para_GoodState_CapitalInt,
        threadsManagement);
    double y_part2_sum_Good_CapitalInt = get<0>(t_Part2_Good_CapitalInt);
    std::vector<double> grad_a_Part2_sum_Good_CapitalInt = get<1>(t_Part2_Good_CapitalInt);
    int Ncount_part2_Good_CapitalInt = get<2>(t_Part2_Good_CapitalInt);
    double Nprob_part2_Good_CapitalInt = get<3>(t_Part2_Good_CapitalInt);
    ArrayXXd logL_part2_Good_CapitalInt = get<4>(t_Part2_Good_CapitalInt);
    ArrayXXd dlogL_d0_part2_Good_CapitalInt = get<5>(t_Part2_Good_CapitalInt);
    ArrayXXd dlogL_d1_part2_Good_CapitalInt = get<6>(t_Part2_Good_CapitalInt);
    ArrayXXd dlogL_d2_part2_Good_CapitalInt = get<7>(t_Part2_Good_CapitalInt);
    // cout << "part2_Good_CapitalInt size = " << dlogL_d0_part2_Good_CapitalInt.rows() << "; " << dlogL_d0_part2_Good_CapitalInt.cols()
    //     << "; " << dlogL_d1_part2_Good_CapitalInt.rows() << "; " << dlogL_d1_part2_Good_CapitalInt.cols()
    //     << "; " << dlogL_d2_part2_Good_CapitalInt.rows() << "; " << dlogL_d2_part2_Good_CapitalInt.cols() << endl;

    GoodState = 0; LaborIntensive = 0;
    ParaEst1 para_est1_Bad_CapitalInt = constructParaEst_Part1(theta_part1,GoodState,LaborIntensive);
    //// estimated productivity
    ArrayXXd lnphi_a_BadState_CapitalInt = Backout_phi_a_PD(sim_data_Para_BadState_CapitalInt.Revenue_P,
        sim_data_Para_BadState_CapitalInt.Capital_P, sim_data_Para_BadState_CapitalInt.Employ_ur_P,
        sim_data_Para_BadState_CapitalInt.Employ_uc_P,
        para_est1_Bad_CapitalInt.alpha_tilde_K, para_est1_Bad_CapitalInt.alpha_tilde_L,
        para_est1_Bad_CapitalInt.alpha_Lr, para_est1_Bad_CapitalInt.alpha_Lc);
    tuple<double,std::vector<double>,int,double,ArrayXXd,ArrayXXd,ArrayXXd,ArrayXXd> t_Part2_Bad_CapitalInt
        = EstimationAuxiliaryModel_part2_ObjFun_FullVersion_Step2(theta_part2,sim_data_Para_BadState_CapitalInt,
        threadsManagement);
    double y_part2_sum_Bad_CapitalInt = get<0>(t_Part2_Bad_CapitalInt);
    std::vector<double> grad_a_Part2_sum_Bad_CapitalInt = get<1>(t_Part2_Bad_CapitalInt);
    int Ncount_part2_Bad_CapitalInt = get<2>(t_Part2_Bad_CapitalInt);
    double Nprob_part2_Bad_CapitalInt = get<3>(t_Part2_Bad_CapitalInt);
    ArrayXXd logL_part2_Bad_CapitalInt = get<4>(t_Part2_Bad_CapitalInt);
    ArrayXXd dlogL_d0_part2_Bad_CapitalInt = get<5>(t_Part2_Bad_CapitalInt);
    ArrayXXd dlogL_d1_part2_Bad_CapitalInt = get<6>(t_Part2_Bad_CapitalInt);
    ArrayXXd dlogL_d2_part2_Bad_CapitalInt = get<7>(t_Part2_Bad_CapitalInt);
    // cout << "part2_Bad_CapitalInt size = " << dlogL_d0_part2_Bad_CapitalInt.rows() << "; " << dlogL_d0_part2_Bad_CapitalInt.cols()
    //     << "; " << dlogL_d1_part2_Bad_CapitalInt.rows() << "; " << dlogL_d1_part2_Bad_CapitalInt.cols()
    //     << "; " << dlogL_d2_part2_Bad_CapitalInt.rows() << "; " << dlogL_d2_part2_Bad_CapitalInt.cols() << endl;

    ArrayXXd logL_part2 = AppendMatrix(logL_part2_Good_LaborInt,logL_part2_Bad_LaborInt,
        logL_part2_Good_CapitalInt,logL_part2_Bad_CapitalInt);
    ArrayXXd dlogL_d0_part2 = AppendMatrix(dlogL_d0_part2_Good_LaborInt,dlogL_d0_part2_Bad_LaborInt,
    dlogL_d0_part2_Good_CapitalInt,dlogL_d0_part2_Bad_CapitalInt);
    ArrayXXd dlogL_d1_part2 = AppendMatrix(dlogL_d1_part2_Good_LaborInt,dlogL_d1_part2_Bad_LaborInt,
    dlogL_d1_part2_Good_CapitalInt,dlogL_d1_part2_Bad_CapitalInt);
    ArrayXXd dlogL_d2_part2 = AppendMatrix(dlogL_d2_part2_Good_LaborInt,dlogL_d2_part2_Bad_LaborInt,
    dlogL_d2_part2_Good_CapitalInt,dlogL_d2_part2_Bad_CapitalInt);
    // cout << "para.EstTbar = " << para.EstTbar << "; dlogL_d0_part2; dlogL_d1_part2; dlogL_d2_part2 done" << endl;
    // cout << "dlogL_d0_part2.rows() = " << dlogL_d0_part2.rows() << "; dlogL_d1_part2.rows() = " << dlogL_d1_part2.rows()
    //     << "; dlogL_d2_part2.rows() = " << dlogL_d2_part2.rows() << endl;

    ArrayXXd dlogL_part2 = ArrayXXd::Zero(dlogL_d0_part2.rows(),para.EstTbar*para.dim2);
    dlogL_part2 << dlogL_d0_part2,dlogL_d1_part2,dlogL_d2_part2;
    // cout << "dlogL_part2 done; size dlogL_part2 = " << dlogL_part2.rows() << "; " << dlogL_part2.cols() << endl;

    double Nprob_part2 = Nprob_part2_Good_LaborInt + Nprob_part2_Bad_LaborInt + Ncount_part2_Good_CapitalInt
        + Ncount_part2_Bad_CapitalInt;
    // cout << "dlogL_part2 done; size dlogL_part2 = " << dlogL_part2.rows() << "; " << dlogL_part2.cols() << endl;
    // cout << "Nprob_part2 = " << Nprob_part2 << endl;

    return tuple<ArrayXXd,ArrayXXd,double>(dlogL_part2,logL_part2,Nprob_part2);
}

tuple<ArrayXXd,ArrayXXd,double> alias::EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
    const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
    const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
    const EquV_Auxiliary & EquV_GoodState_LaborInt, const EquV_Auxiliary & EquV_BadState_LaborInt,
    const EquV_Auxiliary & EquV_GoodState_CapitalInt, const EquV_Auxiliary & EquV_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;

    // cout << "Start EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2;" << endl;
    // *********
    /** Part 3 of the full model estimation **/
    GoodState = 1; LaborIntensive = 1;
    tuple<ArrayXXd,ArrayXXd,double> t_part3_GoodState_LaborInt = Cal_dlogL_Part3_SimData(sim_data_Para_GoodState_LaborInt,
        EquV_GoodState_LaborInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Good_LaborInt = get<0>(t_part3_GoodState_LaborInt);
    ArrayXXd logL_part3_Good_LaborInt = get<1>(t_part3_GoodState_LaborInt);
    double Nprob_part3_Good_LaborInt = get<2>(t_part3_GoodState_LaborInt);

    GoodState = 0; LaborIntensive = 1;
    tuple<ArrayXXd,ArrayXXd,double> t_part3_BadState_LaborInt = Cal_dlogL_Part3_SimData(sim_data_Para_BadState_LaborInt,
        EquV_BadState_LaborInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Bad_LaborInt = get<0>(t_part3_BadState_LaborInt);
    ArrayXXd logL_part3_Bad_LaborInt = get<1>(t_part3_BadState_LaborInt);
    double Nprob_part3_Bad_LaborInt = get<2>(t_part3_BadState_LaborInt);

    GoodState = 1; LaborIntensive = 0;
    tuple<ArrayXXd,ArrayXXd,double> t_part3_GoodState_CapitalInt = Cal_dlogL_Part3_SimData(sim_data_Para_GoodState_CapitalInt,
        EquV_GoodState_CapitalInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Good_CapitalInt = get<0>(t_part3_GoodState_CapitalInt);
    ArrayXXd logL_part3_Good_CapitalInt = get<1>(t_part3_GoodState_CapitalInt);
    double Nprob_part3_Good_CapitalInt = get<2>(t_part3_GoodState_CapitalInt);

    GoodState = 0; LaborIntensive = 0;
    tuple<ArrayXXd,ArrayXXd,double> t_part3_BadState_CapitalInt = Cal_dlogL_Part3_SimData(sim_data_Para_BadState_CapitalInt,
        EquV_BadState_CapitalInt, GoodState, LaborIntensive, threadsManagement);
    ArrayXXd dlogL_part3_Bad_CapitalInt = get<0>(t_part3_BadState_CapitalInt);
    ArrayXXd logL_part3_Bad_CapitalInt = get<1>(t_part3_BadState_CapitalInt);
    double Nprob_part3_Bad_CapitalInt = get<2>(t_part3_BadState_CapitalInt);

    // cout << "grad_a3_sum_Bad_CapitalInt = " << grad_a3_sum_Bad_CapitalInt.transpose() << endl;
    // cout << "part3_Bad_CapitalInt size = " << dlogL_part3_Bad_CapitalInt.rows() << "; " << dlogL_part3_Bad_CapitalInt.cols() << endl;
    // cout << "dlogL_part3_Bad_CapitalInt = " << dlogL_part3_Bad_CapitalInt.colwise().sum() << endl;
    // throw runtime_error("895");
    ArrayXXd logL_part3 = AppendMatrix(logL_part3_Good_LaborInt,logL_part3_Bad_LaborInt,
        logL_part3_Good_CapitalInt,logL_part3_Bad_CapitalInt);
    ArrayXXd dlogL_part3 = AppendMatrix(dlogL_part3_Good_LaborInt,dlogL_part3_Bad_LaborInt,
        dlogL_part3_Good_CapitalInt,dlogL_part3_Bad_CapitalInt);

    double Nprob_part3 = Nprob_part3_Good_LaborInt + Nprob_part3_Bad_LaborInt + Nprob_part3_Good_CapitalInt + Nprob_part3_Bad_CapitalInt;
    // cout << "Nprob_part3 = " << Nprob_part3 << endl;
    // cout << "logL_part3.rows() = " << logL_part3.rows() << endl;
    return tuple<ArrayXXd,ArrayXXd,double>(dlogL_part3,logL_part3,Nprob_part3);
}
//
tuple<ArrayXXd,ArrayXXd,double> alias::Cal_dlogL_Part3_SimData(const SimData & SimData, const EquV_Auxiliary & EquV,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement) {

    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3 =
        EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(SimData,
        EquV.para_est_a_point, EquV.para_vec_a_point, EquV.EquV_a_point, EquV.Evalmat_a_point,
        GoodState, LaborIntensive, threadsManagement);
    double obj_point_sum = get<0>(t_Part3);
    int Ncount_part3 = get<1>(t_Part3);
    double Nprob_part3 = get<2>(t_Part3);
    ArrayXXd logL_part3 = get<3>(t_Part3);
    ArrayXXi K_index_max_P_mat = get<4>(t_Part3);
    ArrayXXi Lur_index_max_P_mat = get<5>(t_Part3);
    ArrayXXi Luc_index_max_P_mat = get<6>(t_Part3);
    ArrayXXd K_max_P_mat = get<7>(t_Part3);
    ArrayXXd Lur_max_P_mat = get<8>(t_Part3);
    ArrayXXd Luc_max_P_mat = get<9>(t_Part3);
    ArrayXXi K_index_max_D_mat = get<10>(t_Part3);
    ArrayXXi Lur_index_max_D_mat = get<11>(t_Part3);
    ArrayXXi Luc_index_max_D_mat = get<12>(t_Part3);
    ArrayXXd K_max_D_mat = get<13>(t_Part3);
    ArrayXXd Lur_max_D_mat = get<14>(t_Part3);
    ArrayXXd Luc_max_D_mat = get<15>(t_Part3);

    ArrayXXd dlogL_part3 = ArrayXXd::Zero(SimData.good.rows(),para.EstTbar * para.dim3);
    ArrayXd grad_a3_sum = ArrayXd::Zero(para.dim3);
    for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
    //        cout << "i_a = " << i_a << endl;
        double obj_plus_sum = 0.0;
        if ( ( (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15) and GoodState == 1 ) or
            ( (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14) and GoodState == 0 ) ) {
            obj_plus_sum = 0.0;
            grad_a3_sum(i_a) = 0.0;
        }
        else {
            tuple<double,ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_Simulation_FullVersion_Diff_Step2(SimData,
                EquV.para_est_a_plus[i_a], EquV.para_vec_a_plus[i_a],
                EquV.EquV_a_plus[i_a], EquV.Evalmat_a_plus[i_a],
                K_index_max_P_mat, Lur_index_max_P_mat, Luc_index_max_P_mat,
                K_max_P_mat, Lur_max_P_mat, Luc_max_P_mat,
                K_index_max_D_mat, Lur_index_max_D_mat, Luc_index_max_D_mat,
                K_max_D_mat, Lur_max_D_mat, Luc_max_D_mat,
                GoodState, LaborIntensive, threadsManagement);
            obj_plus_sum = get<0>(t_plus);
            grad_a3_sum(i_a) = (obj_plus_sum - obj_point_sum)
                / EquV.eps_diff(i_a);
            ArrayXXd logL_part3_plus = get<1>(t_plus);
            dlogL_part3.middleCols(i_a*para.EstTbar,para.EstTbar)
                = ( -logL_part3_plus - (-logL_part3) ) / EquV.eps_diff(i_a);
        }
    }

    return tuple<ArrayXXd,ArrayXXd,double>(dlogL_part3,logL_part3,Nprob_part3);
}

ArrayXXd alias::Cal_db_dtheta_a(const ArrayXd & theta_Est,const ArrayXd & theta_a,
    const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
    const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;

    GoodState = 1; LaborIntensive = 1;
    SimData sim_data_Para_GoodState_LaborInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 1;
    SimData sim_data_Para_BadState_LaborInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);
    GoodState = 1; LaborIntensive = 0;
    SimData sim_data_Para_GoodState_CapitalInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    GoodState = 0; LaborIntensive = 0;
    SimData sim_data_Para_BadState_CapitalInt = SimulationData_Step2Estimation_FullVersion(theta_Est,
        sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    cout << "Simulate the full model" << endl;

    /*** Part 1 ***/
    ArrayXXd db_dtheta_a_part1_temp = ArrayXXd::Zero(para.dim1, para.dim1);
    ArrayXd NumProb_1_Good = ArrayXd::Zero(para.dim1);
    ArrayXd NumProb_1_Bad = ArrayXd::Zero(para.dim1);
    for (size_t k = 0; k < para.dim1; k++) {
        cout << "k = " << k << endl;
        double eps = abs(theta_a(para.dim3 + k)) * 1e-8;

        ArrayXd theta_a_plus = theta_a;
        theta_a_plus(para.dim3 + k) = theta_a(para.dim3 + k) + eps;
        ArrayXd theta_a_minus = theta_a;
        theta_a_minus(para.dim3 + k) = theta_a(para.dim3 + k) - eps;
        cout << theta_a_plus.transpose() << endl;
        cout << theta_a_minus.transpose() << endl;

        tuple<ArrayXXd,ArrayXXd,double> t_part1_plus = EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(theta_a_plus,
            sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part1_plus = get<0>(t_part1_plus);
        double Nprob_part1_plus = get<2>(t_part1_plus);
        tuple<ArrayXXd,ArrayXXd,double> t_part1_minus = EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(theta_a_minus,
            sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part1_minus = get<0>(t_part1_minus);
        double Nprob_part1_minus = get<2>(t_part1_minus);
        ArrayXXd ddL_part1 = ( -dL_part1_plus - (-dL_part1_minus) ) / 2.0 / eps;
        // cout << "dL_part1_plus = " << dL_part1_plus.rows() << "; " << dL_part1_plus.cols() << endl;
        // cout << "dL_part1_minus = " << dL_part1_minus.rows() << "; " << dL_part1_minus.cols() << endl;

        ArrayXXd dL_part1_plus_vec = Vectorize_dl_dtheta_part(dL_part1_plus, para.dim1);
        ArrayXXd dL_part1_minus_vec = Vectorize_dl_dtheta_part(dL_part1_minus, para.dim1);
        // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;
        // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;

        ArrayXi temp_index = ( (dL_part1_plus_vec - dL_part1_minus_vec).abs() > 0 ).cast<int>().colwise().sum() + 1;
        cout << "temp_index = " << temp_index.transpose() << endl;
        cout << "Nprob_part1_plus = " << Nprob_part1_plus << endl;
        cout << "dL_part1_plus = " << dL_part1_plus.rows() << endl;

        db_dtheta_a_part1_temp.row(k) = ( -dL_part1_plus_vec - (-dL_part1_minus_vec) ).colwise().sum() / 2.0 / eps / temp_index.cast<double>().transpose();
        // cout << "double(ddL_part3.rows()) = " << double(ddL_part3.rows()) << endl;
        cout << "db_dtheta_a_part1_temp.row(k) = " << db_dtheta_a_part1_temp.row(k) << endl;
    }
    cout << "db_dtheta_a_part1_temp = " << db_dtheta_a_part1_temp << endl;
    writeToCSVfile("db_dtheta_a_part1_temp.csv", db_dtheta_a_part1_temp.cast<double>().matrix());
    ArrayXXd db_dtheta_a_part1_tran = db_dtheta_a_part1_temp.transpose();
    ArrayXXd db_dtheta_a_part1 = db_dtheta_a_part1_temp*0.5+db_dtheta_a_part1_tran*0.5;
    // throw runtime_error("1451");

    /*** Part 2 ***/
    ArrayXXd db_dtheta_a_part2_temp = ArrayXXd::Zero(para.dim2, para.dim2);
    for (size_t k = 0; k < para.dim2; k++) {
        cout << "k = " << k << endl;
        double eps = abs(theta_a(para.dim3 + para.dim1 + k)) * 1e-8;

        ArrayXd theta_a_plus = theta_a;
        theta_a_plus(para.dim3 + para.dim1 + k) = theta_a(para.dim3 + para.dim1 + k) + eps;
        ArrayXd theta_a_minus = theta_a;
        theta_a_minus(para.dim3 + para.dim1 + k) = theta_a(para.dim3 + para.dim1 + k) - eps;
        cout << theta_a_plus.transpose() << endl;
        cout << theta_a_minus.transpose() << endl;

        tuple<ArrayXXd,ArrayXXd,double> t_part2_plus = EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(theta_a_plus,
            sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part2_plus = get<0>(t_part2_plus);
        double Nprob_part2_plus = get<2>(t_part2_plus);
        tuple<ArrayXXd,ArrayXXd,double> t_part2_minus = EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(theta_a_minus,
            sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
            threadsManagement);
        ArrayXXd dL_part2_minus = get<0>(t_part2_minus);
        double Nprob_part2_minus = get<2>(t_part2_minus);
        ArrayXXd ddL_part2 = ( -dL_part2_plus - (-dL_part2_minus) ) / 2.0 / eps;
        // cout << "dL_part1_plus = " << dL_part1_plus.rows() << "; " << dL_part1_plus.cols() << endl;
        // cout << "dL_part1_minus = " << dL_part1_minus.rows() << "; " << dL_part1_minus.cols() << endl;

        ArrayXXd dL_part2_plus_vec = Vectorize_dl_dtheta_part(dL_part2_plus, para.dim2);
        ArrayXXd dL_part2_minus_vec = Vectorize_dl_dtheta_part(dL_part2_minus, para.dim2);
        // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;
        // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;

        ArrayXi temp_index = ( (dL_part2_plus_vec - dL_part2_minus_vec).abs() > 0 ).cast<int>().colwise().sum() + 1;
        cout << "temp_index = " << temp_index.transpose() << endl;

        db_dtheta_a_part2_temp.row(k) = ( -dL_part2_plus_vec - (-dL_part2_minus_vec) ).colwise().sum() / 2.0 / eps;
        cout << "db_dtheta_a_part2_temp.row(k) = " << db_dtheta_a_part2_temp.row(k) << endl;

        db_dtheta_a_part2_temp.row(k) = ( -dL_part2_plus_vec - (-dL_part2_minus_vec) ).colwise().sum() / 2.0 / eps / temp_index.cast<double>().transpose();
        cout << "db_dtheta_a_part2_temp.row(k) = " << db_dtheta_a_part2_temp.row(k) << endl;

        cout << "dL_part2_plus = " << dL_part2_plus.rows() << endl;
        cout << "Nprob_part2_plus = " << Nprob_part2_plus << endl;
    }
    cout << "db_dtheta_a_part2_temp = " << db_dtheta_a_part2_temp << endl;
    writeToCSVfile("db_dtheta_a_part2_temp.csv", db_dtheta_a_part2_temp.cast<double>().matrix());
    ArrayXXd db_dtheta_a_part2_tran = db_dtheta_a_part2_temp.transpose();
    ArrayXXd db_dtheta_a_part2 = db_dtheta_a_part2_temp*0.5+db_dtheta_a_part2_tran*0.5;

    /*** Part 3 ***/
    // ArrayXXd db_dtheta_a_part3_temp = ArrayXXd::Zero(para.dim3, para.dim3);
    // for (size_t k = 0; k < para.dim3; k++) {
    //     cout << "k = " << k << endl;
    //     double eps = abs(theta_a(k)) * 1e-6;
    //
    //     ArrayXd theta_a_plus = theta_a;
    //     theta_a_plus(k) = theta_a(k) + eps;
    //     ArrayXd theta_a_minus = theta_a;
    //     theta_a_minus(k) = theta_a(k) - eps;
    //     cout << theta_a_plus.transpose() << endl;
    //     cout << theta_a_minus.transpose() << endl;
    //
    //     GoodState = 1; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_GoodState_LaborInt_plus = CalEquV_plus_minus(theta_a_plus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_BadState_LaborInt_plus = CalEquV_plus_minus(theta_a_plus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 1; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_GoodState_CapitalInt_plus = CalEquV_plus_minus(theta_a_plus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_BadState_CapitalInt_plus = CalEquV_plus_minus(theta_a_plus, GoodState, LaborIntensive, threadsManagement);
    //     tuple<ArrayXXd,ArrayXXd,double> t_part3_plus = EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
    //         sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,
    //         sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
    //         EquV_GoodState_LaborInt_plus,EquV_BadState_LaborInt_plus,
    //         EquV_GoodState_CapitalInt_plus,EquV_BadState_CapitalInt_plus,
    //         threadsManagement);
    //     ArrayXXd dL_part3_plus = get<0>(t_part3_plus);
    //     double Nprob_part3_plus = get<2>(t_part3_plus);
    //
    //     GoodState = 1; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_GoodState_LaborInt_minus = CalEquV_plus_minus(theta_a_minus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 1;
    //     EquV_Auxiliary EquV_BadState_LaborInt_minus = CalEquV_plus_minus(theta_a_minus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 1; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_GoodState_CapitalInt_minus = CalEquV_plus_minus(theta_a_minus, GoodState, LaborIntensive, threadsManagement);
    //     GoodState = 0; LaborIntensive = 0;
    //     EquV_Auxiliary EquV_BadState_CapitalInt_minus = CalEquV_plus_minus(theta_a_minus, GoodState, LaborIntensive, threadsManagement);
    //     tuple<ArrayXXd,ArrayXXd,double> t_part3_minus = EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
    //         sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,
    //         sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
    //         EquV_GoodState_LaborInt_minus,EquV_BadState_LaborInt_minus,
    //         EquV_GoodState_CapitalInt_minus,EquV_BadState_CapitalInt_minus,
    //         threadsManagement);
    //     ArrayXXd dL_part3_minus = get<0>(t_part3_minus);
    //     double Nprob_part3_minus = get<2>(t_part3_minus);
    //
    //     ArrayXXd ddL_part3 = ( -dL_part3_plus - (-dL_part3_minus) ) / 2.0 / eps;
    //     // cout << "dL_part1_plus = " << dL_part1_plus.rows() << "; " << dL_part1_plus.cols() << endl;
    //     // cout << "dL_part1_minus = " << dL_part1_minus.rows() << "; " << dL_part1_minus.cols() << endl;
    //
    //     ArrayXXd dL_part3_plus_vec = Vectorize_dl_dtheta_part(dL_part3_plus, para.dim3);
    //     ArrayXXd dL_part3_minus_vec = Vectorize_dl_dtheta_part(dL_part3_minus, para.dim3);
    //     // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;
    //     // cout << "dL_part1_plus_vec = " << dL_part1_plus_vec.rows() << "; " << dL_part1_plus_vec.cols() << endl;
    //
    //     ArrayXi temp_index = ( (dL_part3_plus_vec - dL_part3_minus_vec).abs() > 0 ).cast<int>().colwise().sum() + 1;
    //     cout << "temp_index = " << temp_index.transpose() << endl;
    //
    //     db_dtheta_a_part3_temp.row(k) = ( -dL_part3_plus_vec - (-dL_part3_minus_vec) ).colwise().sum() / 2.0 / eps;
    //     // cout << "double(ddL_part3.rows()) = " << double(ddL_part3.rows()) << endl;
    //     cout << "db_dtheta_a_part3_temp.row(k) = " << db_dtheta_a_part3_temp.row(k) << endl;
    //
    //     db_dtheta_a_part3_temp.row(k) = ( -dL_part3_plus_vec - (-dL_part3_minus_vec) ).colwise().sum() / 2.0 / eps
    //         / temp_index.cast<double>().transpose();
    //     cout << "db_dtheta_a_part3_temp.row(k) = " << db_dtheta_a_part3_temp.row(k) << endl;
    //
    //     cout << "dL_part3_plus.rows() = " << dL_part3_plus.rows() << endl;
    //     cout << "Nprob_part3_plus = " << Nprob_part3_plus << endl;
    //}

    ArrayXXd db_dtheta_a_part3_temp = ArrayXXd::Zero(para.dim3, para.dim3);
    for (size_t k = 0; k < para.dim3; k++) {
        cout << "k = " << k << endl;

        double eps_theta1 = abs(theta_a(k)) * 1e-8;

        ArrayXd theta_a_plus = theta_a;
        theta_a_plus(k) = theta_a(k) + eps_theta1;
        ArrayXd theta_a_minus = theta_a;
        theta_a_minus(k) = theta_a(k) - eps_theta1;
        cout << theta_a_plus.transpose() << endl;
        cout << theta_a_minus.transpose() << endl;

        for (size_t kk = 0; kk < para.dim3; kk++) {
            cout << "kk = " << kk << endl;
            double eps_theta2 = abs(theta_a(kk)) * 1e-8;

            ArrayXd theta_a_plus_plus = theta_a_plus;
            theta_a_plus_plus(kk) = theta_a_plus(kk) + eps_theta2;
            ArrayXd theta_a_plus_minus = theta_a_plus;
            theta_a_plus_minus(kk) = theta_a_plus(kk) - eps_theta2;

            ArrayXd theta_a_minus_plus = theta_a_minus;
            theta_a_minus_plus(kk) = theta_a_minus(kk) + eps_theta2;
            ArrayXd theta_a_minus_minus = theta_a_minus;
            theta_a_minus_minus(kk) = theta_a_minus(kk) - eps_theta2;

            tuple<ArrayXXd,double> t_logL_part3_plus_plus = EstimationIndirectInference_SD_l_Part3_SimData(theta_a_plus_plus,
                sim_data_Para_GoodState_LaborInt, sim_data_Para_BadState_LaborInt,
                sim_data_Para_GoodState_CapitalInt, sim_data_Para_BadState_CapitalInt,
                threadsManagement);
            ArrayXXd logL_part3_plus_plus = get<0>(t_logL_part3_plus_plus);

            tuple<ArrayXXd,double> t_logL_part3_plus_minus = EstimationIndirectInference_SD_l_Part3_SimData(theta_a_plus_minus,
                sim_data_Para_GoodState_LaborInt, sim_data_Para_BadState_LaborInt,
                sim_data_Para_GoodState_CapitalInt, sim_data_Para_BadState_CapitalInt,
                threadsManagement);
            ArrayXXd logL_part3_plus_minus = get<0>(t_logL_part3_plus_minus);

            tuple<ArrayXXd,double> t_logL_part3_minus_plus = EstimationIndirectInference_SD_l_Part3_SimData(theta_a_minus_plus,
                sim_data_Para_GoodState_LaborInt, sim_data_Para_BadState_LaborInt,
                sim_data_Para_GoodState_CapitalInt, sim_data_Para_BadState_CapitalInt,
                threadsManagement);
            ArrayXXd logL_part3_minus_plus = get<0>(t_logL_part3_minus_plus);

            tuple<ArrayXXd,double> t_logL_part3_minus_minus = EstimationIndirectInference_SD_l_Part3_SimData(theta_a_minus_minus,
                sim_data_Para_GoodState_LaborInt, sim_data_Para_BadState_LaborInt,
                sim_data_Para_GoodState_CapitalInt, sim_data_Para_BadState_CapitalInt,
                threadsManagement);
            ArrayXXd logL_part3_minus_minus = get<0>(t_logL_part3_minus_minus);

            ArrayXXd ddL_part3 = - ( (logL_part3_plus_plus - logL_part3_plus_minus) / 2.0 / eps_theta2
                - (logL_part3_minus_plus - logL_part3_minus_minus) / 2.0 / eps_theta2 ) / 2.0 / eps_theta1;

            // int NumObs = (ddL_part3.abs() > 0).cast<int>().sum();
            ArrayXd ddL_part3_sum = ddL_part3.abs().rowwise().sum();
            int NumObs = (ddL_part3_sum.abs() > 0).cast<int>().sum();

            db_dtheta_a_part3_temp(k,kk) = ddL_part3.sum() / ( double(NumObs)+0.001 );
            cout << "k = " << k << "; kk = " << kk << "; db_dtheta_a_part3_temp(k,kk) = " << db_dtheta_a_part3_temp(k,kk)
                << "; NumObs = " << NumObs << endl;
        }
    }

    cout << "db_dtheta_a_part3_temp = " << db_dtheta_a_part3_temp << endl;
    writeToCSVfile("db_dtheta_a_part3_temp.csv", db_dtheta_a_part3_temp.cast<double>().matrix());
    ArrayXXd db_dtheta_a_part3_tran = db_dtheta_a_part3_temp.transpose();
    ArrayXXd db_dtheta_a_part3 = db_dtheta_a_part3_temp*0.5+db_dtheta_a_part3_tran*0.5;

    ArrayXXd db_dtheta_a = ArrayXXd::Zero(para.dim1+para.dim2+para.dim3,para.dim1+para.dim2+para.dim3);
    cout << "db_dtheta_a size = " << db_dtheta_a.cols() << "; " << db_dtheta_a.rows() << endl;
    db_dtheta_a.block(para.dim3,para.dim3,para.dim1,para.dim1) = db_dtheta_a_part1;
    db_dtheta_a.block(para.dim3+para.dim1,para.dim3+para.dim1,para.dim2,para.dim2) = db_dtheta_a_part2;
    db_dtheta_a.block(0,0,para.dim3,para.dim3) = db_dtheta_a_part3;
    cout << "db_dtheta_a done" << endl;
    cout << "db_dtheta_a = " << db_dtheta_a << endl;

    return db_dtheta_a;
}


ArrayXXd alias::Cal_db_dtheta0_v2(const ArrayXd & theta_Est,const ArrayXd & theta_a,
    const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
    const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    int GoodState; int LaborIntensive;

    /*** part 1 ***/
    ArrayXXd db_dtheta0_temp = ArrayXXd::Zero(para.dim1+para.dim2+para.dim3,para.dim1+para.dim2+para.dim3);
    for (size_t k = 0; k < para.dim1+para.dim2+para.dim3; ++k) {
        cout << "k = " << k << endl;
        double eps_theta0 = abs(theta_Est(k)) * 1e-8;

        ArrayXd theta_Est_plus = theta_Est;
        theta_Est_plus(k) = theta_Est(k) + eps_theta0;
        cout << "theta_Est_plus = " << theta_Est_plus.transpose() << endl;

        GoodState = 1; LaborIntensive = 1;
        SimData sim_data_Para_GoodState_LaborInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
            sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
            GoodState, LaborIntensive, threadsManagement);
        GoodState = 0; LaborIntensive = 1;
        SimData sim_data_Para_BadState_LaborInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
            sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
            GoodState, LaborIntensive, threadsManagement);
        GoodState = 1; LaborIntensive = 0;
        SimData sim_data_Para_GoodState_CapitalInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
            sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
            GoodState, LaborIntensive, threadsManagement);
        GoodState = 0; LaborIntensive = 0;
        SimData sim_data_Para_BadState_CapitalInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
            sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
            GoodState, LaborIntensive, threadsManagement);

        ArrayXd theta_Est_minus = theta_Est;
        theta_Est_minus(k) = theta_Est(k) - eps_theta0;
        cout << "theta_Est_minus = " << theta_Est_minus.transpose() << endl;
        cout << "diff theta_Est = " << theta_Est_plus.transpose() - theta_Est_minus.transpose() << endl;

        GoodState = 1; LaborIntensive = 1;
        SimData sim_data_Para_GoodState_LaborInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
            sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
            GoodState, LaborIntensive, threadsManagement);
        GoodState = 0; LaborIntensive = 1;
        SimData sim_data_Para_BadState_LaborInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
            sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
            GoodState, LaborIntensive, threadsManagement);
        GoodState = 1; LaborIntensive = 0;
        SimData sim_data_Para_GoodState_CapitalInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
            sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
            GoodState, LaborIntensive, threadsManagement);
        GoodState = 0; LaborIntensive = 0;
        SimData sim_data_Para_BadState_CapitalInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
            sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
            GoodState, LaborIntensive, threadsManagement);

        for (size_t kk = 0; kk < para.dim1+para.dim2+para.dim3; ++kk) { //
            cout << "kk = " << kk << endl;
            double eps_theta_a = abs(theta_a(kk)) * 1e-8;

            ArrayXd theta_a_plus = theta_a;
            theta_a_plus(kk) = theta_a(kk) + eps_theta_a;
            ArrayXd theta_a_minus = theta_a;
            theta_a_minus(kk) = theta_a(kk) - eps_theta_a;

            tuple<ArrayXXd,ArrayXXd,ArrayXXd> t_logL_dtheta_Est_plus_theta_a_plus = Cal_logL_dtheta_Est_theta_a(
                theta_a_plus,sim_data_Para_GoodState_LaborInt_plus,
                sim_data_Para_BadState_LaborInt_plus,
                sim_data_Para_GoodState_CapitalInt_plus,
                sim_data_Para_BadState_CapitalInt_plus,
                threadsManagement);
            ArrayXXd logL_part1_theta_Est_plus_theta_a_plus = get<0>(t_logL_dtheta_Est_plus_theta_a_plus);
            ArrayXXd logL_part2_theta_Est_plus_theta_a_plus = get<1>(t_logL_dtheta_Est_plus_theta_a_plus);
            ArrayXXd logL_part3_theta_Est_plus_theta_a_plus = get<2>(t_logL_dtheta_Est_plus_theta_a_plus);

            tuple<ArrayXXd,ArrayXXd,ArrayXXd> t_logL_dtheta_Est_minus_theta_a_plus = Cal_logL_dtheta_Est_theta_a(
                theta_a_plus,sim_data_Para_GoodState_LaborInt_minus,
                sim_data_Para_BadState_LaborInt_minus,
                sim_data_Para_GoodState_CapitalInt_minus,
                sim_data_Para_BadState_CapitalInt_minus,
                threadsManagement);
            ArrayXXd logL_part1_theta_Est_minus_theta_a_plus = get<0>(t_logL_dtheta_Est_minus_theta_a_plus);
            ArrayXXd logL_part2_theta_Est_minus_theta_a_plus = get<1>(t_logL_dtheta_Est_minus_theta_a_plus);
            ArrayXXd logL_part3_theta_Est_minus_theta_a_plus = get<2>(t_logL_dtheta_Est_minus_theta_a_plus);

            tuple<ArrayXXd,ArrayXXd,ArrayXXd> t_logL_dtheta_Est_plus_theta_a_minus = Cal_logL_dtheta_Est_theta_a(
                theta_a_minus,sim_data_Para_GoodState_LaborInt_plus,
                sim_data_Para_BadState_LaborInt_plus,
                sim_data_Para_GoodState_CapitalInt_plus,
                sim_data_Para_BadState_CapitalInt_plus,
                threadsManagement);
            ArrayXXd logL_part1_theta_Est_plus_theta_a_minus = get<0>(t_logL_dtheta_Est_plus_theta_a_minus);
            ArrayXXd logL_part2_theta_Est_plus_theta_a_minus = get<1>(t_logL_dtheta_Est_plus_theta_a_minus);
            ArrayXXd logL_part3_theta_Est_plus_theta_a_minus = get<2>(t_logL_dtheta_Est_plus_theta_a_minus);

            tuple<ArrayXXd,ArrayXXd,ArrayXXd> t_logL_dtheta_Est_minus_theta_a_minus = Cal_logL_dtheta_Est_theta_a(
                theta_a_minus,sim_data_Para_GoodState_LaborInt_minus,
                sim_data_Para_BadState_LaborInt_minus,
                sim_data_Para_GoodState_CapitalInt_minus,
                sim_data_Para_BadState_CapitalInt_minus,
                threadsManagement);
            ArrayXXd logL_part1_theta_Est_minus_theta_a_minus = get<0>(t_logL_dtheta_Est_minus_theta_a_minus);
            ArrayXXd logL_part2_theta_Est_minus_theta_a_minus = get<1>(t_logL_dtheta_Est_minus_theta_a_minus);
            ArrayXXd logL_part3_theta_Est_minus_theta_a_minus = get<2>(t_logL_dtheta_Est_minus_theta_a_minus);

            ArrayXXd dlogL_theta_Est_theta_a;
            if (kk < para.dim3) {
                dlogL_theta_Est_theta_a = (logL_part3_theta_Est_plus_theta_a_plus - logL_part3_theta_Est_minus_theta_a_plus)
                    - (logL_part3_theta_Est_plus_theta_a_minus - logL_part3_theta_Est_minus_theta_a_minus);
            }
            else if (kk >= para.dim3 and kk < para.dim3 + para.dim1) {
                dlogL_theta_Est_theta_a = (logL_part1_theta_Est_plus_theta_a_plus - logL_part1_theta_Est_minus_theta_a_plus)
                    - (logL_part1_theta_Est_plus_theta_a_minus - logL_part1_theta_Est_minus_theta_a_minus);
            }
            else {
                dlogL_theta_Est_theta_a = (logL_part2_theta_Est_plus_theta_a_plus - logL_part2_theta_Est_minus_theta_a_plus)
                    - (logL_part2_theta_Est_plus_theta_a_minus - logL_part2_theta_Est_minus_theta_a_minus);
            }

            ArrayXi FirmIndex = (dlogL_theta_Est_theta_a.abs() > 0).cast<int>().rowwise().maxCoeff();
            int Num = FirmIndex.sum() - 1;
            Num = abs(Num);
            cout << "k = " << k << "; kk = " << kk << "; Num = " << Num << endl;
            db_dtheta0_temp(kk,k) = dlogL_theta_Est_theta_a.sum() / double(Num) / eps_theta_a / 2.0 / eps_theta0 / 2.0;
            cout << "db_dtheta0_temp(kk,k) = " << db_dtheta0_temp(kk,k) << endl;
        }
        // throw runtime_error("1761");
        cout << "db_dtheta0_temp.col(k) = " << db_dtheta0_temp.col(k).transpose() << endl;
        // throw runtime_error("1617");
    }
    writeToCSVfile("db_dtheta0_v2_temp.csv", db_dtheta0_temp.cast<double>().matrix());
    ArrayXXd db_dtheta0_tran = db_dtheta0_temp.transpose();
    ArrayXXd db_dtheta0 = db_dtheta0_temp*0.5+db_dtheta0_tran*0.5;
    // ArrayXXd db_dtheta0 = db_dtheta0_temp;
    return db_dtheta0;
}

tuple<ArrayXXd,ArrayXXd,ArrayXXd> alias::Cal_logL_dtheta_Est_theta_a(const ArrayXd & theta_a,
    const SimData & sim_data_Para_GoodState_LaborInt,const SimData & sim_data_Para_BadState_LaborInt,
    const SimData & sim_data_Para_GoodState_CapitalInt,const SimData & sim_data_Para_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

    tuple<ArrayXXd,ArrayXXd,double> t_part1_theta_Est_theta_a = EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(
        theta_a, sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,
        sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd logL_part1_theta_Est_theta_a = get<1>(t_part1_theta_Est_theta_a);

    tuple<ArrayXXd,ArrayXXd,double> t_part2_theta_Est_theta_a = EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(
        theta_a,sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,
        sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd logL_part2_theta_Est_theta_a = get<1>(t_part2_theta_Est_theta_a);

    tuple<ArrayXXd,double> t_part3_theta_Est_plus_theta_a = EstimationIndirectInference_SD_l_Part3_SimData(
        theta_a,sim_data_Para_GoodState_LaborInt,sim_data_Para_BadState_LaborInt,
        sim_data_Para_GoodState_CapitalInt,sim_data_Para_BadState_CapitalInt,
        threadsManagement);
    ArrayXXd logL_part3_theta_Est_theta_a = get<0>(t_part3_theta_Est_plus_theta_a);

    return tuple<ArrayXXd,ArrayXXd,ArrayXXd>(logL_part1_theta_Est_theta_a,logL_part2_theta_Est_theta_a,logL_part3_theta_Est_theta_a);
}

tuple<ArrayXXd,double> alias::EstimationIndirectInference_SD_l_Part3_SimData(const ArrayXd & theta_a,
    const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
    const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
    MultiThreads::Threads_Management &threadsManagement) {

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
    cout << "theta3_a = " << theta_a_point.transpose() << endl;

    int GoodState; int LaborIntensive;
    GoodState = 1; LaborIntensive = 1;
    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_GoodState_LaborInt = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
    ParaEst para_est_GoodState_LaborInt = get<0>(t_Equ_GoodState_LaborInt);
    ParaVec para_vec_GoodState_LaborInt = get<1>(t_Equ_GoodState_LaborInt);
    EquStateV EquV_a_GoodState_LaborInt= get<2>(t_Equ_GoodState_LaborInt);
    EquStateVmat Evalmat_a_GoodState_LaborInt = get<3>(t_Equ_GoodState_LaborInt);

    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3_sim_GoodState_LaborInt =
        EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_GoodState_LaborInt,
        para_est_GoodState_LaborInt, para_vec_GoodState_LaborInt, EquV_a_GoodState_LaborInt, Evalmat_a_GoodState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_part3_Good_LaborInt = get<3>(t_Part3_sim_GoodState_LaborInt);
    double Nprob_part3_Good_LaborInt = get<2>(t_Part3_sim_GoodState_LaborInt);

    GoodState = 0; LaborIntensive = 1;
    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_BadState_LaborInt = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
        ParaEst para_est_BadState_LaborInt = get<0>(t_Equ_BadState_LaborInt);
        ParaVec para_vec_BadState_LaborInt = get<1>(t_Equ_BadState_LaborInt);
        EquStateV EquV_a_BadState_LaborInt= get<2>(t_Equ_BadState_LaborInt);
        EquStateVmat Evalmat_a_BadState_LaborInt = get<3>(t_Equ_BadState_LaborInt);

    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3_sim_BadState_LaborInt =
        EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_BadState_LaborInt,
        para_est_BadState_LaborInt, para_vec_BadState_LaborInt, EquV_a_BadState_LaborInt, Evalmat_a_BadState_LaborInt,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_part3_Bad_LaborInt = get<3>(t_Part3_sim_BadState_LaborInt);
    double Nprob_part3_Bad_LaborInt = get<2>(t_Part3_sim_BadState_LaborInt);

    GoodState = 1; LaborIntensive = 0;
    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_GoodState_CapitalInt = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
    ParaEst para_est_GoodState_CapitalInt = get<0>(t_Equ_GoodState_CapitalInt);
    ParaVec para_vec_GoodState_CapitalInt = get<1>(t_Equ_GoodState_CapitalInt);
    EquStateV EquV_a_GoodState_CapitalInt = get<2>(t_Equ_GoodState_CapitalInt);
    EquStateVmat Evalmat_a_GoodState_CapitalInt = get<3>(t_Equ_GoodState_CapitalInt);

    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3_sim_GoodState_CapitalInt =
        EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_GoodState_CapitalInt,
        para_est_GoodState_CapitalInt, para_vec_GoodState_CapitalInt, EquV_a_GoodState_CapitalInt, Evalmat_a_GoodState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_part3_Good_CapitalInt = get<3>(t_Part3_sim_GoodState_CapitalInt);
    double Nprob_part3_Good_CapitalInt = get<2>(t_Part3_sim_GoodState_CapitalInt);

    GoodState = 0; LaborIntensive = 0;
    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_BadState_CapitalInt = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
        ParaEst para_est_BadState_CapitalInt = get<0>(t_Equ_BadState_CapitalInt);
        ParaVec para_vec_BadState_CapitalInt = get<1>(t_Equ_BadState_CapitalInt);
        EquStateV EquV_a_BadState_CapitalInt = get<2>(t_Equ_BadState_CapitalInt);
        EquStateVmat Evalmat_a_BadState_CapitalInt = get<3>(t_Equ_BadState_CapitalInt);
    tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
        ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3_sim_BadState_CapitalInt =
        EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_BadState_CapitalInt,
        para_est_BadState_CapitalInt, para_vec_BadState_CapitalInt, EquV_a_BadState_CapitalInt, Evalmat_a_BadState_CapitalInt,
        GoodState, LaborIntensive, threadsManagement);
    ArrayXXd logL_part3_Bad_CapitalInt = get<3>(t_Part3_sim_BadState_CapitalInt);
    double Nprob_part3_Bad_CapitalInt = get<2>(t_Part3_sim_BadState_CapitalInt);

    // cout << "grad_a3_sum_Bad_CapitalInt = " << grad_a3_sum_Bad_CapitalInt.transpose() << endl;
    // cout << "part3_Bad_CapitalInt size = " << dlogL_part3_Bad_CapitalInt.rows() << "; " << dlogL_part3_Bad_CapitalInt.cols() << endl;
    // cout << "dlogL_part3_Bad_CapitalInt = " << dlogL_part3_Bad_CapitalInt.colwise().sum() << endl;
    // throw runtime_error("895");
    ArrayXXd logL_part3 = AppendMatrix(logL_part3_Good_LaborInt,logL_part3_Bad_LaborInt,
        logL_part3_Good_CapitalInt,logL_part3_Bad_CapitalInt);

    double Nprob_part3 = Nprob_part3_Good_LaborInt + Nprob_part3_Bad_LaborInt + Nprob_part3_Good_CapitalInt + Nprob_part3_Bad_CapitalInt;
    // cout << "Nprob_part3 = " << Nprob_part3 << endl;

    return tuple<ArrayXXd,double>(logL_part3,Nprob_part3);
}




ArrayXXd alias::Vectorize_dl_dtheta_part(const ArrayXXd & dL_part, const int & dim) {

    ArrayXXd dL = ArrayXXd::Zero(dL_part.rows()*para.EstTbar,dim);

    for (size_t kk = 0; kk < dim; kk++) {
        ArrayXXd temp1 = dL_part.middleCols(kk*para.EstTbar,para.EstTbar);
        ArrayXd temp1_vec = ArrayXd::Zero(dL_part.rows()*para.EstTbar);
        for (size_t tt = 0; tt < para.EstTbar; tt++) {
            temp1_vec.segment(tt*temp1.col(tt).rows(),temp1.col(tt).rows()) = temp1.col(tt);
        }
        dL.col(kk) = temp1_vec;
    }

    return dL;
}

ArrayXXd alias::Vectorize_dl_dtheta(const ArrayXXd & dL_part1,const ArrayXXd & dL_part2,const ArrayXXd & dL_part3) {

    ArrayXXd dL = ArrayXXd::Zero(dL_part1.rows()*para.EstTbar,para.dim1+para.dim2+para.dim3);

    for (size_t kk = 0; kk < para.dim3; kk++) {
        ArrayXXd temp1 = dL_part3.middleCols(kk*para.EstTbar,para.EstTbar);
        ArrayXd temp1_vec = ArrayXd::Zero(dL_part3.rows()*para.EstTbar);
        for (size_t tt = 0; tt < para.EstTbar; tt++) {
            temp1_vec.segment(tt*temp1.col(tt).rows(),temp1.col(tt).rows()) = temp1.col(tt);
        }
        dL.col(kk) = temp1_vec;
    }
    for (size_t kk = 0; kk < para.dim1; kk++) {
        ArrayXXd temp1 = dL_part1.middleCols(kk*para.EstTbar,para.EstTbar);
        ArrayXd temp1_vec = ArrayXd::Zero(dL_part1.rows()*para.EstTbar);
        for (size_t tt = 0; tt < para.EstTbar; tt++) {
            temp1_vec.segment(tt*temp1.col(tt).rows(),temp1.col(tt).rows()) = temp1.col(tt);
        }
        dL.col(para.dim3+kk) = temp1_vec;
    }
    for (size_t kk = 0; kk < para.dim2; kk++) {
        ArrayXXd temp1 = dL_part2.middleCols(kk*para.EstTbar,para.EstTbar);
        ArrayXd temp1_vec = ArrayXd::Zero(dL_part2.rows()*para.EstTbar);
        for (size_t tt = 0; tt < para.EstTbar; tt++) {
            temp1_vec.segment(tt*temp1.col(tt).rows(),temp1.col(tt).rows()) = temp1.col(tt);
        }
        dL.col(para.dim3+para.dim1+kk) = temp1_vec;
    }

    return dL;
}


ArrayXXd alias::AppendMatrix(const ArrayXXd & X1, const ArrayXXd & X2,const ArrayXXd & X3, const ArrayXXd & X4) {

    int N1 = X1.rows();
    int N2 = X2.rows();
    int N3 = X3.rows();
    int N4 = X4.rows();

    int Ncol = X1.cols();

    ArrayXXd X = ArrayXXd::Zero(N1+N2+N3+N4,Ncol);
    X.middleRows(0,N1) = X1;
    X.middleRows(N1,N2) = X2;
    X.middleRows(N1+N2,N3) = X3;
    X.middleRows(N1+N2+N3,N4) = X4;

    return X;
}



ArrayXXd alias::VarianceMatrix_dL(const ArrayXXd & dL_AllPart) {

    ArrayXXd Var_dL = ArrayXXd::Zero(dL_AllPart.cols(),dL_AllPart.cols());
    // for (int i = 0; i < dL_AllPart.cols(); ++i) {
    //     for (int j = 0; j < dL_AllPart.cols(); ++j) {
    //         ArrayXd temp_i = dL_AllPart.col(i);
    //         ArrayXd temp_j = dL_AllPart.col(j);
    //
    //         double temp_i_mean = temp_i.mean();
    //         double temp_j_mean = temp_j.mean();
    //
    //         ArrayXd dtemp_i = (temp_i - temp_i_mean);
    //         ArrayXd dtemp_j = (temp_j - temp_j_mean);
    //
    //         Var_dL(i,j) = ( dtemp_i * dtemp_j ).mean();
    //     }
    // }
    for (int i = 0; i < dL_AllPart.cols(); ++i) {
        for (int j = 0; j < dL_AllPart.cols(); ++j) {

            ArrayXd temp_i = dL_AllPart.col(i);
            ArrayXd temp_j = dL_AllPart.col(j);

            ArrayXi temp_i_nonzero = (temp_i.abs() > 0.0).cast<int>();
            ArrayXi temp_j_nonzero = (temp_j.abs() > 0.0).cast<int>();
            ArrayXd temp_ij_nonzero = (temp_i_nonzero * temp_j_nonzero).cast<double>();

            // cout << "temp_ij_nonzero.sum() = " << temp_ij_nonzero.sum()  << endl;

            if (temp_ij_nonzero.sum() == 0) {
                Var_dL(i,j) = 0.0;
            }
            else
            {
                double temp_i_mean = (temp_i*temp_ij_nonzero).sum() / double( temp_ij_nonzero.sum() );
                double temp_j_mean = (temp_j*temp_ij_nonzero).sum() / double( temp_ij_nonzero.sum() );

                ArrayXd dtemp_i = (temp_i - temp_i_mean) * temp_ij_nonzero.cast<double>();
                ArrayXd dtemp_j = (temp_j - temp_j_mean) * temp_ij_nonzero.cast<double>();

                // Var_dL(i,j) = (dtemp_i * dtemp_j).sum() / double( temp_i_nonzero.sum() );
                Var_dL(i,j) = ( dtemp_i * dtemp_j * temp_ij_nonzero).sum() / double( temp_ij_nonzero.sum() );
                // Var_dL(i,j) = (dtemp_i * dtemp_j * temp_ij_nonzero).mean();

            }
        }
    }
    return Var_dL;
}

//



// //

//

//
// //


//
//

//
//


//


// ArrayXXd alias::Cal_db_dtheta0_v1(const ArrayXd & theta_Est,const ArrayXd & theta_a,
//     const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
//     const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
//     const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
//     const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
//     MultiThreads::Threads_Management &threadsManagement) {
//
//     int GoodState; int LaborIntensive;
//
//     GoodState = 1; LaborIntensive = 1;
//     EquV_Auxiliary EquV_GoodState_LaborInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
//     GoodState = 0; LaborIntensive = 1;
//     EquV_Auxiliary EquV_BadState_LaborInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
//     GoodState = 1; LaborIntensive = 0;
//     EquV_Auxiliary EquV_GoodState_CapitalInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
//     GoodState = 0; LaborIntensive = 0;
//     EquV_Auxiliary EquV_BadState_CapitalInt = CalEquV_plus_minus(theta_a, GoodState, LaborIntensive, threadsManagement);
//
//     /*** part 1 ***/
//     ArrayXXd db_dtheta0_temp = ArrayXXd::Zero(para.dim1+para.dim2+para.dim3,para.dim1+para.dim2+para.dim3);
//     for (size_t k = 0; k < para.dim1+para.dim2+para.dim3; ++k) {
//         cout << "k = " << k << endl;
//         double eps = abs(theta_Est(k)) * 1e-6;
//
//         ArrayXd theta_Est_plus = theta_Est;
//         theta_Est_plus(k) = theta_Est(k) + eps;
//         cout << "theta_Est_plus = " << theta_Est_plus.transpose() << endl;
//
//         GoodState = 1; LaborIntensive = 1;
//         SimData sim_data_Para_GoodState_LaborInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
//             sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
//             GoodState, LaborIntensive, threadsManagement);
//         GoodState = 0; LaborIntensive = 1;
//         SimData sim_data_Para_BadState_LaborInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
//             sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
//             GoodState, LaborIntensive, threadsManagement);
//         GoodState = 1; LaborIntensive = 0;
//         SimData sim_data_Para_GoodState_CapitalInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
//             sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
//             GoodState, LaborIntensive, threadsManagement);
//         GoodState = 0; LaborIntensive = 0;
//         SimData sim_data_Para_BadState_CapitalInt_plus = SimulationData_Step2Estimation_FullVersion(theta_Est_plus,
//             sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
//             GoodState, LaborIntensive, threadsManagement);
//
//         tuple<ArrayXXd,ArrayXXd,double> t_part1_plus = EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(theta_a,
//             sim_data_Para_GoodState_LaborInt_plus,sim_data_Para_BadState_LaborInt_plus,
//             sim_data_Para_GoodState_CapitalInt_plus,sim_data_Para_BadState_CapitalInt_plus,
//             threadsManagement);
//         ArrayXXd dL_part1_plus = get<0>(t_part1_plus);
//         ArrayXXd L_part1_plus = get<0>(t_part1_plus);
//         double Nprob_part1_plus = get<2>(t_part1_plus);
//         tuple<ArrayXXd,ArrayXXd,double> t_part2_plus  = EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(theta_a,
//             sim_data_Para_GoodState_LaborInt_plus,sim_data_Para_BadState_LaborInt_plus,
//             sim_data_Para_GoodState_CapitalInt_plus,sim_data_Para_BadState_CapitalInt_plus,
//             threadsManagement);
//         ArrayXXd dL_part2_plus = get<0>(t_part2_plus);
//         double Nprob_part2_plus = get<2>(t_part2_plus);
//         tuple<ArrayXXd,ArrayXXd,double> t_part3_plus = EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
//             sim_data_Para_GoodState_LaborInt_plus,sim_data_Para_BadState_LaborInt_plus,
//             sim_data_Para_GoodState_CapitalInt_plus,sim_data_Para_BadState_CapitalInt_plus,
//             EquV_GoodState_LaborInt,EquV_BadState_LaborInt,
//             EquV_GoodState_CapitalInt,EquV_BadState_CapitalInt,
//             threadsManagement);
//         ArrayXXd dL_part3_plus = get<0>(t_part3_plus);
//
//         ArrayXd theta_Est_minus = theta_Est;
//         theta_Est_minus(k) = theta_Est(k) - eps;
//         cout << "theta_Est_minus = " << theta_Est_minus.transpose() << endl;
//         cout << "diff theta_Est = " << theta_Est_plus.transpose() - theta_Est_minus.transpose() << endl;
//
//         GoodState = 1; LaborIntensive = 1;
//         SimData sim_data_Para_GoodState_LaborInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
//             sim_var_2step_GoodState_LaborInt,para.p_K_StatusQuo, simdata_panel_GoodState_LaborInt,
//             GoodState, LaborIntensive, threadsManagement);
//         GoodState = 0; LaborIntensive = 1;
//         SimData sim_data_Para_BadState_LaborInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
//             sim_var_2step_BadState_LaborInt,para.p_K_StatusQuo, simdata_panel_BadState_LaborInt,
//             GoodState, LaborIntensive, threadsManagement);
//         GoodState = 1; LaborIntensive = 0;
//         SimData sim_data_Para_GoodState_CapitalInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
//             sim_var_2step_GoodState_CapitalInt,para.p_K_StatusQuo, simdata_panel_GoodState_CapitalInt,
//             GoodState, LaborIntensive, threadsManagement);
//         GoodState = 0; LaborIntensive = 0;
//         SimData sim_data_Para_BadState_CapitalInt_minus = SimulationData_Step2Estimation_FullVersion(theta_Est_minus,
//             sim_var_2step_BadState_CapitalInt,para.p_K_StatusQuo, simdata_panel_BadState_CapitalInt,
//             GoodState, LaborIntensive, threadsManagement);
//
//         tuple<ArrayXXd,ArrayXXd,double> t_part1_minus = EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(theta_a,
//             sim_data_Para_GoodState_LaborInt_minus,sim_data_Para_BadState_LaborInt_minus,
//             sim_data_Para_GoodState_CapitalInt_minus,sim_data_Para_BadState_CapitalInt_minus,
//             threadsManagement);
//         ArrayXXd dL_part1_minus = get<0>(t_part1_minus);
//         ArrayXXd L_part1_minus = get<1>(t_part1_minus);
//         double Nprob_part1_minus = get<2>(t_part1_minus);
//         tuple<ArrayXXd,ArrayXXd,double> t_part2_minus = EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(theta_a,
//             sim_data_Para_GoodState_LaborInt_minus,sim_data_Para_BadState_LaborInt_minus,
//             sim_data_Para_GoodState_CapitalInt_minus,sim_data_Para_BadState_CapitalInt_minus,
//             threadsManagement);
//         ArrayXXd dL_part2_minus = get<0>(t_part2_minus);
//         double Nprob_part2_minus = get<2>(t_part2_minus);
//         tuple<ArrayXXd,ArrayXXd,double> t_part3_minus = EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
//             sim_data_Para_GoodState_LaborInt_minus,sim_data_Para_BadState_LaborInt_minus,
//             sim_data_Para_GoodState_CapitalInt_minus,sim_data_Para_BadState_CapitalInt_minus,
//             EquV_GoodState_LaborInt,EquV_BadState_LaborInt,
//             EquV_GoodState_CapitalInt,EquV_BadState_CapitalInt,
//             threadsManagement);
//         ArrayXXd dL_part3_minus = get<0>(t_part3_minus);
//
//         ArrayXXd dL_plus_vec = Vectorize_dl_dtheta(dL_part1_plus, dL_part2_plus, dL_part3_plus);
//         ArrayXXd dL_minus_vec = Vectorize_dl_dtheta(dL_part1_minus, dL_part2_minus, dL_part3_minus);
//
//         cout << "dL_plus_vec = " << dL_plus_vec.colwise().sum() << endl;
//         cout << "dL_minus_vec = " << dL_minus_vec.colwise().sum() << endl;
//
//         db_dtheta0_temp.row(k) = ( -dL_plus_vec - (-dL_minus_vec) ).colwise().sum() / 2.0 / eps;
//         // cout << "double(ddL_part3.rows()) = " << double(ddL_part3.rows()) << endl;
//         cout << "db_dtheta0_temp.row(k) = " << db_dtheta0_temp.row(k) << endl;
//     }
//     writeToCSVfile("db_dtheta0_v1_temp.csv", db_dtheta0_temp.cast<double>().matrix());
//     ArrayXXd db_dtheta0 = db_dtheta0_temp;
//     return db_dtheta0;
// }
//
//

//


//
// tuple<ArrayXXd,int> alias::Cal_dlogL_Part3_RealData(const SimData & ReadData, const ArrayXd & theta_a,
//     const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement) {
//
//     for (size_t kk = 0; kk < para.dim3; kk++) {
//         double eps = abs(theta_a(kk)) * 1e-6;
//         ArrayXd theta_a_plus = theta_a;
//         theta_a_plus(kk) = theta_a(kk) + eps;
//
//
//     ArrayXd theta1_a(para.dim1);
//     for (size_t i = 0; i < para.dim1; ++i) {theta1_a(i) = theta_a(para.dim3+i);}
//     // cout << "theta1_a = " << theta1_a.transpose() << endl;
//
//
//     ArrayXd theta2_a(para.dim2);
//     for (size_t i = 0; i < para.dim2; ++i) {theta2_a(i) = theta_a(para.dim3+para.dim1+i);}
//     // cout << "theta2_a = " << theta2_a.transpose() << endl;
//
//
//     ArrayXd theta_a_point = theta_a;
//     std::vector<double> theta3_a_point(para.dim3);
//     for (size_t i = 0; i < para.dim3; ++i) {
//         theta3_a_point[i] = theta_a_point(i);
//     }
//     cout << "theta3_a = " << theta_a_point.transpose() << endl;
//
//
//     int GoodState; int LaborIntensive;
//     GoodState = 1; LaborIntensive = 1;
//     tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_GoodState_LaborInt = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
//         theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
//     ParaEst para_est_GoodState_LaborInt = get<0>(t_Equ_GoodState_LaborInt);
//     ParaVec para_vec_GoodState_LaborInt = get<1>(t_Equ_GoodState_LaborInt);
//     EquStateV EquV_a_GoodState_LaborInt= get<2>(t_Equ_GoodState_LaborInt);
//     EquStateVmat Evalmat_a_GoodState_LaborInt = get<3>(t_Equ_GoodState_LaborInt);
//
//
//     tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
//         ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3_sim_GoodState_LaborInt =
//         EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_GoodState_LaborInt,
//         para_est_GoodState_LaborInt, para_vec_GoodState_LaborInt, EquV_a_GoodState_LaborInt, Evalmat_a_GoodState_LaborInt,
//         GoodState, LaborIntensive, threadsManagement);
//     ArrayXXd logL_part3_Good_LaborInt = get<3>(t_Part3_sim_GoodState_LaborInt);
//     double Nprob_part3_Good_LaborInt = get<2>(t_Part3_sim_GoodState_LaborInt);
//
//
//     GoodState = 0; LaborIntensive = 1;
//     tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ_BadState_LaborInt = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
//         theta3_a_point, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
//         ParaEst para_est_BadState_LaborInt = get<0>(t_Equ_BadState_LaborInt);
//         ParaVec para_vec_BadState_LaborInt = get<1>(t_Equ_BadState_LaborInt);
//         EquStateV EquV_a_BadState_LaborInt= get<2>(t_Equ_BadState_LaborInt);
//         EquStateVmat Evalmat_a_BadState_LaborInt = get<3>(t_Equ_BadState_LaborInt);
//
//
//     tuple<double,int,double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd,
//         ArrayXXi,ArrayXXi,ArrayXXi,ArrayXXd,ArrayXXd,ArrayXXd> t_Part3_sim_BadState_LaborInt =
//         EstimationAuxiliaryModel_part3_Simulation_FullVersion_Step2(sim_data_Para_BadState_LaborInt,
//         para_est_BadState_LaborInt, para_vec_BadState_LaborInt, EquV_a_BadState_LaborInt, Evalmat_a_BadState_LaborInt,
//         GoodState, LaborIntensive, threadsManagement);
//     ArrayXXd logL_part3_Bad_LaborInt = get<3>(t_Part3_sim_BadState_LaborInt);
//     double Nprob_part3_Bad_LaborInt = get<2>(t_Part3_sim_BadState_LaborInt);
//
//         ArrayXd theta_a_minus = theta_a;
//         theta_a_minus(kk) = theta_a(kk) - eps;
//     }
//     tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi> t_Part3 = EstimationAuxiliaryModel_part3_likelihood(
//         ReadData, EquV.para_est_a_point, EquV.para_vec_a_point, EquV.EquV_a_point, EquV.Evalmat_a_point,
//         GoodState, LaborIntensive, threadsManagement);
//     double obj_point_sum = get<0>(t_Part3);
//     int Ncount_part3 = get<1>(t_Part3);
//     ArrayXXd logL_part3 = get<2>(t_Part3);
//     ArrayXXi K_index_max_mat = get<3>(t_Part3);
//     ArrayXXi Lur_index_max_mat = get<4>(t_Part3);
//     ArrayXXi Luc_index_max_mat = get<5>(t_Part3);
//
//     ArrayXXd dlogL_part3 = ArrayXXd::Zero(ReadData.good.rows(),para.EstTbar * para.dim3);
//     ArrayXd grad_a3_sum = ArrayXd::Zero(para.dim3);
//     for (size_t i_a = 0; i_a < para.dim3; ++i_a) {
//     //        cout << "i_a = " << i_a << endl;
//         double obj_plus_sum = 0.0;
//         if ( ( (i_a == 6 or i_a == 7 or i_a == 10 or i_a == 15) and GoodState == 1 ) or
//             ( (i_a == 4 or i_a == 5 or i_a == 9 or i_a == 14) and GoodState == 0 ) ) {
//             obj_plus_sum = 0.0;
//             grad_a3_sum(i_a) = 0.0;
//         }
//         else {
//             tuple<double,int,ArrayXXd> t_plus = EstimationAuxiliaryModel_part3_likelihood_Diff(ReadData,
//                 EquV.para_est_a_plus[i_a], EquV.para_vec_a_plus[i_a], EquV.EquV_a_plus[i_a], EquV.Evalmat_a_plus[i_a],
//                 K_index_max_mat, Lur_index_max_mat, Luc_index_max_mat,
//                 GoodState, LaborIntensive, threadsManagement);
//             obj_plus_sum = get<0>(t_plus);
//             grad_a3_sum(i_a) = (obj_plus_sum - obj_point_sum) / EquV.eps_diff(i_a);
//             ArrayXXd logL_part3_plus = get<2>(t_plus);
//             dlogL_part3.middleCols(i_a*para.EstTbar,para.EstTbar) = (-logL_part3_plus - (-logL_part3)) / EquV.eps_diff(i_a);
//         }
//     }
//
//     return tuple<ArrayXXd,int>(dlogL_part3,Ncount_part3);
// }
