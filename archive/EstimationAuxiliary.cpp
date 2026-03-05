#include "EstimationAuxiliary.h"
using namespace Eigen;
using namespace std;
using namespace alias;
//
using Ipopt_Wrapper::Ipopt_Solver;
using Ipopt_Wrapper::Option;
using Ipopt_Wrapper::Ipopt_Status;
using Ipopt_Wrapper::Simple_Opt_Problem;

/**************************************************************
* Version 2: Choose a specific inital value and see the estimation
**************************************************************/
tuple<ArrayXd,ArrayXd,ArrayXd> alias::EstimationAuxiliaryModel_Output_Version2(const SimData & ReadData_GoodState_LaborInt,
    const SimData & ReadData_BadState_LaborInt, const SimData & ReadData_GoodState_CapitalInt,
    const SimData & ReadData_BadState_CapitalInt, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd theta_Initial = ArrayXd::Zero(para.dim3);
    /**** Para Est 3: only use the previous year data ****/
    theta_Initial << 106.387, 120.898, 1661.64, 106.151, 236.549, 169.232, 291.717, 200.49,
        2.57276, 41.85, 37.7904, 0.495642, 0.379922, 0.804774, -423.02, -628.051, 536.874;

    // theta_Initial << 106.387, 120.898, 1661.64, 106.151, 236.549, 169.232, 291.717, 200.49,
    //     2.57276, 41.85, 37.7904, 0.495642, 0.379922, 0.804774, -423.02, -628.051, 536.874;
    // f = 1.77780256173, error = 0.00720603768382
    // --------------------------------------
    // optimal_x = 64.3591932411; 168.407327531; 1662.57579102; 115.745612122; 247.061467265; 178.294113248; 273.716457993; 195.673380953; 5.36315306194; 31.7249969942; 30.3631367811; 0.380090876926; 0.377281880635; 0.683080451129; -414.517906307; -621.987707083; 508.86704831; ;
    // optimal f = 1.77780256173
    // Part3 optimal value: theta3_a =  64.3591932411  168.407327531  1662.57579102  115.745612122  247.061467265  178.294113248  273.716457993  195.673380953  5.36315306194  31.7249969942  30.3631367811 0.380090876926 0.377281880635 0.683080451129 -414.517906307 -621.987707083   508.86704831

    //
    // theta_Initial <<  61.6331, 179.402, 1662.86,103.318,218.532,178.482,256.496,200.002,
    //     5.19303,38.2064,36.8367,0.381529,0.3766,0.680907,-413.768,-626.425,528.122;
    // f = 1.77775, error = 0.00709515
    // --------------------------------------
    // optimal_x = 56.9062; 167.793; 1663.14; 115.166; 243.303; 180.599; 251.446; 194.367; 5.36631; 31.593; 30.2341; 0.379918; 0.377281; 0.68109; -409.068; -622.348; 502.87; ;
    // optimal f = 1.77775
    // Part3 optimal value: theta3_a =  56.9062  167.793  1663.14  115.166  243.303  180.599  251.446  194.367  5.36631   31.593  30.2341 0.379918 0.377281  0.68109 -409.068 -622.348   502.87

    tuple<ArrayXd,ArrayXd,ArrayXd> t_AuxPara = SetupInitialParaGuess_forEstimation(theta_Initial);
    ArrayXd thetaData_1_ini = get<0>(t_AuxPara);
    ArrayXd thetaData_2_ini = get<1>(t_AuxPara);
    ArrayXd thetaData_3_ini = get<2>(t_AuxPara);
    ArrayXd theta1_a = thetaData_1_ini;
    ArrayXd theta2_a = thetaData_2_ini;
    ArrayXd theta3_a = thetaData_3_ini;
    double obj3_a = 0.0;
//
    tuple<ArrayXd,ArrayXd,ArrayXd,double> t_AuxillaryEst = EstimationAuxiliaryModel(
        ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,
        ReadData_GoodState_CapitalInt, ReadData_BadState_CapitalInt,
        thetaData_1_ini,thetaData_2_ini,thetaData_3_ini,threadsManagement);
    theta1_a = get<0>(t_AuxillaryEst);
    theta2_a = get<1>(t_AuxillaryEst);
    theta3_a = get<2>(t_AuxillaryEst);
    obj3_a = get<3>(t_AuxillaryEst);

    return tuple<ArrayXd,ArrayXd,ArrayXd>(theta1_a,theta2_a,theta3_a);
}

/**************************************************************
////// Initial Guess for Parameters in the Auxilliary model
**************************************************************/
tuple<ArrayXd,ArrayXd,ArrayXd> alias::SetupInitialParaGuess_forEstimation(const ArrayXd & theta_RandomInitial_select) {

    // deflated
    //// gamma0;gamma1;sigma_phi_eps; alpha_tilde_K;alpha_Lr
    ArrayXd thetaData_1(para.dim1);
    thetaData_1 << 0.159361,0.165752,0.902305,0.787875; // use all periods

    //// F_P; sigma_PD
    ArrayXd thetaData_2(para.dim2);
    thetaData_2 << -14.2557, 14.3706, 11.6135;  // good states --labor intensive sector

    //// H_K(0);c_HK(1);F_K(2);c_FK(3);    H_ur(4);c_H_ur(5);F_L_G_ur(6);F_H_G_ur(7);c_F_ur(8);
    //// H_uc(9);c_H_uc(10);F_G_uc(11);c_F_uc(12);     sigma_Kerror(13);sigma_Lerror_ur(14);sigma_Lerror_uc(15);
    //// F_G_E(16);simga_SE(17)
    ArrayXd thetaData_3(para.dim3);
    thetaData_3 = theta_RandomInitial_select;

    return tuple<ArrayXd,ArrayXd,ArrayXd>(thetaData_1,thetaData_2,thetaData_3);
}

/**************************************************************
* Estimating the auxilliary model
**************************************************************/
tuple<ArrayXd,ArrayXd,ArrayXd,double> alias::EstimationAuxiliaryModel(const SimData & RealData_GoodState_LaborInt,
    const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
    const SimData & RealData_BadState_CapitalInt,
    const ArrayXd & theta1_a_ini, const ArrayXd & theta2_a_ini, const ArrayXd & theta3_a_ini,
    MultiThreads::Threads_Management & threadsManagement) {

    //// Part 1: The production function and productivity process
    cout << "Part1 initial value: theta1_a_ini = " << theta1_a_ini.transpose() << endl;
    ArrayXd theta1_a(para.dim1);
    theta1_a = EstimationAuxiliaryModel_part1_theta1(RealData_GoodState_LaborInt,
        RealData_BadState_LaborInt, RealData_GoodState_CapitalInt, RealData_BadState_CapitalInt,
        theta1_a_ini,threadsManagement);
    cout << "Part1 optimal value: theta1_a = " << std::setprecision(12) << theta1_a.transpose() << endl;
//    throw runtime_error("89");

    cout << "Part2 initial value: theta2_a_ini = " << theta2_a_ini.transpose() << endl;
    ArrayXd theta2_a(para.dim2);
    theta2_a = EstimationAuxiliaryModel_part2_theta2(RealData_GoodState_LaborInt,
        RealData_BadState_LaborInt, RealData_GoodState_CapitalInt, RealData_BadState_CapitalInt,
        theta1_a,theta2_a_ini,threadsManagement);
    cout << "Part2 optimal value: theta2_a = " << std::setprecision(12) << theta2_a.transpose() << endl;
//    throw runtime_error("96");
//
//    cout << "Part3 initial value: theta3_a_ini = " << theta3_a_ini.transpose() << endl;
    ArrayXd theta3_a(para.dim3);
    tuple<ArrayXd,double> t = EstimationAuxiliaryModel_part3_theta3(RealData_GoodState_LaborInt,
        RealData_BadState_LaborInt, RealData_GoodState_CapitalInt, RealData_BadState_CapitalInt,
        theta1_a,theta2_a,theta3_a_ini,threadsManagement);
    theta3_a = get<0>(t);
    double obj3_a = get<1>(t);
//    theta3_a = theta3_a_ini;
//    double obj3_a = 0.0;
    cout << "Part3 optimal value: theta3_a = " << std::setprecision(12) << theta3_a.transpose() << endl;
    throw runtime_error("1179");
//////    EstimationAuxiliaryModel_part3_theta3_checking(theta3_a,RealData,theta1_a,theta2_a,threadsManagement);

    return tuple<ArrayXd,ArrayXd,ArrayXd,double>(theta1_a,theta2_a,theta3_a,obj3_a);
}
//
/******************************************************************
 * Auxilliary Model Estimation: Part 1
 *******************************************************************/
//// Given the objective function, find the estimates by maximize the objective function
ArrayXd alias::EstimationAuxiliaryModel_part1_theta1(const SimData & RealData_GoodState_LaborInt,
    const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
    const SimData & RealData_BadState_CapitalInt, const ArrayXd & theta1_a_ini,
    MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd theta1_a(para.dim1);

    std::vector<double> x0_vec(para.dim1);
    for (size_t k = 0; k < para.dim1; ++k) { x0_vec[k] = theta1_a_ini(k); }

    auto f_without_hessian = [&RealData_GoodState_LaborInt,&RealData_BadState_LaborInt,
        &RealData_GoodState_CapitalInt,&RealData_BadState_CapitalInt,&threadsManagement]
        (std::vector<double> x, std::vector<double> & grad){

        int GoodState = 1; int LaborIntensive = 1;
        tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Labor = EstimationAuxiliaryModel_part1_ObjFun(x,
            RealData_GoodState_LaborInt,GoodState,LaborIntensive,threadsManagement);
        double obj_sum_Good_Labor = get<0>(t_Good_Labor);
        std::vector<double> grad_sum_Good_Labor = get<1>(t_Good_Labor);
        int Ncount_Good_Labor = get<2>(t_Good_Labor);;

        GoodState = 0; LaborIntensive = 1;
        tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Labor = EstimationAuxiliaryModel_part1_ObjFun(x,
            RealData_BadState_LaborInt,GoodState,LaborIntensive,threadsManagement);
        double obj_sum_Bad_Labor = get<0>(t_Bad_Labor);
        std::vector<double> grad_sum_Bad_Labor = get<1>(t_Bad_Labor);
        int Ncount_Bad_Labor = get<2>(t_Bad_Labor);
////
//        /* only use labor intensive sector */
//        double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor) / double(Ncount_Good_Labor+Ncount_Bad_Labor);
//        grad[0] = (grad_sum_Good_Labor[0]) / double(Ncount_Good_Labor);
//        grad[1] = (grad_sum_Bad_Labor[0]) / double(Ncount_Bad_Labor);
//        grad[2] = (grad_sum_Good_Labor[1] + grad_sum_Bad_Labor[1]) / double(Ncount_Good_Labor+Ncount_Bad_Labor);
//        grad[3] = (grad_sum_Good_Labor[2] + grad_sum_Bad_Labor[2]) / double(Ncount_Good_Labor+Ncount_Bad_Labor);

        /* use both labor and capital intensive sector */
        GoodState = 1; LaborIntensive = 0;
        tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Capital = EstimationAuxiliaryModel_part1_ObjFun(x,
            RealData_GoodState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
        double obj_sum_Good_Capital = get<0>(t_Good_Capital);
        std::vector<double> grad_sum_Good_Capital = get<1>(t_Good_Capital);
        int Ncount_Good_Capital = get<2>(t_Good_Capital);;

        GoodState = 0; LaborIntensive = 0;
        tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Capital = EstimationAuxiliaryModel_part1_ObjFun(x,
            RealData_BadState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
        double obj_sum_Bad_Capital = get<0>(t_Bad_Capital);
        std::vector<double> grad_sum_Bad_Capital = get<1>(t_Bad_Capital);
        int Ncount_Bad_Capital = get<2>(t_Bad_Capital);
//
//        /* only use labor intensive sector */
//        double obj = (obj_sum_Good_Capital + obj_sum_Bad_Capital) / double(Ncount_Good_Capital+Ncount_Bad_Capital);
//        grad[0] = (grad_sum_Good_Capital[0]) / double(Ncount_Good_Capital);
//        grad[1] = (grad_sum_Bad_Capital[0]) / double(Ncount_Bad_Capital);
//        grad[2] = (grad_sum_Good_Capital[1] + grad_sum_Bad_Capital[1]) / double(Ncount_Good_Capital+Ncount_Bad_Capital);
//        grad[3] = (grad_sum_Good_Capital[2] + grad_sum_Bad_Capital[2]) / double(Ncount_Good_Capital+Ncount_Bad_Capital);
//

        double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital)
            / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
        grad[0] = (grad_sum_Good_Labor[0] + grad_sum_Good_Capital[0])
            / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
        grad[1] = (grad_sum_Bad_Labor[0] + grad_sum_Bad_Capital[0])
            / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);

        grad[2] = (grad_sum_Good_Labor[1] + grad_sum_Bad_Labor[1] + grad_sum_Good_Capital[1] + grad_sum_Bad_Capital[1])
            / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
        grad[3] = (grad_sum_Good_Labor[2] + grad_sum_Bad_Labor[2] + grad_sum_Good_Capital[2] + grad_sum_Bad_Capital[2])
            / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
//
//        cout << "obj = " << obj << "; grad[0] = " << grad[0] << "; grad[1] = " << grad[1] << "; grad[2] = " << grad[2]
//            << "; grad[3] = " << grad[3] << endl;
//        throw runtime_error("145");
        return obj;
    };
    //// gamma0;gamma1;sigma_phi_eps; alpha_tilde_L;alpha_tilde_K;alpha_Lc;mu
    // lower bound
    std::vector<double> x_lb(para.dim1);
    x_lb[0] = -5.0; x_lb[1] = -5.0; x_lb[2] = 0.0; x_lb[3] = 0.1;
    // upper bound
    std::vector<double> x_ub(para.dim1);
    x_ub[0] = 5.0; x_ub[1] = 5.0; x_ub[2] = 1.0; x_ub[3] = 2;

    Simple_Opt_Problem problem(f_without_hessian, para.dim1);
    problem.set_x_lb(x_lb);
    problem.set_x_ub(x_ub);

    Option option;
    option.display_level = 5;
    option.test_derivative = 0;
    Ipopt_Solver solver(option);

    auto result = solver.optimize(problem, x0_vec);
    if (result.status == Ipopt_Wrapper::success || result.status == Ipopt_Wrapper::acceptable){
//        std::cout << "opt_f = " << result.opt_obj << ", opt_x = " << endl;;
        for (size_t k = 0; k < para.dim1; ++k) { theta1_a(k) = result.opt_x[k]; }
    }

    cout << "****************** checking the gradients ******************" << endl;
    std::vector<double> theta1_a_vec(para.dim1);
    for (size_t k = 0; k < para.dim1; ++k) { theta1_a_vec[k] = theta1_a(k); }

    std::vector<double> grad(para.dim1);

    int GoodState = 1; int LaborIntensive = 1;
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Labor = EstimationAuxiliaryModel_part1_ObjFun(theta1_a_vec,
        RealData_GoodState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double obj_sum_Good_Labor = get<0>(t_Good_Labor);
    std::vector<double> grad_sum_Good_Labor = get<1>(t_Good_Labor);
    int Ncount_Good_Labor = get<2>(t_Good_Labor);;

    GoodState = 0; LaborIntensive = 1;
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Labor = EstimationAuxiliaryModel_part1_ObjFun(theta1_a_vec,
        RealData_BadState_LaborInt,GoodState,LaborIntensive,threadsManagement);
    double obj_sum_Bad_Labor = get<0>(t_Bad_Labor);
    std::vector<double> grad_sum_Bad_Labor = get<1>(t_Bad_Labor);
    int Ncount_Bad_Labor = get<2>(t_Bad_Labor);

//    /* use only labor intensive sector */
//    double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor) / double(Ncount_Good_Labor + Ncount_Bad_Labor);
//    grad[0] = (grad_sum_Good_Labor[0]) / double(Ncount_Good_Labor);
//    grad[1] = (grad_sum_Bad_Labor[0]) / double(Ncount_Bad_Labor);
//    grad[2] = (grad_sum_Good_Labor[1] + grad_sum_Bad_Labor[1]) / double(Ncount_Good_Labor+Ncount_Bad_Labor);
//    grad[3] = (grad_sum_Good_Labor[2] + grad_sum_Bad_Labor[2]) / double(Ncount_Good_Labor+Ncount_Bad_Labor);

    GoodState = 1; LaborIntensive = 0;
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Capital = EstimationAuxiliaryModel_part1_ObjFun(theta1_a_vec,
        RealData_GoodState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
    double obj_sum_Good_Capital = get<0>(t_Good_Capital);
    std::vector<double> grad_sum_Good_Capital = get<1>(t_Good_Capital);
    int Ncount_Good_Capital = get<2>(t_Good_Capital);;

    GoodState = 0; LaborIntensive = 0;
    tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Capital = EstimationAuxiliaryModel_part1_ObjFun(theta1_a_vec,
        RealData_BadState_CapitalInt,GoodState,LaborIntensive,threadsManagement);
    double obj_sum_Bad_Capital = get<0>(t_Bad_Capital);
    std::vector<double> grad_sum_Bad_Capital = get<1>(t_Bad_Capital);
    int Ncount_Bad_Capital = get<2>(t_Bad_Capital);

//    double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital)
//                 / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
//    grad[0] = (grad_sum_Good_Labor[0] + grad_sum_Good_Capital[0])
//              / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
//    grad[1] = (grad_sum_Bad_Labor[0] + grad_sum_Bad_Capital[0])
//              / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
//
//    grad[2] = (grad_sum_Good_Labor[1] + grad_sum_Bad_Labor[1] + grad_sum_Good_Capital[1] + grad_sum_Bad_Capital[1])
//              / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);
//    grad[3] = (grad_sum_Good_Labor[2] + grad_sum_Bad_Labor[2] + grad_sum_Good_Capital[2] + grad_sum_Bad_Capital[2])
//              / double(Ncount_Good_Labor+Ncount_Bad_Labor+Ncount_Good_Capital+Ncount_Bad_Capital);

    double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital);
    grad[0] = (grad_sum_Good_Labor[0] + grad_sum_Good_Capital[0]);
    grad[1] = (grad_sum_Bad_Labor[0] + grad_sum_Bad_Capital[0]);

    grad[2] = (grad_sum_Good_Labor[1] + grad_sum_Bad_Labor[1] + grad_sum_Good_Capital[1] + grad_sum_Bad_Capital[1]);
    grad[3] = (grad_sum_Good_Labor[2] + grad_sum_Bad_Labor[2] + grad_sum_Good_Capital[2] + grad_sum_Bad_Capital[2]);
    cout << "theta1_a = ";
    for (size_t k = 0; k < para.dim1; ++k) { cout << theta1_a(k) << ", "; }
    cout << ";" << endl;

    cout << "grad = ";
    cout << grad[0] << ", "; cout << grad[1] << ", ";
    for (size_t k = 2; k < para.dim1; ++k) { cout << grad[k] << ", "; }
    cout << ";" << endl;

    cout << "x_ub = ";
    for (size_t k = 0; k < para.dim1; ++k) { cout << x_ub[k] << ", "; }
    cout << ";" << endl;

    cout << "x_lb = ";
    for (size_t k = 0; k < para.dim1; ++k) { cout << x_lb[k] << ", "; }
    cout << ";" << endl;

    return theta1_a;
}

//// Compute the objective function
tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd> alias::EstimationAuxiliaryModel_part1_ObjFun(
    const std::vector<double> & theta1,const SimData & RealData, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    double pi = 3.1415926;
    int N_RealData = RealData.good.rows();

    ArrayXi Period_vec = ArrayXi::LinSpaced(para.EstTbar,0,para.EstTbar-1);
    ArrayXXi Period_mat(N_RealData,para.EstTbar);
    Period_mat.rowwise() = Period_vec.transpose();
    ArrayXd Ones = ArrayXd::Ones(N_RealData);

    ParaEst1 para_est = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    ArrayXXd logResR = ArrayXXd::Zero(N_RealData,para.EstTbar);
    logResR = log(RealData.Revenue)
        - para_est.alpha_tilde_L*para_est.alpha_Lc*log(RealData.Employ_uc)
        - para_est.alpha_tilde_L*para_est.alpha_Lr*log(RealData.Employ_ur)
        - para_est.alpha_tilde_K*log(RealData.Capital);

    ArrayXXd lnLtemp = ArrayXXd::Zero(N_RealData,para.EstTbar);
    lnLtemp = para_est.alpha_Lc * log(RealData.Employ_uc) + para_est.alpha_Lr * log(RealData.Employ_ur);

    ArrayXXd logL = ArrayXXd::Zero(N_RealData,para.EstTbar);
    ArrayXXd dlogL_d0 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // gamma0
    ArrayXXd dlogL_d1 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // gamma1
    ArrayXXd dlogL_d2 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // sigma_phi_eps
//    ArrayXXd dlogL_d3 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // alpha_tilde_L
//    ArrayXXd dlogL_d4 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // alpha_Lc
//    ArrayXXd dlogL_d5 = ArrayXXd::Zero(para.SimN,para.EstTbar);  // mu

    ArrayXXi count = ArrayXXi::Zero(N_RealData,para.EstTbar);
    ArrayXi count_firms = ArrayXi::Zero(N_RealData);

    //// Other period
    for (size_t n = 0; n < N_RealData; ++n) {
        for (size_t t = RealData.first_year_in_Data(n) + 1; t < para.EstTbar; ++t) {
//            cout << "n = " << n << "; t = " << t << endl;
            if (RealData.State_P_nonMiss(n, t) == 1 & RealData.Revenue_Miss(n, t) == 0
                & RealData.Employ_urMiss(n, t) == 0 & RealData.Employ_ucMiss(n, t) == 0
                & RealData.CapitalMiss(n, t) == 0) {

                ArrayXi period_vec_t = Period_mat(n, seqN(0, t))
                    * RealData.State_P_nonMiss(n, seqN(0, t))
                    * (1 - RealData.Revenue_Miss(n, seqN(0, t)))
                    * (1 - RealData.Employ_urMiss(n, seqN(0, t)))
                    * (1 - RealData.Employ_ucMiss(n, seqN(0, t)))
                    * (1 - RealData.CapitalMiss(n, seqN(0, t)));

                int max_period = period_vec_t.maxCoeff();
                int diff_period = t - max_period;
                if ( diff_period <= 1 ) { //// use only the previous year productivity
                    if ( RealData.State_P_nonMiss(n, max_period) == 1 & RealData.Revenue_Miss(n, max_period) == 0
                        & RealData.Employ_urMiss(n, max_period) == 0 & RealData.Employ_ucMiss(n, max_period) == 0
                        & RealData.CapitalMiss(n, max_period) == 0) {

                        double mu1 = -para_est.gamma0 * (1.0 - pow(para_est.gamma1, diff_period))
                            / (1.0 - para_est.gamma1);
                        double mu2 = logResR(n, t)
                            - pow(para_est.gamma1, diff_period) * logResR(n, max_period);
                        double mu = mu1 + mu2;

                        double var = pow(para_est.sigma_phi_eps, 2.0) * (1.0 - pow(para_est.gamma1, 2.0 * diff_period))
                                / (1.0 - pow(para_est.gamma1, 2.0));
                        double sigma = sqrt(var);

                        double dmu_dgamma0 = -(1.0 - pow(para_est.gamma1, diff_period)) / (1.0 - para_est.gamma1);
                        double dmu_dgamma1 = -diff_period * pow(para_est.gamma1, diff_period - 1) * logResR(n, max_period)
                            + para_est.gamma0 * diff_period * pow(para_est.gamma1, diff_period - 1) / (1.0 - para_est.gamma1)
                            - para_est.gamma0 * (1 - pow(para_est.gamma1, diff_period)) / pow(1.0 - para_est.gamma1, 2);

                        double dsigma_dgamma1 = -para_est.sigma_phi_eps * diff_period * pow(para_est.gamma1, 2.0 * diff_period - 1)
                            * pow(1 - pow(para_est.gamma1, 2 * diff_period), -0.5)
                            * pow(1 - pow(para_est.gamma1, 2), -0.5)
                            + para_est.gamma1 * para_est.sigma_phi_eps
                            * pow(1 - pow(para_est.gamma1, 2 * diff_period), 0.5)
                            / pow(1 - pow(para_est.gamma1, 2.0), 1.5);
                        double dsigma_dsigmaphi = sqrt( (1 - pow(para_est.gamma1, 2.0 * diff_period)) / (1 - pow(para_est.gamma1, 2.0)) );

                        logL(n, t) = -0.5 * log(2.0 * pi) - log(sigma) - 0.5 * pow(mu / sigma, 2);

                        dlogL_d0(n, t) = -mu / pow(sigma, 2) * dmu_dgamma0; // gamma0
                        //                    cout << "dlogL_d0(n,t) = " << dlogL_d0(n,t) << endl;
                        dlogL_d1(n, t) = -1.0 / sigma * dsigma_dgamma1 - mu / pow(sigma, 2) * dmu_dgamma1
                            + pow(mu, 2) / pow(sigma, 3) * dsigma_dgamma1; // gamma1
                        //                    cout << "dlogL_d1(n,t) = " << dlogL_d1(n,t) << endl;
                        dlogL_d2(n, t) = -1.0 / sigma * dsigma_dsigmaphi
                            + pow(mu, 2) / pow(sigma, 3) * dsigma_dsigmaphi; // sigma_phi_eps
                        //                    cout << "dlogL_d2(n,t) = " << dlogL_d2(n,t) << endl;

                        double dy_dalpha_tilde_L_1 = mu / pow(sigma, 2) * (lnLtemp(n, t)
                            - pow(para_est.gamma1, diff_period) * lnLtemp(n,max_period));// alpha_tilde_L
                        double dy_dalpha_tilde_L_2 = mu / pow(sigma, 2) * (log(RealData.Capital(n, t))
                            - pow(para_est.gamma1, diff_period) * log(RealData.Capital(n,max_period)));// alpha_tilde_K

//                        dlogL_d3(n, t) = -dy_dalpha_tilde_L_1 + dy_dalpha_tilde_L_2;
//                        //                    cout << "dlogL_d3(n,t) = " << dlogL_d3(n,t) << endl;

                        double temp_alpha_Lc = mu / pow(sigma, 2) * para_est.alpha_tilde_L
                            * (log(RealData.Employ_uc(n, t)) - pow(para_est.gamma1, diff_period)
                            * log(RealData.Employ_uc(n,max_period))); // alpha_uc

                        double temp_alpha_Lr = mu / pow(sigma, 2) * para_est.alpha_tilde_L
                            * (log(RealData.Employ_ur(n, t)) - pow(para_est.gamma1, diff_period)
                            * log(RealData.Employ_ur(n,max_period))); // alpha_ur

//                        dlogL_d4(n, t) = -temp_alpha_Lc + temp_alpha_Lr;
//                        //                    cout << "dlogL_d4(n,t) = " << dlogL_d4(n,t) << endl;

                        if (isfinite(logL(n, t)) == 0) {
                            cout << "n = " << n << "; t = " << t << "; para_est.sigma_phi_eps = "
                                 << para_est.sigma_phi_eps
                                 << "; para_est.gamma1 = " << para_est.gamma1 << "; sigma = " << sigma << endl;
                            cout << "period_vec_t = " << period_vec_t.transpose() << endl;
                            cout << "RealData.State_P_nonMiss(n, seqN(0, t)) = "
                                 << RealData.State_P_nonMiss(n, seqN(0, t)) << endl;
                            cout << "RealData.Revenue_Miss(n, seqN(0, t)) = "
                                 << RealData.Revenue_Miss(n, seqN(0, t))
                                 << endl;
                            cout << "RealData.Employ_urMiss(n, seqN(0, t)) = "
                                 << RealData.Employ_urMiss(n, seqN(0, t))
                                 << endl;
                            cout << "RealData.Employ_ucMiss(n, seqN(0, t)) = "
                                 << RealData.Employ_ucMiss(n, seqN(0, t))
                                 << endl;
                            cout << "RealData.CapitalMiss(n, seqN(0, t)) = " << RealData.CapitalMiss(n, seqN(0, t))
                                 << endl;
                            cout << "diff_period = " << diff_period << "; max_period = " << max_period
                                 << "; RealData.first_year_in_Data(n) = " << RealData.first_year_in_Data(n)
                                 << "; logL(n, t) = " << logL(n, t) << endl;
                            cout << "; log(sim_data.Revenue) = " << RealData.Revenue(n, max_period)
                                 << "; log(sim_data.Employ_uc) = " << RealData.Employ_uc(n, max_period)
                                 << "; log(sim_data.Employ_ur) = " << RealData.Employ_ur(n, max_period) << endl;
                            cout << "; sim_data.Capital = " << RealData.Capital(n, max_period) << endl;
                            cout << "; sim_data.State_P_nonMiss(n,max_period) = "
                                 << RealData.State_P_nonMiss(n, max_period) << endl;
                            throw runtime_error("100");
                        }
                        count(n, t) = 1;
                        count_firms(n) = 1;
                    }
                }
            }
        }
    }
    int Ncount = count.sum();
    int Ncount_firms = count_firms.sum();
//    logL = (logL.isFinite()).select(logL, 0.0);
    double y_sum = -logL.sum(); // value of objective function
//    cout << "y_sum = " << y_sum << endl;

    std::vector<double> grad_sum(para.dim1);
    dlogL_d0 = (dlogL_d0.isFinite()).select(dlogL_d0, 0.0); // gamma0
    grad_sum[0] = -dlogL_d0.sum();
//    cout << "dlogL_d0 = " << grad_sum[0] << endl;

    dlogL_d1 = (dlogL_d1.isFinite()).select(dlogL_d1, 0.0); // gamma1
    grad_sum[1] = -dlogL_d1.sum();
//    cout << "dlogL_d1 = " << grad_sum[1] << endl;

    dlogL_d2 = (dlogL_d2.isFinite()).select(dlogL_d2, 0.0); // sigma_phi_eps
    grad_sum[2] = -dlogL_d2.sum();
//    cout << "dlogL_d2 = " << grad_sum[2] << endl;

    return tuple<double,std::vector<double>,int,int,ArrayXXd,ArrayXXd,ArrayXXd>(y_sum,grad_sum,Ncount,Ncount_firms,
        dlogL_d0,dlogL_d1,dlogL_d2);
}

//// Backout phi based on theta1_a_ini
ArrayXXd alias::Backout_phi_a(const SimData & RealData,const double & alpha_tilde_K,const double & alpha_tilde_L,
                              const double & alpha_Lr,const double & alpha_Lc) {

//    ArrayXXd L = pow( para_est.alpha_Lc*pow(sim_data.Employ_uc,(para_est.mu-1)/para_est.mu)
//        + para_est.alpha_Lr*pow(sim_data.Employ_ur,(para_est.mu-1)/para_est.mu), para_est.mu/(para_est.mu-1));

    ArrayXXd logL = alpha_Lc * log(RealData.Employ_uc) + alpha_Lr * log(RealData.Employ_ur);
    ArrayXXd lnphi_a = log(RealData.Revenue) - alpha_tilde_K*log(RealData.Capital)-alpha_tilde_L*logL;

    return lnphi_a;
}

/********************************************************************************************************************
* Auxilliary Model Estimation: Part 2
********************************************************************************************************************/
//// Given the objective function, find the estimates by maximize the objective function
ArrayXd alias::EstimationAuxiliaryModel_part2_theta2(const SimData & RealData_GoodState_LaborInt,
    const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
    const SimData & RealData_BadState_CapitalInt, const ArrayXd & theta1_a, const ArrayXd & theta2_a_ini,
    MultiThreads::Threads_Management & threadsManagement) {

    std::vector<double> theta1(para.dim1);
    for (size_t k = 0; k < para.dim1; ++k) { theta1[k] = theta1_a(k); }
    int GoodState = 1; int LaborIntensive = 1;
    ParaEst1 para_est1_Good_Labor = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    GoodState = 0; LaborIntensive = 1;
    ParaEst1 para_est1_Bad_Labor = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    GoodState = 1; LaborIntensive = 0;
    ParaEst1 para_est1_Good_Capital = constructParaEst_Part1(theta1,GoodState,LaborIntensive);
    GoodState = 0; LaborIntensive = 0;
    ParaEst1 para_est1_Bad_Capital = constructParaEst_Part1(theta1,GoodState,LaborIntensive);

    //// estimated productivity
    ArrayXXd lnphi_a_GoodState_Labor = Backout_phi_a(RealData_GoodState_LaborInt,
        para_est1_Good_Labor.alpha_tilde_K,para_est1_Good_Labor.alpha_tilde_L,
        para_est1_Good_Labor.alpha_Lr,para_est1_Good_Labor.alpha_Lc);
    ArrayXXd lnphi_a_BadState_Labor = Backout_phi_a(RealData_BadState_LaborInt,
        para_est1_Bad_Labor.alpha_tilde_K,para_est1_Bad_Labor.alpha_tilde_L,
        para_est1_Bad_Labor.alpha_Lr,para_est1_Bad_Labor.alpha_Lc);
    ArrayXXd lnphi_a_GoodState_Capital = Backout_phi_a(RealData_GoodState_CapitalInt,
        para_est1_Good_Capital.alpha_tilde_K,para_est1_Good_Capital.alpha_tilde_L,
        para_est1_Good_Capital.alpha_Lr,para_est1_Good_Capital.alpha_Lc);
    ArrayXXd lnphi_a_BadState_Capital = Backout_phi_a(RealData_BadState_CapitalInt,
        para_est1_Bad_Capital.alpha_tilde_K,para_est1_Bad_Capital.alpha_tilde_L,
        para_est1_Bad_Capital.alpha_Lr,para_est1_Bad_Capital.alpha_Lc);

    ArrayXd theta2_a(para.dim2);
    std::vector<double> x0_vec(para.dim2);
    for (size_t k = 0; k < para.dim2; ++k) { x0_vec[k] = theta2_a_ini(k); }
    auto f_without_hessian= [&RealData_GoodState_LaborInt,&RealData_BadState_LaborInt,
            &RealData_GoodState_CapitalInt,&RealData_BadState_CapitalInt,&lnphi_a_GoodState_Labor,&lnphi_a_BadState_Labor,
            &lnphi_a_GoodState_Capital,&lnphi_a_BadState_Capital,&para_est1_Good_Labor,&para_est1_Bad_Labor,
            &para_est1_Good_Capital,&para_est1_Bad_Capital] (std::vector<double> x, std::vector<double> & grad){

        tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Labor = EstimationAuxiliaryModel_part2_ObjFun(x,
            para_est1_Good_Labor,RealData_GoodState_LaborInt,lnphi_a_GoodState_Labor);
        double obj_sum_Good_Labor = get<0>(t_Good_Labor);
        std::vector<double> grad_sum_Good_Labor= get<1>(t_Good_Labor);
        int Ncount_Good_Labor = get<2>(t_Good_Labor);

        tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Labor = EstimationAuxiliaryModel_part2_ObjFun(x,
            para_est1_Bad_Labor,RealData_BadState_LaborInt,lnphi_a_BadState_Labor);
        double obj_sum_Bad_Labor = get<0>(t_Bad_Labor);
        std::vector<double> grad_sum_Bad_Labor = get<1>(t_Bad_Labor);
        int Ncount_Bad_Labor = get<2>(t_Bad_Labor);

//        /* use only labor intensive sector */
//        double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor) / double(Ncount_Good_Labor + Ncount_Bad_Labor);
//        for (size_t k = 0; k < para.dim2; ++k) {
//            grad[k] = (grad_sum_Good_Labor[k] + grad_sum_Bad_Labor[k]) / double(Ncount_Good_Labor + Ncount_Bad_Labor);
//        }

        tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Capital = EstimationAuxiliaryModel_part2_ObjFun(x,
            para_est1_Good_Capital,RealData_GoodState_CapitalInt,lnphi_a_GoodState_Capital);
        double obj_sum_Good_Capital = get<0>(t_Good_Capital);
        std::vector<double> grad_sum_Good_Capital = get<1>(t_Good_Capital);
        int Ncount_Good_Capital = get<2>(t_Good_Capital);

        tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Capital = EstimationAuxiliaryModel_part2_ObjFun(x,
            para_est1_Bad_Capital,RealData_BadState_CapitalInt,lnphi_a_BadState_Capital);
        double obj_sum_Bad_Capital = get<0>(t_Bad_Capital);
        std::vector<double> grad_sum_Bad_Capital = get<1>(t_Bad_Capital);
        int Ncount_Bad_Capital = get<2>(t_Bad_Capital);

        double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital)
                     / double(Ncount_Good_Labor + Ncount_Bad_Labor + Ncount_Good_Capital + Ncount_Bad_Capital);
        for (size_t k = 0; k < para.dim2; ++k) {
            grad[k] = (grad_sum_Good_Labor[k] + grad_sum_Bad_Labor[k] + grad_sum_Good_Capital[k] + grad_sum_Bad_Capital[k])
                      / double(Ncount_Good_Labor + Ncount_Bad_Labor + Ncount_Good_Capital + Ncount_Bad_Capital);
        }
        return obj;
    };

    // lower bound
    std::vector<double> x_lb(para.dim2);
    x_lb[0] = -20; x_lb[1] = -20; x_lb[2] = 0.1;

    // upper bound
    std::vector<double> x_ub(para.dim2);
    x_ub[0] = 100; x_ub[1] = 100; x_ub[2] = 100;

    Simple_Opt_Problem problem(f_without_hessian, para.dim2);
    problem.set_x_lb(x_lb);
    problem.set_x_ub(x_ub);

    Option option;
    option.display_level = 5;
    option.test_derivative = 0;
    Ipopt_Solver solver(option);

    auto result = solver.optimize(problem, x0_vec);
    if (result.status == Ipopt_Wrapper::success || result.status == Ipopt_Wrapper::acceptable){
//        std::cout << "opt_f = " << result.opt_obj << ", opt_x = " << endl;;
        theta2_a(0) = result.opt_x[0]; theta2_a(1) = result.opt_x[1]; theta2_a(2) = result.opt_x[2];
    }

    cout << "****************** checking the gradients ******************" << endl;
    std::vector<double> theta2_a_vec(para.dim2);
    for (size_t k = 0; k < para.dim2; ++k) {theta2_a_vec[k] = theta2_a(k);}
    std::vector<double> grad(para.dim2);

    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Labor = EstimationAuxiliaryModel_part2_ObjFun(
        theta2_a_vec,para_est1_Good_Labor,RealData_GoodState_LaborInt,lnphi_a_GoodState_Labor);
    double obj_sum_Good_Labor = get<0>(t_Good_Labor);
    std::vector<double> grad_sum_Good_Labor= get<1>(t_Good_Labor);
    int Ncount_Good_Labor = get<2>(t_Good_Labor);

    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Labor = EstimationAuxiliaryModel_part2_ObjFun(
        theta2_a_vec,para_est1_Bad_Labor,RealData_BadState_LaborInt,lnphi_a_BadState_Labor);
    double obj_sum_Bad_Labor = get<0>(t_Bad_Labor);
    std::vector<double> grad_sum_Bad_Labor = get<1>(t_Bad_Labor);
    int Ncount_Bad_Labor = get<2>(t_Bad_Labor);

//    /* use only labor intensive sector */
//    double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor) / double(Ncount_Good_Labor + Ncount_Bad_Labor);
//    for (size_t k = 0; k < para.dim2; ++k) {
//        grad[k] = (grad_sum_Good_Labor[k] + grad_sum_Bad_Labor[k]) / double(Ncount_Good_Labor + Ncount_Bad_Labor);
//    }

    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Good_Capital = EstimationAuxiliaryModel_part2_ObjFun(
        theta2_a_vec,para_est1_Good_Capital,RealData_GoodState_CapitalInt,lnphi_a_GoodState_Capital);
    double obj_sum_Good_Capital = get<0>(t_Good_Capital);
    std::vector<double> grad_sum_Good_Capital = get<1>(t_Good_Capital);
    int Ncount_Good_Capital = get<2>(t_Good_Capital);

    tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> t_Bad_Capital = EstimationAuxiliaryModel_part2_ObjFun(
        theta2_a_vec,para_est1_Bad_Capital,RealData_BadState_CapitalInt,lnphi_a_BadState_Capital);
    double obj_sum_Bad_Capital = get<0>(t_Bad_Capital);
    std::vector<double> grad_sum_Bad_Capital = get<1>(t_Bad_Capital);
    int Ncount_Bad_Capital = get<2>(t_Bad_Capital);

    double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital)
                 / double(Ncount_Good_Labor + Ncount_Bad_Labor + Ncount_Good_Capital + Ncount_Bad_Capital);
    for (size_t k = 0; k < para.dim2; ++k) {
        grad[k] = (grad_sum_Good_Labor[k] + grad_sum_Bad_Labor[k] + grad_sum_Good_Capital[k] + grad_sum_Bad_Capital[k])
                  / double(Ncount_Good_Labor + Ncount_Bad_Labor + Ncount_Good_Capital + Ncount_Bad_Capital);
    }

    cout << "theta2_a_vec = " << theta2_a_vec[0] << "; " << theta2_a_vec[1] << "; " << theta2_a_vec[2] << endl;
    cout << "grad = " << grad[0] << ", " << grad[1] << ", " << grad[2] << endl;
    cout << "x_ub = " << x_ub[0] << ", " << x_ub[1] << ", " << x_ub[2] << endl;
    cout << "x_lb = " << x_lb[0] << ", " << x_lb[1] << ", " << x_lb[2] << endl;

    return theta2_a;
}

//// Compute the objective function
tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd> alias::EstimationAuxiliaryModel_part2_ObjFun(
    const std::vector<double> & theta2, const ParaEst1 & para_est1, const SimData & RealData,const ArrayXXd & lnphi_a) {

    int N_RealData = RealData.good.rows();

    ArrayXi Period_vec = ArrayXi::LinSpaced(para.EstTbar,0,para.EstTbar-1);
    ArrayXXi Period_mat(N_RealData,para.EstTbar);
    Period_mat.rowwise() = Period_vec.transpose();

    ArrayXd Ones = ArrayXd::Ones(N_RealData);
    ArrayXXd ValueAdd = RealData.Revenue;
    ParaEst2 para_est2 = constructParaEst_Part2(theta2);

    ArrayXXd logL = ArrayXXd::Zero(N_RealData,para.EstTbar);
    ArrayXXd dlogL_d0 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // F_P
    ArrayXXd dlogL_d1 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // F_D
    ArrayXXd dlogL_d2 = ArrayXXd::Zero(N_RealData,para.EstTbar);  // sigma_PD

    ArrayXXi count = ArrayXXi::Zero(N_RealData,para.EstTbar);
    ArrayXi countP_firms = ArrayXi::Zero(N_RealData);
    ArrayXi countD_firms = ArrayXi::Zero(N_RealData);
    ArrayXi count_firms = ArrayXi::Zero(N_RealData);
    //// Other period
    for (size_t n = 0; n < N_RealData; ++n) {
        for (size_t t = RealData.first_year_in_Data(n)+1; t < para.EstTbar; ++t) {
            ArrayXi period_vec_t = Period_mat(n, seqN(0, t+1))
                * RealData.State_P_nonMiss(n, seqN(0, t+1))
                * (1 - RealData.Revenue_Miss(n, seqN(0, t+1)))
                * (1 - RealData.Employ_urMiss(n, seqN(0, t+1)))
                * (1 - RealData.Employ_ucMiss(n, seqN(0, t+1)))
                * (1 - RealData.CapitalMiss(n, seqN(0, t+1)));

            int max_period = period_vec_t.maxCoeff();
            int t_lag = max_period;

            if ( RealData.State_P_nonMiss(n, t-1) == 1 ) {

                double VV_PD;
                double VV_PD_temp = (log(ValueAdd(n, t_lag)) - para_est2.F_P_PP) / para_est2.sigma_PD;
                VV_PD = min(VV_PD_temp, 7.5);
                VV_PD = max(VV_PD, -7.5);

                if (RealData.State_P_nonMiss(n, t) == 1) {
                    //// choose Production
                    double lnProb_P = log(normCDF(VV_PD, 0, 1));
                    double dlnProb_P_d0 = normPDF(VV_PD, 0, 1) / normCDF(VV_PD, 0, 1)
                                          * ( -1.0 / para_est2.sigma_PD );
                    double dlnProb_P_d2 = normPDF(VV_PD, 0, 1) / normCDF(VV_PD, 0, 1)
                                          * VV_PD * (-1.0 / para_est2.sigma_PD);

                    logL(n, t) = lnProb_P;
                    dlogL_d0(n, t) = dlnProb_P_d0;
                    dlogL_d1(n, t) = 0.0;
                    dlogL_d2(n, t) = dlnProb_P_d2;

                    count(n, t) = 1;
                    countP_firms(n) = 1;
                    count_firms(n) = 1;
                }
                else if (RealData.State_D_nonMiss(n, t) == 1) {
                    //// choose Dormancy
                    double lnProb_D = log(normCDF(-VV_PD, 0, 1));
                    double dlnProb_D_d0 = normPDF(-VV_PD, 0, 1) / normCDF(-VV_PD, 0, 1)
                        * ( 1.0 / para_est2.sigma_PD);
                    double dlnProb_D_d2 = normPDF(-VV_PD, 0, 1) / normCDF(-VV_PD, 0, 1)
                        * (-VV_PD) * (-1.0 / para_est2.sigma_PD);

                    logL(n, t) = lnProb_D;
                    dlogL_d0(n, t) = dlnProb_D_d0;
                    dlogL_d1(n, t) = 0.0;
                    dlogL_d2(n, t) = dlnProb_D_d2;
                    count(n, t) = 1;
                    countP_firms(n) = 1;
                    count_firms(n) = 1;
                }
            }
            else if ( RealData.State_D_nonMiss(n, t-1) == 1 ) {

                double VV_PD;
                double VV_PD_temp = (log(ValueAdd(n, t_lag)) - para_est2.F_P_DP) / para_est2.sigma_PD;
                VV_PD = min(VV_PD_temp, 7.5);
                VV_PD = max(VV_PD, -7.5);

                if (RealData.State_P_nonMiss(n, t) == 1) {
                    //// choose Production
                    double lnProb_P = log(normCDF(VV_PD, 0, 1));
                    double dlnProb_P_d1 = normPDF(VV_PD, 0, 1) / normCDF(VV_PD, 0, 1)
                                          * (-1.0 / para_est2.sigma_PD);
                    double dlnProb_P_d2 = normPDF(VV_PD, 0, 1) / normCDF(VV_PD, 0, 1)
                                          * VV_PD * (-1.0 / para_est2.sigma_PD);

                    logL(n, t) = lnProb_P;
                    dlogL_d0(n, t) = 0.0;
                    dlogL_d1(n, t) = dlnProb_P_d1;
                    dlogL_d2(n, t) = dlnProb_P_d2;
                    count(n, t) = 1;
                    countD_firms(n) = 1;
                    count_firms(n) = 1;
                }
                else if (RealData.State_D_nonMiss(n, t) == 1) {
                    //// choose Dormancy
                    double lnProb_D = log(normCDF(-VV_PD, 0, 1));
                    double dlnProb_D_d1 = normPDF(-VV_PD, 0, 1) / normCDF(-VV_PD, 0, 1)
                                          * (1.0 / para_est2.sigma_PD);
                    double dlnProb_D_d2 = normPDF(-VV_PD, 0, 1) / normCDF(-VV_PD, 0, 1)
                                          * (-VV_PD) * (-1.0 / para_est2.sigma_PD);

                    logL(n, t) = lnProb_D;
                    dlogL_d0(n, t) = 0.0;
                    dlogL_d1(n, t) = dlnProb_D_d1;
                    dlogL_d2(n, t) = dlnProb_D_d2;
                    count(n, t) = 1;
                    countD_firms(n) = 1;
                    count_firms(n) = 1;
                }
            }
        }
    }

    int Ncount = count.sum();
    int NcountP_firms = countP_firms.sum();
    int NcountD_firms = countD_firms.sum();
    int Ncount_firms = count_firms.sum();

//    logL = (logL.isFinite()).select(logL, 0.0);
    double y_sum = -logL.sum();

    std::vector<double> grad_sum(para.dim2);
//    dlogL_d0 = (dlogL_d0.isFinite()).select(dlogL_d0, 0.0);
    grad_sum[0] = -dlogL_d0.sum();

//    dlogL_d1 = (dlogL_d1.isFinite()).select(dlogL_d1, 0.0);
    grad_sum[1] = -dlogL_d1.sum();

    grad_sum[2] = -dlogL_d2.sum();

    return tuple<double,std::vector<double>,int,int,int,int,ArrayXXd,ArrayXXd,ArrayXXd>(y_sum,grad_sum,
        Ncount,NcountP_firms,NcountD_firms,Ncount_firms,dlogL_d0,dlogL_d1,dlogL_d2);
}

/********************************************************************************************************************
* Auxilliary Model Estimation: Part 3
********************************************************************************************************************/
//// Given the objective function, find the estimates by maximize the objective function
tuple<ArrayXd,double> alias::EstimationAuxiliaryModel_part3_theta3(const SimData & RealData_GoodState_LaborInt,
    const SimData & RealData_BadState_LaborInt, const SimData & RealData_GoodState_CapitalInt,
    const SimData & RealData_BadState_CapitalInt,
    const ArrayXd & theta1_a, const ArrayXd & theta2_a, const ArrayXd & theta3_a_ini,
    MultiThreads::Threads_Management & threadsManagement) {

//////////////
////    cout << "theta3_a_ini = " << theta3_a_ini.transpose() << endl;
////    ArrayXd temp_vec = ArrayXd::LinSpaced(11,0,10);
////    for (size_t n = 0; n<11; ++n) {
////        std::vector<double> theta3_a_temp(para.dim3); // initial values
////        for (size_t k = 0; k < para.dim3; ++k) { theta3_a_temp[k] = theta3_a_ini(k); }
////        theta3_a_temp[0] = theta3_a_ini(0) + temp_vec(n);
////        cout << " theta3_a_temp = " <<  theta3_a_temp[0] << endl;
////
//////        cout << "theta3_a_temp = " << theta3_a_temp[0] << "; " << theta3_a_temp[1] << "; " << theta3_a_temp[2] << "; "
//////                << theta3_a_temp[3] << "; " << theta3_a_temp[4] << "; " << theta3_a_temp[5] << "; "
//////                << theta3_a_temp[6] << "; " << theta3_a_temp[7] << "; " << theta3_a_temp[8] << "; "
//////                << theta3_a_temp[9] << "; " << theta3_a_temp[10] << "; " << theta3_a_temp[11] << "; "
//////                << theta3_a_temp[12] << "; " << theta3_a_temp[13] << "; " << theta3_a_temp[14] << endl;
////
////        int GoodState = 1; int LaborIntensive = 1;
////        tuple<double,int,EquStateV> t_obj_Good_Labor = EstimationAuxiliaryModel_part3_ObjFun(theta3_a_temp,
////            RealData_GoodState_LaborInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
////        double obj_sum_Good_Labor = get<0>(t_obj_Good_Labor);
////        int Ncount_Good_Labor = get<1>(t_obj_Good_Labor);
////        EquStateV EquV_Good_Labor = get<2>(t_obj_Good_Labor);
////
////        GoodState = 1; LaborIntensive = 0;
////        tuple<double,int,EquStateV> t_obj_Good_Capital = EstimationAuxiliaryModel_part3_ObjFun(theta3_a_temp,
////             RealData_GoodState_CapitalInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
////        double obj_sum_Good_Capital = get<0>(t_obj_Good_Capital);
////        int Ncount_Good_Capital = get<1>(t_obj_Good_Capital);
////        EquStateV EquV_Good_Capital = get<2>(t_obj_Good_Capital);
////
////        GoodState = 0; LaborIntensive = 1;
////        tuple<double,int,EquStateV> t_obj_Bad_Labor = EstimationAuxiliaryModel_part3_ObjFun(theta3_a_temp,
////            RealData_BadState_LaborInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
////        double obj_sum_Bad_Labor = get<0>(t_obj_Bad_Labor);
////        int Ncount_Bad_Labor = get<1>(t_obj_Bad_Labor);
////        EquStateV EquV_Bad_Labor = get<2>(t_obj_Bad_Labor);
////
////        GoodState = 0; LaborIntensive = 0;
////        tuple<double,int,EquStateV> t_obj_Bad_Capital = EstimationAuxiliaryModel_part3_ObjFun(theta3_a_temp,
////            RealData_BadState_CapitalInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
////        double obj_sum_Bad_Capital = get<0>(t_obj_Bad_Capital);
////        int Ncount_Bad_Capital = get<1>(t_obj_Bad_Capital);
////        EquStateV EquV_Bad_Capital = get<2>(t_obj_Bad_Capital);
//////        throw runtime_error("811");
////
////        cout << setprecision(10) << "obj_sum_Good_Labor = " << obj_sum_Good_Labor << "; obj_sum_Good_Capital = " << obj_sum_Good_Capital
////            << "; obj_sum_Bad_Labor = " << obj_sum_Bad_Labor << "; obj_sum_Bad_Capital = " << obj_sum_Bad_Capital
////            << "; obj_sum_Good_Labor + obj_sum_Good_Capital + obj_sum_Bad_Labor + obj_sum_Bad_Capital = "
////            << obj_sum_Good_Labor + obj_sum_Good_Capital + obj_sum_Bad_Labor + obj_sum_Bad_Capital << endl;
////    }
////    throw runtime_error("800");
//
    //// objective function
    auto f_derivative_free
        = [&RealData_GoodState_LaborInt,&RealData_BadState_LaborInt,&RealData_GoodState_CapitalInt,
           &RealData_BadState_CapitalInt,&theta1_a,&theta2_a,&threadsManagement](std::vector<double> x){

            int GoodState = 1; int LaborIntensive = 1;
            tuple<double,int,EquStateV> t_obj_Good_Labor = EstimationAuxiliaryModel_part3_ObjFun(x,
                RealData_GoodState_LaborInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
            double obj_sum_Good_Labor = get<0>(t_obj_Good_Labor);
            int Ncount_Good_Labor = get<1>(t_obj_Good_Labor);
            EquStateV EquV_Good_Labor = get<2>(t_obj_Good_Labor);

            GoodState = 0; LaborIntensive = 1;
            tuple<double,int,EquStateV> t_obj_Bad_Labor = EstimationAuxiliaryModel_part3_ObjFun(x,
                RealData_BadState_LaborInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
            double obj_sum_Bad_Labor = get<0>(t_obj_Bad_Labor);
            int Ncount_Bad_Labor = get<1>(t_obj_Bad_Labor);
            EquStateV EquV_Bad_Labor = get<2>(t_obj_Bad_Labor);
//
            /* use both labor and capital intensive sector */
            GoodState = 1; LaborIntensive = 0;
            tuple<double,int,EquStateV> t_obj_Good_Capital = EstimationAuxiliaryModel_part3_ObjFun(x,
                RealData_GoodState_CapitalInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
            double obj_sum_Good_Capital = get<0>(t_obj_Good_Capital);
            int Ncount_Good_Capital = get<1>(t_obj_Good_Capital);
            EquStateV EquV_Good_Capital = get<2>(t_obj_Good_Capital);
////////
            GoodState = 0; LaborIntensive = 0;
            tuple<double,int,EquStateV> t_obj_Bad_Capital = EstimationAuxiliaryModel_part3_ObjFun(x,
                RealData_BadState_CapitalInt,theta1_a,theta2_a,GoodState,LaborIntensive,threadsManagement);
            double obj_sum_Bad_Capital = get<0>(t_obj_Bad_Capital);
            int Ncount_Bad_Capital = get<1>(t_obj_Bad_Capital);
            EquStateV EquV_Bad_Capital = get<2>(t_obj_Bad_Capital);
////
//            /* use only labor intensive sector */
//            double obj = (obj_sum_Good_Labor + obj_sum_Good_Capital) / double(Ncount_Good_Labor + Ncount_Good_Capital);
////                double obj = (obj_sum_Good_Labor + obj_sum_Good_Capital);
//            {
//                std::scoped_lock<std::mutex> lock(threadsManagement.print_mutex);
//                cout << "para_est = ";
//                for (size_t k = 0; k < para.dim3; ++k) { cout << x[k] << ", "; }
//                cout << ";" << endl;
//                cout << "******* Good State and Labor Intensive Sector " << endl;
//                PrintResultEstimationAuxiliaryModel(x,EquV_Good_Labor);
//                cout << "******* Bad State and Labor Intensive Sector " << endl;
//                PrintResultEstimationAuxiliaryModel(x,EquV_Good_Capital);
//                cout << "obj_sum_Good_Labor + obj_sum_Good_Capital = " << (obj_sum_Good_Labor + obj_sum_Good_Capital)
//                    << "; obj = " << setprecision(10) << -obj << endl;
//            }
//
            {
                std::scoped_lock<std::mutex> lock(threadsManagement.print_mutex);
                cout << "para_est = ";
                for (size_t k = 0; k < para.dim3; ++k) { cout << x[k] << ", "; }
                cout << ";" << endl;
                cout << "******* Good State and Labor Intensive Sector " << endl;
                PrintResultEstimationAuxiliaryModel(x, EquV_Good_Labor);
                cout << "******* Bad State and Labor Intensive Sector " << endl;
                PrintResultEstimationAuxiliaryModel(x, EquV_Bad_Labor);
                cout << "******* Good State and Capital Intensive Sector " << endl;
                PrintResultEstimationAuxiliaryModel(x, EquV_Good_Capital);
                cout << "******* Bad State and Capital Intensive Sector " << endl;
                PrintResultEstimationAuxiliaryModel(x, EquV_Bad_Capital);
            }
            double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital)
                / (Ncount_Good_Labor + Ncount_Bad_Labor + Ncount_Good_Capital + Ncount_Bad_Capital);
//////            double obj = (obj_sum_Good_Labor + obj_sum_Bad_Labor + obj_sum_Good_Capital + obj_sum_Bad_Capital);
            cout << "-obj = " << -obj << endl;
            return -obj;
        };

    std::vector<double> x0_vec(para.dim3); // initial values
    for (size_t k = 0; k < para.dim3; ++k) { x0_vec[k] = theta3_a_ini(k); }

    //// H_K(0);c_HK(1);F_K(2);c_FK(3);    H_ur(4);c_H_ur(5);F_L_G_ur(6);F_H_G_ur(7);c_F_ur(8);
    //// H_uc(9);c_H_uc(10);F_G_uc(11);c_F_uc(12);     sigma_Kerror(13);sigma_Lerror_ur(14);sigma_Lerror_uc(15);
    //// F_G_E(16);simga_SE(17)
    std::vector<double> lb(para.dim3);
    int nn = 1;
    lb[0] = 0; lb[1] = 0; lb[2] = 0;
    for (size_t k = nn+0; k < nn+para.dim3-6; ++k) { lb[k] = 0.0; }
    lb[para.dim3-6] = 0.01; lb[para.dim3-5] = 0.01; lb[para.dim3-4] = 0.01;
    lb[para.dim3-3] = -5000; lb[para.dim3-2] = -5000;
    lb[para.dim3-1] = 0.1;

    std::vector<double> ub(para.dim3);
    ub[0] = 5000; ub[1] = 5000; ub[2] = 5000;
    for (size_t k = nn+0; k < nn+para.dim3-6; ++k) { ub[k] = 5000; }
    ub[para.dim3-6] = 2; ub[para.dim3-5] = 2; ub[para.dim3-4] = 2;
    ub[para.dim3-3] = 5000; ub[para.dim3-2] = 5000;
    ub[para.dim3-1] = 5000;

//
////    std::vector<double> lb(para.dim3);
////    int nn = 1;
////    lb[0] = 0.0;
////    for (size_t k = nn+0; k < nn+para.dim3-5; ++k) { lb[k] = 0; }
////    lb[para.dim3-5] = 0.01; lb[para.dim3-4] = 0.01; lb[para.dim3-3] = 0.01;
////    lb[para.dim3-2] = -5000;
////    lb[para.dim3-1] = 0.1;
//////
////    std::vector<double> ub(para.dim3);
////    ub[0] = 5000.0;
////    for (size_t k = nn+0; k < nn+para.dim3-5; ++k) { ub[k] = 5000; }
////    ub[para.dim3-5] = 2; ub[para.dim3-4] = 2; ub[para.dim3-3] = 2;
////    ub[para.dim3-2] = 5000;
////    ub[para.dim3-1] = 5000;
////    cout << "lb[0] = " << lb[0] << "," << lb[1] << ","<< lb[2] << ","<< lb[3] << ","<< lb[4] << ","<< lb[5] << ","
////        << lb[6] << ","<< lb[7] << ","<< lb[8] << "," << lb[9] << ","<< lb[10] << ","<< lb[11] << ","<< lb[12] << endl;
////    cout << "ub[0] = " << ub[0] << ","<< ub[1] << ","<< ub[2] << ","<< ub[3] << ","<< ub[4] << ","<< ub[5] << ","
////        << ub[6] << ","<< ub[7] << ","<< ub[8] << ","<< ub[9] << ub[10] << ","<< ub[11] << ","<< ub[12] << endl;
////    cout << "918" << endl;

    DFS::Optimization_Problem problem(para.dim3, f_derivative_free);
    problem.set_x_lb(lb);
    problem.set_x_ub(ub);

    DFS::Direct_Search_Option option;
    option.display_level = 2;

//    auto result = DFS::Coordinate_Direct_Parallel_Search(problem,x0_vec,option,threadsManagement);
    auto result = DFS::Coordinate_Direct_Search(problem, x0_vec, option);
    ArrayXd theta3_a(para.dim3);
    for (size_t k = 0; k < para.dim3; ++k) { theta3_a(k) = result.optimal_x[k]; }
    double obj3_a = result.optimal_f;

    cout << "optimal_x = " ;
    for (size_t k = 0; k < para.dim3; ++k) { cout << result.optimal_x[k] << "; ";}
    cout << ";" << endl;
    std::cout << "optimal f = " << result.optimal_f << std::endl;

    return tuple<ArrayXd,double>(theta3_a,obj3_a);
}

//// Compute the objective function
tuple<double,int,EquStateV> alias::EstimationAuxiliaryModel_part3_ObjFun(const std::vector<double> & theta3,
    const SimData & RealData, const ArrayXd & theta1_a, const ArrayXd & theta2_a,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {

    //// Solve the value function and policy functions given parameter values
    tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
        theta3, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
    ParaEst para_est = get<0>(t_Equ);
    ParaVec para_vec = get<1>(t_Equ);
    EquStateV EquV = get<2>(t_Equ);
    EquStateVmat Evalmat = get<3>(t_Equ);

    tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi> t_obj = EstimationAuxiliaryModel_part3_likelihood(
        RealData,para_est,para_vec,EquV,Evalmat,GoodState,LaborIntensive,threadsManagement);
    double y_sum = get<0>(t_obj);
    int Ncount = get<1>(t_obj);
////    ArrayXXd likelihood = get<2>(t_obj);
////    ArrayXd obj_grad = ArrayXd::Zero(para.dim3);
////    cout << "y_sum = " << y_sum << "; Ncount = " << Ncount << endl;
////    throw runtime_error("634");
    return tuple<double,int,EquStateV>(y_sum,Ncount,EquV);
}

//// solve the value function/policy function for a given set of parameters
tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> alias::EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
    const std::vector<double> & theta3, const ArrayXd & theta1_a, const ArrayXd & theta2_a,
    const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management & threadsManagement) {

    ArrayXd theta3_a_vec(para.dim3);
    for (size_t k = 0; k < para.dim3; ++k) { theta3_a_vec(k) = theta3[k]; }
    ArrayXd theta_a(para.dim1 + para.dim2 +para.dim3);
    theta_a << theta3_a_vec,theta1_a,theta2_a;
//    cout << "***NEW*** theta_a = " << theta_a.transpose() << endl;

    //// construct the set of parameters
    ParaEst para_est = constructParaEst(theta_a,GoodState,LaborIntensive);
    para_est.PI = 1.0;
//    cout << "para_est = " << "F_P = " << para_est.F_P_PP << ", " << para_est.F_P_DP << ", " << para_est.sigma_PD << ", "
//        << "Capital adj: " << para_est.H_K_Cons << ", " << para_est.H_K << ", "
//        << para_est.c_HK << "," << para_est.F_K_Cons << "," << para_est.F_K << "," << para_est.c_FK << ","
//        << "; RegWorker adj: " << para_est.H_ur_Cons << "," << para_est.H_ur << "," << para_est.c_H_ur << ","
//        << para_est.F_high_ur_Cons << "," << para_est.F_high_ur << "," << para_est.c_high_F_ur << ","
//        << para_est.F_low_ur_Cons << "," << para_est.F_low_ur << "," << para_est.c_low_F_ur << ","
//        << "; NonRegWorker adj: " << para_est.H_uc_Cons << "," << para_est.H_uc << "," << para_est.c_H_uc << ","
//        << para_est.F_uc_Cons << "," << para_est.F_uc << "," << para_est.c_F_uc << ","
//        << "; Sigma adj: " << para_est.sigma_Kerror << "," << para_est.sigma_Lerror_ur << "," << para_est.sigma_Lerror_uc << ","
//        << "; Exit Cost: " << para_est.F_E << "," << para_est.F_E_c_FK << "," << para_est.F_E_c_F_uc << ","
//        << para_est.F_E_c_F_ur << "," << para_est.sigma_SE << endl;
//    throw runtime_error("661");
//
    //// construct the grid of state space
    ParaVec para_vec = constructParaFull(para_est,para.p_K_StatusQuo);

    //// solve the value function policy fucntion
    ArrayXd VEnd_PD = ArrayXd::Zero(para.N_phi * para.N_K * para.N_Lur * para.N_Luc * 2);
    tuple<EquStateV,EquStateVmat> t = solveV(para_est,para_vec,VEnd_PD,threadsManagement);
    EquStateV EquV = get<0>(t);
    EquStateVmat Evalmat = get<1>(t);

    return tuple<ParaEst,ParaVec,EquStateV,EquStateVmat>(para_est,para_vec,EquV,Evalmat);
}

//// Calculate the likelihood for the data
tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi> alias::EstimationAuxiliaryModel_part3_likelihood(
    const SimData & RealData, const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV,
    const EquStateVmat & Evalmat, const int & GoodState, const int & LaborIntensive,
    MultiThreads::Threads_Management & threadsManagement) {

    int N_RealData = RealData.good.rows();

    ArrayXi Period_vec = ArrayXi::LinSpaced(para.EstTbar+1, 0, para.EstTbar);
    ArrayXXi Period_mat(N_RealData, para.EstTbar);
    Period_mat.rowwise() = Period_vec.transpose();

    ArrayXd Ones = ArrayXd::Ones(N_RealData);
    //// estimated productivity
    ArrayXXd lnphi_a = Backout_phi_a(RealData, para_est.alpha_tilde_K, para_est.alpha_tilde_L,
        para_est.alpha_Lr, para_est.alpha_Lc);

    ArrayXXi count = ArrayXXi::Zero(N_RealData, para.EstTbar);
    ArrayXXd logLikelihood = ArrayXXd::Zero(N_RealData, para.EstTbar);

    ArrayXXi K_index_max_mat = ArrayXXi::Zero(N_RealData, para.EstTbar);
    ArrayXXi Lur_index_max_mat = ArrayXXi::Zero(N_RealData, para.EstTbar);
    ArrayXXi Luc_index_max_mat = ArrayXXi::Zero(N_RealData, para.EstTbar);
//    EstArray est_array = initializeEstArray(N_RealData);

    int ExitCount = 0;
    int EventCount = 0;
    auto worker = [&](size_t n, unsigned thread_id) {
//        size_t thread_id = 0;
//    for (size_t n = 10242; n < N_RealData; ++n) {
        for (size_t t = RealData.first_year_in_Data(n) + 1; t <= min(RealData.last_year_in_Data(n), para.EstTbar - 1); ++t) {

            ArrayXi period_vec_t = Period_mat(n, seqN(0, t+1))
                * RealData.State_P_nonMiss(n, seqN(0, t+1))
                * (1 - RealData.Revenue_Miss(n, seqN(0, t+1)))
                * (1 - RealData.Employ_urMiss(n, seqN(0, t+1)))
                * (1 - RealData.Employ_ucMiss(n, seqN(0, t+1)))
                * (1 - RealData.CapitalMiss(n, seqN(0, t+1)));
            int t_lag = period_vec_t.maxCoeff();
            double lnphi_a_lag = lnphi_a(n, t_lag);

            period_vec_t = Period_mat(n, seqN(0, t))
                * (1 - RealData.CapitalMiss(n, seqN(0, t)));
            t_lag = period_vec_t.maxCoeff();
            double Capital_lag = RealData.Capital(n, t_lag);
//            cout << "Capital_lag = " << Capital_lag << "; RealData.Capital(n, t) = " << RealData.Capital(n, t) << endl;

            if (RealData.State_E_Miss(n, t) == 1 and t_lag<t-1) {
                Capital_lag = 0.0/(double(t-t_lag+1)) * RealData.Capital(n, t_lag);
            }

            period_vec_t = Period_mat(n, seqN(0, t))
                * (1 - RealData.Employ_urMiss(n, seqN(0, t)));
            t_lag = period_vec_t.maxCoeff();
            double Employ_ur_lag = RealData.Employ_ur(n, t_lag);

            if (RealData.State_E_Miss(n, t) == 1 and t_lag<t-1) {
                double temp = min( 99.0, double(RealData.Employ_ur(n, t_lag)) ) ;
                Employ_ur_lag = 1.0/(double(t-t_lag+1)) * temp;
////                Employ_ur_lag = 0.2*RealData.Employ_ur(n, t_lag);
//                double temp = min( 99.0, double(RealData.Employ_ur(n, t_lag)) ) ;
////                Employ_ur_lag = 0.5*temp;
            }

            period_vec_t = Period_mat(n, seqN(0, t))
                           * (1 - RealData.Employ_ucMiss(n, seqN(0, t)));
            t_lag = period_vec_t.maxCoeff();
            double Employ_uc_lag = RealData.Employ_uc(n, t_lag);
            if (RealData.State_E_Miss(n, t) == 1 and t_lag<t-1) {
                Employ_uc_lag = 1.0/(double(t-t_lag+1)) * RealData.Employ_uc(n, t_lag);
            }

            if ( (RealData.State_P_nonMiss(n,t-1) == 1 or RealData.State_D_nonMiss(n,t-1) == 1))  {
                ExitCount = ExitCount+1;
//
                int PD_lag = 0;
                if (RealData.State_P_nonMiss(n,t-1) == 1) {PD_lag = 0;}
                else if (RealData.State_D_nonMiss(n,t-1) == 1) {PD_lag = 1;}

                PhiW phi_w = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_lag));
                LKState K_state = defineLKState(para_vec.vec_K, Capital_lag);
                LKState Lur_state = defineLKState(para_vec.vec_Lur, Employ_ur_lag);
                LKState Luc_state = defineLKState(para_vec.vec_Luc, Employ_uc_lag);
//                cout << "t = " << t << "; Employ_uc_lag = " << Employ_uc_lag << "; RealData.Employ_uc(n, t) = " << RealData.Employ_uc(n, t) << endl;
                // cout << "para_est.F_E = " << para_est.F_E << "; para_est.sigma = " << para_est.sigma_SE << endl;
                tuple<double, double, double, double> t_SE = calStayExit_Prob_logit(para_est, para_vec, Evalmat,
                    phi_w, K_state, Lur_state, Luc_state, PD_lag);
                double Prob_S = get<0>(t_SE);
                double Prob_E = get<1>(t_SE);
                double lnProb_S = get<2>(t_SE);
                double lnProb_E = get<3>(t_SE);
//
//                est_array.Est_ProbS(n, t) = Prob_S;
//                est_array.Est_ProbE(n, t) = Prob_E;
//                est_array.Est_lnProbS(n, t) = lnProb_S;
//                est_array.Est_lnProbE(n, t) = lnProb_E;
//                cout << "Prob_S = " << Prob_S << "; Prob_E = " << Prob_E << "; lnProb_S = " << lnProb_S << "; lnProb_E = " << lnProb_E << endl;

                tuple<double, double, double, double, double, double, double, double, double, double, int, int, int> t_Opt
                    = cal_prob_KLurLuc_Data(para_est, para_vec, Evalmat, phi_w, K_state, Lur_state, Luc_state, PD_lag,
                    RealData.Capital(n, t), RealData.Employ_ur(n, t),
                    RealData.Employ_uc(n, t),
                    RealData.CapitalMiss(n,t), RealData.Employ_urMiss(n,t),
                    RealData.Employ_ucMiss(n,t));
                double Prob_P = get<0>(t_Opt);
                double Prob_D = get<1>(t_Opt);
                double lnProb_P = get<2>(t_Opt);
                double lnProb_D = get<3>(t_Opt);
                double lnpdf_State_K = get<4>(t_Opt);
                double lnpdf_State_Lur = get<5>(t_Opt);
                double lnpdf_State_Luc = get<6>(t_Opt);

                double K_opt = get<7>(t_Opt);
                double Lur_opt = get<8>(t_Opt);
                double Luc_opt = get<9>(t_Opt);
                int K_index_max = get<10>(t_Opt);
                int Lur_index_max = get<11>(t_Opt);
                int Luc_index_max = get<12>(t_Opt);

                K_index_max_mat(n, t) = K_index_max;
                Lur_index_max_mat(n, t) = Lur_index_max;
                Luc_index_max_mat(n, t) = Luc_index_max;

                double Prob_Missing = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                    RealData.Revenue(n, t - 1),
                    0,0,0,
                    RealData.State_D_nonMiss(n, t - 1),
                    RealData.State_DP_Miss(n, t - 1),
                    RealData.State_P_nonMiss(n, t - 1),
                    RealData.yr_DP(n,t-1),
                    RealData.state_code(n), RealData.industry(n), t,
                    0, GoodState, LaborIntensive);

//                if (1==1) {
//                    double Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
//                        RealData.Revenue(n, t - 1),
//                        0,0,0,
//                        RealData.State_D_nonMiss(n, t - 1),
//                        RealData.State_DP_Miss(n, t - 1),
//                        RealData.State_P_nonMiss(n, t - 1),
//                        RealData.state_code(n), RealData.industry(n), t,
//                        0, GoodState, LaborIntensive);
//
//                    // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
//                    double Prob_S_Missing = Prob_S * Prob_Missing;
//                    double Prob_Missing_cum = Prob_E;
//                    for (size_t tt = t; tt < para.SimTbar; ++tt) {
//                        Prob_Missing_cum = Prob_Missing_cum + Prob_S_Missing * Prob_E;
//                        Prob_S_Missing = Prob_S_Missing * Prob_S * Prob_Miss_DP;
//                    }
//                    est_array.Est_Prob_E_Miss(n, t) = Prob_Missing_cum;
//                    est_array.Prob_Miss_DP(n, t) = Prob_Miss_DP;
//                }
//
//                if (0==0) {
//                    double Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
//                        RealData.Revenue(n, t - 1),
//                        0,0,0,
//                        RealData.State_D_nonMiss(n, t - 1),
//                        RealData.State_DP_Miss(n, t - 1),
//                        RealData.State_P_nonMiss(n, t - 1),
//                        RealData.state_code(n), RealData.industry(n), t,
//                        0, GoodState, LaborIntensive);
//
//                    // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
//                    double Prob_S_Missing = Prob_S * Prob_Missing;
//
//                    double Prob_S_Missing_next = Prob_S * Prob_Miss_DP;
//                    double Prob_Missing_next_cum = Prob_E;
//                    for (size_t tt = t+1; tt < para.SimTbar; ++tt) {
//                        Prob_Missing_next_cum = Prob_Missing_next_cum + Prob_S_Missing_next * Prob_E;
//                        Prob_S_Missing_next = Prob_S_Missing_next * Prob_S * Prob_Miss_DP;
//                    }
//
//                    est_array.Est_Prob_DP_Miss(n, t) = Prob_S_Missing*(1-Prob_Missing_next_cum);
////                    logLikelihood(n, t) = 0.0;
//                }


                if (RealData.State_P_nonMiss(n, t) == 1) {
                    if (RealData.Employ_ur_Extreme(n,t) == 1) {lnpdf_State_Lur = 0.0;}
                    if (RealData.Employ_uc_Extreme(n,t) == 1) {lnpdf_State_Luc = 0.0;}
                    if (RealData.Capital_Extreme(n,t) == 1) {lnpdf_State_K = 0.0;}

                    logLikelihood(n, t) = lnpdf_State_K + lnpdf_State_Lur + lnpdf_State_Luc
                        + lnProb_P + lnProb_S + log(1 - Prob_Missing);
//                    logLikelihood(n, t) = lnProb_S;
                    count(n, t) = 1;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "lnpdf_State_K = " << lnpdf_State_K << "; lnpdf_State_Lur = " << lnpdf_State_Lur
                             << "; lnpdf_State_Luc = " << lnpdf_State_Luc << "; lnProb_S = " << lnProb_S
                             << "; Prob_Missing = " << Prob_Missing << endl;
                        cout << "lnProb_P = " << lnProb_P << "; Prob_P = " << Prob_P << endl;
                        throw runtime_error("863");
                    }

//                    est_array.Est_lnpdf_KOpt(n, t) = lnpdf_State_K;
//                    est_array.Est_lnpdf_LurOpt(n, t) = lnpdf_State_Lur;
//                    est_array.Est_lnpdf_LucOpt(n, t) = lnpdf_State_Luc;
                }
                else if (RealData.State_D_nonMiss(n, t) == 1) {
                    if (RealData.Employ_ur_Extreme(n,t) == 1) {lnpdf_State_Lur = 0.0;}
                    if (RealData.Employ_uc_Extreme(n,t) == 1) {lnpdf_State_Luc = 0.0;}
                    if (RealData.Capital_Extreme(n,t) == 1) {lnpdf_State_K = 0.0;}

                    logLikelihood(n, t) = lnpdf_State_K + lnpdf_State_Lur + lnpdf_State_Luc
                        + lnProb_D + lnProb_S + log(1 - Prob_Missing);
                    count(n, t) = 1;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "lnpdf_State_K = " << lnpdf_State_K << "; lnpdf_State_Lur = " << lnpdf_State_Lur
                             << "; lnpdf_State_Luc = " << lnpdf_State_Luc << "; lnProb_P = " << lnProb_P
                             << "; lnProb_S = " << lnProb_S << "; Prob_Missing = " << Prob_Missing << endl;
                        throw runtime_error("874");
                    }
                    count(n, t) = 1;

//                    est_array.Est_lnpdf_KOpt(n, t) = lnpdf_State_K;
//                    est_array.Est_lnpdf_LurOpt(n, t) = lnpdf_State_Lur;
//                    est_array.Est_lnpdf_LucOpt(n, t) = lnpdf_State_Luc;
                }
                else if (RealData.State_E_Miss(n, t) == 1) {
                    EventCount = EventCount + 1;
                    double Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                        RealData.Revenue(n, t - 1),
                        0,0,0,
                        RealData.State_D_nonMiss(n, t - 1),
                        RealData.State_DP_Miss(n, t - 1),
                        RealData.State_P_nonMiss(n, t - 1),
                        0,
                        RealData.state_code(n), RealData.industry(n), t,
                        0, GoodState, LaborIntensive);

                    // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
                    double Prob_S_Missing = Prob_S * Prob_Missing;
                    double Prob_Missing_cum = Prob_E;
                    for (size_t tt = t; tt < para.SimTbar; ++tt) {
                            Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                                RealData.Revenue(n, t - 1),
                                0,0,0,
                                RealData.State_D_nonMiss(n, t - 1),
                                RealData.State_DP_Miss(n, t - 1),
                                RealData.State_P_nonMiss(n, t - 1),
                                tt-t+1,
                                RealData.state_code(n), RealData.industry(n), t,
                                0, GoodState, LaborIntensive);
                            Prob_Missing_cum = Prob_Missing_cum + Prob_S_Missing * Prob_E;
                            Prob_S_Missing = Prob_S_Missing * Prob_S * Prob_Miss_DP;
                    }
                    Prob_Missing_cum = Prob_Missing_cum + Prob_S_Missing;
                    logLikelihood(n, t) = log(Prob_Missing_cum);
//                    logLikelihood(n, t) = 0.0;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "Prob_E = " << Prob_E << "; Prob_S = " << Prob_S << endl;
                        throw runtime_error("954");
                    }
                    count(n, t) = 1;
                }
                else if (RealData.State_DP_Miss(n, t) == 1) {
                    double Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                        RealData.Revenue(n, t - 1),
                        0,0,0,
                        RealData.State_D_nonMiss(n, t - 1),
                        RealData.State_DP_Miss(n, t - 1),
                        RealData.State_P_nonMiss(n, t - 1),
                        0,
                        RealData.state_code(n), RealData.industry(n), t,
                        0, GoodState, LaborIntensive);

                    // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
                    double Prob_S_Missing = Prob_S * Prob_Missing;

                    double Prob_S_Missing_next = Prob_S * Prob_Miss_DP;
                    double Prob_Missing_next_cum = Prob_E;
                    for (size_t tt = t+1; tt < para.SimTbar; ++tt) {
                        Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                            RealData.Revenue(n, t - 1),
                            0,0,0,
                            RealData.State_D_nonMiss(n, t - 1),
                            RealData.State_DP_Miss(n, t - 1),
                            RealData.State_P_nonMiss(n, t - 1),
                            tt-t+1,
                            RealData.state_code(n), RealData.industry(n), t,
                            0, GoodState, LaborIntensive);
                        Prob_Missing_next_cum = Prob_Missing_next_cum + Prob_S_Missing_next * Prob_E;
                        Prob_S_Missing_next = Prob_S_Missing_next * Prob_S * Prob_Miss_DP;
                    }
                    Prob_Missing_next_cum = Prob_Missing_next_cum + Prob_S_Missing_next;

                    logLikelihood(n, t) = log(Prob_S_Missing*(1-Prob_Missing_next_cum));
//                    logLikelihood(n, t) = 0.0;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "para_est.F_E" << para_est.F_E << "; para_est.sigma = " << para_est.sigma_SE << endl;
                        cout << "RealData.State_D_nonMiss(n, t - 1) = " << RealData.State_D_nonMiss(n, t - 1)
                             << "; RealData.State_DP_Miss(n, t - 1) = " << RealData.State_DP_Miss(n, t - 1)
                             << "; RealData.State_P_nonMiss(n, t - 1) = " << RealData.State_P_nonMiss(n, t - 1) << endl;
                        cout << "Prob_Missing = " << Prob_Missing << "; Prob_S = " << Prob_S
                             << "; Prob_S_Missing = " << Prob_S_Missing << "; Prob_Missing_next_cum = " << Prob_Missing_next_cum << endl;
                        throw runtime_error("881");
                    }
                    count(n, t) = 1;
                }
            }
        }
//    }
    };
    MultiThreads::simple_parallel_for(worker,N_RealData, threadsManagement);

    int Ncount = count.sum();
    double y_sum = logLikelihood.sum();
////////////
//    writeToCSVfile("lnphi_a.csv", lnphi_a.cast<double>().matrix());
//    writeToCSVfile("Est_KOpt.csv", est_array.Est_KOpt.cast<double>().matrix());
//    writeToCSVfile("Est_LurOpt.csv", est_array.Est_LurOpt.cast<double>().matrix());
//    writeToCSVfile("Est_LucOpt.csv", est_array.Est_LucOpt.cast<double>().matrix());
//
//    writeToCSVfile("K_index_max_mat.csv", K_index_max_mat.cast<double>().matrix());
//    writeToCSVfile("Lur_index_max_mat.csv", Lur_index_max_mat.cast<double>().matrix());
//    writeToCSVfile("Luc_index_max_mat.csv", Luc_index_max_mat.cast<double>().matrix());
//
//    writeToCSVfile("Est_lnpdf_KOpt.csv", est_array.Est_lnpdf_KOpt.cast<double>().matrix());
//    writeToCSVfile("Est_lnpdf_LurOpt.csv", est_array.Est_lnpdf_LurOpt.cast<double>().matrix());
//    writeToCSVfile("Est_lnpdf_LucOpt.csv", est_array.Est_lnpdf_LucOpt.cast<double>().matrix());
////
//    writeToCSVfile("Est_ProbS.csv", est_array.Est_ProbS.cast<double>().matrix());
//    writeToCSVfile("Est_ProbE.csv", est_array.Est_ProbE.cast<double>().matrix());
//    writeToCSVfile("Est_Prob_E_Miss.csv", est_array.Est_Prob_E_Miss.cast<double>().matrix());
//    writeToCSVfile("Est_Prob_DP_Miss.csv", est_array.Est_Prob_DP_Miss.cast<double>().matrix());
//    writeToCSVfile("Est_Prob_Miss_DP.csv", est_array.Prob_Miss_DP.cast<double>().matrix());
//
//
//    writeToCSVfile("Capital.csv", RealData.Capital.cast<double>().matrix());
//    writeToCSVfile("Employ_ur.csv", RealData.Employ_ur.cast<double>().matrix());
//    writeToCSVfile("Employ_uc.csv", RealData.Employ_uc.cast<double>().matrix());
//    writeToCSVfile("Revenue.csv", RealData.Revenue.cast<double>().matrix());
//    writeToCSVfile("State_P_nonMiss.csv", RealData.State_P_nonMiss.cast<double>().matrix());
//    writeToCSVfile("State_D_nonMiss.csv", RealData.State_D_nonMiss.cast<double>().matrix());
//    writeToCSVfile("State_DP_Miss.csv", RealData.State_DP_Miss.cast<double>().matrix());
//    writeToCSVfile("State_E_Miss.csv", RealData.State_E_Miss.cast<double>().matrix());
//////
//    cout << "average Exit = " << double(EventCount) / double(ExitCount) << endl;
//    cout << "average Prob_E_Miss = " << est_array.Est_Prob_E_Miss.sum() / (est_array.Est_Prob_E_Miss > 0).cast<double>().sum()
//        << "; (est_array.Est_Prob_E_Miss > 0).cast<double>().sum() = " << (est_array.Est_Prob_E_Miss > 0).cast<double>().sum()
//        << "; est_array.Est_Prob_E_Miss.sum() = " << est_array.Est_Prob_E_Miss.sum() << endl;
////    throw runtime_error("1349");
////////    throw runtime_error("898");
    return tuple<double,int,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi>(y_sum,Ncount,logLikelihood,
        K_index_max_mat,Lur_index_max_mat,Luc_index_max_mat);
}

tuple<double,int,ArrayXXd> alias::EstimationAuxiliaryModel_part3_likelihood_Diff(const SimData & RealData,
    const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV, const EquStateVmat & Evalmat,
    const ArrayXXi & K_index_max_mat, const ArrayXXi & Lur_index_max_mat, const ArrayXXi & Luc_index_max_mat,
    const int & GoodState, const int & LaborIntensive,MultiThreads::Threads_Management & threadsManagement) {

    int N_RealData = RealData.good.rows();

    ArrayXi Period_vec = ArrayXi::LinSpaced(para.EstTbar+1, 0, para.EstTbar);
    ArrayXXi Period_mat(N_RealData, para.EstTbar);
    Period_mat.rowwise() = Period_vec.transpose();

    ArrayXd Ones = ArrayXd::Ones(N_RealData);
    //// estimated productivity
    ArrayXXd lnphi_a = Backout_phi_a(RealData, para_est.alpha_tilde_K, para_est.alpha_tilde_L,
        para_est.alpha_Lr, para_est.alpha_Lc);

    ArrayXXi count = ArrayXXi::Zero(N_RealData, para.EstTbar);
    ArrayXXd logLikelihood = ArrayXXd::Zero(N_RealData, para.EstTbar);

    int ExitCount = 0;
    int EventCount = 0;
    auto worker = [&](size_t n, unsigned thread_id) {
//        size_t thread_id = 0;
//    for (size_t n = 10242; n < N_RealData; ++n) {
        for (size_t t = RealData.first_year_in_Data(n) + 1; t <= min(RealData.last_year_in_Data(n), para.EstTbar - 1); ++t) {

            ArrayXi period_vec_t = Period_mat(n, seqN(0, t+1))
                * RealData.State_P_nonMiss(n, seqN(0, t+1))
                * (1 - RealData.Revenue_Miss(n, seqN(0, t+1)))
                * (1 - RealData.Employ_urMiss(n, seqN(0, t+1)))
                * (1 - RealData.Employ_ucMiss(n, seqN(0, t+1)))
                * (1 - RealData.CapitalMiss(n, seqN(0, t+1)));
            int t_lag = period_vec_t.maxCoeff();
            double lnphi_a_lag = lnphi_a(n, t_lag);

            period_vec_t = Period_mat(n, seqN(0, t))
                * (1 - RealData.CapitalMiss(n, seqN(0, t)));
            t_lag = period_vec_t.maxCoeff();
            double Capital_lag = RealData.Capital(n, t_lag);
//            cout << "Capital_lag = " << Capital_lag << "; RealData.Capital(n, t) = " << RealData.Capital(n, t) << endl;

            if (RealData.State_E_Miss(n, t) == 1 and t_lag<t-1) {
                Capital_lag = 0.0/(double(t-t_lag+1)) * RealData.Capital(n, t_lag);
            }

            period_vec_t = Period_mat(n, seqN(0, t))
                * (1 - RealData.Employ_urMiss(n, seqN(0, t)));
            t_lag = period_vec_t.maxCoeff();
            double Employ_ur_lag = RealData.Employ_ur(n, t_lag);

            if (RealData.State_E_Miss(n, t) == 1 and t_lag<t-1) {
                double temp = min( 99.0, double(RealData.Employ_ur(n, t_lag)) ) ;
                Employ_ur_lag = 1.0/(double(t-t_lag+1)) * temp;
////                Employ_ur_lag = 0.2*RealData.Employ_ur(n, t_lag);
//                double temp = min( 99.0, double(RealData.Employ_ur(n, t_lag)) ) ;
////                Employ_ur_lag = 0.5*temp;
            }

            period_vec_t = Period_mat(n, seqN(0, t))
                           * (1 - RealData.Employ_ucMiss(n, seqN(0, t)));
            t_lag = period_vec_t.maxCoeff();
            double Employ_uc_lag = RealData.Employ_uc(n, t_lag);
            if (RealData.State_E_Miss(n, t) == 1 and t_lag<t-1) {
                Employ_uc_lag = 1.0/(double(t-t_lag+1)) * RealData.Employ_uc(n, t_lag);
            }

            if ( (RealData.State_P_nonMiss(n,t-1) == 1 or RealData.State_D_nonMiss(n,t-1) == 1))  {
                ExitCount = ExitCount+1;
//
                int PD_lag = 0;
                if (RealData.State_P_nonMiss(n,t-1) == 1) {PD_lag = 0;}
                else if (RealData.State_D_nonMiss(n,t-1) == 1) {PD_lag = 1;}

                PhiW phi_w = definePhiWeights(para_vec.vec_phi, exp(lnphi_a_lag));
                LKState K_state = defineLKState(para_vec.vec_K, Capital_lag);
                LKState Lur_state = defineLKState(para_vec.vec_Lur, Employ_ur_lag);
                LKState Luc_state = defineLKState(para_vec.vec_Luc, Employ_uc_lag);
//                cout << "t = " << t << "; Employ_uc_lag = " << Employ_uc_lag << "; RealData.Employ_uc(n, t) = " << RealData.Employ_uc(n, t) << endl;
                // cout << "para_est.F_E = " << para_est.F_E << "; para_est.sigma = " << para_est.sigma_SE << endl;
                tuple<double, double, double, double> t_SE = calStayExit_Prob_logit(para_est, para_vec, Evalmat,
                    phi_w, K_state, Lur_state, Luc_state, PD_lag);
                double Prob_S = get<0>(t_SE);
                double Prob_E = get<1>(t_SE);
                double lnProb_S = get<2>(t_SE);
                double lnProb_E = get<3>(t_SE);
//
//                est_array.Est_ProbS(n, t) = Prob_S;
//                est_array.Est_ProbE(n, t) = Prob_E;
//                est_array.Est_lnProbS(n, t) = lnProb_S;
//                est_array.Est_lnProbE(n, t) = lnProb_E;
//                cout << "Prob_S = " << Prob_S << "; Prob_E = " << Prob_E << "; lnProb_S = " << lnProb_S << "; lnProb_E = " << lnProb_E << endl;

                tuple<double, double, double, double, double, double, double> t_Opt
                    = cal_prob_KLurLuc_Data_Diff(para_est, para_vec, Evalmat, phi_w, K_state, Lur_state, Luc_state, PD_lag,
                    RealData.Capital(n, t), RealData.Employ_ur(n, t), RealData.Employ_uc(n, t),
                    RealData.CapitalMiss(n,t), RealData.Employ_urMiss(n,t), RealData.Employ_ucMiss(n,t),
                    K_index_max_mat(n, t),Lur_index_max_mat(n, t),Luc_index_max_mat(n, t));
                double Prob_P = get<0>(t_Opt);
                double Prob_D = get<1>(t_Opt);
                double lnProb_P = get<2>(t_Opt);
                double lnProb_D = get<3>(t_Opt);
                double lnpdf_State_K = get<4>(t_Opt);
                double lnpdf_State_Lur = get<5>(t_Opt);
                double lnpdf_State_Luc = get<6>(t_Opt);

                double Prob_Missing = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                    RealData.Revenue(n, t - 1),
                    0,0,0,
                    RealData.State_D_nonMiss(n, t - 1),
                    RealData.State_DP_Miss(n, t - 1),
                    RealData.State_P_nonMiss(n, t - 1),
                    RealData.yr_DP(n,t-1),
                    RealData.state_code(n), RealData.industry(n), t,
                    0, GoodState, LaborIntensive);


                if (RealData.State_P_nonMiss(n, t) == 1) {
                    if (RealData.Employ_ur_Extreme(n,t) == 1) {lnpdf_State_Lur = 0.0;}
                    if (RealData.Employ_uc_Extreme(n,t) == 1) {lnpdf_State_Luc = 0.0;}
                    if (RealData.Capital_Extreme(n,t) == 1) {lnpdf_State_K = 0.0;}

                    logLikelihood(n, t) = lnpdf_State_K + lnpdf_State_Lur + lnpdf_State_Luc
                        + lnProb_P + lnProb_S + log(1 - Prob_Missing);
//                    logLikelihood(n, t) = lnProb_S;
                    count(n, t) = 1;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "lnpdf_State_K = " << lnpdf_State_K << "; lnpdf_State_Lur = " << lnpdf_State_Lur
                             << "; lnpdf_State_Luc = " << lnpdf_State_Luc << "; lnProb_S = " << lnProb_S
                             << "; Prob_Missing = " << Prob_Missing << endl;
                        cout << "lnProb_P = " << lnProb_P << "; Prob_P = " << Prob_P << endl;
                        throw runtime_error("863");
                    }

//                    est_array.Est_lnpdf_KOpt(n, t) = lnpdf_State_K;
//                    est_array.Est_lnpdf_LurOpt(n, t) = lnpdf_State_Lur;
//                    est_array.Est_lnpdf_LucOpt(n, t) = lnpdf_State_Luc;
                }
                else if (RealData.State_D_nonMiss(n, t) == 1) {
                    if (RealData.Employ_ur_Extreme(n,t) == 1) {lnpdf_State_Lur = 0.0;}
                    if (RealData.Employ_uc_Extreme(n,t) == 1) {lnpdf_State_Luc = 0.0;}
                    if (RealData.Capital_Extreme(n,t) == 1) {lnpdf_State_K = 0.0;}

                    logLikelihood(n, t) = lnpdf_State_K + lnpdf_State_Lur + lnpdf_State_Luc
                        + lnProb_D + lnProb_S + log(1 - Prob_Missing);
                    count(n, t) = 1;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "lnpdf_State_K = " << lnpdf_State_K << "; lnpdf_State_Lur = " << lnpdf_State_Lur
                             << "; lnpdf_State_Luc = " << lnpdf_State_Luc << "; lnProb_P = " << lnProb_P
                             << "; lnProb_S = " << lnProb_S << "; Prob_Missing = " << Prob_Missing << endl;
                        throw runtime_error("874");
                    }
                    count(n, t) = 1;

//                    est_array.Est_lnpdf_KOpt(n, t) = lnpdf_State_K;
//                    est_array.Est_lnpdf_LurOpt(n, t) = lnpdf_State_Lur;
//                    est_array.Est_lnpdf_LucOpt(n, t) = lnpdf_State_Luc;
                }
                else if (RealData.State_E_Miss(n, t) == 1) {
                    EventCount = EventCount + 1;
                    double Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                        RealData.Revenue(n, t - 1),
                        0,0,0,
                        RealData.State_D_nonMiss(n, t - 1),
                        RealData.State_DP_Miss(n, t - 1),
                        RealData.State_P_nonMiss(n, t - 1),
                        0,
                        RealData.state_code(n), RealData.industry(n), t,
                        0, GoodState, LaborIntensive);

                    // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
                    double Prob_S_Missing = Prob_S * Prob_Missing;
                    double Prob_Missing_cum = Prob_E;
                    for (size_t tt = t; tt < para.SimTbar; ++tt) {
                            Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                                RealData.Revenue(n, t - 1),
                                0,0,0,
                                RealData.State_D_nonMiss(n, t - 1),
                                RealData.State_DP_Miss(n, t - 1),
                                RealData.State_P_nonMiss(n, t - 1),
                                tt-t+1,
                                RealData.state_code(n), RealData.industry(n), t,
                                0, GoodState, LaborIntensive);
                            Prob_Missing_cum = Prob_Missing_cum + Prob_S_Missing * Prob_E;
                            Prob_S_Missing = Prob_S_Missing * Prob_S * Prob_Miss_DP;
                    }
                    Prob_Missing_cum = Prob_Missing_cum + Prob_S_Missing;
                    logLikelihood(n, t) = log(Prob_Missing_cum);
//                    logLikelihood(n, t) = 0.0;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "Prob_E = " << Prob_E << "; Prob_S = " << Prob_S << endl;
                        throw runtime_error("954");
                    }
                    count(n, t) = 1;
                }
                else if (RealData.State_DP_Miss(n, t) == 1) {
                    double Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                        RealData.Revenue(n, t - 1),
                        0,0,0,
                        RealData.State_D_nonMiss(n, t - 1),
                        RealData.State_DP_Miss(n, t - 1),
                        RealData.State_P_nonMiss(n, t - 1),
                        0,
                        RealData.state_code(n), RealData.industry(n), t,
                        0, GoodState, LaborIntensive);

                    // Prob_Missing_cum is a cummulative prob. The firm will never appear in the data, either exit or disappear
                    double Prob_S_Missing = Prob_S * Prob_Missing;

                    double Prob_S_Missing_next = Prob_S * Prob_Miss_DP;
                    double Prob_Missing_next_cum = Prob_E;
                    for (size_t tt = t+1; tt < para.SimTbar; ++tt) {
                        Prob_Miss_DP = ProbMissing(Capital_lag,Employ_ur_lag,Employ_uc_lag,
                            RealData.Revenue(n, t - 1),
                            0,0,0,
                            RealData.State_D_nonMiss(n, t - 1),
                            RealData.State_DP_Miss(n, t - 1),
                            RealData.State_P_nonMiss(n, t - 1),
                            tt-t+1,
                            RealData.state_code(n), RealData.industry(n), t,
                            0, GoodState, LaborIntensive);
                        Prob_Missing_next_cum = Prob_Missing_next_cum + Prob_S_Missing_next * Prob_E;
                        Prob_S_Missing_next = Prob_S_Missing_next * Prob_S * Prob_Miss_DP;
                    }
                    Prob_Missing_next_cum = Prob_Missing_next_cum + Prob_S_Missing_next;

                    logLikelihood(n, t) = log(Prob_S_Missing*(1-Prob_Missing_next_cum));
//                    logLikelihood(n, t) = 0.0;

                    if (isfinite(logLikelihood(n, t)) == 0) {
                        cout << "para_est.F_E" << para_est.F_E << "; para_est.sigma = " << para_est.sigma_SE << endl;
                        cout << "RealData.State_D_nonMiss(n, t - 1) = " << RealData.State_D_nonMiss(n, t - 1)
                             << "; RealData.State_DP_Miss(n, t - 1) = " << RealData.State_DP_Miss(n, t - 1)
                             << "; RealData.State_P_nonMiss(n, t - 1) = " << RealData.State_P_nonMiss(n, t - 1) << endl;
                        cout << "Prob_Missing = " << Prob_Missing << "; Prob_S = " << Prob_S
                             << "; Prob_S_Missing = " << Prob_S_Missing << "; Prob_Missing_next_cum = " << Prob_Missing_next_cum << endl;
                        throw runtime_error("881");
                    }
                    count(n, t) = 1;
                }
            }
        }
//    }
    };
    MultiThreads::simple_parallel_for(worker,N_RealData, threadsManagement);

    int Ncount = count.sum();
    double y_sum = logLikelihood.sum();

    return tuple<double,int,ArrayXXd>(y_sum,Ncount,logLikelihood);
}

/********************************************************************************************************************
* Auxilliary Model Estimation: Part 3: functions used in EstimationAuxiliaryModel_part3_likelihood
********************************************************************************************************************/
////**** Define the state of productivity ****////
PhiW alias::definePhiWeights(const ArrayXd & vec_phi, const double & phi_state) {

    PhiW phi_w;

    phi_w.phi = phi_state;
    if (phi_state <= vec_phi(0)) {
        phi_w.w1 = 1; phi_w.w2 = 0; phi_w.phi_index = 0; phi_w.phi_grid_choice = 1;
    }
    else if (phi_state >= vec_phi(para.N_phi-1)) {
        phi_w.w1 = 0; phi_w.w2 = 1; phi_w.phi_index = para.N_phi-1; phi_w.phi_grid_choice = 1;
    }
    else {
        phi_w.phi_index = find_lower_bound_scale(vec_phi, phi_state);
        phi_w.w2 = (log(phi_state) - log(vec_phi(phi_w.phi_index)))
                   / (log(vec_phi(phi_w.phi_index + 1)) - log(vec_phi(phi_w.phi_index)));
        phi_w.w1 = (log(vec_phi(phi_w.phi_index + 1)) - log(phi_state))
                   / (log(vec_phi(phi_w.phi_index + 1)) - log(vec_phi(phi_w.phi_index)));
        phi_w.phi_grid_choice = 0;
        phi_w.phi = phi_state;
    }
    return phi_w;
}

////**** Define the state of labor and capital ****////
LKState alias::defineLKState(const ArrayXd & vec_LK,const double & LK_val) {

    LKState LK_state;
    int N = vec_LK.size();

    if (LK_val <= vec_LK(0)) {
        LK_state.i_val_state = 0;
        LK_state.w_1 = 1;
        LK_state.w_2 = 0;
        LK_state.val_state = vec_LK(0);
    }
    else if (LK_val >= vec_LK(N-1)) {
        LK_state.i_val_state = N-2;
        LK_state.w_1 = 0;
        LK_state.w_2 = 1;
        LK_state.val_state = vec_LK(N-1);
    }
    else {
        int i_temp = find_lower_bound_scale(vec_LK, LK_val);
        LK_state.i_val_state = min(i_temp,N-2);
        LK_state.w_2 = (log(LK_val) - log(vec_LK(LK_state.i_val_state)))
            / (log(vec_LK(LK_state.i_val_state+1)) - log(vec_LK(LK_state.i_val_state)));
        LK_state.w_1 = (log(vec_LK(LK_state.i_val_state+1)) - log(LK_val))
            / (log(vec_LK(LK_state.i_val_state+1)) - log(vec_LK(LK_state.i_val_state)));
        LK_state.val_state = LK_val;
    }

    return LK_state;
}

////**** Calculate the exit/stay probability for each data point  ****////
//// with logit distribution assumption
tuple<double,double,double,double> alias::calStayExit_Prob_logit(const ParaEst & para_est, const ParaVec & para_vec,
    const EquStateVmat & Evalmat, const PhiW & phi_w, const LKState & K_state, const LKState & Lur_state,
    const LKState & Luc_state, const int & PD_state) {

    //// choose Stay or Exit State
    ArrayXXd EVal_SE_mat = Evalmat.EVal_PD_SE_mat(seqN(PD_state*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
    ArrayXd EVal_K_Interp = cal_EVal_K_phi_interpolation(EVal_SE_mat, phi_w,
        Lur_state, Luc_state, para.N_K);

    double beta = (EVal_K_Interp(K_state.i_val_state+1) - EVal_K_Interp(K_state.i_val_state))
        / (para_vec.vec_K(K_state.i_val_state+1) - para_vec.vec_K(K_state.i_val_state));
    double EVal_S = beta * (K_state.val_state - para_vec.vec_K(K_state.i_val_state))
        + EVal_K_Interp(K_state.i_val_state);
    double ResV_E = calResidual_Value_Exit_double(para_est, para_vec.p_K, K_state.val_state,
        Lur_state.val_state,Luc_state.val_state);
//    cout << "*** ResV_E = " << ResV_E << "; EVal_S = " << EVal_S << "; beta = " << beta
//         << "; K_state.val_state = " << K_state.val_state
//         << "; para_vec.vec_K(K_state.i_val_state) = " << para_vec.vec_K(K_state.i_val_state)
//         << "; EVal_K_Interp(K_state.i_val_state+1) = " << EVal_K_Interp(K_state.i_val_state+1) << endl;
////    cout << "Evalmat.EVal_SE_mat = " << Evalmat.EVal_SE_mat.transpose() << endl;

    double size_firm_ur = para_vec.vec_Lur(Lur_state.i_val_state);
    double size_firm_K = para_vec.vec_K(K_state.i_val_state);
    double dEVal_SE_sigma = (EVal_S - ResV_E - para_est.F_E) / para_est.sigma_SE;
    double Prob_S; double Prob_E;
    double lnProb_S; double lnProb_E;

    double dEVal_SE_sigma_temp = max(dEVal_SE_sigma,-700.0);
    dEVal_SE_sigma_temp = min(dEVal_SE_sigma_temp,700.0);
    Prob_S = exp(dEVal_SE_sigma_temp - log(exp(dEVal_SE_sigma_temp) + 1.0));
    Prob_E = exp(-dEVal_SE_sigma_temp - log(exp(-dEVal_SE_sigma_temp) + 1.0));

    lnProb_S = dEVal_SE_sigma_temp - log(exp(dEVal_SE_sigma_temp) + 1.0);
    lnProb_E = -dEVal_SE_sigma_temp - log(exp(-dEVal_SE_sigma_temp) + 1.0);

    return tuple<double,double,double,double>(Prob_S,Prob_E,lnProb_S,lnProb_E);
}

////**** Interpolation the value function (EVal_K_mat; EVal_Lur_mat; EVal_Luc_mat)  ****////
//// EVal_K_mat
ArrayXd alias::cal_EVal_K_phi_interpolation(const ArrayXXd & EVal_K_mat, const PhiW & phi_w,
    const LKState & Lur_state, const LKState & Luc_state, const int & NK) {

    // Eigen::Map<const Eigen::ArrayXXd> EVprime_EK_mat(EVprime_EK.data(),para.N_Luc*para.N_Lur,para.N_phi*para.N_K);
    ArrayXXd EVal_K_phi_mat;
    if (phi_w.phi_grid_choice == 1) {
        EVal_K_phi_mat = EVal_K_mat(all,seqN(phi_w.phi_index*NK,NK));
    }
    else {
        EVal_K_phi_mat = phi_w.w1*EVal_K_mat(all,seqN(phi_w.phi_index*NK,NK))
            + phi_w.w2*EVal_K_mat(all,seqN((phi_w.phi_index+1)*NK,NK));
    }

    int Start_id_Lur1Luc1 = Lur_state.i_val_state * para.N_Luc + Luc_state.i_val_state;
    int Start_id_Lur2Luc1 = (Lur_state.i_val_state+1) * para.N_Luc + Luc_state.i_val_state;
    int Start_id_Lur1Luc2 = Lur_state.i_val_state * para.N_Luc + Luc_state.i_val_state+1;
    int Start_id_Lur2Luc2 = (Lur_state.i_val_state+1) * para.N_Luc + Luc_state.i_val_state+1;

    ArrayXd EVal_K_Lur1Luc1 = EVal_K_phi_mat.row(Start_id_Lur1Luc1);
    ArrayXd EVal_K_Lur2Luc1 = EVal_K_phi_mat.row(Start_id_Lur2Luc1);
    ArrayXd EVal_K_Lur1Luc2 = EVal_K_phi_mat.row(Start_id_Lur1Luc2);
    ArrayXd EVal_K_Lur2Luc2 = EVal_K_phi_mat.row(Start_id_Lur2Luc2);

    ArrayXd temp2 = EVal_K_Lur1Luc1 * Lur_state.w_2 * Luc_state.w_2;
    ArrayXd temp3 = EVal_K_Lur2Luc1 * Lur_state.w_1 * Luc_state.w_2;
    ArrayXd temp4 = EVal_K_Lur1Luc2 * Lur_state.w_2 * Luc_state.w_1;
    ArrayXd temp5 = EVal_K_Lur2Luc2 * Lur_state.w_1 * Luc_state.w_1;

    ArrayXd EVal_K_vec = temp2+temp3+temp4+temp5;
    return EVal_K_vec;
}
//// EVal_Lur_mat
ArrayXd alias::cal_EVal_Lur_phi_interpolation(const ArrayXXd & EVal_Lur_mat, const PhiW & phi_w,
    const LKState & K_state, const LKState & Luc_state, const int & NLur) {

    // Eigen::Map<const Eigen::ArrayXXd> EVprime_ELur_mat(EVprime_ELur.data(),para.N_Luc,para.N_phi*para.N_K*para.N_Lur);
    ArrayXXd EVal_Lur_phi_mat;
    if (phi_w.phi_grid_choice == 1) {
        EVal_Lur_phi_mat = EVal_Lur_mat(all,
            seqN(phi_w.phi_index*para.N_K*NLur,para.N_K*NLur));
    }
    else {
        EVal_Lur_phi_mat = phi_w.w1*EVal_Lur_mat(all,
            seqN(phi_w.phi_index*para.N_K*NLur,para.N_K*NLur))
            + phi_w.w2*EVal_Lur_mat(all,
            seqN((phi_w.phi_index+1)*para.N_K*NLur,para.N_K*NLur));
    }

    ArrayXd EVal_Lur_K1Luc1 = EVal_Lur_phi_mat(Luc_state.i_val_state,
        seqN(K_state.i_val_state*NLur,NLur));
    ArrayXd EVal_Lur_K2Luc1 = EVal_Lur_phi_mat(Luc_state.i_val_state,
        seqN((K_state.i_val_state+1)*NLur,NLur));
    ArrayXd EVal_Lur_K1Luc2 = EVal_Lur_phi_mat(Luc_state.i_val_state+1,
        seqN(K_state.i_val_state*NLur,NLur));
    ArrayXd EVal_Lur_K2Luc2 = EVal_Lur_phi_mat(Luc_state.i_val_state+1,
        seqN((K_state.i_val_state+1)*NLur,NLur));

    ArrayXd temp2 = EVal_Lur_K1Luc1 * K_state.w_2 * Luc_state.w_2;
    ArrayXd temp3 = EVal_Lur_K2Luc1 * K_state.w_1 * Luc_state.w_2;
    ArrayXd temp4 = EVal_Lur_K1Luc2 * K_state.w_2 * Luc_state.w_1;
    ArrayXd temp5 = EVal_Lur_K2Luc2 * K_state.w_1 * Luc_state.w_1;

    ArrayXd EVal_Lur_vec = temp2+temp3+temp4+temp5;
    return EVal_Lur_vec;
}
//// EVal_Luc_mat
ArrayXd alias::cal_EVal_Luc_phi_interpolation(const ArrayXXd & EVal_Luc_mat, const PhiW & phi_w,
    const LKState & K_state, const LKState & Lur_state, const int & NLuc) {

//    cout << "EVal_Luc_mat.size() = " << EVal_Luc_mat.rows() << "; " << EVal_Luc_mat.cols() << endl;
    // Eigen::Map<const Eigen::ArrayXXd> EVprime_ELuc_mat(EVprime_ELuc.data(),para.N_Luc,para.N_phi*para.N_K*para.N_Lur);
    ArrayXXd EVal_Luc_phi_mat;
    if (phi_w.phi_grid_choice == 1) {
        EVal_Luc_phi_mat = EVal_Luc_mat(all,
            seqN(phi_w.phi_index*para.N_K*para.N_Lur,para.N_K*para.N_Lur));
    }
    else {
        EVal_Luc_phi_mat = phi_w.w1*EVal_Luc_mat(all,
            seqN(phi_w.phi_index*para.N_K*para.N_Lur,para.N_K*para.N_Lur))
            + phi_w.w2*EVal_Luc_mat(all,
            seqN((phi_w.phi_index+1)*para.N_K*para.N_Lur,para.N_K*para.N_Lur));
    }

    int Start_id_K1Lur1 = K_state.i_val_state * para.N_Lur + Lur_state.i_val_state;
    int Start_id_K2Lur1 = (K_state.i_val_state+1) * para.N_Lur + Lur_state.i_val_state;
    int Start_id_K1Lur2 = K_state.i_val_state * para.N_Lur + (Lur_state.i_val_state+1);
    int Start_id_K2Lur2 = (K_state.i_val_state+1) * para.N_Lur + (Lur_state.i_val_state+1);

    ArrayXd EVal_Luc_K1Lur1 = EVal_Luc_phi_mat.col(Start_id_K1Lur1);
    ArrayXd EVal_Luc_K2Lur1 = EVal_Luc_phi_mat.col(Start_id_K2Lur1);
    ArrayXd EVal_Luc_K1Lur2 = EVal_Luc_phi_mat.col(Start_id_K1Lur2);
    ArrayXd EVal_Luc_K2Lur2 = EVal_Luc_phi_mat.col(Start_id_K2Lur2);

    ArrayXd temp2 = EVal_Luc_K1Lur1 * K_state.w_2 * Lur_state.w_2;
    ArrayXd temp3 = EVal_Luc_K2Lur1 * K_state.w_1 * Lur_state.w_2;
    ArrayXd temp4 = EVal_Luc_K1Lur2 * K_state.w_2 * Lur_state.w_1;
    ArrayXd temp5 = EVal_Luc_K2Lur2 * K_state.w_1 * Lur_state.w_1;

    ArrayXd EVal_Luc_vec = temp2+temp3+temp4+temp5;
    return EVal_Luc_vec;
}

////**** calculation for a specific capital / labor employment  ****////
double alias::calResidual_Value_Exit_double(const ParaEst & para_est, const double & p_K, const double & K,
    const double & Lur, const double & Luc) {

    double F_Kbase = para_est.F_E_c_FK * log(1.0+p_K*K*para.deltaK);

//    double F_Lbase_ur = para_est.F_E_c_F_ur * log(1.0+Lur);

    double F_Lbase_ur = 0.0;
//    F_Lbase_ur = para_est.c_low_F_ur * log(1.0 + Lur);
    if (Lur < para.Lur_cutoff1) {
        F_Lbase_ur = para_est.c_low_F_ur * log(1.0 + Lur);
    }
    else {
        F_Lbase_ur = para_est.c_high_F_ur * log(1.0 + Lur);
    }

//    double F_Lbase_uc = para_est.F_E_c_F_uc * log(1.0+Luc);
    double F_Lbase_uc = para_est.c_F_uc * log(1.0 + Luc);

    double ResVal_E = - F_Lbase_ur - F_Lbase_uc - F_Kbase + p_K*para.deltaK*K;
//    cout << "ResVal_E  = " << ResVal_E << endl;
    return ResVal_E;
}

////**** Calculate the probability of Lur/Luc/K for each data point  ****////
//// the main function to calculate the probability of Lur/Luc/K
tuple<double,double,double,double,double,double,double,double,double,double,int,int,int> alias::cal_prob_KLurLuc_Data(
    const ParaEst & para_est, const ParaVec & para_vec, const EquStateVmat & Evalmat, const PhiW & phi_w,
    const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state, const int & PD_lag,
    const double & K_state_val, const double & Lur_state_val, const double & Luc_state_val,
    const int & K_missing, const int & Lur_missing, const int & Luc_missing) {

    double K_opt = 0.0;
    int K_index_max = 0;
    double lnpdf_State_K = 0.0;
    if (K_missing == 0 & PD_lag == 0) {
        ArrayXXd EVal_K_mat = Evalmat.EVal_PD_K_mat(seqN(PD_lag*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
        tuple<double,int> t_K = sol_opt_K(para_est,para_vec,phi_w,EVal_K_mat,Lur_state,Luc_state,K_state.val_state,0);
        K_opt = get<0>(t_K);
        K_index_max = get<1>(t_K);
        lnpdf_State_K = - log(para_est.sigma_Kerror) - 0.5*log(2.0*3.1415926)
            - 0.5 * pow( (log(K_state_val) - log(K_opt))/para_est.sigma_Kerror,2 );
//    cout << "lnpdf_State_K = " << lnpdf_State_K << endl;
    }

    double Lur_opt = 0.0;
    int Lur_index_max = 0;
    double lnpdf_State_Lur = 0.0;
    if (Lur_missing == 0 & PD_lag == 0) {
        ArrayXXd EVal_Lur_mat = Evalmat.EVal_PD_Lur_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
        tuple<double, int> t_Lur = sol_opt_Lur(para_est, para_vec, phi_w, EVal_Lur_mat, K_state, Luc_state,
            Lur_state.val_state,0);
        Lur_opt = get<0>(t_Lur);
        Lur_index_max = get<1>(t_Lur);
        lnpdf_State_Lur = -log(para_est.sigma_Lerror_ur) - 0.5 * log(2.0 * 3.1415926)
            - 0.5 * pow((log(Lur_state_val) - log(Lur_opt)) / para_est.sigma_Lerror_ur, 2);
        //    cout << "lnpdf_State_Lur = " << lnpdf_State_Lur << endl;
    }

    double Luc_opt = 0.0;
    int Luc_index_max = 0;
    double lnpdf_State_Luc = 0.0;
    if (Luc_missing == 0 & PD_lag == 0) {
        ArrayXXd EVal_Luc_mat = Evalmat.EVal_PD_Luc_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
        tuple<double, int> t_Luc = sol_opt_Luc(para_est, para_vec, phi_w, EVal_Luc_mat, K_state, Lur_state,
        Luc_state.val_state,0);
        Luc_opt = get<0>(t_Luc);
        Luc_index_max = get<1>(t_Luc);
        lnpdf_State_Luc = -log(para_est.sigma_Lerror_uc) - 0.5 * log(2.0 * 3.1415926)
            - 0.5 * pow((log(Luc_state_val) - log(Luc_opt)) / para_est.sigma_Lerror_uc, 2);
    }

//    cout << "K_state = " << K_state.val_state << "; Lur_state = " << Lur_state.val_state << "; Luc_state = " << Luc_state.val_state << endl;
    tuple<double,double,double,double,double,double> t_PD = calProductionDormancy_Prob_lognormal(para_est,para_vec,Evalmat,
        PD_lag,phi_w.phi,K_state.val_state,Lur_state.val_state,Luc_state.val_state);
    double Prob_P = get<0>(t_PD);
    double Prob_D = get<1>(t_PD);
    double lnProb_P = get<2>(t_PD);
    double lnProb_D = get<3>(t_PD);

    return tuple<double,double,double,double,double,double,double,double,double,double,int,int,int> (Prob_P,Prob_D,
        lnProb_P,lnProb_D,lnpdf_State_K,lnpdf_State_Lur,lnpdf_State_Luc,K_opt,Lur_opt,Luc_opt,
        K_index_max,Lur_index_max,Luc_index_max);
}

tuple<double,double,double,double,double,double,double> alias::cal_prob_KLurLuc_Data_Diff(
    const ParaEst & para_est, const ParaVec & para_vec, const EquStateVmat & Evalmat, const PhiW & phi_w,
    const LKState & K_state, const LKState & Lur_state, const LKState & Luc_state, const int & PD_lag,
    const double & K_state_val, const double & Lur_state_val, const double & Luc_state_val,
    const int & K_missing, const int & Lur_missing, const int & Luc_missing,
    const int & K_index_max, const int & Lur_index_max, const int & Luc_index_max) {

    double K_opt = 0.0;
    double lnpdf_State_K = 0.0;
    if (K_missing == 0 & PD_lag == 0) {
        ArrayXXd EVal_K_mat = Evalmat.EVal_PD_K_mat(seqN(PD_lag*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
        K_opt = sol_opt_K_Diff(para_est,para_vec,phi_w,EVal_K_mat,Lur_state,Luc_state,K_state.val_state,0,K_index_max);
        lnpdf_State_K = - log(para_est.sigma_Kerror) - 0.5*log(2.0*3.1415926)
            - 0.5 * pow( (log(K_state_val) - log(K_opt))/para_est.sigma_Kerror,2 );
//    cout << "lnpdf_State_K = " << lnpdf_State_K << endl;
    }

    double Lur_opt = 0.0;
    double lnpdf_State_Lur = 0.0;
    if (Lur_missing == 0 & PD_lag == 0) {
        ArrayXXd EVal_Lur_mat = Evalmat.EVal_PD_Lur_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
        Lur_opt = sol_opt_Lur_Diff(para_est,para_vec,phi_w,EVal_Lur_mat,K_state,Luc_state,Lur_state.val_state,0,Lur_index_max);
        lnpdf_State_Lur = -log(para_est.sigma_Lerror_ur) - 0.5 * log(2.0 * 3.1415926)
            - 0.5 * pow((log(Lur_state_val) - log(Lur_opt)) / para_est.sigma_Lerror_ur, 2);
        //    cout << "lnpdf_State_Lur = " << lnpdf_State_Lur << endl;
    }

    double Luc_opt = 0.0;
    double lnpdf_State_Luc = 0.0;
    if (Luc_missing == 0 & PD_lag == 0) {
        ArrayXXd EVal_Luc_mat = Evalmat.EVal_PD_Luc_mat(seqN(PD_lag * para.N_Luc, para.N_Luc), all);
        Luc_opt = sol_opt_Luc_Diff(para_est,para_vec,phi_w,EVal_Luc_mat,K_state,Lur_state,Luc_state.val_state,0,Luc_index_max);
        lnpdf_State_Luc = -log(para_est.sigma_Lerror_uc) - 0.5 * log(2.0 * 3.1415926)
            - 0.5 * pow((log(Luc_state_val) - log(Luc_opt)) / para_est.sigma_Lerror_uc, 2);
    }

//    cout << "K_state = " << K_state.val_state << "; Lur_state = " << Lur_state.val_state << "; Luc_state = " << Luc_state.val_state << endl;
    tuple<double,double,double,double,double,double> t_PD = calProductionDormancy_Prob_lognormal(para_est,para_vec,Evalmat,
        PD_lag,phi_w.phi,K_state.val_state,Lur_state.val_state,Luc_state.val_state);
    double Prob_P = get<0>(t_PD);
    double Prob_D = get<1>(t_PD);
    double lnProb_P = get<2>(t_PD);
    double lnProb_D = get<3>(t_PD);

    return tuple<double,double,double,double,double,double,double> (Prob_P,Prob_D,
        lnProb_P,lnProb_D,lnpdf_State_K,lnpdf_State_Lur,lnpdf_State_Luc);
}


//// solve the optimal Capital target given the states in the data
tuple<double,int> alias::sol_opt_K(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd EVal_K_mat, const LKState & Lur_state, const LKState & Luc_state, const double & K_state_val,
    const int & InitialPeriod) {

    ArrayXd EVal_K_vec = cal_EVal_K_phi_interpolation(EVal_K_mat, phi_w, Lur_state, Luc_state, para.N_K);

    double rstar = exp(0.5*pow(para_est.sigma_Kerror,2)) * para_vec.p_K * para_est.r_K;
    double Kbase = para.deltaK * K_state_val;

    double H_Kbase = 0.0;
    double F_Kbase = 0.0;
    if (InitialPeriod == 0) {
        H_Kbase = para_est.H_K / pow(Kbase, 2) + para_est.c_HK * log(1.0+Kbase) / pow(Kbase, 2)
            + para_est.H_K_Cons * para_vec.p_K * para_vec.p_K;
        F_Kbase = para_est.F_K / pow(Kbase, 2) + para_est.c_FK * log(1.0+Kbase) / pow(Kbase, 2)
            + para_est.F_K_Cons * para_vec.p_K * para_vec.p_K;
    }
    double c_K = 1.0;

    tuple<double,double,int> t_K = solve1D_K_LinearSpline(H_Kbase,F_Kbase,
        c_K,rstar,para_vec.p_K,Kbase,para_vec.vec_K,EVal_K_vec);
    double K = get<0>(t_K);
    double OptV = get<1>(t_K);
    int K_index_max = get<2>(t_K);

    return tuple<double,int>(K,K_index_max);
}

//// solve the optimal Capital target given the states in the data
double alias::sol_opt_K_Diff(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd EVal_K_mat, const LKState & Lur_state, const LKState & Luc_state, const double & K_state_val,
    const int & InitialPeriod, const int & K_index_max) {

    ArrayXd EVal_K_vec = cal_EVal_K_phi_interpolation(EVal_K_mat, phi_w, Lur_state, Luc_state, para.N_K);

    double rstar = exp(0.5*pow(para_est.sigma_Kerror,2)) * para_vec.p_K * para_est.r_K;
    double Kbase = para.deltaK * K_state_val;

    double H_Kbase = 0.0;
    double F_Kbase = 0.0;
    if (InitialPeriod == 0) {
        H_Kbase = para_est.H_K / pow(Kbase, 2) + para_est.c_HK * log(1.0+Kbase) / pow(Kbase, 2)
            + para_est.H_K_Cons * para_vec.p_K * para_vec.p_K;
        F_Kbase = para_est.F_K / pow(Kbase, 2) + para_est.c_FK * log(1.0+Kbase) / pow(Kbase, 2)
            + para_est.F_K_Cons * para_vec.p_K * para_vec.p_K;
    }
    double c_K = 1.0;

    double K = solve1D_K_LinearSpline_Diff(H_Kbase,F_Kbase,
        c_K,rstar,para_vec.p_K,Kbase,para_vec.vec_K,EVal_K_vec,K_index_max);

    return K;
}

//// solve the optimal Lur target given the states in the data
tuple<double,int> alias::sol_opt_Lur(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd & EVal_Lur_mat, const LKState & K_state, const LKState & Luc_state, const double & Lur_state_val,
    const int & InitialPeriod) {

    ArrayXd EVal_Lur_vec = cal_EVal_Lur_phi_interpolation(EVal_Lur_mat, phi_w, K_state, Luc_state, para.N_Lur);

    //// Value of production
    double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
    double Lbase = Lur_state_val;

    double H_Lbase_ur = 0.0;
    double F_Lbase_ur = 0.0;
    if (InitialPeriod == 0) {
        H_Lbase_ur =
                para_est.H_ur / pow(Lbase, 2) + para_est.c_H_ur * log(1.0+Lbase) / pow(Lbase, 2) + para_est.H_ur_Cons;
        F_Lbase_ur = 0;
        if (Lbase < para.Lur_cutoff1) {
            F_Lbase_ur = para_est.F_low_ur / pow(Lbase, 2) + para_est.c_low_F_ur * log(1.0+Lbase) / pow(Lbase, 2) +
                         para_est.F_low_ur_Cons;
        } else {
            F_Lbase_ur = para_est.F_high_ur / pow(Lbase, 2) + para_est.c_high_F_ur * log(1.0+Lbase) / pow(Lbase, 2) +
                         para_est.F_high_ur_Cons;
        }
    }

    tuple<double,double,int> t_Lur = solve1D_L_LinearSpline(H_Lbase_ur, F_Lbase_ur, wstar, Lbase,
        para_vec.vec_Lur, EVal_Lur_vec);
    double Lur = get<0>(t_Lur);
    double TotVOpt = get<1>(t_Lur);
    int Lur_index_max = get<2>(t_Lur);

    return tuple<double,int>(Lur,Lur_index_max);
}

//// solve the optimal Lur target given the states in the data
double alias::sol_opt_Lur_Diff(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd & EVal_Lur_mat, const LKState & K_state, const LKState & Luc_state, const double & Lur_state_val,
    const int & InitialPeriod, const int & Lur_index_max) {

    ArrayXd EVal_Lur_vec = cal_EVal_Lur_phi_interpolation(EVal_Lur_mat, phi_w, K_state, Luc_state, para.N_Lur);

    //// Value of production
    double wstar = para_est.w_ur * exp(0.5 * pow(para_est.sigma_Lerror_ur, 2));
    double Lbase = Lur_state_val;

    double H_Lbase_ur = 0.0;
    double F_Lbase_ur = 0.0;
    if (InitialPeriod == 0) {
        H_Lbase_ur =
                para_est.H_ur / pow(Lbase, 2) + para_est.c_H_ur * log(1.0+Lbase) / pow(Lbase, 2) + para_est.H_ur_Cons;
        F_Lbase_ur = 0;
        if (Lbase < para.Lur_cutoff1) {
            F_Lbase_ur = para_est.F_low_ur / pow(Lbase, 2) + para_est.c_low_F_ur * log(1.0+Lbase) / pow(Lbase, 2) +
                         para_est.F_low_ur_Cons;
        } else {
            F_Lbase_ur = para_est.F_high_ur / pow(Lbase, 2) + para_est.c_high_F_ur * log(1.0+Lbase) / pow(Lbase, 2) +
                         para_est.F_high_ur_Cons;
        }
    }

    double Lur = solve1D_L_LinearSpline_Diff(H_Lbase_ur, F_Lbase_ur, wstar, Lbase,
        para_vec.vec_Lur, EVal_Lur_vec, Lur_index_max);

    return Lur;
}


//// solve the optimal Luc target given the states in the data
tuple<double,int> alias::sol_opt_Luc(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd & EVal_Luc_mat, const LKState & K_state, const LKState & Lur_state, const double & Luc_state_val,
    const int & InitialPeriod) {

    //// solve Luc
//    cout << "solve Luc" << endl;
    ArrayXd EVal_Luc_vec= cal_EVal_Luc_phi_interpolation(EVal_Luc_mat, phi_w, K_state, Lur_state, para.N_Luc);

    double wstar = para_est.w_uc * exp(0.5 * pow(para_est.sigma_Lerror_uc, 2));
    double Lbase = Luc_state_val;

    double H_Lbase_uc = 0.0;
    double F_Lbase_uc = 0.0;
    if (InitialPeriod == 0) {
        H_Lbase_uc = para_est.H_uc / pow(Lbase, 2) + para_est.c_H_uc * log(Lbase+1.0) / pow(Lbase, 2) + para_est.H_uc_Cons;
        F_Lbase_uc = para_est.F_uc / pow(Lbase, 2) + para_est.c_F_uc * log(Lbase+1.0) / pow(Lbase, 2) + para_est.F_uc_Cons;
    }

    tuple<double,double,int> t_Luc = solve1D_L_LinearSpline(H_Lbase_uc, F_Lbase_uc, wstar, Lbase,
        para_vec.vec_Luc, EVal_Luc_vec);
    double Luc = get<0>(t_Luc);
    double TotVOpt = get<1>(t_Luc);
    int Luc_index_max = get<2>(t_Luc);

    return tuple<double,int>(Luc,Luc_index_max);
}

//// solve the optimal Luc target given the states in the data
double alias::sol_opt_Luc_Diff(const ParaEst & para_est,const ParaVec & para_vec,const PhiW & phi_w,
    const ArrayXXd & EVal_Luc_mat, const LKState & K_state, const LKState & Lur_state, const double & Luc_state_val,
    const int & InitialPeriod, const int & Luc_index_max) {

    //// solve Luc
    //    cout << "solve Luc" << endl;
    ArrayXd EVal_Luc_vec= cal_EVal_Luc_phi_interpolation(EVal_Luc_mat, phi_w, K_state, Lur_state, para.N_Luc);

    double wstar = para_est.w_uc * exp(0.5 * pow(para_est.sigma_Lerror_uc, 2));
    double Lbase = Luc_state_val;

    double H_Lbase_uc = 0.0;
    double F_Lbase_uc = 0.0;
    if (InitialPeriod == 0) {
        H_Lbase_uc = para_est.H_uc / pow(Lbase, 2) + para_est.c_H_uc * log(Lbase+1.0) / pow(Lbase, 2) + para_est.H_uc_Cons;
        F_Lbase_uc = para_est.F_uc / pow(Lbase, 2) + para_est.c_F_uc * log(Lbase+1.0) / pow(Lbase, 2) + para_est.F_uc_Cons;
    }

    double Luc = solve1D_L_LinearSpline_Diff(H_Lbase_uc, F_Lbase_uc, wstar, Lbase,
        para_vec.vec_Luc, EVal_Luc_vec,Luc_index_max);

    return Luc;
}

//// Calculate the exit/stay probability for each data point ////
tuple<double,double,double,double,double,double> alias::calProductionDormancy_Prob_lognormal(const ParaEst & para_est,
    const ParaVec & para_vec, const EquStateVmat & Evalmat, const int & PD_lag, const double & phi,
    const double & K_val_state, const double & Lur_val_state, const double & Luc_val_state) {

    PhiW phi_w = definePhiWeights(para_vec.vec_phi, phi);
    LKState K_state = defineLKState(para_vec.vec_K, K_val_state);
    LKState Lur_state = defineLKState(para_vec.vec_Lur, Lur_val_state);
    LKState Luc_state = defineLKState(para_vec.vec_Luc, Luc_val_state);

    int PD = 0;
    ArrayXXd EVal_SE_mat_P = Evalmat.EVal_PD_SE_mat(seqN(PD*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
    ArrayXd EVal_K_Interp_P = cal_EVal_K_phi_interpolation(EVal_SE_mat_P, phi_w,
        Lur_state, Luc_state, para.N_K);
    double beta_P = (EVal_K_Interp_P(K_state.i_val_state+1) - EVal_K_Interp_P(K_state.i_val_state))
                  / (para_vec.vec_K(K_state.i_val_state+1) - para_vec.vec_K(K_state.i_val_state));
    double Vprime_P = beta_P * (K_state.val_state - para_vec.vec_K(K_state.i_val_state))
                  + EVal_K_Interp_P(K_state.i_val_state);

    PD = 1;
    ArrayXXd EVal_SE_mat_D = Evalmat.EVal_PD_SE_mat(seqN(PD*para.N_Luc*para.N_Lur,para.N_Luc*para.N_Lur),all);
    ArrayXd EVal_K_Interp_D = cal_EVal_K_phi_interpolation(EVal_SE_mat_D, phi_w,
        Lur_state, Luc_state, para.N_K);
    double beta_D = (EVal_K_Interp_D(K_state.i_val_state+1) - EVal_K_Interp_D(K_state.i_val_state))
                    / (para_vec.vec_K(K_state.i_val_state+1) - para_vec.vec_K(K_state.i_val_state));
    double Vprime_D = beta_D * (K_state.val_state - para_vec.vec_K(K_state.i_val_state))
                    + EVal_K_Interp_D(K_state.i_val_state);

    double Revenue = RevOpt1_sim(para_est, phi, K_val_state, Lur_val_state,Luc_val_state);

    double F_P; double sigma_PD;
    if (PD_lag == 0) {
        F_P = para_est.F_P_PP;
        sigma_PD = para_est.sigma_PD;
    } else {
        F_P = para_est.F_P_DP;
        sigma_PD = para_est.sigma_PD;
    }

    double diff_PD = Vprime_P - Vprime_D;
    if (Vprime_P < Vprime_D) {diff_PD = 0.0;}
    double dEVal_PD_sigma = (log(Revenue + diff_PD) - F_P) / sigma_PD;

    dEVal_PD_sigma = min(dEVal_PD_sigma,7.5);
    dEVal_PD_sigma = max(dEVal_PD_sigma,-7.5);

    double Prob_P = normCDF(dEVal_PD_sigma, 0, 1);
    double Prob_D = normCDF(-dEVal_PD_sigma, 0, 1);

    double lnProb_P = log(Prob_P);
    double lnProb_D = log(Prob_D);

    if (isfinite(lnProb_P) == 0) {
        cout << "Revenue = " << Revenue << "; Vprime_P = " << Vprime_P << "; Vprime_D = " << Vprime_D
            << "; Vprime_P - Vprime_D = " << Vprime_P - Vprime_D
            << "; Revenue + Vprime_P - Vprime_D = " << Revenue + Vprime_P - Vprime_D << endl;
        cout << "lnProb_P = " << lnProb_P << "; Prob_P = " << Prob_P << endl;
        throw runtime_error("1338");
    }

    double CutoffVal = Revenue + diff_PD;
    return tuple<double,double,double,double,double,double>(Prob_P,Prob_D,lnProb_P,lnProb_D,Revenue,CutoffVal);
}

//// Revenue given the states
double alias::RevOpt1_sim(const ParaEst & para_est, const double & phi, const double & K, const double & L_ur,
    const double & L_uc) {

    /*** CES production function: Cobb Douglas ***/
    double L = pow(L_uc, para_est.alpha_Lc) * pow(L_ur, para_est.alpha_Lr);
    double KL = pow(K, para_est.alpha_K) * pow(L, para_est.alpha_L);

    double temp2_phi = exp(para_est.gamma0 + pow(para_est.sigma_phi_eps,2)/2.0);
    double temp3_phi = pow(phi,para_est.gamma1);

    double Rev = temp2_phi*temp3_phi * pow(KL, (para.sigma-1.0)/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) )
        * pow( para_est.PI,1.0/para.sigma / (1.0 - para_est.alpha_M*(para.sigma-1.0)/para.sigma) );
    //
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


void alias::PrintResultEstimationAuxiliaryModel(const std::vector<double> & x,const EquStateV & EquV) {

    cout << "******* Prob_E = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV.ProbPD_E(i_phi*para.N_KL+0*para.N_Lur*para.N_Luc + 0*para.N_Lur + 0) << ", "; }
    cout << ";" << endl;
    cout << "******* Prob_D = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV.ProbPD_D(i_phi*para.N_KL+0*para.N_Lur*para.N_Luc + 0*para.N_Lur + 0) << ", "; }
    cout << ";" << endl;
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV.ProbPD_D(para.N + i_phi*para.N_KL+0*para.N_Lur*para.N_Luc + 0*para.N_Lur + 0) << ", "; }
    cout << ";" << endl;
    cout << "******* Eval = ";
    for (size_t i_phi = 0; i_phi < para.N_phi; ++i_phi) { cout << EquV.EVal_PD(i_phi*para.N_KL) << ", "; }
    cout << ";" << endl;
//
//    cout << "******* K max = " << est_array.Est_KOpt.maxCoeff()
//         << "; min = " << est_array.Est_KOpt.minCoeff()
//         << "; mean = " << ( est_array.Est_KOpt * RealData.State_P_nonMiss(all,seqN(0,para.EstTbar)).cast<double>() ).sum()
//         / double(RealData.State_P_nonMiss.sum()) << endl;
//    cout << "******* States Lur max = " << est_array.Est_LurOpt.maxCoeff()
//         << "; min = " << est_array.Est_LurOpt.minCoeff()
//         << "; mean = " << ( est_array.Est_LurOpt * RealData.State_P_nonMiss(all,seqN(0,para.EstTbar)).cast<double>() ).sum()
//         / double(RealData.State_P_nonMiss.sum()) << endl;
//    cout << "******* Luc max = " << est_array.Est_LucOpt.maxCoeff()
//         << "; min = " << est_array.Est_LucOpt.minCoeff()
//         << "; mean = " <<  ( est_array.Est_LucOpt * RealData.State_P_nonMiss(all,seqN(0,para.EstTbar)).cast<double>() ).sum()
//        / double(RealData.State_P_nonMiss.sum()) << endl;
}


//// Store the model predicted optimal capital/labor for each data point
//EstArray alias::initializeEstArray(const size_t & N) {
//
//    EstArray est_array;
//    est_array.Est_KOpt = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_LurOpt = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_LucOpt = ArrayXXd::Zero(N,para.EstTbar);
//
//    est_array.Est_lnpdf_KOpt = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_lnpdf_LurOpt = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_lnpdf_LucOpt = ArrayXXd::Zero(N,para.EstTbar);
//
//    est_array.Est_lnProbP = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_lnProbD = ArrayXXd::Zero(N,para.EstTbar);
//
//    est_array.Est_lnProbS = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_lnProbE = ArrayXXd::Zero(N,para.EstTbar);
//
//    est_array.Est_ProbP = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_ProbD = ArrayXXd::Zero(N,para.EstTbar);
//
//    est_array.Est_ProbS = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_ProbE = ArrayXXd::Zero(N,para.EstTbar);
//
//    est_array.Est_Prob_E_Miss = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Est_Prob_DP_Miss = ArrayXXd::Zero(N,para.EstTbar);
//    est_array.Prob_Miss_DP = ArrayXXd::Zero(N,para.EstTbar);
//
//    return est_array;
//}
//
///**************************************************************
//* Version 1: Randomly choose inital value and see the estimation
//**************************************************************/
//tuple<ArrayXd,ArrayXd,ArrayXd> alias::EstimationAuxiliaryModel_Output_Version1(const SimData & ReadData_GoodState_LaborInt,
//    const SimData & ReadData_BadState_LaborInt, const SimData & ReadData_GoodState_CapitalInt,
//    const SimData & ReadData_BadState_CapitalInt, MultiThreads::Threads_Management & threadsManagement) {
//
//    ArrayXd theta1_a; ArrayXd theta2_a; ArrayXd theta3_a;
//    // Generate 10 initial guesses
//    ArrayXXd theta_RandomInitial_select = GenerateInitialGuessParaRandom();
//    // Estimate the model with 10 different initial guess
//    for (size_t col = 0; col < 10; ++col) {
//        // set up the parameters: group them into three parts
//        tuple<ArrayXd,ArrayXd,ArrayXd> t_AuxPara = SetupInitialParaGuess_forEstimation(
//                theta_RandomInitial_select.col(col));
//        ArrayXd thetaData_1_ini = get<0>(t_AuxPara);
//        ArrayXd thetaData_2_ini = get<1>(t_AuxPara);
//        ArrayXd thetaData_3_ini = get<2>(t_AuxPara);
//        ArrayXd thetaData_ini(para.dim1 + para.dim2 + para.dim3);
//        thetaData_ini << thetaData_1_ini, thetaData_3_ini;
//        cout << "col = " << col << "; The Initial Parameter Guess for the auxilliary model:" << endl;
//        cout << "Group 1 parameters: thetaData_1_ini = " << thetaData_1_ini.transpose() << endl;
//        cout << "Group 2&3 parameters: thetaData_3_ini = " << thetaData_3_ini.transpose() << endl;
//
//        // Start the estimation of the auxiliary model: theta1_a and theta3_a are the estimates
//        tuple<ArrayXd,ArrayXd,ArrayXd,double> t_AuxillaryEst = EstimationAuxiliaryModel(
//            ReadData_GoodState_LaborInt,ReadData_BadState_LaborInt,
//            ReadData_GoodState_CapitalInt, ReadData_BadState_CapitalInt,
//            thetaData_1_ini, thetaData_2_ini, thetaData_3_ini,threadsManagement);
//        theta1_a = get<0>(t_AuxillaryEst); theta2_a = get<1>(t_AuxillaryEst); theta3_a = get<2>(t_AuxillaryEst);
//        double obj3_a = get<3>(t_AuxillaryEst);
//        ArrayXd theta_a(para.dim1 + para.dim2 + para.dim3);
//        theta_a << theta1_a, theta2_a, theta3_a; // the full set of parameter estimates from the auxiliary model
//        cout << "The auxilliary model is estimated:" << endl;
//        cout << "Group 1 parameters: theta1_a = " << theta1_a.transpose() << endl;
//        cout << "Group 2 parameters: theta3_a = " << theta2_a.transpose() << endl;
//        cout << "Group 3 parameters: theta3_a = " << theta3_a.transpose() << endl;
////
////        // print relevant esimation results to file
////        std::vector<double> theta3_a_a(para.dim3);
////        for (size_t k = 0; k < para.dim3; ++k) { theta3_a_a[k] = theta3_a(k); }
////        int GoodState = 1; int LaborIntensive = 1;
////        tuple<ParaEst,ParaVec,EquStateV,EquStateVmat> t_Equ = EstimationAuxiliaryModel_part3_ObjFun_solveEquV(
////            theta3_a_a, theta1_a, theta2_a, GoodState, LaborIntensive, threadsManagement);
////        EquStateV EquV_a = get<2>(t_Equ);
////        writeToCSVfile("theta_a_Prob_E_State_" + to_string(GoodState) + "_Industry_" + to_string(LaborIntensive)
////            + "_col_" + to_string(col) + ".csv", EquV_a.Prob_E.cast<double>().matrix());
////
////        writeToCSVfile("theta_a_ini_State_" + to_string(GoodState) + "_Industry_" + to_string(LaborIntensive)
////            + "_col_" + to_string(col) + ".csv", thetaData_ini.cast<double>().matrix());
////        writeToCSVfile("theta_a_State_" + to_string(GoodState) + "_Industry_" + to_string(LaborIntensive)
////            + "_col_" + to_string(col) + ".csv", theta_a.cast<double>().matrix());
////        ArrayXd obj3_a_vec(1);
////        obj3_a_vec(0) = obj3_a;
////        writeToCSVfile("obj3_a_State_" + to_string(GoodState) + "_Industry_" + to_string(LaborIntensive)
////                       + "_col_" + to_string(col) + ".csv", obj3_a_vec.cast<double>().matrix());
//    }
//    return tuple<ArrayXd,ArrayXd,ArrayXd>(theta1_a,theta2_a,theta3_a);
//}
//

//
//////void alias::EstimationAuxiliaryModel_part2_theta2_checking(const std::vector<double> & theta2,const ArrayXd & theta1,
//////    const SimData & sim_data, const int & GoodState, const int & LaborIntensive) {
//////
//////    double alpha_tilde_L = theta1(3);
//////    double alpha_tilde_K = (1-para.alpha_M)*(para.sigma-1.0)/para.sigma / (1.0-para.alpha_M*(para.sigma-1)/para.sigma)
//////        - alpha_tilde_L;
//////    double alpha_Lr = theta1(4);
//////    double alpha_Lc = 1 - alpha_Lr;
//////
//////    //// estimated productivity
//////    ArrayXXd logL = alpha_Lc * log(sim_data.Employ_uc) + alpha_Lr * log(sim_data.Employ_ur);
//////    ArrayXXd lnphi_a = log(sim_data.Revenue) - alpha_tilde_K*log(sim_data.Capital)-alpha_tilde_L*logL;
//////
//////    std::vector<double> x_lb(para.dim2);
//////    x_lb[0] = -1; x_lb[1] = -0.99;
//////
//////    std::vector<double> x_ub(para.dim2);
//////    x_ub[0] = 1; x_ub[1] = 0.999;
//////
//////    int Ntemp = 11;
//////    ArrayXXd Obj_vec(Ntemp,para.dim2);
//////    ArrayXXd dObj_dx_vec(Ntemp,para.dim2);
//////    ArrayXXd dObj_dx_vec_sim(Ntemp,para.dim2);
//////    ArrayXXd ParaSpace_vec(Ntemp,para.dim2);
//////    for (size_t i = 0; i < para.dim2; ++i) {
//////        ArrayXd tempvec0 = ArrayXd::LinSpaced(Ntemp,theta2[i]*0.5,theta2[i]*5);
//////        for (size_t n = 0; n < Ntemp; ++n) {
//////            std::vector<double> x1_vec = theta2;
//////            x1_vec[i] = tempvec0(n);
//////            cout << "i = " << i << "; n = " << n << "; x1_vec= " << x1_vec[0] << "; " << x1_vec[1] << endl;
//////            ParaSpace_vec(n,i) = tempvec0(n);
//////
//////            tuple<double,std::vector<double>> t = EstimationAuxiliaryModel_part2_ObjFun(x1_vec, theta1,
//////                                                                                        sim_data,lnphi_a,GoodState,LaborIntensive);
//////            Obj_vec(n,i) = get<0>(t);
//////            std::vector<double> grad = get<1>(t);
//////            dObj_dx_vec(n,i) = grad[i];
//////
//////            std::vector<double> x1_vec_plus = x1_vec;
//////            x1_vec_plus[i] = x1_vec[i] + 1e-6;
//////            tuple<double,std::vector<double>> t_plus = EstimationAuxiliaryModel_part2_ObjFun(x1_vec_plus,theta1,
//////                                                                                             sim_data,lnphi_a,GoodState,LaborIntensive);
//////            double Obj_vec_plus = get<0>(t_plus);
//////
//////            std::vector<double> x1_vec_minus = x1_vec;
//////            x1_vec_minus[i] = x1_vec[i] - 1e-6;
//////            tuple<double,std::vector<double>> t_minus = EstimationAuxiliaryModel_part2_ObjFun(x1_vec_minus,theta1,
//////                                                                                              sim_data,lnphi_a,GoodState,LaborIntensive);
//////            double Obj_vec_minus = get<0>(t_minus);
//////
//////            dObj_dx_vec_sim(n,i) = (Obj_vec_plus - Obj_vec_minus) / (x1_vec_plus[i] - x1_vec_minus[i]);
//////
//////            cout << "***** Obj_vec(n,i) = " << Obj_vec(n,i) << "; dObj_vec(n,i) = " << dObj_dx_vec(n,i)
//////                 << "; dObj_dx_vec_sim(n,i) = " << dObj_dx_vec_sim(n,i) << endl;
////////            throw runtime_error("1091");
//////        }
//////    }
//////
//////    writeToCSVfile("Obj_vec.csv", Obj_vec.matrix());
//////    writeToCSVfile("dObj_dx_vec.csv", dObj_dx_vec.matrix());
//////    writeToCSVfile("ParaSpace_vec.csv", ParaSpace_vec.matrix());
//////    throw runtime_error("843");
//////}
//////
//////void alias::EstimationAuxiliaryModel_part3_theta3_checking(const ArrayXd & theta3, const SimData & RealData,
//////    const ArrayXd & theta1_a, const ArrayXd & theta2_a, const int & GoodState, const int & LaborIntensive,
//////    MultiThreads::Threads_Management & threadsManagement) {
//////
//////    cout << "theta3 = " << theta3.transpose() << endl;
//////    std::vector<double> x0_vec(para.dim3);
//////    x0_vec[0] = theta3(0); x0_vec[1] = theta3(1); x0_vec[2] = theta3(2); x0_vec[3] = theta3(3);
//////    x0_vec[4] = theta3(4); x0_vec[5] = theta3(5); x0_vec[6] = theta3(6); x0_vec[7] = theta3(7);
//////    x0_vec[8] = theta3(8); x0_vec[9] = theta3(9); x0_vec[10] = theta3(10); x0_vec[11] = theta3(11);
//////    x0_vec[12] = theta3(12); x0_vec[13] = theta3(13); x0_vec[14] = theta3(14); x0_vec[15] = theta3(15);
//////    x0_vec[16] = theta3(16); x0_vec[17] = theta3(17);
//////
//////    int Ntemp = 5;
//////    ArrayXXd Obj_vec(Ntemp,para.dim3);
//////    ArrayXXd dObj_dx_vec_sim(Ntemp,para.dim3);
//////    ArrayXXd dObj_dx_vec(Ntemp,para.dim3);
//////    ArrayXXd ParaSpace_vec(Ntemp,para.dim3);
//////
//////    for (size_t i = 0; i < para.dim3; ++i) {
//////        ArrayXd tempvec0 = ArrayXd::LinSpaced(Ntemp,x0_vec[i]*0.5,x0_vec[i]*1.5);
//////        for (size_t n = 0; n < Ntemp; ++n) {
//////            std::vector<double> x1_vec = x0_vec;
//////            x1_vec[i] = tempvec0(n);
//////            cout << "**** i = " << i << "; n = " << n << "; x1_vec= " << x1_vec[0] << "; " << x1_vec[1]
//////                 << "; " << x1_vec[2] << "; " << x1_vec[3] << "; " << x1_vec[4] << "; " << x1_vec[5]
//////                 << "; " << x1_vec[6] << "; " << x1_vec[7] << "; " << x1_vec[8] << "; " << x1_vec[9]
//////                 << "; " << x1_vec[10] << "; " << x1_vec[11] << "; " << x1_vec[12] << "; " << x1_vec[13]
//////                 << "; " << x1_vec[14] << "; " << x1_vec[15] << "; " << x1_vec[16] << "; " << x1_vec[17] << endl;
//////
//////            tuple<double,ArrayXXd,ArrayXXi,ArrayXXi,ArrayXXi,EstArray,EquStateV,DiffEVal,ArrayXd>
//////                    t_obj = EstimationAuxiliaryModel_part3_ObjFun(x1_vec,RealData, theta1_a, theta2_a, GoodState,
//////                                                                  LaborIntensive, threadsManagement);
//////            double Obj = get<0>(t_obj);
//////            Obj_vec(n,i) = Obj;
//////            ArrayXXd likelihood = get<1>(t_obj);
//////            ArrayXXi K_index_max_mat = get<2>(t_obj);
//////            ArrayXXi Lur_index_max_mat = get<3>(t_obj);
//////            ArrayXXi Luc_index_max_mat = get<4>(t_obj);
//////            EstArray est_array = get<5>(t_obj);
//////            EquStateV EquV = get<6>(t_obj);
//////            DiffEVal diff_Eval = get<7>(t_obj);
//////            ArrayXd obj_grad = get<8>(t_obj);
//////            dObj_dx_vec_sim(n,i) = obj_grad(i);
//////
//////            cout << "max = " << est_array.Est_KOpt.maxCoeff() << "; min = " << est_array.Est_KOpt.minCoeff()
//////                 << "; mean = " << ( est_array.Est_KOpt * RealData.State_P_nonMiss.cast<double>() ).sum()
//////                                   / double(RealData.State_P_nonMiss.sum()) << endl;
//////            cout << "max = " << est_array.Est_LurOpt.maxCoeff() << "; min = " << est_array.Est_LurOpt.minCoeff()
//////                 << "; mean = " << ( est_array.Est_LurOpt * RealData.State_P_nonMiss.cast<double>() ).sum()
//////                                   / double(RealData.State_P_nonMiss.sum()) << endl;
//////            cout << "max = " << est_array.Est_LucOpt.maxCoeff() << "; min = " << est_array.Est_LucOpt.minCoeff()
//////                 << "; mean = " <<  ( est_array.Est_LucOpt * RealData.State_P_nonMiss.cast<double>() ).sum()
//////                                    / double(RealData.State_P_nonMiss.sum()) << endl;
////////////////////
//////            std::vector<double> x1_plus = x1_vec;
//////            x1_plus[i] = x1_vec[i] + 1e-8;
//////            tuple<double,ArrayXXd,EstArray,EquStateV,ArrayXd> t_obj_plus = EstimationAuxiliaryModel_part3_ObjFun_Diff(
//////                    x1_plus, RealData, theta1_a, theta2_a, EquV, diff_Eval,
//////                    K_index_max_mat, Lur_index_max_mat, Luc_index_max_mat, GoodState, LaborIntensive, threadsManagement);
//////            double Obj_plus = get<0>(t_obj_plus);
//////            ArrayXXd likelihood_plus = get<1>(t_obj_plus);
//////            EstArray est_array_plus = get<2>(t_obj_plus);
//////            EquStateV EquV_plus = get<3>(t_obj_plus);
//////
////////////////////////
//////            std::vector<double> x1_minus = x1_vec;
//////            double Obj_minus = Obj;
//////            EstArray est_array_minus = est_array;
//////            EquStateV EquV_minus = EquV;
//////
//////////////////////
//////            cout << "est_array_plus.Est_KOpt = " << (abs(est_array_plus.Est_KOpt - est_array_minus.Est_KOpt)).maxCoeff() << endl;
//////            cout << "est_array_plus.Est_LurOpt = " << (abs(est_array_plus.Est_LurOpt - est_array_minus.Est_LurOpt)).maxCoeff() << endl;
//////            cout << "est_array_plus.Est_LucOpt = " << (abs(est_array_plus.Est_LucOpt - est_array_minus.Est_LucOpt)).maxCoeff() << endl;
//////
//////            cout << "est_array_plus.Est_lnpdf_KOpt = " << abs(est_array_plus.Est_lnpdf_KOpt - est_array_minus.Est_lnpdf_KOpt).maxCoeff() << endl;
//////            cout << "est_array_plus.Est_lnpdf_LurOpt = " << abs(est_array_plus.Est_lnpdf_LurOpt - est_array_minus.Est_lnpdf_LurOpt).maxCoeff() << endl;
//////            cout << "est_array_plus.Est_lnpdf_LucOpt = " << abs(est_array_plus.Est_lnpdf_LucOpt - est_array_minus.Est_lnpdf_LucOpt).maxCoeff() << endl;
//////
//////            cout << "est_array_plus.Est_lnProbS = " << abs(est_array_plus.Est_lnProbS - est_array_minus.Est_lnProbS).maxCoeff() << endl;
//////            cout << "est_array_plus.Est_lnProbE = " << abs(est_array_plus.Est_lnProbE - est_array_minus.Est_lnProbE).maxCoeff() << endl;
//////
//////            cout << "est_array_plus.Est_lnProbP = " << abs(est_array_plus.Est_lnProbP - est_array_minus.Est_lnProbP).maxCoeff() << endl;
//////            cout << "est_array_plus.Est_lnProbD = " << abs(est_array_plus.Est_lnProbD - est_array_minus.Est_lnProbD).maxCoeff() << endl;
//////
//////            cout << "diff in EquV = " << (abs(EquV_plus.EVal - EquV_minus.EVal)/abs(EquV_minus.EVal)).maxCoeff() << endl;
//////            cout << "diff in EquV.OptK = " << (abs(EquV_plus.OptK - EquV_minus.OptK)/abs(EquV_minus.OptK)).maxCoeff() << endl;
//////            cout << "diff in EquV.OptLur = " << (abs(EquV_plus.OptLur - EquV_minus.OptLur)/abs(EquV_minus.OptLur)).maxCoeff() << endl;
//////            cout << "diff in EquV.OptLuc = " << (abs(EquV_plus.OptLuc - EquV_minus.OptLuc)/abs(EquV_minus.OptLuc)).maxCoeff() << endl;
//////
//////            dObj_dx_vec(n,i) = (Obj_plus - Obj_minus)/(x1_plus[i] - x1_minus[i]);
////////            dObj_dx_vec(n,i) = (Obj_plus - Obj_minus);
//////
//////            cout << "i = " << i << "; n = " << n << "; Obj_vec(n,i) = " << Obj_vec(n,i)
//////                 << "; dObj_dx_vec(n,i) = " << dObj_dx_vec(n,i)
//////                 << "; dObj_dx_vec_sim(n,i) = " << dObj_dx_vec_sim(n,i)<< endl;
////////            throw runtime_error("1358");
//////        }
//////        throw runtime_error("1358");
//////        ParaSpace_vec.col(i) = tempvec0;
//////    }
////////    throw runtime_error("1358");
//////
//////    writeToCSVfile("Obj_vec.csv", Obj_vec.matrix());
//////    writeToCSVfile("dObj_dx_vec.csv", dObj_dx_vec.matrix());
//////    writeToCSVfile("ParaSpace_vec.csv", ParaSpace_vec.matrix());
//////    throw runtime_error("886");
//////}
//////
//
////// with lognormal distribution assumption
////tuple<double,double,double,double> alias::calStayExit_Prob_lognormal(const ParaEst & para_est, const ParaVec & para_vec,
////                                                                     const EquStateVmat & Evalmat, const PhiW & phi_w, const LKState & K_state, const LKState & Lur_state,
////                                                                     const LKState & Luc_state) {
////
////    //// choose Stay or Exit State
////    ArrayXd EVal_K_Interp = cal_EVal_K_phi_interpolation(Evalmat.EVal_SE_mat, phi_w, Lur_state, Luc_state, para.N_K);
////
////    double beta = (EVal_K_Interp(K_state.i_val_state+1) - EVal_K_Interp(K_state.i_val_state))
////                  / (para_vec.vec_K(K_state.i_val_state+1) - para_vec.vec_K(K_state.i_val_state));
////    double EVal_S = beta * (K_state.val_state - para_vec.vec_K(K_state.i_val_state)) + EVal_K_Interp(K_state.i_val_state);
////
////    double ResV_E = calResidual_Value_Exit_double(para_est, para_vec.p_K, K_state.val_state,
////                                                  Lur_state.val_state,Luc_state.val_state);
////
////    double dEVal_SE_sigma = (EVal_S - ResV_E);
////    if (isfinite(dEVal_SE_sigma) == 0) {
////        cout << "EVal_S = " << EVal_S << "; ResV_E = " << ResV_E << "; dEVal_SE_sigma = " << dEVal_SE_sigma
////             << " ; para_est.F_E = " << para_est.F_E << endl;
////        throw runtime_error("462");
////    }
////    double Prob_S; double Prob_E;
////    if (dEVal_SE_sigma <= 0) {
////        Prob_S = 0.0; Prob_E = 1.0;
////    }
////    else {
////        Prob_S = normCDF( (log(dEVal_SE_sigma) - para_est.F_E)/para_est.sigma_SE,0.0,1.0 );
////        Prob_E = normCDF( -(log(dEVal_SE_sigma) - para_est.F_E)/para_est.sigma_SE,0.0,1.0 );
////    }
////    double lnProb_S;
////    if (Prob_S == 0) { lnProb_S = log(normCDF( -7.5,0.0,1.0 )); }
////    else { lnProb_S = log(Prob_S); }
////
////    double lnProb_E;
////    if (Prob_E == 0) { lnProb_E = log(normCDF( -7.5,0.0,1.0 )); }
////    else { lnProb_E = log(Prob_E); }
////
////    return tuple<double,double,double,double>(Prob_S,Prob_E,lnProb_S,lnProb_E);
////}