#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <map>
#include "multithread_loop.h"

#include "Basics.h"
#include "SimulationData.h"
#include "EstimationAuxiliary.h"
#include "Estimation_IndirectInference.h"

namespace alias {
    using namespace std;
    using namespace alias;
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;
    using MatrixXd = Eigen::MatrixXd;

    ArrayXd EstimationIndirectInference_StandardDeviation(const ArrayXd & theta_Est, const ArrayXd & theta_a,
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);


    /*** Calculate the JMatrix ***/
    ArrayXXd EstimationIndirectInference_SD_JMatrix_RealData(const ArrayXd & theta,
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    tuple<ArrayXXd,int,int,int,int,int> EstimationIndirectInference_SD_dl_dtheta_Part1_RealData(const ArrayXd & theta,
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);
    tuple<ArrayXXd,int,int,int,int> EstimationIndirectInference_SD_dl_dtheta_Part2_RealData(const ArrayXd & theta,
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    // ArrayXd EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(const ArrayXd & theta_a,
    //     const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    //     const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    //     MultiThreads::Threads_Management &threadsManagement);
    tuple<ArrayXXd,int,int> EstimationIndirectInference_SD_dl_dtheta_Part3_RealData(
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        const EquV_Auxiliary & EquV_GoodState_LaborInt, const EquV_Auxiliary & EquV_BadState_LaborInt,
        const EquV_Auxiliary & EquV_GoodState_CapitalInt, const EquV_Auxiliary & EquV_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    ArrayXXd EstimationIndirectInference_SD_l_Part3_RealData(const ArrayXd & theta_a,
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    tuple<ArrayXXd,int> Cal_logL_Part3_RealData(const SimData & ReadData, const ArrayXd & theta_a,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement);
    tuple<ArrayXXd,int> Cal_dlogL_Part3_RealData(const SimData & ReadData, const EquV_Auxiliary & EquV,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement);

    ArrayXXd EstimationIndirectInference_SD_IMatrix_RealData(const ArrayXd & theta,
        const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
        const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    ArrayXXd EstimationIndirectInference_SD_IStarMatrix_SimData(const ArrayXd & theta_Est,const ArrayXd & theta_a,
        const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
        const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
        const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
        const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    tuple<ArrayXXd,ArrayXXd,double> EstimationIndirectInference_SD_dl_dtheta_Part1_SimData(const ArrayXd & theta,
        const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
        const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);
    tuple<ArrayXXd,ArrayXXd,double> EstimationIndirectInference_SD_dl_dtheta_Part2_SimData(const ArrayXd & theta,
        const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
        const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);
    tuple<ArrayXXd,ArrayXXd,double> EstimationIndirectInference_SD_dl_dtheta_Part3_SimData(
        const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
        const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
        const EquV_Auxiliary & EquV_GoodState_LaborInt, const EquV_Auxiliary & EquV_BadState_LaborInt,
        const EquV_Auxiliary & EquV_GoodState_CapitalInt, const EquV_Auxiliary & EquV_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);
    tuple<ArrayXXd,ArrayXXd,double> Cal_dlogL_Part3_SimData(const SimData & SimData, const EquV_Auxiliary & EquV,
        const int & GoodState, const int & LaborIntensive, MultiThreads::Threads_Management &threadsManagement);



    ArrayXXd Cal_db_dtheta_a(const ArrayXd & theta_Est,const ArrayXd & theta_a,
        const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
        const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
        const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
        const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);


    ArrayXXd Cal_db_dtheta0_v2(const ArrayXd & theta_Est,const ArrayXd & theta_a,
        const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
        const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
        const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
        const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    tuple<ArrayXXd,ArrayXXd,ArrayXXd> Cal_logL_dtheta_Est_theta_a(const ArrayXd & theta_a,
        const SimData & sim_data_Para_GoodState_LaborInt,const SimData & sim_data_Para_BadState_LaborInt,
        const SimData & sim_data_Para_GoodState_CapitalInt,const SimData & sim_data_Para_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    tuple<ArrayXXd,double> EstimationIndirectInference_SD_l_Part3_SimData(const ArrayXd & theta_a,
        const SimData & sim_data_Para_GoodState_LaborInt, const SimData & sim_data_Para_BadState_LaborInt,
        const SimData & sim_data_Para_GoodState_CapitalInt, const SimData & sim_data_Para_BadState_CapitalInt,
        MultiThreads::Threads_Management &threadsManagement);

    ArrayXXd Vectorize_dl_dtheta_part(const ArrayXXd & dL_part, const int & dim);
    ArrayXXd Vectorize_dl_dtheta(const ArrayXXd & dL_part1,const ArrayXXd & dL_part2,const ArrayXXd & dL_part3);

    ArrayXXd AppendMatrix(const ArrayXXd & X1, const ArrayXXd & X2, const ArrayXXd & X3, const ArrayXXd & X4);

    ArrayXXd VarianceMatrix_dL(const ArrayXXd & dL_AllPart);





    //
    //

    //
    //

    //
    //
    //
    // tuple<ArrayXXd,double> EstimationIndirectInference_SD_l_Part3_RealData(const ArrayXd & theta_a,
    //     const SimData & ReadData_GoodState_LaborInt, const SimData & ReadData_BadState_LaborInt,
    //     const SimData & ReadData_GoodState_CapitalInt, const SimData & ReadData_BadState_CapitalInt,
    //     MultiThreads::Threads_Management &threadsManagement) {
    //

    //

    // tuple<ArrayXXd,double> Cal_logL_Part3_SimData(const SimData & SimData,
    //     const ParaEst & para_est, const ParaVec & para_vec, const EquStateV & EquV,
    //     const EquStateVmat & Evalmat, const int & GoodState, const int & LaborIntensive,
    //     MultiThreads::Threads_Management &threadsManagement);
    //
    //
    //

    //

    //

    //

    //
    // MatrixXd EstimationIndirectInference_SD_dtheta_a_dtheta0Matrix(const ArrayXd & theta_Est,const ArrayXd & theta_a,
    //     const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    //     const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    //     const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
    //     const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
    //     MultiThreads::Threads_Management &threadsManagement);
    //
    // ArrayXXd Cal_db_dtheta0_v1(const ArrayXd & theta_Est,const ArrayXd & theta_a,
    //     const SimData & simdata_panel_GoodState_LaborInt, const SimData & simdata_panel_BadState_LaborInt,
    //     const SimData & simdata_panel_GoodState_CapitalInt, const SimData & simdata_panel_BadState_CapitalInt,
    //     const SimVar & sim_var_2step_GoodState_LaborInt, const SimVar & sim_var_2step_BadState_LaborInt,
    //     const SimVar & sim_var_2step_GoodState_CapitalInt, const SimVar & sim_var_2step_BadState_CapitalInt,
    //     MultiThreads::Threads_Management &threadsManagement);
    //

    //


}