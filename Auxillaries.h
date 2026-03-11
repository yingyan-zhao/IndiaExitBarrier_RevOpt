#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/distributions/normal.hpp>
//#include <iostream>
//#include <fstream>
//#include "Ipopt_wrapper.h"
//#include <unsupported/Eigen/SpecialFunctions>
//
namespace alias
{
    using namespace std;
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;
    using ArrayXXi = Eigen::ArrayXXi;
    using ArrayXi = Eigen::ArrayXi;
    using MatrixXd = Eigen::MatrixXd;
    using MatrixXi = Eigen::MatrixXi;

    //// ////////////////////////////////////////////////////////////////////////////////////
    /// Generate Normal Grid
    //// //////////////////////////////////////////////////////////////////////////////////
    ArrayXd GenerateNormalGrid(const size_t N);
    ArrayXd NormalGridProbabilities(const ArrayXd& x, double mu, double sigma);
    MatrixXd GenerateTransitionMatrixAR1(const ArrayXd& x_grid, double gamma0, double gamma1, double sigma_eps);


    // ////////////////////////////////////////////////////////////////////////////////////
    /// CDF of Normal distribution
    // ////////////////////////////////////////////////////////////////////////////////////
    ArrayXd normCDF_vec(const ArrayXd& X);

    /***********************************************************************************************
    * log normal distribution Interval CDF
    ***********************************************************************************************/
    ArrayXd lognormalEY_lessX(const ArrayXd& X, double mu, double sigma);

    /***********************************************************************************************
    * Linear Spine : find the coefficient
    * xvec = [x1,x2]; yvec = [y1,y2]
    ***********************************************************************************************/
    ArrayXd LinearInterpolation1D_Coeff(const double & x1, const double & x2, const double & y1, const double & y2);

    /***********************************************************************************************
    // Solve a*x^alpha + b*x + c = 0, with x > 0.
    // Uses Newton in log-space: x = exp(y), so positivity is guaranteed.
    ***********************************************************************************************/
    struct NewtonResult {
        bool converged;
        double x;   // positive root estimate
        int iterations;
    };
    NewtonResult SolvePowerEqPositiveNewton(const double & a, const double & alpha, const double & b, const double & c,
        const double & x_ini);


    double solve_axalpha_plus_c_eq0(const double & a, const double & alpha, const double & c);

    struct RootResult {
        bool converged;
        bool bracket_found;
        double k;         // root estimate (k>0)
        int iterations;
    };
    RootResult solve_k_equation_hybrid(const double & d1, const double & d0, const double & d2, const double & a);

}




//
//    /***********************************************************************************************
//    *   log normal distribution Interval CDF
//    ***********************************************************************************************/
//    double lognormConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma);
//    double lognormConditionalEpectation_right_tail(const double & x1, const double & mu, const double & sigma);
//    double lognormCDF(const double & x1, const double & x2, const double & mu, const double & sigma);
//
//    // ////////////////////////////////////////////////////////////////////////////////////
//    /// CDF of Normal distribution
//    // ////////////////////////////////////////////////////////////////////////////////////
//    double normPDF(const double & x, const double & mu, const double & sigma);
//    ArrayXd normPDF_vec(const ArrayXd & x, const double & mu, const double & sigma);
//
//    double normCDF(const double & x, double mu, double sigma); // Phi(-∞, x) aka N(x)
//    double normCDF_interval(const double & x1, const double & x2, double mu, double sigma);
//    double normConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma);
//    ArrayXd normCDF_vec(const ArrayXd & x, double mu, double sigma);
//    ArrayXXd normCDF_mat(const ArrayXd & x, const ArrayXd & mu, const double & sigma);
//
//
//
//    int find_lower_bound_scale(const ArrayXd &vect, double xx);
//    ArrayXi find_lower_bound_vec(const ArrayXd &vect, ArrayXd xx);
//    ArrayXXi find_lower_bound_mat(const ArrayXd &vect, ArrayXXd xx);
//
//    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//
//    ArrayXi mat_max_index(ArrayXXd mat);
//
//    tuple<double, double, double> SolveMax_cubicspline2D(const ArrayXd Coef_spline, const ArrayXd x,
//        const ArrayXd y, const ArrayXXd Y);
//
//    void writeToCSVfile(string name, MatrixXd matrix);
//    ArrayXXd readCSV(std::string file, int rows, int cols);
//
//    ArrayXXd Cwisepow(const double a, const ArrayXXd & mat);
//
//    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    void load_csv (Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & M,
//                   const std::string & path) {
//        std::ifstream indata;
//        indata.open(path);
//        if (!indata.good()){
//            throw std::runtime_error("error in opening files");
//        }
//        std::string line;
//        std::vector<_Scalar> values;
//        uint rows = 0;
//        while (std::getline(indata, line)) {
//            std::stringstream lineStream(line);
//            std::string cell;
//            while (std::getline(lineStream, cell, ',')) {
//                values.push_back(std::stod(cell));
//            }
//            ++rows;
//        }
//        Eigen::Map<const Eigen::Matrix<_Scalar, Eigen::Dynamic,Eigen::Dynamic,
//                Eigen::RowMajor>> data_map(values.data(), rows, values.size()/rows);
//        M = data_map;
//        indata.close();
//    }
//    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    void load_csv (Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & M,
//                   const std::string & path) {
//        std::ifstream indata;
//        indata.open(path);
//        if (!indata.good()){
//            throw std::runtime_error("error in opening files");
//        }
//        std::string line;
//        std::vector<_Scalar> values;
//        uint rows = 0;
//        while (std::getline(indata, line)) {
//            std::stringstream lineStream(line);
//            std::string cell;
//            while (std::getline(lineStream, cell, ',')) {
//                values.push_back(std::stod(cell));
//            }
//            ++rows;
//        }
//        Eigen::Map<const Eigen::Array<_Scalar, Eigen::Dynamic,Eigen::Dynamic,
//                Eigen::RowMajor>> data_map(values.data(), rows, values.size()/rows);
//        M = data_map;
//        indata.close();
//    }
////    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
////    void save_csv(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
////                  & M, std::string file_name){
////        const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
////        std::ofstream file(file_name.c_str());
////        file << M.format(CSVFormat);
////
////        file.close();
////    }
////    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
////    void save_csv(const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
////                  & M, std::string file_name){
////        const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
////        std::ofstream file(file_name.c_str());
////        file << M.format(CSVFormat);
////        file.close();
////    }
//}