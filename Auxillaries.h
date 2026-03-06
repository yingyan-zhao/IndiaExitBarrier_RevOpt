#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <fstream>
#include "Ipopt_wrapper.h"
#include <unsupported/Eigen/SpecialFunctions>

namespace alias {
    using namespace std;
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;
    using ArrayXXi = Eigen::ArrayXXi;
    using ArrayXi = Eigen::ArrayXi;
    using MatrixXd = Eigen::MatrixXd;
    using MatrixXi = Eigen::MatrixXi;

    /***********************************************************************************************
    *   log normal distribution Interval CDF
    ***********************************************************************************************/
    double lognormConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma);
    double lognormConditionalEpectation_right_tail(const double & x1, const double & mu, const double & sigma);
    double lognormCDF(const double & x1, const double & x2, const double & mu, const double & sigma);

    // ////////////////////////////////////////////////////////////////////////////////////
    /// CDF of Normal distribution
    // ////////////////////////////////////////////////////////////////////////////////////
    double normPDF(const double & x, const double & mu, const double & sigma);
    ArrayXd normPDF_vec(const ArrayXd & x, const double & mu, const double & sigma);

    double normCDF(const double & x, double mu, double sigma); // Phi(-∞, x) aka N(x)
    double normCDF_interval(const double & x1, const double & x2, double mu, double sigma);
    double normConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma);
    ArrayXd normCDF_vec(const ArrayXd & x, double mu, double sigma);
    ArrayXXd normCDF_mat(const ArrayXd & x, const ArrayXd & mu, const double & sigma);
    // ////////////////////////////////////////////////////////////////////////////////////
    /// Generate Normal Grid
    // ////////////////////////////////////////////////////////////////////////////////////
    tuple<ArrayXd,ArrayXd> GenerateNormalGrid(const size_t N);

    ArrayXi find_matlab(const ArrayXi A);

    int find_lower_bound_scale(const ArrayXd &vect, double xx);
    ArrayXi find_lower_bound_vec(const ArrayXd &vect, ArrayXd xx);
    ArrayXXi find_lower_bound_mat(const ArrayXd &vect, ArrayXXd xx);

    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

    ArrayXi mat_max_index(ArrayXXd mat);

    tuple<double, double, double> SolveMax_cubicspline2D(const ArrayXd Coef_spline, const ArrayXd x,
        const ArrayXd y, const ArrayXXd Y);

    void writeToCSVfile(string name, MatrixXd matrix);
    ArrayXXd readCSV(std::string file, int rows, int cols);

    ArrayXXd Cwisepow(const double a, const ArrayXXd & mat);

    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    void load_csv (Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & M,
                   const std::string & path) {
        std::ifstream indata;
        indata.open(path);
        if (!indata.good()){
            throw std::runtime_error("error in opening files");
        }
        std::string line;
        std::vector<_Scalar> values;
        uint rows = 0;
        while (std::getline(indata, line)) {
            std::stringstream lineStream(line);
            std::string cell;
            while (std::getline(lineStream, cell, ',')) {
                values.push_back(std::stod(cell));
            }
            ++rows;
        }
        Eigen::Map<const Eigen::Matrix<_Scalar, Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>> data_map(values.data(), rows, values.size()/rows);
        M = data_map;
        indata.close();
    }
    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    void load_csv (Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & M,
                   const std::string & path) {
        std::ifstream indata;
        indata.open(path);
        if (!indata.good()){
            throw std::runtime_error("error in opening files");
        }
        std::string line;
        std::vector<_Scalar> values;
        uint rows = 0;
        while (std::getline(indata, line)) {
            std::stringstream lineStream(line);
            std::string cell;
            while (std::getline(lineStream, cell, ',')) {
                values.push_back(std::stod(cell));
            }
            ++rows;
        }
        Eigen::Map<const Eigen::Array<_Scalar, Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>> data_map(values.data(), rows, values.size()/rows);
        M = data_map;
        indata.close();
    }
//    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    void save_csv(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
//                  & M, std::string file_name){
//        const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
//        std::ofstream file(file_name.c_str());
//        file << M.format(CSVFormat);
//
//        file.close();
//    }
//    template <class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    void save_csv(const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
//                  & M, std::string file_name){
//        const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
//        std::ofstream file(file_name.c_str());
//        file << M.format(CSVFormat);
//        file.close();
//    }
}