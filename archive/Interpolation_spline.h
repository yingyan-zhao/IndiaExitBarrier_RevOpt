#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include "Auxillaries.h"
#include <unsupported/Eigen/KroneckerProduct>
//#include "Basics.h"


using namespace Eigen;
//using namespace Ipopt_Wrapper;

namespace alias {
    ArrayXd construct_dv_linear(const ArrayXd & x, const ArrayXd & v);
    ArrayXd construct_dv_hermite_Schumaker(const ArrayXd & x, const ArrayXd & v);

    /******************************************************************************************************************/
    ArrayXXd LinearInterpolation1D_Coeff(const ArrayXd & xvec, const ArrayXd & yvec);

    tuple<ArrayXd, int> SolveZeroQuadratic(const ArrayXd & coef_quadratic);

    /******************************************************************************************************************/
    ArrayXd CalLinearspline1D_Coeff(const double & x1, const double & x2, const double & Y1, const double & Y2);
    ArrayXd CalLinearspline2D_Coeff(const ArrayXd & x, const ArrayXd & y, const ArrayXd & v);
    ArrayXd CalLinearspline3D_Coeff(const ArrayXd & x, const ArrayXd & y, const ArrayXd & z, const ArrayXd & v);

    /******************************************************************************************************************/

    tuple<double, double> SolveMax_Linear1D_K_part(const ArrayXd & coef, const double & Adj, const double & rstar,
        const double & p_K, const double & Kbase);
    tuple<double, double> SolveMax_Linear1D_K(const ArrayXd & coef,const double & x1, const double & x2,
        const double & TotV1, const double TotV2, const double & H, const double & F, const double & c_K,
        const double & rstar, const double & p_K, const double & Kbase);
    tuple<double, double> SolveMax_Linear1D_K_temp(const ArrayXd & coef,const double & x1, const double & x2,
                                              const double & TotV1, const double TotV2, const double & H, const double & F, const double & c_K,
                                              const double & rstar, const double & p_K, const double & Kbase);

    tuple<double, double> SolveMax_Linear1D_L_part(const ArrayXd & coef,const double & Adj, const double & wstar,
        const double & Lbase);
    tuple<double, double> SolveMax_Linear1D_L(const ArrayXd & coef,const double & x1, const double & x2,
        const double & TotV1, const double TotV2, const double & H, const double & F, const double & wstar,
        const double & Lbase);

    /******************************************************************************************************************/
    tuple<ArrayXXd,ArrayXXd> construct_dv_dx_cubicspline2D(const ArrayXd & x, const ArrayXd & y, const ArrayXXd & Y);
    ArrayXXd construct_dv_cross_cubicspline2D(const ArrayXd & x, const ArrayXd & y, const ArrayXXd & dY_dx, const ArrayXXd & dY_dy);
    ArrayXd Calcubicspline2D_Coeff(const ArrayXd & x, const ArrayXd & y, const ArrayXXd & Y,
        const ArrayXXd & dY_dx, const ArrayXXd & dY_dy, const ArrayXXd & d2Y_dxdy);

    /******************************************************************************************************************/
    double solveZero_test(const double & x1, const double & x2, const double & dy1, const double & dy2,
        const double & H, const double & F, const double & w, const double & sigma_Lerror, const double & Lbase);
    double Vmax_solveZero_test(const double & x_opt, const double & x1, const double & x2, const double & y1,
        const double & y2, const double & H, const double & F, const double & w, const double & sigma_Lerror,
        const double & Lbase);
}