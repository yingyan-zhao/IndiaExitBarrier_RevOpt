
#include "Interpolation_spline.h"

using namespace Eigen;
using namespace std;
using namespace alias;
//using Ipopt_Wrapper::Ipopt_Solver;
//using Ipopt_Wrapper::Option;
//using Ipopt_Wrapper::Ipopt_Status;
//using Ipopt_Wrapper::Simple_Opt_Problem;

//// Dynamic Programming with Hermite Interpolation : KEN JUDD AND CAI (2011); Alogrithm 3. Revised Schumaker Shape-Preserving Interpolation
//// (x,v) is a decreasing function

/************************************************************************************************************
**** calculate first diff
************************************************************************************************************/
ArrayXd alias::construct_dv_linear(const ArrayXd & x, const ArrayXd & v) {
    int n = x.size();

    ArrayXd delta = (v.segment(1,n-1) - v.segment(0,n-1))
                    / (x.segment(1,n-1) - x.segment(0,n-1));

    return delta;
}

ArrayXd alias::construct_dv_hermite_Schumaker(const ArrayXd & x, const ArrayXd & v) {
    int n = x.size();
    ArrayXd dv(n);

    ArrayXd delta = (v.segment(1,n-1) - v.segment(0,n-1))
                    / (x.segment(1,n-1) - x.segment(0,n-1));

    ArrayXd L1 = ( (x.segment(1,n-2) - x.segment(0,n-2)) )
                 / (x.segment(2,n-2) - x.segment(0,n-2));
    ArrayXd L2 = ( (x.segment(2,n-2) - x.segment(1,n-2)) )
                 / (x.segment(2,n-2) - x.segment(0,n-2));

    ArrayXd dv_temp = L1*delta.segment(0,n-2) + L2*delta.segment(1,n-2);

    double dv_Start = delta(0) - (dv_temp(0) - delta(0));
    double dv_End = delta(n-2) + (delta(n-2) - dv_temp(n-3));

    dv << dv_Start,dv_temp,dv_End;
    return dv;
}

ArrayXXd alias::LinearInterpolation1D_Coeff(const ArrayXd & xvec, const ArrayXd & yvec) {
    int N = xvec.size();
    ArrayXXd coef(N+1,2);

    ArrayXd temp_coef1 = (yvec.segment(1,N-1) - yvec.segment(0,N-1))
                         / (xvec.segment(1,N-1) - xvec.segment(0,N-1));
    double temp_coef1_Start = temp_coef1(0);
    double temp_coef1_End = temp_coef1(N-2);

    coef.col(1) << temp_coef1_Start,temp_coef1,temp_coef1_End;

    ArrayXd temp_coef0 = yvec.segment(0,N-1) - temp_coef1*xvec.segment(0,N-1);
    double temp_coef0_Start = yvec(0) - temp_coef1_Start*xvec(0);
    double temp_coef0_End = yvec(N-1) - temp_coef1_End*xvec(N-1);

    coef.col(0) << temp_coef0_Start,temp_coef0,temp_coef0_End;

    return coef;
}

tuple<ArrayXd, int> alias::SolveZeroQuadratic(const ArrayXd & coef_quadratic) {
    double a = coef_quadratic(2);
    double b = coef_quadratic(1);
    double c = coef_quadratic(0);

    double Delta = pow(b,2) - 4.0*a*c;
//    cout << "Delta = " << Delta << endl;

    ArrayXd x_opt(2); int Indicator;
    if (Delta < 0) {
        Indicator = 0;
    }
    else {
        Indicator = 1;
        x_opt(0) = (-b + sqrt(Delta)) / 2.0 / a;
        x_opt(1) = (-b - sqrt(Delta)) / 2.0 / a;
    }
    return tuple<ArrayXd, int>(x_opt,Indicator);
}


/************************************************************************************************************
**** Linear Interpolation
************************************************************************************************************/
ArrayXd alias::CalLinearspline1D_Coeff(const double & x1, const double & x2, const double & Y1, const double & Y2) {

    double b = (Y2-Y1)/(x2-x1);
    double a = Y1 - b * x1;
    ArrayXd coef(2); coef << a,b;

    return coef;
}

ArrayXd alias::CalLinearspline2D_Coeff(const ArrayXd & x, const ArrayXd & y, const ArrayXd & v) {

    double temp1 = 1.0 / (x(1) - x(0)) / (y(1) - y(0));
    MatrixXd A(4,4);
    A(0,0) = x(1) * y(1); A(0,1) = - x(1) * y(0);
    A(0,2) = - x(0) * y(1); A(0,3) = x(0) * y(0);

    A(1,0) = - y(1); A(1,1) = y(0);
    A(1,2) = y(1); A(1,3) = - y(0);

    A(2,0) = - x(1); A(2,1) = x(1);
    A(2,2) = x(0); A(2,3) = - x(0);

    A(3,0) = 1.0; A(3,1) = - 1.0; A(3,2) = - 1.0; A(3,3) = 1.0;

    MatrixXd B = v.matrix();

    ArrayXd coef = temp1 * A * B;
    return coef;
}

ArrayXd alias::CalLinearspline3D_Coeff(const ArrayXd & x, const ArrayXd & y, const ArrayXd & z, const ArrayXd & v) {

    double denominator = (x(0) - x(1)) * (y(0) - y(1)) * (z(0) - z(1));

    double a0 = -v(0)*x(1)*y(1)*z(1) + v(4)*x(1)*y(1)*z(0)
                + v(2)*x(1)*y(0)*z(1) - v(6)*x(1)*y(0)*z(0)
                + v(1)*x(0)*y(1)*z(1) - v(5)*x(0)*y(1)*z(0)
                - v(3)*x(0)*y(0)*z(1) + v(7)*x(0)*y(0)*z(0);

    double a1 = v(0)*y(1)*z(1) - v(4)*y(1)*z(0)
                - v(2)*y(0)*z(1) + v(6)*y(0)*z(0)
                - v(1)*y(1)*z(1) + v(5)*y(1)*z(0)
                + v(3)*y(0)*z(1) - v(7)*y(0)*z(0);

    double a2 = v(0)*x(1)*z(1) - v(4)*x(1)*z(0)
                - v(2)*x(1)*z(1) + v(6)*x(1)*z(0)
                - v(1)*x(0)*z(1) + v(5)*x(0)*z(0)
                + v(3)*x(0)*z(1) - v(7)*x(0)*z(0);

    double a3 = v(0)*x(1)*y(1) - v(4)*x(1)*y(1)
                - v(2)*x(1)*y(0) + v(6)*x(1)*y(0)
                - v(1)*x(0)*y(1) + v(5)*x(0)*y(1)
                + v(3)*x(0)*y(0) - v(7)*x(0)*y(0);

    double a4 = -v(0)*z(1) + v(4)*z(0) + v(2)*z(1) - v(6)*z(0)
                + v(1)*z(1) - v(5)*z(0) - v(3)*z(1) + v(7)*z(0);

    double a5 = -v(0)*y(1) + v(4)*y(1) + v(2)*y(0) - v(6)*y(0)
                + v(1)*y(1) - v(5)*y(1)
                - v(3)*y(0) + v(7)*y(0);

    double a6 = -v(0)*x(1) + v(4)*x(1) + v(2)*x(1) - v(6)*x(1)
                + v(1)*x(0) - v(5)*x(0) - v(3)*x(0) + v(7)*x(0);

    double a7 = v(0) - v(4) - v(2) + v(6)
                - v(1) + v(5) + v(3) - v(7);

    ArrayXd coef_temp(8); coef_temp << a0,a1,a2,a3,a4,a5,a6,a7;
    ArrayXd coef = coef_temp / denominator;

    return coef;
}

/************************************************************************************
**** solve optimal K
************************************************************************************/
tuple<double, double> alias::SolveMax_Linear1D_K_part(const ArrayXd & coef, const double & Adj, const double & rstar,
                                                      const double & p_K, const double & Kbase) {

    double x_opt = (coef(1) - rstar - p_K)/2.0/Adj + Kbase;
    double v_max = coef(0) + coef(1)*x_opt - Adj*pow(x_opt - Kbase,2) - rstar*x_opt - p_K*(x_opt-Kbase);

    return tuple<double,double>(x_opt,v_max);
}

tuple<double,double> alias::SolveMax_Linear1D_K(const ArrayXd & coef,const double & x1, const double & x2,
    const double & TotV1, const double TotV2, const double & H, const double & F, const double & c_K,
    const double & rstar, const double & p_K, const double & Kbase) {

    double x_opt; double v_max;
    if (Kbase <= x1) {
        double Adj = H;
        double p_K_adj = p_K;

        tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef,Adj,rstar,p_K_adj,Kbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);
        if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
        else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else if (Kbase >= x2) {
        double Adj = F;
        double p_K_adj = p_K * c_K;

        tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef,Adj,rstar,p_K_adj,Kbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);
        if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
        else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else {
        double VKbase = coef(0) + coef(1) * Kbase - rstar * Kbase;
        double dVKbase = coef(1) - rstar - p_K;
        if (dVKbase >= 0) {
            double Adj = H;
            double p_K_adj = p_K;

            tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef, Adj, rstar, p_K_adj, Kbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < Kbase) { x_opt = Kbase; v_max = VKbase; }
            else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        }
        else {
            double Adj = F;
            double p_K_adj = p_K * c_K;

            tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef, Adj, rstar, p_K_adj, Kbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
            else if (x_opt_temp > Kbase) { x_opt = Kbase; v_max = VKbase; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        }
    }

    return tuple<double, double>(x_opt,v_max);
}


tuple<double,double> alias::SolveMax_Linear1D_K_temp(const ArrayXd & coef,const double & x1, const double & x2,
    const double & TotV1, const double TotV2, const double & H, const double & F, const double & c_K,
    const double & rstar, const double & p_K, const double & Kbase) {
//    cout << "???? rstar = " << rstar << "; p_K = " << p_K << "; Kbase = " << Kbase << "; x1 = " << x1 << "; x2 = " << x2
//        << "; coef = " << coef.transpose() << endl;
    double x_opt; double v_max;
    if (Kbase <= x1) {
        double Adj = H;
        double p_K_adj = p_K;

        tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef,Adj,rstar,p_K_adj,Kbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);
        if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
        else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else if (Kbase >= x2) {
        double Adj = F;
        double p_K_adj = p_K * c_K;
//        cout << "Adj = " << Adj << endl;

        tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef,Adj,rstar,p_K_adj,Kbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);

//        cout << "x_opt_temp = " << x_opt_temp << endl;
        if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
        else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else {
        double VKbase = coef(0) + coef(1) * Kbase - rstar * Kbase;
        double dVKbase = coef(1) - rstar - p_K;
        if (dVKbase >= 0) {
            double Adj = H;
            double p_K_adj = p_K;

            tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef, Adj, rstar, p_K_adj, Kbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < Kbase) { x_opt = Kbase; v_max = VKbase; }
            else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        }
        else {
            double Adj = F;
            double p_K_adj = p_K * c_K;

            tuple<double, double> t_max = SolveMax_Linear1D_K_part(coef, Adj, rstar, p_K_adj, Kbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
            else if (x_opt_temp > Kbase) { x_opt = Kbase; v_max = VKbase; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        }
    }

    return tuple<double, double>(x_opt,v_max);
}

/************************************************************************************
**** solve optimal L
************************************************************************************/
tuple<double, double> alias::SolveMax_Linear1D_L_part(const ArrayXd & coef,const double & Adj, const double & wstar,
    const double & Lbase) {

    double x_opt = (coef(1) - wstar)/2.0/Adj + Lbase ;
    double v_max = coef(0) + coef(1)*x_opt - Adj*pow(x_opt - Lbase,2) - wstar*x_opt;

    return tuple<double,double>(x_opt,v_max);
}

tuple<double, double> alias::SolveMax_Linear1D_L(const ArrayXd & coef,const double & x1, const double & x2,
    const double & TotV1, const double TotV2, const double & H, const double & F, const double & wstar,
    const double & Lbase) {

    double x_opt; double v_max;
    if (Lbase <= x1) {
        double Adj = H;

        tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef,Adj,wstar,Lbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);
        if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
        else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else if (Lbase >= x2) {
        double Adj = F;

        tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef,Adj,wstar,Lbase);
        double x_opt_temp = get<0>(t_max);
        double v_max_temp = get<1>(t_max);
        if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
        else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
        else { x_opt = x_opt_temp; v_max = v_max_temp; }
    }
    else {
        double VLbase = coef(0) + coef(1) * Lbase - wstar * Lbase;
        double dVLbase = coef(1) - wstar;
        if (dVLbase >= 0) {
            double Adj = H;

            tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef, Adj, wstar, Lbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < Lbase) { x_opt = Lbase; v_max = VLbase; }
            else if (x_opt_temp > x2) { x_opt = x2; v_max = TotV2; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        } else {
            double Adj = F;

            tuple<double, double> t_max = SolveMax_Linear1D_L_part(coef, Adj, wstar, Lbase);
            double x_opt_temp = get<0>(t_max);
            double v_max_temp = get<1>(t_max);
            if (x_opt_temp < x1) { x_opt = x1; v_max = TotV1; }
            else if (x_opt_temp > Lbase) { x_opt = Lbase; v_max = VLbase; }
            else { x_opt = x_opt_temp; v_max = v_max_temp; }
        }
    }

    return tuple<double, double>(x_opt,v_max);
}
/************************************************************************************************************/
tuple<ArrayXXd,ArrayXXd> alias::construct_dv_dx_cubicspline2D(const ArrayXd & x, const ArrayXd & y, const ArrayXXd & Y) {

    int n_x = x.size();
    int n_y = y.size();

    ArrayXXd dY_dy(n_x,n_y);
    for (size_t i_x = 0; i_x < n_x; ++i_x) {
        dY_dy.row(i_x) = construct_dv_hermite_Schumaker(y, Y.row(i_x).transpose()).transpose();
    }
    ArrayXXd dY_dx(n_x,n_y);
    for (size_t i_y = 0; i_y < n_y; ++i_y) {
        dY_dx.col(i_y) = construct_dv_hermite_Schumaker(x, Y.col(i_y));
    }

    return tuple<ArrayXXd,ArrayXXd>(dY_dx,dY_dy);
}

ArrayXXd alias::construct_dv_cross_cubicspline2D(const ArrayXd & x, const ArrayXd & y, const ArrayXXd & dY_dx, const ArrayXXd & dY_dy) {

    int n_x = x.size();
    int n_y = y.size();

    ArrayXXd d2Y_dxdy1(n_x,n_y);
    for (size_t i_x = 0; i_x < n_x; ++i_x) {
        d2Y_dxdy1.row(i_x) = construct_dv_hermite_Schumaker(y, dY_dx.row(i_x).transpose()).transpose();
    }
    ArrayXXd d2Y_dxdy2(n_x,n_y);
    for (size_t i_y = 0; i_y < n_y; ++i_y) {
        d2Y_dxdy2.col(i_y) = construct_dv_hermite_Schumaker(x, dY_dy.col(i_y));
    }

    ArrayXXd d2Y_dxdy = 0.5*d2Y_dxdy1 + 0.5*d2Y_dxdy2;
    return d2Y_dxdy;
}

ArrayXd alias::Calcubicspline2D_Coeff(const ArrayXd & x, const ArrayXd & y, const ArrayXXd & Y,
                                      const ArrayXXd & dY_dx, const ArrayXXd & dY_dy, const ArrayXXd & d2Y_dxdy) {

    int n_row = 2;
    int n_col = 2;

    int x_digit = 0;
    if (x(0) > 0) {
        x_digit = log10(x(0));
    }
    int y_digit = 0;
    if (y(0) > 0) {
        y_digit = log10(y(0));
    }

    ArrayXd xx = x/pow(10,x_digit);
    ArrayXd yy = y/pow(10,y_digit);

    ArrayXXd dY_dxx = dY_dx*pow(10,x_digit);
    ArrayXXd dY_dyy = dY_dy*pow(10,y_digit);
    ArrayXXd d2Y_dxxdyy = d2Y_dxdy*pow(10,x_digit)*pow(10,y_digit);

    ArrayXXd temp_x(n_row,4);
    temp_x.col(0) = ArrayXd::Ones(n_row);
    temp_x.col(1) = xx;
    temp_x.col(2) = xx.pow(2);
    temp_x.col(3) = xx.pow(3);

    ArrayXXd temp1_x(n_row,4);
    temp1_x.col(0) = ArrayXd::Zero(n_row);
    temp1_x.col(1) = ArrayXd::Ones(n_row);
    temp1_x.col(2) = 2*xx.pow(1);
    temp1_x.col(3) = 3*xx.pow(2);

//    ArrayXXd temp2_x(n_row,4);
//    temp2_x.col(0) = ArrayXd::Zero(n_row);
//    temp2_x.col(1) = ArrayXd::Zero(n_row);
//    temp2_x.col(2) = 2*ArrayXd::Ones(n_row);
//    temp2_x.col(3) = 6*x;

    ArrayXXd temp_y(n_col,4);
    temp_y.col(0) = ArrayXd::Ones(n_col);
    temp_y.col(1) = yy;
    temp_y.col(2) = yy.pow(2);
    temp_y.col(3) = yy.pow(3);

    ArrayXXd temp1_y(n_col,4);
    temp1_y.col(0) = ArrayXd::Zero(n_col);
    temp1_y.col(1) = ArrayXd::Ones(n_col);
    temp1_y.col(2) = 2*yy.pow(1);
    temp1_y.col(3) = 3*yy.pow(2);

//    ArrayXXd temp2_y(n_row,4);
//    temp2_y.col(0) = ArrayXd::Zero(n_row);
//    temp2_y.col(1) = ArrayXd::Zero(n_row);
//    temp2_y.col(2) = 2*ArrayXd::Ones(n_col);
//    temp2_y.col(3) = 6*y;

    // set 1
    ArrayXXd A1 = KroneckerProduct(temp_x,temp_y);
    // set 2: dy_dx
    ArrayXXd A2 = KroneckerProduct(temp1_x,temp_y);
    // set 3: dy_dy
    ArrayXXd A3 = KroneckerProduct(temp_x,temp1_y);
    // set 3: d2y_dxdy
    ArrayXXd A4 = KroneckerProduct(temp1_x,temp1_y);
    ArrayXXd A = ArrayXXd::Zero(16,16);
    A(seqN(0,4),all) = A1;
    A(seqN(4,4),all) = A2;
    A(seqN(8,4),all) = A3;
    A(seqN(12,4),all) = A4;
//    A.row(12) = KroneckerProduct(temp_x.row(0),temp2_y.row(1));
//    A.row(13) = KroneckerProduct(temp2_x.row(1),temp_y.row(0));
//    A.row(14) = KroneckerProduct(temp_x.row(1),temp2_y.row(1));
//    A.row(15) = KroneckerProduct(temp2_x.row(1),temp_y.row(1));

    ArrayXd B(16);
    B.segment(0,2) = Y.row(0).transpose();
    B.segment(2,2) = Y.row(1).transpose();
    B.segment(4,2) = dY_dxx.row(0).transpose();
    B.segment(6,2) = dY_dxx.row(1).transpose();
    B.segment(8,2) = dY_dyy.row(0).transpose();
    B.segment(10,2) = dY_dyy.row(1).transpose();
    B.segment(12,2) = d2Y_dxxdyy.row(0).transpose();
    B.segment(14,2) = d2Y_dxxdyy.row(1).transpose();
//    B(12) = d2Y_d2Lur(0,1);
//    B(13) = d2Y_d2Luc(1,0);
//    B(14) = d2Y_d2Lur(1,1);
//    B(15) = d2Y_d2Luc(1,1);
//    cout << "B = " << B.transpose() << endl;

    MatrixXd coef_Y = A.matrix().colPivHouseholderQr().solve(B.matrix());

    ArrayXd coefx(4); coefx << 1,1.0/pow(10,x_digit),1.0/pow(10,x_digit*2),
            1.0/pow(10,x_digit*3);
    ArrayXd coefy(4); coefy << 1,1.0/pow(10,y_digit),1.0/pow(10,y_digit*2),
            1.0/pow(10,y_digit*3);
    ArrayXd coefxy = KroneckerProduct(coefx,coefy);

    ArrayXd coef_Y_array = coef_Y.array() * coefxy;

    return coef_Y_array;
}


//tuple<double, double> alias::SolveMax_rationalspline1D_K_part(const ArrayXd & coef,const double & coef_x1,
//    const double & coef_x2, const double & x1, const double & x2, const double & y1, const double y2,
//    const double & Adj, const double & rstar, const double & p_K, const double & Kbase,
//    Ipopt_Solver & solver) {
//
//    double Vdefault; double xdefault;
//    if (y1 > y2) {Vdefault = y1; xdefault = x1;}
//    else {Vdefault = y2; xdefault = x2;}
//
//    auto f = [&coef,&coef_x1,&coef_x2,&Adj,&rstar,&p_K,&Kbase]
//        (std::vector<double> x, std::vector<double> & grad){
//
//        double obj_temp = SolveVal_rationalspline1D(coef,coef_x1,coef_x2,x[0])
//            - Adj * pow(x[0] - Kbase,2) - rstar*x[0] - p_K*(x[0] - Kbase);
//
//        double grad_temp = SolvediffVal_rationalspline1D(coef,coef_x1,coef_x2,x[0])
//            - 2.0*Adj*(x[0] - Kbase) - rstar - p_K;
//
//        double obj = -obj_temp;
//        grad[0] = -grad_temp;
//        return obj;
//    };
//
//    Simple_Opt_Problem problem(f, 1);
//
//    std::vector<double> x_lb(1); x_lb[0] = x1;
//    std::vector<double> x_ub(1); x_ub[0] = x2;
//    problem.set_x_lb(x_lb);
//    problem.set_x_ub(x_ub);
//
//    std::vector<double> x_ini(1); x_ini[0] = xdefault;
//
//    auto result = solver.optimize(problem, x_ini);
//
//    double x_opt; double v_max;
//    if (result.status == Ipopt_Wrapper::success){
////        std::cout << "opt_f = " << result.opt_obj << ", opt_x = " << result.opt_x[0] << endl;
//        x_opt = result.opt_x[0];
//        v_max = result.opt_obj;
//    }
//
//    return tuple<double, double>(x_opt,v_max);
//}
//
//tuple<double, double> alias::SolveMax_rationalspline1D_K_part(const ArrayXd & coef,const double & coef_x1,
//    const double & coef_x2, const double & x1, const double & x2, const double & y1, const double y2,
//    const double & AdjC, const double & rstar, const double & p_K, const double & Kbase,Ipopt_Solver & solver) {
//
//    double xdefault; double Vdefault;
//    if (y1 <= y2) {Vdefault = y2; xdefault = x2;}
//    else {Vdefault = y1; xdefault = x1;}
//
//    cout << "coef = " << coef.transpose() << endl;
//    cout << "coef_x1 = " << coef_x1 << "; coef_x2 = " << coef_x2 << "; x1 = " << x1 << "; x2 = " << x2
//         << "; y1 = " << y1 << "; y2 = " << y2 << endl;
//    double c1 = coef(0); double c2 = coef(1); double c3 = coef(2); double c4 = coef(3);
//    cout << "rstar = " << rstar << "; p_K = " << p_K << "; AdjC = " << AdjC << endl;
////    ArrayXd x_temp = ArrayXd::LinSpaced(100,x2-0.05,x2);
////    ArrayXd y_temp(100);
////    for (size_t i = 0; i < 100; ++i) {
////        cout << "i = " << i << endl;
////        y_temp(i) = SolvediffVal_rationalspline1D(coef, coef_x1, coef_x2, x_temp(i))
////                - rstar - p_K - 2.0*AdjC*(x_temp(i) - Kbase);
////    }
////    cout << "y_temp = " << y_temp << endl;
////    throw runtime_error("310");
////// Part1
//    ArrayXd coef_part1(4);
//    coef_part1(0) = (c2 - rstar - p_K) * pow(c4, 2) * pow(coef_x1 - coef_x2, 2);
//    coef_part1(1) = 2.0 * (c2 - rstar - p_K) * (c3 + c4) * c4 * (coef_x1 - coef_x2);
//    coef_part1(2) = (c2 - rstar - p_K) * pow((c3 + c4), 2);
//    coef_part1(3) = 0.0;
//
//    cout << "coef_part1 = " << coef_part1.transpose() << endl;
////// Part2
//    ArrayXd coef_part2(4);
//    coef_part2(0) = c3 * pow(c4, 2) * pow(coef_x1 - coef_x2, 2);
//    coef_part2(1) = 2.0 * c3 * pow(c4, 2) * (coef_x1 - coef_x2);
//    coef_part2(2) = c3 * pow(c4, 2) + pow(c3, 2) * c4;
//    coef_part2(3) = 0.0;
//
//    cout << "coef_part2 = " << coef_part2.transpose() << endl;
////// Part3
//
//    ArrayXd coef_part3(4);
//    cout << "AdjC = " << AdjC << "; coef_x1 = " << coef_x1 << "; Kbase = " << Kbase << "; c4 = " << c4 << "; coef_x2= " << coef_x2 << endl;
//    coef_part3(0) = (-2*AdjC*(coef_x1-Kbase)) * pow(c4, 2) * pow(coef_x1 - coef_x2, 2);
//    coef_part3(1) = 2.0 * (-2*AdjC*(coef_x1-Kbase)) * (c3 + c4) * c4 * (coef_x1 - coef_x2);
//    coef_part3(2) = (-2*AdjC*(coef_x1-Kbase)) * pow((c3 + c4), 2);
//    coef_part3(3) = 0.0;
//    cout << " coef_part3(0) = " <<  coef_part3(0) << endl;
//
//    ArrayXd coef_part4(4);
//    coef_part4(0) = 0.0;
//    coef_part4(1) = -2*AdjC * pow(c4, 2) * pow(coef_x1 - coef_x2, 2);
//    coef_part4(2) = 2.0 * (-2*AdjC) * (c3 + c4) * c4 * (coef_x1 - coef_x2);
//    coef_part4(3) = (-2*AdjC) * pow((c3 + c4), 2);
//
//    cout << "coef_part3+coef_part4 = " << (coef_part3+coef_part4).transpose() << endl;
//
//    ArrayXd coef_cubic = coef_part1 + coef_part2 + coef_part3 + coef_part4;
//
//    cout << "coef_cubic = " << coef_cubic.transpose() << endl;
//
//    tuple<ArrayXd, int> t_max = SolveZeroCubic(coef_cubic);
//    ArrayXd x_opt_vec = get<0>(t_max);
//    cout << "x_opt_vec = " << x_opt_vec.transpose() << endl;
//    x_opt_vec = x_opt_vec + coef_x1;
//    cout << "x_opt_vec = " << x_opt_vec.transpose() << endl;
//    int Indicator = get<1>(t_max);
//    cout << "Indicator = " << Indicator << endl;
//    cout << "xdefault = " << xdefault << "; Vdefault = " << Vdefault << endl;
//    double x_opt = xdefault; double v_max = Vdefault;
//    if (Indicator == 0) {
//        if (x_opt_vec(0) < x2 & x_opt_vec(0) > x1) {
//            double y = SolveVal_rationalspline1D(coef_cubic,x1,x2,x_opt_vec(0))
//                       - AdjC * pow(x_opt_vec(0) - Kbase,2) - rstar*x_opt_vec(0)
//                       - p_K * (x_opt_vec(0) - Kbase);
//            if (y > Vdefault) {
//                x_opt = x_opt_vec(0);
//                v_max = y;
//            }
//        }
//    } else {
//        double yopt_test = Vdefault;
//        double xopt_test = xdefault;
//        if (x_opt_vec(0) < x2 & x_opt_vec(0) > x1) {
//            double y = SolveVal_rationalspline1D(coef_cubic,x1,x2,x_opt_vec(0))
//                       - AdjC * pow(x_opt_vec(0) - Kbase,2) - rstar*x_opt_vec(0)
//                       - p_K * (x_opt_vec(0) - Kbase);
//            if (y > yopt_test) {
//                yopt_test = y;
//                xopt_test = x_opt_vec(0);
//            }
//        }
//        else if (x_opt_vec(1) < x2 & x_opt_vec(1) > x1) {
//            double y = SolveVal_rationalspline1D(coef_cubic,x1,x2,x_opt_vec(1))
//                       - AdjC * pow(x_opt_vec(1) - Kbase,2) - rstar*x_opt_vec(1)
//                       - p_K * (x_opt_vec(1) - Kbase);
//            if (y > yopt_test) {
//                yopt_test = y;
//                xopt_test = x_opt_vec(1);
//            }
//        } else if (x_opt_vec(2) < x2 & x_opt_vec(2) > x1) {
//            double y = SolveVal_rationalspline1D(coef_cubic,x1,x2,x_opt_vec(2))
//                       - AdjC * pow(x_opt_vec(2) - Kbase,2) - rstar*x_opt_vec(2)
//                       - p_K * (x_opt_vec(2) - Kbase);
//            if (y > yopt_test) {
//                yopt_test = y;
//                xopt_test = x_opt_vec(2);
//            }
//        }
//
//        if (yopt_test > Vdefault) {
//            x_opt = xopt_test;
//            v_max = yopt_test;
//        }
//    }
//
//    cout << "x_opt = " << x_opt << "; v_max = " << v_max << endl;
//    throw runtime_error("365");
//    return tuple<double, double>(x_opt,v_max);
//}
