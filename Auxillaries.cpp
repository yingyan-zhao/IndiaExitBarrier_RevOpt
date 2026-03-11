#include "Auxillaries.h"

using namespace alias;
using namespace Eigen;
//using Ipopt_Wrapper::Ipopt_Solver;
//using Ipopt_Wrapper::Option;
//using Ipopt_Wrapper::Ipopt_Status;
//using Ipopt_Wrapper::Simple_Opt_Problem;
//
//
//
/***********************************************************************************************
* Generate grid for X with each probability. X follows a normal distribution
***********************************************************************************************/
ArrayXd alias::GenerateNormalGrid(const size_t N)
{
    if (N == 0) throw std::invalid_argument("N must be > 0");

    boost::math::normal dist(0.0, 1.0); // mu=0, sigma=1
    ArrayXd X(static_cast<Eigen::Index>(N));

    const double invN = 1.0 / static_cast<double>(N);
    for (Eigen::Index i = 0; i < X.size(); ++i) {
        const double p = (static_cast<double>(i) + 0.5) * invN; // equal-prob bins
        X(i) = boost::math::quantile(dist, p);
    }

    return X;
}

/***********************************************************************************************
* I have a grid X. I want to assign Bin probability around each grid point. X follows normal distribution.
***********************************************************************************************/
ArrayXd alias::NormalGridProbabilities(const ArrayXd& x, double mu, double sigma) {
    if (x.size() == 0) throw std::invalid_argument("x is empty");
    if (sigma <= 0.0) throw std::invalid_argument("sigma must be > 0");

    const Eigen::Index n = x.size();
    ArrayXd p = ArrayXd::Zero(n);
    boost::math::normal dist(mu, sigma);

    if (n == 1) { p(0) = 1.0; return p; }

    // midpoints between grid points
    ArrayXd mid(n - 1);
    for (Eigen::Index i = 0; i < n - 1; ++i) {
        mid(i) = 0.5 * (x(i) + x(i + 1));
    }

    p(0) = boost::math::cdf(dist, mid(0));
    for (Eigen::Index i = 1; i < n - 1; ++i) {
        p(i) = boost::math::cdf(dist, mid(i)) - boost::math::cdf(dist, mid(i - 1));
    }
    p(n - 1) = 1.0 - boost::math::cdf(dist, mid(n - 2));

    return p; // sums to ~1
}

/***********************************************************************************************
* I have an random variable X. X_t = \gamma_0 + \gamma_1 X_{t-1} + \varepsilon. I have a grid for X. Please generate a transition probability matrix for me
***********************************************************************************************/
MatrixXd alias::GenerateTransitionMatrixAR1(const ArrayXd& x_grid, double gamma0, double gamma1, double sigma_eps)
{
    const Eigen::Index n = x_grid.size();
    if (n == 0) throw std::invalid_argument("x_grid is empty");
    if (sigma_eps <= 0.0) throw std::invalid_argument("sigma_eps must be > 0");

    for (Eigen::Index k = 1; k < n; ++k) {
        if (x_grid(k) <= x_grid(k - 1)) {
            throw std::invalid_argument("x_grid must be strictly increasing");
        }
    }

    MatrixXd P = MatrixXd::Zero(n, n);
    if (n == 1) {
        P(0, 0) = 1.0;
        return P;
    }

    // Bin boundaries: midpoints between adjacent grid points
    ArrayXd b(n - 1);
    for (Eigen::Index j = 0; j < n - 1; ++j) {
        b(j) = 0.5 * (x_grid(j) + x_grid(j + 1));
    }

    for (Eigen::Index i = 0; i < n; ++i) {
        const double mu_cond = gamma0 + gamma1 * x_grid(i);
        boost::math::normal dist(mu_cond, sigma_eps);

        P(i, 0) = boost::math::cdf(dist, b(0));
        for (Eigen::Index j = 1; j < n - 1; ++j) {
            P(i, j) = boost::math::cdf(dist, b(j)) - boost::math::cdf(dist, b(j - 1));
        }
        P(i, n - 1) = 1.0 - boost::math::cdf(dist, b(n - 2));

        // Optional numeric cleanup
        double rowsum = P.row(i).sum();
        if (rowsum > 0.0) P.row(i) /= rowsum;
    }

    return P;
}

/***********************************************************************************************
* Normal distribution
***********************************************************************************************/
ArrayXd alias::normCDF_vec(const ArrayXd& X) {
    static const boost::math::normal stdn(0.0, 1.0);
    return X.unaryExpr([](double x) { return boost::math::cdf(stdn, x); });
}

/***********************************************************************************************
* log normal distribution Interval CDF
***********************************************************************************************/
ArrayXd alias::lognormalEY_lessX(const ArrayXd& X, double mu, double sigma) {
    boost::math::normal stdn(0.0, 1.0);
    const double scale = std::exp(mu + 0.5 * sigma * sigma);

    return X.unaryExpr([&](double x) {
        if (x <= 0.0) return std::numeric_limits<double>::quiet_NaN(); // P(Y<x)=0
        double z1 = (std::log(x) - mu - sigma * sigma) / sigma;
        double z0 = (std::log(x) - mu) / sigma;
        double p0 = boost::math::cdf(stdn, z0);
        if (p0 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return scale * boost::math::cdf(stdn, z1) / p0;
    });
}

/***********************************************************************************************
* Linear Spine : find the coefficient
* x1,y1,x2,y2
***********************************************************************************************/
ArrayXd alias::LinearInterpolation1D_Coeff(const double & x1, const double & x2, const double & y1, const double & y2) {
    if (x2 == x1) {
        throw std::invalid_argument("x1 and x2 must be different");
    }

    ArrayXd coef(2);
    const double b = (y2 - y1) / (x2 - x1);
    coef(0) = y1 - b * x1;   // intercept a
    coef(1) = b;             // slope b
    return coef;
}

/***********************************************************************************************
// Solve a*x^alpha + b*x + c = 0, with x > 0.
// Uses Newton in log-space: x = exp(y), so positivity is guaranteed.
***********************************************************************************************/
NewtonResult alias::SolvePowerEqPositiveNewton(const double & a, const double & alpha, const double & b, const double & c,
    const double & x_ini) {

    double tol = 1e-5;
    int max_iter = 100;
    if (x_ini <= 0.0) return {false, x_ini, 0};

    double y = std::log(x_ini);
    const double tiny = 1e-14;

    auto g = [&](double yy) {
        return a * std::exp(alpha * yy) + b * std::exp(yy) + c;
    };
    auto gp = [&](double yy) {
        return a * alpha * std::exp(alpha * yy) + b * std::exp(yy);
    };

    for (int it = 1; it <= max_iter; ++it) {
        double val = g(y);
        if (!std::isfinite(val)) return {false, std::exp(y), it - 1};
        if (std::abs(val) < tol) return {true, std::exp(y), it};

        double der = gp(y);
        if (!std::isfinite(der) || std::abs(der) < tiny) {
            return {false, std::exp(y), it - 1};
        }

        double step = -val / der;
        double y_new = y + step;

        // Damping/backtracking for stability
        bool accepted = false;
        for (int bt = 0; bt < 20; ++bt) {
            double vnew = g(y_new);
            if (std::isfinite(vnew) && std::abs(vnew) < std::abs(val)) {
                accepted = true;
                break;
            }
            step *= 0.5;
            y_new = y + step;
        }
        if (!accepted) return {false, std::exp(y), it - 1};

        if (std::abs(y_new - y) <= tol * (1.0 + std::abs(y_new))) {
            return {true, std::exp(y_new), it};
        }
        y = y_new;
    }

    return {false, std::exp(y), max_iter};
}

/***********************************************************************************************
// Solve a*x^alpha + c = 0, with x > 0.
***********************************************************************************************/
double alias::solve_axalpha_plus_c_eq0(const double & a, const double & alpha, const double & c) {
    if (a == 0.0) throw std::invalid_argument("a must be nonzero");
    if (alpha == 0.0) throw std::invalid_argument("alpha must be nonzero");

    const double rhs = -c / a; // x^alpha = rhs

    // If alpha is not an integer, require rhs > 0 for real solution
    const bool alpha_is_integer = (std::floor(alpha) == alpha);
    if (!alpha_is_integer && rhs <= 0.0) {
        return -1.0;
    }

    // For integer alpha, negative rhs is fine only when 1/alpha leads to real odd root.
    // Here we use pow for real-valued principal branch.
    if (rhs < 0.0) {
        // Conservative: reject negative rhs in this real-solver.
        return -1.0;
    }

    return std::pow(rhs, 1.0 / alpha);
}

/***********************************************************************************************
// Solve: d1*k^a + d0*k^(a-1) + d2 = 0, with k > 0
***********************************************************************************************/
RootResult alias::solve_k_equation_hybrid(const double & d1, const double & d0, const double & d2, const double & a) {
    double k0 = 1000.0;
    double tol = 1e-5;
    int max_iter = 100;

    auto f = [&](double k) {
        return d1 * std::pow(k, a) + d0 * std::pow(k, a - 1.0) + d2;
    };

    auto fp = [&](double k) {
        return d1 * a * std::pow(k, a - 1.0)
             + d0 * (a - 1.0) * std::pow(k, a - 2.0);
    };

    // 1) Find positive bracket [kl, ku] with sign change
    double kl = k0, ku = k0;
    double fl = f(kl), fu = fl;

    bool bracket = false;
    double growth = 2.0;

    for (int i = 0; i < 80; ++i) {
        kl = k0 / std::pow(growth, i + 1);   // shrink left toward 0+
        ku = k0 * std::pow(growth, i + 1);   // expand right
        fl = f(kl);
        fu = f(ku);

        if (std::isfinite(fl) && std::isfinite(fu) && (fl == 0.0 || fu == 0.0 || fl * fu < 0.0)) {
            bracket = true;
            break;
        }
    }

    if (!bracket) {
        return {false, false, k0, 0};
    }
    if (fl == 0.0) return {true, true, kl, 0};
    if (fu == 0.0) return {true, true, ku, 0};

    // 2) Hybrid Newton + bisection inside bracket
    double k = std::sqrt(kl * ku); // geometric midpoint works well on positive domain

    for (int it = 1; it <= max_iter; ++it) {
        double fk = f(k);
        if (!std::isfinite(fk)) {
            k = std::sqrt(kl * ku);
            fk = f(k);
        }

        if (std::abs(fk) < tol) {
            return {true, true, k, it};
        }

        // Try Newton step
        double knew = k;
        double d = fp(k);
        if (std::isfinite(d) && std::abs(d) > 1e-14) {
            double kn = k - fk / d;
            // keep step inside bracket and positive
            if (kn > kl && kn < ku && kn > 0.0 && std::isfinite(kn)) {
                knew = kn;
            } else {
                knew = std::sqrt(kl * ku); // bisection fallback (geometric)
            }
        } else {
            knew = std::sqrt(kl * ku);
        }

        double fnew = f(knew);

        // Update bracket using sign
        if (fl * fnew <= 0.0) {
            ku = knew;
            fu = fnew;
        } else {
            kl = knew;
            fl = fnew;
        }

        k = knew;

        // Bracket-width stopping check
        if (std::abs(ku - kl) <= tol * (1.0 + std::abs(k))) {
            return {true, true, k, it};
        }
    }

    return {false, true, k, max_iter};
}
///***********************************************************************************************
//* normal pdf distribution
//***********************************************************************************************/
//double alias::normPDF(const double & x, const double & mu, const double & sigma) {
//
//    double pi = 3.1415926;
//    double y = 1.0/sigma/sqrt(2.0*pi) * exp(-0.5 * pow( (x-mu)/sigma,2 ) );
//    return y;
//}
//
//ArrayXd alias::normPDF_vec(const ArrayXd & x, const double & mu, const double & sigma) {
//
//    double pi = 3.1415926;
//    int N = x.size();
//    ArrayXd mu_sq_vec = ( (x-mu*ArrayXd::Ones(N))/sigma ).pow(2);
//    ArrayXd y = 1.0/sqrt(2.0*pi)/sigma * exp( -0.5*mu_sq_vec );
//
//    return y;
//}
//
///***********************************************************************************************
//* normal cdf distribution
//***********************************************************************************************/
//double alias::normCDF(const double & x, double mu, double sigma) // Phi(-∞, x) aka N(x)
//{
//    return 1 - erfc((x-mu)/sqrt(2)/sigma)/2;
//}
//
//ArrayXd alias::normCDF_vec(const ArrayXd & x, double mu, double sigma) // Phi(-∞, x) aka N(x)
//{
//    int n = x.size();
//    ArrayXd xvec = (x - mu*ArrayXd::Ones(n))/sigma/sqrt(2);
//
//    ArrayXd y = ArrayXd::Ones(n) - xvec.erfc()/2.0;
//    return y;
//}
//
//ArrayXXd alias::normCDF_mat(const ArrayXd & x, const ArrayXd & mu, const double & sigma) {
//
//    int N_x = x.size();
//    int N_mu = mu.size();
//
//    ArrayXXd mu_mat(N_mu,N_x);
//    mu_mat.rowwise() = x.transpose();
//    mu_mat.colwise() -= mu;
//    mu_mat = mu_mat/sigma/sqrt(2);
//
//    ArrayXXd y = ArrayXXd::Ones(N_mu,N_x) - mu_mat.erfc()/2.0;
//    return y;
//}
//
//
//double alias::normCDF_interval(const double & x1, const double & x2, double mu, double sigma) // Phi(-∞, x) aka N(x)
//{
//    double y;
//    if (x1<1e-20) {
//        y = normCDF(x2, mu, sigma);
//    }
//    else if (x2>1e20) {
//        y = 1 - normCDF(x1, mu, sigma);
//    }
//    else {
//        y = normCDF(x2, mu, sigma) - normCDF(x1, mu, sigma);
//    }
//    return y;
//}
//
//
//double alias::normConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma) {
//
//    double y;
//    if (x1<1e-20) {
//        double temp1 = mu * normCDF_interval(1e-20, x2, mu, sigma);
//        double temp2 = pow(sigma,2) * (normPDF(x2, mu, sigma) - normPDF(1e-20, mu, sigma));
//
//        y = temp1 - temp2;
//    }
//    else if (x2>1e20) {
//        double temp1 = mu * (1.0 - normCDF(x1, mu, sigma));
//        double temp2 = pow(sigma,2) * (0 - normPDF(x1, mu, sigma));
//
//        y = temp1 - temp2;
//    }
//    else {
//        double temp1 = mu * normCDF_interval(x1, x2, mu, sigma);
//        double temp2 = pow(sigma,2) * (normPDF(x2, mu, sigma) - normPDF(x1, mu, sigma));
//
//        y = temp1 - temp2;
//    }
//
//    return y;
//}
//
///***********************************************************************************************
//* log normal distribution Interval CDF
//***********************************************************************************************/
//double alias::lognormConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma) {
//
//    double temp1 = exp(mu+pow(sigma,2)/2.0);
//    double y;
//    if (x1 < 1e-16) {
//        double temp2 = normCDF((log(x2)-mu-pow(sigma,2))/sigma,0,1);
//        y = temp1 * temp2;
//    }
//    else if (x2 > exp(mu + 10*sigma)) {
//        double temp2 = normCDF((mu+pow(sigma,2)-log(x1))/sigma,0,1);
//        y = temp1 * temp2;
//    }
//    else {
//        double temp2 = normCDF((log(x2)-mu-pow(sigma,2))/sigma,0,1)
//                       - normCDF((log(x1)-mu-pow(sigma,2))/sigma,0,1);
//        y = temp1 * temp2;
//    }
//
//    return y;
//}
//
///***********************************************************************************************
//* log normal distribution Interval CDF
//***********************************************************************************************/
//double alias::lognormConditionalEpectation_right_tail(const double & x1, const double & mu,
//    const double & sigma) {
//
//    double temp1 = exp(mu+pow(sigma,2)/2.0);
//    double temp2 = normCDF((mu+pow(sigma,2)-log(x1))/sigma,0,1);
//    double y = temp1 * temp2;
//
//    return y;
//}
//
//double alias::lognormCDF(const double & x1, const double & x2, const double & mu, const double & sigma) {
//
//    double y;
//    if (x1<1e-16) {
//        y = normCDF((log(x2)-mu)/sigma,0,1);
//    }
//    else if (x2 > exp(mu + 10*sigma)) {
//        y = 1.0 - normCDF((log(x1)-mu)/sigma,0,1);
//    }
//    else {
//        y = normCDF((log(x2)-mu)/sigma,0,1) - normCDF((log(x1)-mu)/sigma,0,1);
//    }
//
//    return y;
//}
//

//
///***********************************************************************************************
//* Mis functions
//***********************************************************************************************/
//// find in matlab
//ArrayXi alias::find_matlab(const ArrayXi A) {
//
//    std::vector<int> idxs;
//    for(size_t i = 0; i < A.size(); ++i) {
//        if(A(i)) {
//            idxs.push_back(i);
//        }
//    }
//
//    int* ptr_data = &idxs[0];
//    ArrayXi index = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(idxs.data(), idxs.size());
//
//    return index;
//}
//
//int alias::find_lower_bound_scale(const ArrayXd &vect, double xx) {
//    int yy;
//
//    auto it = lower_bound(vect.cbegin(), vect.cend(), xx);
//    if (it == vect.cend()){
//        yy = ((it-1) - vect.cbegin());
//    }
//    else if (it == vect.cbegin()){
//        yy = 0;
//    }
//    else {
//        yy = it - vect.cbegin()-1;
//    }
//    return yy;
//}
//
//void alias::writeToCSVfile(string name, MatrixXd matrix)
//{
//    ofstream file(name.c_str());
//    file << matrix.format(CSVFormat);
//}
//
//ArrayXXd alias::readCSV(std::string file, int rows, int cols) {
//
//    std::ifstream in(file);
//
//    std::string line;
//
//    int row = 0;
//    int col = 0;
//
//    ArrayXXd res(rows, cols);
//
//    if (in.is_open()) {
//
//        while (std::getline(in, line)) {
//
//            char *ptr = (char *) line.c_str();
//            int len = line.length();
//
//            col = 0;
//
//            char *start = ptr;
//            for (int i = 0; i < len; i++) {
//
//                if (ptr[i] == ',') {
//                    res(row, col++) = atof(start);
//                    start = ptr + i + 1;
//                }
//            }
//            res(row, col) = atof(start);
//
//            row++;
//        }
//
//        in.close();
//    }
//    return res;
//}
