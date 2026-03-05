#include "Auxillaries.h"

using namespace alias;
using namespace Eigen;
using Ipopt_Wrapper::Ipopt_Solver;
using Ipopt_Wrapper::Option;
using Ipopt_Wrapper::Ipopt_Status;
using Ipopt_Wrapper::Simple_Opt_Problem;



/***********************************************************************************************
* normal pdf distribution
***********************************************************************************************/
double alias::normPDF(const double & x, const double & mu, const double & sigma) {

    double pi = 3.1415926;
    double y = 1.0/sigma/sqrt(2.0*pi) * exp(-0.5 * pow( (x-mu)/sigma,2 ) );
    return y;
}

ArrayXd alias::normPDF_vec(const ArrayXd & x, const double & mu, const double & sigma) {

    double pi = 3.1415926;
    int N = x.size();
    ArrayXd mu_sq_vec = ( (x-mu*ArrayXd::Ones(N))/sigma ).pow(2);
    ArrayXd y = 1.0/sqrt(2.0*pi)/sigma * exp( -0.5*mu_sq_vec );

    return y;
}

/***********************************************************************************************
* normal cdf distribution
***********************************************************************************************/
double alias::normCDF(const double & x, double mu, double sigma) // Phi(-∞, x) aka N(x)
{
    return 1 - erfc((x-mu)/sqrt(2)/sigma)/2;
}

ArrayXd alias::normCDF_vec(const ArrayXd & x, double mu, double sigma) // Phi(-∞, x) aka N(x)
{
    int n = x.size();
    ArrayXd xvec = (x - mu*ArrayXd::Ones(n))/sigma/sqrt(2);

    ArrayXd y = ArrayXd::Ones(n) - xvec.erfc()/2.0;
    return y;
}

ArrayXXd alias::normCDF_mat(const ArrayXd & x, const ArrayXd & mu, const double & sigma) {

    int N_x = x.size();
    int N_mu = mu.size();

    ArrayXXd mu_mat(N_mu,N_x);
    mu_mat.rowwise() = x.transpose();
    mu_mat.colwise() -= mu;
    mu_mat = mu_mat/sigma/sqrt(2);

    ArrayXXd y = ArrayXXd::Ones(N_mu,N_x) - mu_mat.erfc()/2.0;
    return y;
}


double alias::normCDF_interval(const double & x1, const double & x2, double mu, double sigma) // Phi(-∞, x) aka N(x)
{
    double y;
    if (x1<1e-20) {
        y = normCDF(x2, mu, sigma);
    }
    else if (x2>1e20) {
        y = 1 - normCDF(x1, mu, sigma);
    }
    else {
        y = normCDF(x2, mu, sigma) - normCDF(x1, mu, sigma);
    }
    return y;
}


double alias::normConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma) {

    double y;
    if (x1<1e-20) {
        double temp1 = mu * normCDF_interval(1e-20, x2, mu, sigma);
        double temp2 = pow(sigma,2) * (normPDF(x2, mu, sigma) - normPDF(1e-20, mu, sigma));

        y = temp1 - temp2;
    }
    else if (x2>1e20) {
        double temp1 = mu * (1.0 - normCDF(x1, mu, sigma));
        double temp2 = pow(sigma,2) * (0 - normPDF(x1, mu, sigma));

        y = temp1 - temp2;
    }
    else {
        double temp1 = mu * normCDF_interval(x1, x2, mu, sigma);
        double temp2 = pow(sigma,2) * (normPDF(x2, mu, sigma) - normPDF(x1, mu, sigma));

        y = temp1 - temp2;
    }

    return y;
}

/***********************************************************************************************
* log normal distribution Interval CDF
***********************************************************************************************/
double alias::lognormConditionalEpectation(const double & x1, const double & x2, const double & mu, const double & sigma) {

    double temp1 = exp(mu+pow(sigma,2)/2.0);
    double y;
    if (x1 < 1e-16) {
        double temp2 = normCDF((log(x2)-mu-pow(sigma,2))/sigma,0,1);
        y = temp1 * temp2;
    }
    else if (x2 > exp(mu + 10*sigma)) {
        double temp2 = normCDF((mu+pow(sigma,2)-log(x1))/sigma,0,1);
        y = temp1 * temp2;
    }
    else {
        double temp2 = normCDF((log(x2)-mu-pow(sigma,2))/sigma,0,1)
                       - normCDF((log(x1)-mu-pow(sigma,2))/sigma,0,1);
        y = temp1 * temp2;
    }

    return y;
}

/***********************************************************************************************
* log normal distribution Interval CDF
***********************************************************************************************/
double alias::lognormConditionalEpectation_right_tail(const double & x1, const double & mu,
    const double & sigma) {

    double temp1 = exp(mu+pow(sigma,2)/2.0);
    double temp2 = normCDF((mu+pow(sigma,2)-log(x1))/sigma,0,1);
    double y = temp1 * temp2;

    return y;
}

double alias::lognormCDF(const double & x1, const double & x2, const double & mu, const double & sigma) {

    double y;
    if (x1<1e-16) {
        y = normCDF((log(x2)-mu)/sigma,0,1);
    }
    else if (x2 > exp(mu + 10*sigma)) {
        y = 1.0 - normCDF((log(x1)-mu)/sigma,0,1);
    }
    else {
        y = normCDF((log(x2)-mu)/sigma,0,1) - normCDF((log(x1)-mu)/sigma,0,1);
    }

    return y;
}

/***********************************************************************************************
* Generate normal grid
***********************************************************************************************/
tuple<ArrayXd,ArrayXd> alias::GenerateNormalGrid(const size_t N)
{

    boost::math::normal normal(0, 1);

    ArrayXd vec_norm_phi(N);
    ArrayXd vec_norm_phi_cut(N-1);

    double interval = 1/float(N);

    vec_norm_phi(0) = boost::math::quantile(normal, interval/2);
    double cut = interval/2;
    double cutcut = 0;
    for (long j = 1; j < N; ++j)
    {
        cutcut = cutcut + interval;
        vec_norm_phi_cut(j-1) = boost::math::quantile(normal, cutcut);

        cut = cut + interval;
        vec_norm_phi(j) = boost::math::quantile(normal, cut);
    }
    return tuple<ArrayXd,ArrayXd>(vec_norm_phi,vec_norm_phi_cut);
}

/***********************************************************************************************
* Mis functions
***********************************************************************************************/
// find in matlab
ArrayXi alias::find_matlab(const ArrayXi A) {

    std::vector<int> idxs;
    for(size_t i = 0; i < A.size(); ++i) {
        if(A(i)) {
            idxs.push_back(i);
        }
    }

    int* ptr_data = &idxs[0];
    ArrayXi index = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(idxs.data(), idxs.size());

    return index;
}

int alias::find_lower_bound_scale(const ArrayXd &vect, double xx) {
    int yy;

    auto it = lower_bound(vect.cbegin(), vect.cend(), xx);
    if (it == vect.cend()){
        yy = ((it-1) - vect.cbegin());
    }
    else if (it == vect.cbegin()){
        yy = 0;
    }
    else {
        yy = it - vect.cbegin()-1;
    }
    return yy;
}

void alias::writeToCSVfile(string name, MatrixXd matrix)
{
    ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
}

ArrayXXd alias::readCSV(std::string file, int rows, int cols) {

    std::ifstream in(file);

    std::string line;

    int row = 0;
    int col = 0;

    ArrayXXd res(rows, cols);

    if (in.is_open()) {

        while (std::getline(in, line)) {

            char *ptr = (char *) line.c_str();
            int len = line.length();

            col = 0;

            char *start = ptr;
            for (int i = 0; i < len; i++) {

                if (ptr[i] == ',') {
                    res(row, col++) = atof(start);
                    start = ptr + i + 1;
                }
            }
            res(row, col) = atof(start);

            row++;
        }

        in.close();
    }
    return res;
}