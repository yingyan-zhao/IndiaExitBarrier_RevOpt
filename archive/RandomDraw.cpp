#include "RandomDraw.h"

using namespace Eigen;
using namespace std;
using namespace alias;

// Return a 'dim' x 'n_draws' matrix simulated from quasi-random normal distribution
ArrayXXd alias::quasi_random_normal(long dim, long n_draws, double mean, double stderr)
{
    boost::random::sobol rng(dim);
    boost::random::uniform_01<double> uniform_dist;
    boost::math::normal normal(mean, stderr);

    ArrayXXd M(dim, n_draws);
    for (long i = 0; i < n_draws; ++i)
    {
        ArrayXd v(dim);
        do
        {
            for (long j = 0; j < dim; ++j)
            {
                v(j) = uniform_dist(rng);
            }
        }
        while (v.minCoeff() < 1e-30);

        for (long j = 0; j < dim; ++j)
        {
            M(j, i) = boost::math::quantile(normal, v[j]);
        }
    }

    return M;
}


ArrayXXd alias::quasi_random_uniform(long dim, long n_draws)
{
    boost::random::sobol rng(dim);
    boost::random::uniform_01<double> uniform_dist;

    ArrayXXd M(n_draws, dim);
    for (long i = 0; i < n_draws; ++i)
    {
        for (long j = 0; j < dim; ++j)
        {
            M(i, j) = uniform_dist(rng);
        }
    }

    return M;
}


// Return a 'dim' x 'n_draws' matrix simulated from quasi-random normal distribution
ArrayXXd alias::quasi_random_uniform_to_normal(const Eigen::Ref<const ArrayXXd>& M,
                                               double mean, double stderr) {
    const Eigen::Index n_rows = M.rows();
    const Eigen::Index n_cols = M.cols();

    boost::math::normal normal(mean, stderr);
    ArrayXXd MM(n_rows, n_cols);

    for (Eigen::Index i = 0; i < n_rows; ++i) {
        for (Eigen::Index j = 0; j < n_cols; ++j) {
            MM(i, j) = boost::math::quantile(normal, M(i, j));
        }
    }
    return MM;
}
