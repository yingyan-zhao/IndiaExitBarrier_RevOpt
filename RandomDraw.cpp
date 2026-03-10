#include "RandomDraw.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

using namespace Eigen;
using namespace std;
using namespace alias;

namespace {
inline double clamp_unit_interval(double u) {
    constexpr double eps = 1e-12;
    return std::min(1.0 - eps, std::max(eps, u));
}
} // namespace

// Return a 'dim' x 'n_draws' matrix simulated from quasi-random normal distribution
ArrayXXd alias::quasi_random_normal(long dim, long n_draws, double mean, double stderr)
{
    if (dim <= 0) {
        throw std::invalid_argument("quasi_random_normal: dim must be > 0");
    }
    if (n_draws < 0) {
        throw std::invalid_argument("quasi_random_normal: n_draws must be >= 0");
    }
    if (stderr <= 0.0) {
        throw std::invalid_argument("quasi_random_normal: stderr must be > 0");
    }

    boost::random::sobol rng(dim);
    boost::random::uniform_01<double> uniform_dist;
    boost::math::normal normal(mean, stderr);

    ArrayXXd M(dim, n_draws);
    for (long i = 0; i < n_draws; ++i)
    {
        for (long j = 0; j < dim; ++j)
        {
            const double u = clamp_unit_interval(uniform_dist(rng));
            M(j, i) = boost::math::quantile(normal, u);
        }
    }

    return M;
}


ArrayXXd alias::quasi_random_uniform(long dim, long n_draws)
{
    if (dim <= 0) {
        throw std::invalid_argument("quasi_random_uniform: dim must be > 0");
    }
    if (n_draws < 0) {
        throw std::invalid_argument("quasi_random_uniform: n_draws must be >= 0");
    }

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
    if (stderr <= 0.0) {
        throw std::invalid_argument("quasi_random_uniform_to_normal: stderr must be > 0");
    }

    const Eigen::Index n_rows = M.rows();
    const Eigen::Index n_cols = M.cols();

    boost::math::normal normal(mean, stderr);
    ArrayXXd MM(n_rows, n_cols);

    for (Eigen::Index i = 0; i < n_rows; ++i) {
        for (Eigen::Index j = 0; j < n_cols; ++j) {
            const double u = clamp_unit_interval(M(i, j));
            MM(i, j) = boost::math::quantile(normal, u);
        }
    }
    return MM;
}
