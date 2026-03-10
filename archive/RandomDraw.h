#ifndef RandomDraw_Guard_H
#define RandomDraw_Guard_H

#include <iostream>
#include <Eigen/Dense>
#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/distributions/normal.hpp>


namespace alias {
    using ArrayXd = Eigen::ArrayXd;
    using ArrayXXd = Eigen::ArrayXXd;

    // Return a 'dim' x 'n_draws' matrix simulated from quasi-random normal distribution
    ArrayXXd quasi_random_normal(long dim, long n_draws, double mean, double stderr);

    // Return a 'dim' x 'n_draws' matrix simulated from quasi-random uniform distribution
    ArrayXXd quasi_random_uniform(long dim, long n_draws);

    ArrayXXd alias::quasi_random_uniform_to_normal(const Eigen::Ref<const ArrayXXd>& M, double mean, double stderr);
}

#endif