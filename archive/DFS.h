/**
 * This library has the following structure:
 *  1. The problem should first be described as "Optimization Problem" (defined in 'DFS.h')
 *  2. Then, choose solver along with its option to solve the optimization problem.
 *     Right now, algorithms that have been implemented are:
 *     a. direct search algorithm ('direct_search.h') --- local optimization problem
 *     b. SOO algorithm ('soo.h')                     --- global optimization problem
 *  3. The solver will then return the optimization result as 'Result' (defined in 'DFS.h')
 */

#pragma once
#ifndef DFS_DFS_H
#define DFS_DFS_H

#include <vector>
#include <typeinfo>
#include <string>

namespace DFS {

    /**
     * The following structure provides a trivial inequality restriction
     */
    struct Trivial_G {
        std::vector<double> operator()(const std::vector<double> & x){return std::vector<double>();};
    };
    static Trivial_G trivial_g;

    /**
     * The following structure describes the following optimization problem:
     *           max / min  f(x)
     *                s.t.  g(x) <= 0
     *                      lb <= x <= ub
     *  where x is a vector of dimension dim_x.
     *
     *  It's implementation is based on two templates:
     * @tparam F        type of objective function. If f is of type F, then the following must be defined:
     *                  double f(const std::vector<double> & x)
     *
     * @tparam G        type of constrains. If g is of type G, then the following must be defined:
     *                  std::vector<double> g(const std::vector<double> & x)
     *                  where the return of g is of 'dim_g' dimension. If 'dim_g = 0', it means that there are no
     *                  effective constrains in 'g'.
     *
     *  The whole structure has four functions to be used by the user:
     *  1. constructor based on (dim_x, f, g, dim_g)
     *  2. a copy constructor
     *  3. 'set_x_lb' which sets the lower bound of x. If this function is not called, x is assumed to have no
     *      lower bound.
     *  4. 'set_x_ub' which sets the upper bound of x. If this function is not called, x is assumed to have no
     *      upper bound.
     *
     *  We allow for lb[i] = ub[i] for some i, then this constrain on x[i] is considered as x[i] = lb[i] = ub[i].
     *  This is allowed by the framework.
     *
     */
    template<typename F, typename G = decltype(trivial_g)>
    struct Optimization_Problem {
        Optimization_Problem(unsigned dim_x, F &f, G &g = trivial_g, unsigned dim_g = 0);
        Optimization_Problem(Optimization_Problem & problem);

        void set_x_lb(const std::vector<double> & lb);
        void set_x_ub(const std::vector<double> & ub);

        /* ************************************************************************************
         *  Implementation Details:
         * ************************************************************************************/

        void to_be_minimized();
        void to_be_maximized();

        F &f_;
        G &g_;

        bool minimize;
        unsigned dim_effective_x;
        std::vector<double> effective_x_lb;
        std::vector<double> effective_x_ub;
        unsigned dim_g;

        static inline constexpr bool has_nonlinear_constrains = !std::is_same<G, Trivial_G>::value;

        double obj_for_minimization(const std::vector<double> & effective_x);
        double obj_for_maximization(const std::vector<double> & effective_x);


        std::vector<double> g(const std::vector<double> & x);

        std::vector<double> from_effective_x_to_x(const std::vector<double> & effective_x);
        std::vector<double> from_x_to_effective_x(const std::vector<double> & x);

    private:
        unsigned dim_x;
        std::vector<double> x_lb;
        std::vector<double> x_ub;
    };


    /**
     * The optimization result
     */

    enum Flag {converged, max_iter_reached, max_eval_reached, target_reached, other};
    struct Result{
        std::vector<double> optimal_x;
        double optimal_f;
        unsigned n_eval = 0;
        Flag flag;
    };

    template<typename F, typename G>
    void project_box_constrain(std::vector<double> & effective_x,  const Optimization_Problem<F, G> & problem);

/************************************************************************************************/
/*                           Implementation                                                     */
/************************************************************************************************/

template<typename F, typename G>
Optimization_Problem<F, G>::Optimization_Problem(unsigned int dim_x, F &f, G &g, unsigned int dim_g)
: f_(f)
, g_(g)
, minimize(true)
, dim_effective_x(dim_x)
, effective_x_lb(dim_x, std::numeric_limits<double>::lowest())
, effective_x_ub(dim_x, std::numeric_limits<double>::max())
, dim_g(dim_g)
, dim_x(dim_x)
, x_lb(dim_x, std::numeric_limits<double>::lowest())
, x_ub(dim_x, std::numeric_limits<double>::max())
{
    if ((dim_g == 0) != std::is_same<G, Trivial_G>::value){
        throw std::runtime_error("Error! g is provided but dim_g = 0. ");
    }
}

template<typename F, typename G>
Optimization_Problem<F, G>::Optimization_Problem(Optimization_Problem<F, G> &problem)
: f_(problem.f_)
, g_(problem.g_)
, minimize(problem.minimize)
, dim_effective_x(problem.dim_effective_x)
, effective_x_lb(problem.effective_x_lb)
, effective_x_ub(problem.effective_x_ub)
, dim_g(problem.dim_g)
, dim_x(problem.dim_x)
, x_lb(problem.x_lb)
, x_ub(problem.x_ub)
{}



template<typename F, typename G>
void Optimization_Problem<F, G>::to_be_maximized() {
    minimize = false;
}

template<typename F, typename G>
void Optimization_Problem<F, G>::to_be_minimized() {
    minimize = true;
}

template<typename F, typename G>
double Optimization_Problem<F, G>::obj_for_minimization(const std::vector<double> &x) {
    if (x.size() != dim_effective_x){
        throw std::runtime_error("Wrong dimension!");
    }
    double offset = minimize ? 1 : -1;
    return  dim_x == dim_effective_x ? offset * f_(x) : offset * f_(from_effective_x_to_x(x));
}

    template<typename F, typename G>
    double Optimization_Problem<F, G>::obj_for_maximization(const std::vector<double> &x) {
        if (x.size() != dim_effective_x){
            throw std::runtime_error("Wrong dimension!");
        }
        double offset = minimize ? -1 : 1;
        return  dim_x == dim_effective_x ? offset * f_(x) : offset * f_(from_effective_x_to_x(x));
    }



template<typename F, typename G>
std::vector<double> Optimization_Problem<F, G>::g(const std::vector<double> &x) {
    if (x.size() != dim_effective_x){
        throw std::runtime_error("Wrong dimension!");
    }
    return dim_x == dim_effective_x ? g_(x) : g_(from_effective_x_to_x(x));
}

template<typename F, typename G>
void Optimization_Problem<F, G>::set_x_lb(const std::vector<double> &lb) {
    if (lb.size() != dim_x){
        throw std::runtime_error("wrong dimension in lb");
    }
    x_lb = lb;

    dim_effective_x = 0;
    effective_x_lb.clear();
    effective_x_ub.clear();
    for (unsigned i = 0; i < dim_x; ++i){
        if (x_lb[i] < x_ub[i]){
            effective_x_lb.push_back(x_lb[i]);
            effective_x_ub.push_back(x_ub[i]);
            ++dim_effective_x;
        } else if (x_lb[i] > x_ub[i]){
            throw std::runtime_error("lower bound cannot be larger than upper bound!");
        }
    }
}

template<typename F, typename G>
void Optimization_Problem<F, G>::set_x_ub(const std::vector<double> &ub) {
    if (ub.size() != dim_x){
        throw std::runtime_error("wrong dimension in lb");
    }
    x_ub = ub;
    dim_effective_x = 0;
    effective_x_lb.clear();
    effective_x_ub.clear();
    for (unsigned i = 0; i < dim_x; ++i){
        if (x_lb[i] < x_ub[i]){
            effective_x_lb.push_back(x_lb[i]);
            effective_x_ub.push_back(x_ub[i]);
            ++dim_effective_x;
        } else if (x_lb[i] > x_ub[i]){
            throw std::runtime_error("lower bound cannot be larger than upper bound!");
        }
    }
}

template<typename F, typename G>
std::vector<double> Optimization_Problem<F, G>::from_effective_x_to_x(const std::vector<double> &effective_x) {
    if (effective_x.size() != dim_effective_x){
        throw std::runtime_error("Error! Wrong dimension in effective_x");
    }
    std::vector<double> x(dim_x);
    unsigned count = 0;
    for (unsigned i = 0; i < dim_x; ++i){
        x[i] = x_lb[i] < x_ub[i] ? effective_x[count++] : x_lb[i];
    }
    if (count != dim_effective_x){
        throw std::runtime_error("Error!");
    }
    return x;
}

template<typename F, typename G>
std::vector<double> Optimization_Problem<F, G>::from_x_to_effective_x(const std::vector<double> &x) {
    if (x.size() != dim_x){
        throw std::runtime_error("Error! Wrong dimension in x. dim(x) = " + std::to_string(x.size()) + ". it should be  " + std::to_string(dim_x));
    }
    std::vector<double> effective_x;
    effective_x.reserve(dim_effective_x);
    for (unsigned i = 0; i < dim_x; ++i){
        if (x_lb[i] < x_ub[i]){
            effective_x.push_back(x[i]);
        }
    }
    return effective_x;
}

    template<typename F, typename G>
    void project_box_constrain(std::vector<double> & effective_x,  const Optimization_Problem<F, G> & problem){
        if (effective_x.size() != problem.dim_effective_x){
            throw std::runtime_error("Error!");
        }

        for (unsigned i = 0; i < problem.dim_effective_x; ++i){
            if (effective_x[i] < problem.effective_x_lb[i]){
                effective_x[i] = problem.effective_x_lb[i];
            } else if (effective_x[i] > problem.effective_x_ub[i]){
                effective_x[i] = problem.effective_x_ub[i];
            }
        }
    }

    inline
    std::vector<double> product(double a, const std::vector<double> & x){
        std::vector<double> y = x;
        for (auto & yy : y ){
            yy *= a;
        }
        return y;
    }

    inline
    std::vector<double> product(const std::vector<double> & x, double a){
        return product(a, x);
    }

    inline
    std::vector<double> plus(const std::vector<double> & x, const std::vector<double> & y){
        std::vector<double> z(x.size());
        for (unsigned i = 0; i < x.size(); ++i){
            z[i] = x[i] + y[i];
        }
        return z;
    }

    inline
    std::vector<std::vector<double>> find_orthogonal_basis(std::vector<double> & v){
        std::vector<double> w = v;
        if (v[0] > 0){
            w[0] += 1;
        } else {
            w[0] -= 1;
        }

        double w_product = 0;
        for (auto & a : w){
            w_product += a * a;
        }

        const double s = -2 / w_product;
        std::vector<double> Q(v.size() * v.size());
        for (unsigned i = 0; i < v.size(); ++i){
            for (unsigned j = 0; j < v.size(); ++j){
                Q[j + i * v.size()] = s * w[j] * w[i];
            }
            Q[i + i * v.size()] += 1;
        }

        std::vector<std::vector<double>> basis;
        basis.reserve(v.size() - 1);
        for (unsigned i = 1; i < v.size(); ++i){
            basis.emplace_back(Q.begin() + i * v.size(), Q.begin() + (i + 1) * v.size());
        }
        return basis;
    }

    inline void
    check(bool should_be_true, const std::string& message) {
        if (not should_be_true)
        {
            throw std::runtime_error(message);
        }
    }
}


#endif //DFS_DFS_H
