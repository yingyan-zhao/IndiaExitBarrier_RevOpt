/**
 * This header implements two direct search algorithms:
 *  - 'Direct_Search_Simple'
 *  - 'Coordinate_Direct_Search'
 *
 *  Both algorithms have good theoretical fundation. When the objective function is smooth,
 *  'Coordinate_Direct_Search' should perform better than 'Direct_Search_Simple'.
 *
 *  Both algorithms are based on the following paper:
 *
 *  A LINESEARCH-BASED DERIVATIVE-FREE APPROACH FOR NONSMOOTH CONSTRAINED OPTIMIZATION
 *                 by Fasano, Liuzzi, Lucidi and Rinaldi (2014), SIAM J. OPTIM
 */

#ifndef DFS_DIRECT_SEARCH_H
#define DFS_DIRECT_SEARCH_H

#include <cmath>
#include "DFS.h"
#include "multithread_loop.h"

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <numeric>

namespace DFS{

    struct Direct_Search_Option{
        unsigned max_iter = std::numeric_limits<unsigned>::max();
        unsigned max_eval = std::numeric_limits<unsigned>::max();

        double tol = 1e-2;    // if error < tol, return converged
        int display_level = 1; // 0 --- no display
                               // 1 --- display some progress
                               // 2 --- display details

        const double* f_target = nullptr;   // stop when optimal f value is below *f_target (if minimizing)
                                      // stop when optimal f value is above *f_target (if maximizing)

        double theta = 0.5;
        double gamma = 1e-6;
        double eta = 1e-3;
        double delta = 0.5;
    };

    /**
     * Solve the optimization problem using the DFN_simple algorithm in
     *             A LINESEARCH-BASED DERIVATIVE-FREE APPROACH FOR NONSMOOTH CONSTRAINED OPTIMIZATION
     *                     by Fasano, Liuzzi, Lucidi and Rinaldi (2014), SIAM J. OPTIM
     * @tparam F
     * @tparam G
     * @param problem       the optimization problem
     * @param initial_x     initial starting point
     * @param option        option
     * @return  Result (defined in 'DFS.h')
     */
    template<typename F, typename G>
    Result Direct_Search_Simple(Optimization_Problem<F, G> & problem,
                                const std::vector<double> & initial_x,
                                const Direct_Search_Option & option = Direct_Search_Option());

    /**
     * Solve the optimization problem using the CS-DFN algorithm in
     *             A LINESEARCH-BASED DERIVATIVE-FREE APPROACH FOR NONSMOOTH CONSTRAINED OPTIMIZATION
     *                     by Fasano, Liuzzi, Lucidi and Rinaldi (2014), SIAM J. OPTIM
     * @tparam F
     * @tparam G
     * @param problem       the optimization problem
     * @param initial_x     initial starting point
     * @param option        option
     * @return  Result (defined in 'DFS.h')
     */
    template<typename F, typename G>
    Result Coordinate_Direct_Search(Optimization_Problem<F, G> & problem,
                                    const std::vector<double> & initial_x,
                                    const Direct_Search_Option & option = Direct_Search_Option());

    /**
     * Solve the optimization problem using the
     *  Parallelized hybrid optimization methods for nonsmooth problems using NOMAD and linesearch
     *                     by Liuzzi and Truemper, Comp. Appl. Math
     * @tparam F
     * @tparam G
     * @param problem               the optimization problem
     * @param initial_x             initial starting point
     * @param option                option
     * @param threadsManagement     thread pool
     * @param n_threads             maximum threads to be used
     * @return  Result (defined in 'DFS.h')
     */

    template<typename F, typename G>
    Result Coordinate_Direct_Parallel_Search(Optimization_Problem<F, G> & problem,
                                             const std::vector<double> & initial_x,
                                             const Direct_Search_Option & option,
                                             MultiThreads::Threads_Management& threadsManagement,
                                             unsigned n_threads = std::thread::hardware_concurrency());


    /********************************************************************************************************/
    /*                                      Implementation                                                  */
    /********************************************************************************************************/
    struct Quasi_Random_Direction{
        Quasi_Random_Direction(unsigned dim_x);
        std::vector<double> operator ()();

        const unsigned dim_x;
        boost::random::sobol sobol_engine;
        boost::random::uniform_01<double> uniform_dist;

        using sobol_random_gen = boost::variate_generator<boost::random::sobol, boost::random::uniform_01<double>>;
        sobol_random_gen uniform;
    };

    inline Quasi_Random_Direction::Quasi_Random_Direction(unsigned int dim_x)
    : dim_x(dim_x)
    , sobol_engine(dim_x - 1)
    , uniform(sobol_engine, uniform_dist)
    {}

    inline
    std::vector<double> Quasi_Random_Direction::operator()() {
        constexpr double pi = 3.1415926535897932384626433;
        constexpr double twopi = 2 * pi;
        std::vector<double> z(dim_x - 1);
        for (int i = 0; i < dim_x - 1; ++i){
            z[i] = uniform();
        }

        double s = 1;
        std::vector<double> d(dim_x);
        for (uint i = 0; i + 2 < dim_x; ++i){
            double t = pi * z[i];
            d[i] = s * std::cos(t);
            s *= std::sin(t);
        }
        if (dim_x >= 2){
            d[dim_x - 2] = s * std::cos(twopi * z[dim_x - 2]);
            s *= std::sin(twopi * z[dim_x - 2]);
        }
        if (dim_x >= 1){
            d[dim_x - 1] = s;
        }
        return d;
    }

    template<typename F, typename G>
    unsigned projected_continuous_search(Optimization_Problem<F, G> & problem,
                                     const Direct_Search_Option & option,
                                     double f_y,
                                     double alpha_tilde,
                                     const std::vector<double> & y,
                                     const std::vector<double> & p,
                                     double & alpha,
                                     std::vector<double> & x,
                                     double & f_x){
        unsigned n_eval = 0;
        if (y.size() != p.size()){
            throw std::runtime_error("Error!");
        }

        bool go_to_step_4 = false;
        double sign;

        // Step 0:
        alpha = alpha_tilde;

        // Step 1:
        x = plus(product(alpha, p), y);
        project_box_constrain(x, problem);
        f_x = problem.obj_for_minimization(x);
        ++n_eval;
        if (f_x - f_y <=  - option.gamma * alpha * alpha){
            sign = 1;
            go_to_step_4 = true;
        } else {
            // Step 2:
            x = plus(product(-alpha, p), y);
            project_box_constrain(x, problem);
            f_x = problem.obj_for_minimization(x);
            ++n_eval;

            if (f_x - f_y <= - option.gamma * alpha * alpha){
                sign = -1;
                go_to_step_4 = true;
            }
        }

        // Step 3:
        if (!go_to_step_4){
            alpha = 0;
            x = y;
            f_x = f_y;
        } else {
            while (true){
                // Step 4
                double beta = alpha / option.delta;
                std::vector<double> x_temp = plus(product(sign * beta, p), y);
                project_box_constrain(x_temp, problem);
                double f_x_temp = problem.obj_for_minimization(x_temp);
                ++n_eval;

                // Step 5:
                if (f_x_temp - f_y > -option.gamma * beta * beta){
                    break;
                } else {
                    // Step 6:
                    alpha = beta;
                    x = x_temp;
                    f_x = f_x_temp;
                }
            }
        }

        return n_eval;
    }

    template<typename F, typename G>
    Result Direct_Search_Simple(Optimization_Problem<F, G> & problem,
                                const std::vector<double> & initial_x,
                                const Direct_Search_Option & option){
        const unsigned & dim_x = problem.dim_effective_x;
        const std::vector<double> & x_lb = problem.effective_x_lb;
        const std::vector<double> & x_ub = problem.effective_x_ub;

        // Step 0: check if the inputs are valid
        std::vector<double> ini_x = problem.from_x_to_effective_x(initial_x);
        if (ini_x.size() != dim_x){
            throw std::runtime_error("Initial value is of the wrong dimension");
        }
        if (problem.has_nonlinear_constrains){
            throw std::runtime_error("this algorithm can only deal with box constraint.");
        }

        project_box_constrain(ini_x, problem);

        // deal with trivial case
        if (dim_x == 0){
            Result result;
            result.optimal_x = initial_x;
            result.optimal_f = problem.obj_for_minimization(ini_x);
            if (!problem.minimize){
                result.optimal_f *= -1;
            }
            result.flag = converged;
            return result;
        }

        // Step 1: preparation
        double alpha_tilde = 0;
        for (unsigned i = 0; i < dim_x; ++i){
            alpha_tilde += std::max(1e-3, std::min(1.0, std::fabs(ini_x[i]))) / dim_x;
        }
        double f_target = std::numeric_limits<double>::lowest();
        if (option.f_target != nullptr && problem.minimize){
            f_target = *option.f_target;
        }
        if (option.f_target != nullptr && !problem.minimize){
            f_target = - *(option.f_target);
        }

        // Step 2: prepare for quasi_random_direction
        if (dim_x > 3667){
            throw std::runtime_error("the sobol sequence in boost has only been implemented for dim <= 3667!");
        }
        Quasi_Random_Direction quasi_random_direction(dim_x);

        // Step 3: prepare record
        double f_x = problem.obj_for_minimization(ini_x);
        std::vector<double> x = ini_x;
        Flag flag = max_iter_reached;

        std::vector<double> opt_x = ini_x;
        double opt_f = f_x;

        unsigned n_eval = 1;

        // Step 4: start iteration
        for (unsigned k = 0; k < option.max_iter; ++k){
            if (option.display_level >= 1){
                std::cout << "at iteration "<< k << std::endl;
            }
            std::vector<double> d = quasi_random_direction();

            if (n_eval >= option.max_eval){
                flag = max_eval_reached;
                break;
            }

            double alpha;
            double f_x_temp;
            std::vector<double> x_temp(dim_x);
            n_eval += projected_continuous_search(problem, option, f_x, alpha_tilde, x, d, alpha, x_temp, f_x_temp);

            double alpha_tilde_old = alpha_tilde;
            if (alpha == 0){
                alpha_tilde *= option.theta;
                if (option.display_level >= 2){
                    std::cout << "search fail with search direction: " << d[0] << ", " << d[1] << std::endl;
                }
            } else {
                if (option.display_level >= 2){
                    std::cout << "search success with search direction: " << d[0] << ", " << d[1] << std::endl;
                }

                alpha_tilde = alpha;
                x = x_temp;
                f_x = f_x_temp;

                if (f_x < opt_f){
                    opt_f = f_x;
                    opt_x = x;
                }
            }

            bool stop = std::max(alpha_tilde_old, alpha) < option.tol;
            if (stop){
                flag = converged;
                break;
            }

            // further refinement
            auto basis = find_orthogonal_basis(d);
            for (const auto & b : basis){
                double alpha;
                double f_x_temp;
                std::vector<double> x_temp(dim_x);
                n_eval += projected_continuous_search(problem, option, f_x, alpha_tilde, x, d, alpha, x_temp, f_x_temp);

                if (f_x_temp < f_x){
                    f_x = f_x_temp;
                    x = x_temp;
                    if (f_x < opt_f){
                        opt_f = f_x;
                        opt_x = x;
                    }
                }
            }

            if (option.display_level >= 1){
                std::cout << "f = " << f_x << ", error = " << std::max(alpha_tilde_old, alpha) << std::endl;
                std::cout << "--------------------------------------" << std::endl;
            }
            if (opt_f < f_target){
                flag = target_reached;
                break;
            }
        }

        Result result;
        result.optimal_x = problem.from_effective_x_to_x(opt_x);
        result.optimal_f = problem.minimize ? opt_f : -opt_f;
        result.flag = flag;
        result.n_eval = n_eval;


        return result;
    }


    template<typename F, typename G>
    unsigned continuous_search(Optimization_Problem<F, G> & problem,
                                    double alpha_tilde,
                                    const std::vector<double> & y,
                                    const double f_y,
                                    int p,
                                    const Direct_Search_Option & option,
                                    double & alpha,
                                    int & p_plus,
                                    double & f_val){
        unsigned n_eval = 0;
        const auto & dim_x = problem.dim_effective_x;
        const unsigned k = p > 0 ? p - 1 : (-p) - 1; // along which coordinate are we searching on ?
        if (k >= dim_x){
            throw std::runtime_error("Error!");
        }
        const double & x_lb = problem.effective_x_lb[k];
        const double & x_ub = problem.effective_x_ub[k];
        const int p_old = p;

        double alpha_bar;

        // Step 1
        if (p > 0 && x_ub == std::numeric_limits<double>::max()){
            alpha_bar = std::numeric_limits<double>::max();
        } else if (p < 0 && x_lb == std::numeric_limits<double>::min()){
            alpha_bar = std::numeric_limits<double>::max();
        } else {
            alpha_bar = p > 0 ? x_ub - y[k] : y[k] - x_lb;
        }
        alpha = std::min(alpha_bar, alpha_tilde);

        // Step 2
        bool go_to_step6 = false;
        if (alpha > 0){
            auto x = y;
            x[k] += p > 0 ? alpha : -alpha;
            f_val = problem.obj_for_minimization(x);
            ++n_eval;
            if (f_val - f_y <= -option.gamma * alpha * alpha){
                go_to_step6 = true;
            }
        }

        if (!go_to_step6){
            // Step 3:
            p *= -1;
            if (p > 0 && x_ub == std::numeric_limits<double>::max()){
                alpha_bar = std::numeric_limits<double>::max();
            } else if (p < 0 && x_lb == std::numeric_limits<double>::min()){
                alpha_bar = std::numeric_limits<double>::max();
            } else {
                alpha_bar = p > 0 ? x_ub - y[k] : y[k] - x_lb;
            }
            alpha = std::min(alpha_bar, alpha_tilde);
            // Step 4
            if (alpha > 0){
                auto x = y;
                x[k] += p > 0 ? alpha : -alpha;
                f_val = problem.obj_for_minimization(x);
                ++n_eval;
                if (f_val - f_y <= -option.gamma * alpha * alpha){
                    go_to_step6 = true;
                }
            }
        }

        // Step 5
        if (!go_to_step6){
            alpha = 0;
            p_plus = p_old;
            f_val = f_y;
        } else {
            // Step 6
            p_plus = p;
            while (true){
                if (alpha == alpha_bar){
                    break;
                } else {
                    double beta = std::min(alpha_bar, alpha / option.delta);
                    auto x = y;
                    x[k] += p > 0 ? beta : -beta;
                    double f_val_t = problem.obj_for_minimization(x);
                    ++n_eval;
                    if (f_val_t - f_y > -option.gamma * beta * beta){
                        break;
                    }
                    alpha = beta;
                    f_val = f_val_t;
                }
            }
        }

        return n_eval;
    }

    template<typename F, typename G>
    Result Coordinate_Direct_Search(Optimization_Problem<F, G> & problem,
                                    const std::vector<double> & initial_x,
                                    const Direct_Search_Option & option){
        const unsigned & dim_x = problem.dim_effective_x;
        const std::vector<double> & x_lb = problem.effective_x_lb;
        const std::vector<double> & x_ub = problem.effective_x_ub;

        // Step 0: check if the inputs are valid
        std::vector<double> ini_x = problem.from_x_to_effective_x(initial_x);
        if (ini_x.size() != dim_x){
            throw std::runtime_error("Initial value is of the wrong dimension");
        }
        if (problem.has_nonlinear_constrains){
            throw std::runtime_error("this algorithm can only deal with box constraint.");
        }
        project_box_constrain(ini_x, problem);

        // deal with trivial case
        if (dim_x == 0){
            Result result;
            result.optimal_x = initial_x;
            result.optimal_f = problem.obj_for_minimization(ini_x);
            if (!problem.minimize){
                result.optimal_f *= -1;
            }
            result.flag = converged;
            return result;
        }

        // Step 1: preparation
        std::vector<double> alpha_tilde(dim_x + 1);
        alpha_tilde.back() = 0;
        for (unsigned i = 0; i < dim_x; ++i){
            alpha_tilde[i] = std::max(1e-3, std::min(1.0, std::fabs(ini_x[i])));
            alpha_tilde.back() += alpha_tilde[i] / dim_x;
        }

        double f_target = std::numeric_limits<double>::lowest();
        if (option.f_target != nullptr && problem.minimize){
            f_target = *option.f_target;
        }
        if (option.f_target != nullptr && !problem.minimize){
            f_target = - *(option.f_target);
        }

        // Step 2: prepare for quasi_random_direction
        if (dim_x > 3667){
            throw std::runtime_error("the sobol sequence in boost has only been implemented for dim <= 3667!");
        }
        Quasi_Random_Direction quasi_random_direction(dim_x);

        // Step 3: prepare coordinate search
        std::vector<int> search_direction_index(dim_x);
        for (unsigned i = 0; i < dim_x; ++i){
            search_direction_index[i] = i + 1;
        }

        // Step 4: prepare record
        double f_x = problem.obj_for_minimization(ini_x);
        std::vector<double> x = ini_x;
        Flag flag = max_iter_reached;
        unsigned n_eval = 1;

        std::vector<double> opt_x = ini_x;
        double opt_f = f_x;
        auto update_optimal_solution = [&](){
            if (f_x < opt_f){
                opt_f = f_x;
                opt_x = x;
            }
        };

        // Step 5: start iteration
        for (unsigned iter = 0; iter < option.max_iter; ++iter){
            auto alpha_tilde_old = alpha_tilde;
            std::vector<double> alpha(dim_x + 1, 0);

            // search along each coordinates
            for (unsigned i = 0; i < dim_x; ++i){
                int p_t;
                double f_y_t;
                n_eval += continuous_search(problem, alpha_tilde[i], x, f_x, search_direction_index[i], option, alpha[i], p_t, f_y_t);
                if (alpha[i] == 0){
                    alpha_tilde[i] *= option.theta;
                } else {
                    alpha_tilde[i] = alpha[i];
                    search_direction_index[i] = p_t;

                    f_x = f_y_t;
                    x[i] += p_t > 0 ? alpha[i] : -alpha[i];
                }
            }
            update_optimal_solution();
            if (opt_f < f_target){
                flag = target_reached;
                break;
            }

            double max_alpha = *std::max_element(alpha_tilde_old.begin(), alpha_tilde_old.begin() + dim_x);
            max_alpha = std::max(max_alpha, *std::max_element(alpha.begin(), alpha.begin() + dim_x));
            if (max_alpha <= option.eta){
                // conduct quasi_random_search
                std::vector<double> d = quasi_random_direction();
                std::vector<double> y_t(dim_x);
                double f_y_t;
                n_eval += projected_continuous_search(problem, option, f_x, alpha_tilde.back(), x, d, alpha.back(), y_t, f_y_t);
                if (alpha.back() == 0){
                    alpha_tilde.back() *= option.theta;
                } else {
                    alpha_tilde.back() = alpha.back();
                    x = y_t;
                    f_x = f_y_t;

                    update_optimal_solution();
                    if (opt_f < f_target){
                        flag = target_reached;
                        break;
                    }
                }

                // further refinement
                auto basis = find_orthogonal_basis(d);
                for (const auto & b : basis){
                    double alpha_temp;
                    double f_x_temp;
                    std::vector<double> x_temp(dim_x);
                    n_eval += projected_continuous_search(problem, option, f_x, alpha_tilde.back(), x, d, alpha_temp, x_temp, f_x_temp);

                    if (f_x_temp < f_x){
                        f_x = f_x_temp;
                        x = x_temp;
                    }
                }
            }
            update_optimal_solution();

            max_alpha = *std::max_element(alpha.begin(), alpha.end());
            max_alpha = std::max(max_alpha, *std::max_element(alpha_tilde_old.begin(), alpha_tilde_old.end()));
            if (option.display_level >= 1){
                std::cout << "f = " << f_x << ", error = " << max_alpha << std::endl;
                std::cout << "--------------------------------------" << std::endl;
            }
            if (max_alpha < option.tol){
                flag = converged;
                break;
            }
            if (n_eval >= option.max_eval){
                flag = max_eval_reached;
            }

            if (opt_f < f_target){
                flag = target_reached;
                break;
            }
        }

        Result result;
        result.optimal_x = problem.from_effective_x_to_x(opt_x);
        result.optimal_f = problem.minimize ? opt_f : -opt_f;
        result.flag = flag;
        result.n_eval = n_eval;

        return result;
    }


    template<typename F, typename G>
    unsigned projected_continuous_parallel_search(Optimization_Problem<F, G> & problem,
                                         const Direct_Search_Option & option,
                                         double f_y,
                                         double alpha_tilde,
                                         const std::vector<double> & y,
                                         const std::vector<double> & p,
                                         double & alpha,
                                         std::vector<double> & x,
                                         double & f_x,
                                         MultiThreads::Threads_Management & threadsManagement){
        unsigned n_eval = 0;
        if (y.size() != p.size()){
            throw std::runtime_error("Error!");
        }


        // Step 0:
        alpha = alpha_tilde;

        // Step 1:
        std::vector<double> x_1, x_2;
        double f_x_1, f_x_2;

        // search along and against one direction
        auto search_direction_1 = [&](){
            x_1 = plus(product(alpha, p), y);
            project_box_constrain(x_1, problem);
            f_x_1 = problem.obj_for_minimization(x_1);
        };

        auto search_direction_2 = [&](){
            x_2 = plus(product(-alpha, p), y);
            project_box_constrain(x_2, problem);
            f_x_2 = problem.obj_for_minimization(x_2);
        };

        double sign;
        bool go_to_step_4 = false;

        {
            auto dir_2 = threadsManagement.submit_task_only_when_idle(search_direction_2);
            search_direction_1();
            ++n_eval;

            x = x_1;
            f_x = f_x_1;
            if (f_x - f_y <=  - option.gamma * alpha * alpha){
                sign = 1;
                go_to_step_4 = true;
            } else {
                if (dir_2 != nullptr){
                    dir_2->wait_until_done();
                } else {
                    search_direction_2();
                }
                ++n_eval;
                x = x_2;
                f_x = f_x_2;
                if (f_x - f_y <= - option.gamma * alpha * alpha){
                    sign = -1;
                    go_to_step_4 = true;
                }
            }
        }

        // Step 3:
        if (!go_to_step_4){
            alpha = 0;
            x = y;
            f_x = f_y;
        } else {
            while (true){
                // Step 4
                const double beta_1 = alpha / option.delta;
                const double beta_2 = beta_1 / option.delta;

                std::vector<double> xx_1, xx_2;
                double ff_1, ff_2;

                auto extend_1 = [&](){
                    xx_1 = plus(product(sign * beta_1, p), y);
                    project_box_constrain(xx_1, problem);
                    ff_1 = problem.obj_for_minimization(xx_1);
                };

                auto extend_2 = [&](){
                    xx_2 = plus(product(sign * beta_2, p), y);
                    project_box_constrain(xx_2, problem);
                    ff_2 = problem.obj_for_minimization(xx_2);
                };

                auto ext_2 = threadsManagement.submit_task_only_when_idle(extend_2);
                extend_1();
                ++n_eval;

                if (ff_1 - f_y > -option.gamma * beta_1 * beta_1){
                    break;
                } else {
                    if (ext_2 != nullptr){
                        ext_2->wait_until_done();
                    } else {
                        extend_2();
                    }
                    ++n_eval;

                    // Step 5:
                    if (ff_2 - f_y > -option.gamma * beta_2 * beta_2){
                        alpha = beta_1;
                        x = xx_1;
                        f_x = ff_1;
                        break;
                    } else {
                        alpha = beta_2;
                        x = xx_2;
                        f_x = ff_2;
                    }
                }
            }
        }

        return n_eval;
    }

    template<typename F, typename G>
    Result Coordinate_Direct_Parallel_Search(Optimization_Problem<F, G> & problem,
                                            const std::vector<double> & initial_x,
                                            const Direct_Search_Option & option,
                                            MultiThreads::Threads_Management& threadsManagement,
                                            unsigned n_threads){
        const unsigned & dim_x = problem.dim_effective_x;
        const std::vector<double> & x_lb = problem.effective_x_lb;
        const std::vector<double> & x_ub = problem.effective_x_ub;

        // Step 0: check if the inputs are valid
        std::vector<double> ini_x = problem.from_x_to_effective_x(initial_x);
        if (ini_x.size() != dim_x){
            throw std::runtime_error("Initial value is of the wrong dimension");
        }
        if (problem.has_nonlinear_constrains){
            throw std::runtime_error("this algorithm can only deal with box constraint.");
        }
        project_box_constrain(ini_x, problem);

        // deal with trivial case
        if (dim_x == 0){
            Result result;
            result.optimal_x = initial_x;
            result.optimal_f = problem.obj_for_minimization(ini_x);
            if (!problem.minimize){
                result.optimal_f *= -1;
            }
            result.flag = converged;
            return result;
        }

        // Step 1: preparation
        double alpha_tilde_ = 0;
        for (unsigned i = 0; i < dim_x; ++i){
            alpha_tilde_ += std::max(1e-3, std::min(1.0, std::fabs(ini_x[i]))) / dim_x;
        }
        std::vector<double> alpha_tilde(dim_x, alpha_tilde_);

        double f_target = std::numeric_limits<double>::lowest();
        if (option.f_target != nullptr && problem.minimize){
            f_target = *option.f_target;
        }
        if (option.f_target != nullptr && !problem.minimize){
            f_target = - *(option.f_target);
        }

        // Step 2: prepare for quasi_random_direction
        if (dim_x > 3667){
            throw std::runtime_error("the sobol sequence in boost has only been implemented for dim <= 3667!");
        }
        Quasi_Random_Direction quasi_random_direction(dim_x);

        // Step 3: prepare record
        double f_x = problem.obj_for_minimization(ini_x);
        std::vector<double> x = ini_x;
        Flag flag = max_iter_reached;
        unsigned n_eval = 1;

        std::vector<double> opt_x = ini_x;
        double opt_f = f_x;
        auto update_optimal_solution = [&](){
            if (f_x < opt_f){
                opt_f = f_x;
                opt_x = x;
            }
        };

        // Step 5: start iteration
        for (unsigned iter = 0; iter < option.max_iter; ++iter){
            std::vector<double> direction = quasi_random_direction();
            auto basis = find_orthogonal_basis(direction);
            basis.push_back(direction);

            std::vector<std::vector<double>> x_temp(dim_x, std::vector<double>(dim_x));
            std::vector<double> f_x_temp(dim_x);
            std::vector<double> alpha(dim_x);
            std::vector<unsigned> n_evals_temp(dim_x);

            auto search_along_direction = [&](unsigned i){
                n_evals_temp[i] = projected_continuous_search(problem, option, f_x, alpha_tilde[i], x, basis[i], alpha[i], x_temp[i], f_x_temp[i]);
                if (alpha[i] == 0){
                    alpha_tilde[i] *= option.theta;
                } else {
                    alpha_tilde[i] = alpha[i];
                }
            };
            MultiThreads::simple_parallel_for(search_along_direction, dim_x, threadsManagement, n_threads);

            n_eval = std::accumulate(n_evals_temp.begin(), n_evals_temp.end(), n_eval);

            double alpha_sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
            std::vector<double> xx(dim_x, 0);

            if (alpha_sum > 0){
                for (unsigned i = 0; i < basis.size(); ++i){
                    if (alpha[i] > 0){
                        xx = plus(product(alpha[i], x_temp[i]), xx);
                    }
                }
                xx = product(1.0 / alpha_sum, xx);
                double yy = problem.obj_for_minimization(xx);
                x_temp.push_back(xx);
                f_x_temp.push_back(yy);

                auto it = std::min_element(f_x_temp.begin(), f_x_temp.end());
                f_x = *it;
                x = x_temp[it - f_x_temp.begin()];
                ++n_eval;

                update_optimal_solution();
            }

            double max_alpha = *std::max_element(alpha_tilde.begin(), alpha_tilde.end());
            if (max_alpha < option.tol){
                flag = converged;
                break;
            }

            if (n_eval >= option.max_eval){
                flag = max_eval_reached;
                break;
            }

            if (opt_f < f_target){
                flag = target_reached;
                break;
            }

            if (option.display_level >= 1){
                std::cout << "f = " << f_x << ", error = " << max_alpha << std::endl;
                std::cout << "--------------------------------------" << std::endl;
            }
        }

        Result result;
        result.optimal_x = problem.from_effective_x_to_x(opt_x);
        result.optimal_f = problem.minimize ? opt_f : -opt_f;
        result.flag = flag;
        result.n_eval = n_eval;

        return result;
    }






}

#endif //DFS_DIRECT_SEARCH_H
