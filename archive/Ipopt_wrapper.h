#pragma once
#ifndef EXAMPLE_IPOPT_WRAPPER_H
#define EXAMPLE_IPOPT_WRAPPER_H

#include <vector>
#include <limits>
#include <stdexcept>

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"


namespace Ipopt_Wrapper{
    using Ipopt::Index;
    using Ipopt::Number;
    using Ipopt::SolverReturn;
    using Ipopt::IpoptData;
    using Ipopt::IpoptCalculatedQuantities;

/* **************************************************************************
 *  Part 1: Optimization Problem
 * **************************************************************************/

    enum Ipopt_Status {success, acceptable, max_iter_reached, max_cpu_time_reached, other_failure};

    struct Option{
        double tol = 1e-8;

        // if the real error < 'acceptable_tol' in at least 'acceptable_iters' iterations,
        // then the solver will return as 'acceptable'
        int acceptable_iters = 100;
        double acceptable_tol = 1e-6;

        double max_cpu_time = 100000; // in seconds
        int max_iter = 1000;

        int display_level = 5; // choose from 0 to 12, the higher the more detailed print

        int test_derivative = 1;    // 0 --- no test,
                                    // 1 --- check graident at the initial value
                                    // 2 --- check hessian at the initial value

    };

    struct Result{
        double opt_obj;
        std::vector<double> opt_x;
        std::vector<double> opt_grad;
        Ipopt_Status status;

        std::vector<double> multiplier_x_lb;
        std::vector<double> multiplier_x_ub;
    };

    template<typename F>
    class Simple_Opt_Problem {
    public:
        Simple_Opt_Problem(F & f, int dim_x);
        void set_x_lb(const std::vector<double> & lb);
        void set_x_ub(const std::vector<double> & ub);

        double obj_val(const std::vector<double> & x);
        const std::vector<double> & gradient(const std::vector<double> & x);
        const std::vector<double> & hessian(const std::vector<double> & x);

        double last_obj_val() const;
        const std::vector<double> & last_gradient() const;
        const std::vector<double> & last_hessian() const;

        bool same_as_last_x(const std::vector<double> & x);
        void eval(const std::vector<double> & x);

        F & f;
        int dim_x;
        std::vector<double>  x_lb;
        std::vector<double>  x_ub;
        const bool use_hassian;

    private:
        std::vector<double> last_x;
        double last_obj;
        std::vector<double> last_grad;
        std::vector<double> last_hassian;
    };

    class Ipopt_Solver{
    public:
        inline Ipopt_Solver(const Option & option = Option());

        template<typename F>
        Result  optimize(Simple_Opt_Problem<F> & opt, const std::vector<double> & initial_value = std::vector<double>(0));

        Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
    };

/* **************************************************************************
 *  Part 2: Ipopt Object
 * **************************************************************************/
    template<typename F>
    class Simple_Ipopt : public Ipopt::TNLP{
    public:
        Simple_Ipopt(Simple_Opt_Problem<F> & problem, Result & result, const std::vector<double> & starting_point = std::vector<double>(0));

        virtual bool get_nlp_info(
                Index&          n,
                Index&          m,
                Index&          nnz_jac_g,
                Index&          nnz_h_lag,
                IndexStyleEnum& index_style
        );

        virtual bool get_bounds_info(
                Index   n,
                Number* x_l,
                Number* x_u,
                Index   m,
                Number* g_l,
                Number* g_u
        );

        virtual bool get_starting_point(
                Index   n,
                bool    init_x,
                Number* x,
                bool    init_z,
                Number* z_L,
                Number* z_U,
                Index   m,
                bool    init_lambda,
                Number* lambda
        );

        /** Method to return the objective value */
        virtual bool eval_f(
                Index         n,
                const Number* x,
                bool          new_x,
                Number&       obj_value
        );

        virtual bool eval_grad_f(
                Index         n,
                const Number* x,
                bool          new_x,
                Number*       grad_f
        );

        virtual bool eval_g(
                Index         n,
                const Number* x,
                bool          new_x,
                Index         m,
                Number*       g
        );

        virtual bool eval_jac_g(
                Index         n,
                const Number* x,
                bool          new_x,
                Index         m,
                Index         nele_jac,
                Index*        iRow,
                Index*        jCol,
                Number*       values
        );

        virtual bool eval_h(
                Index         n,
                const Number* x,
                bool          new_x,
                Number        obj_factor,
                Index         m,
                const Number* lambda,
                bool          new_lambda,
                Index         nele_hess,
                Index*        iRow,
                Index*        jCol,
                Number*       values
        );


        virtual void finalize_solution(
                SolverReturn               status,
                Index                      n,
                const Number*              x,
                const Number*              z_L,
                const Number*              z_U,
                Index                      m,
                const Number*              g,
                const Number*              lambda,
                Number                     obj_value,
                const IpoptData*           ip_data,
                IpoptCalculatedQuantities* ip_cq
        );

        /** default destructor */
        virtual ~Simple_Ipopt();

    private:
        Simple_Opt_Problem<F> & problem;
        std::vector<double> starting_point;

    public:
        Result & result;
    };

/* **************************************************************************
 *  Part 3: Implementation
 * **************************************************************************/

inline
Ipopt_Solver::Ipopt_Solver(const Option &option)
: app(new Ipopt::IpoptApplication(true, false))
{
        // Set Ipopt option
        app->Options()->SetNumericValue("tol", option.tol);
        app->Options()->SetNumericValue("max_cpu_time", option.max_cpu_time);
        app->Options()->SetIntegerValue("max_iter", option.max_iter);
        app->Options()->SetNumericValue("acceptable_tol", option.acceptable_tol);
        app->Options()->SetIntegerValue("acceptable_iter", option.acceptable_iters);
        app->Options()->SetIntegerValue("print_level", option.display_level);

        if (option.test_derivative == 1){
            app->Options()->SetStringValue("derivative_test", "first-order");
        } else if (option.test_derivative == 2){
            app->Options()->SetStringValue("derivative_test", "second-order");
        } else {
            app->Options()->SetStringValue("derivative_test", "none");
        }
        app->Options()->SetNumericValue("derivative_test_tol", 1e-4);
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-10);

        // start solver
        Ipopt::ApplicationReturnStatus status;
        status = app->Initialize();
        if( status != Ipopt::Solve_Succeeded )
        {
            std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
            throw std::runtime_error("fail to start solver!");
        }
}



    template <typename T>
    struct get_arity : get_arity<decltype(&T::operator())> {};
    template <typename R, typename... Args>
    struct get_arity<R(*)(Args...)> : std::integral_constant<unsigned, sizeof...(Args)> {};
    template <typename R, typename C, typename... Args>
    struct get_arity<R(C::*)(Args...)> : std::integral_constant<unsigned, sizeof...(Args)> {};
    template <typename R, typename C, typename... Args>
    struct get_arity<R(C::*)(Args...) const> : std::integral_constant<unsigned, sizeof...(Args)> {};
// All combinations of variadic/non-variadic, cv-qualifiers and ref-qualifiers

    template<typename F>
    Simple_Opt_Problem<F>::Simple_Opt_Problem(F &f, int dim_x)
            : f(f)
            , dim_x(dim_x)
            , x_lb(dim_x, std::numeric_limits<double>::lowest())
            , x_ub(dim_x, std::numeric_limits<double>::max())
            , use_hassian(get_arity<F>{} == 3)
            , last_x(dim_x, std::numeric_limits<double>::lowest())
            , last_obj(0)
            , last_grad(dim_x, 0)
            , last_hassian(0)
    {
        constexpr size_t n_f_args = get_arity<F>{}; // get number of arguments taken by the f
        if (n_f_args < 2){
            throw std::runtime_error("f return gradient through the 2nd argument");
        }
        static_assert(n_f_args == 2 || n_f_args == 3, "f should have 2 or 3 arguments");
        if (use_hassian){
            last_hassian = std::vector<double>(dim_x * dim_x, 0);
        }
    }

    template<typename F>
    void Simple_Opt_Problem<F>::set_x_lb(const std::vector<double> &lb) {
        if (lb.size() != dim_x){
            throw std::runtime_error("inconsistent dim");
        }
        std::copy(lb.begin(), lb.end(), x_lb.begin());
    }

    template<typename F>
    void Simple_Opt_Problem<F>::set_x_ub(const std::vector<double> &ub) {
        if (ub.size() != dim_x){
            throw std::runtime_error("inconsistent dim");
        }
        std::copy(ub.begin(), ub.end(), x_ub.begin());
    }

    template<typename F>
    const std::vector<double> & Simple_Opt_Problem<F>::last_gradient() const {
        return last_grad;
    }

    template<typename F>
    void Simple_Opt_Problem<F>::eval(const std::vector<double> &x) {
        constexpr size_t n_f_args = get_arity<F>{}; // get number of arguments taken by the f
        last_x = x;
        if constexpr(n_f_args == 2){
            last_obj = f(x, last_grad);
        }
        if constexpr(n_f_args == 3){
            last_obj = f(x, last_grad, last_hassian);
        }
    }

    template<typename F>
    double Simple_Opt_Problem<F>::obj_val(const std::vector<double> &x) {
        if (x.size() != dim_x){
            throw std::runtime_error("error in dim");
        }
        if (!same_as_last_x(x)){
            eval(x);
        }
        return last_obj;
    }

    template<typename F>
    const std::vector<double> & Simple_Opt_Problem<F>::gradient(const std::vector<double> &x) {
        if (x.size() != dim_x){
            throw std::runtime_error("error in dim");
        }
        if (!same_as_last_x(x)){
            eval(x);
        }
        return last_grad;
    }

    template<typename F>
    const std::vector<double> &Simple_Opt_Problem<F>::hessian(const std::vector<double> &x) {
        if (x.size() != dim_x){
            throw std::runtime_error("error in dim");
        }
        if (!use_hassian){
            throw std::runtime_error("hassian is not provided");
        }
        if (!same_as_last_x(x)){
            eval(x);
        }
        return last_hassian;
    }

    template<typename F>
    double Simple_Opt_Problem<F>::last_obj_val() const {
        return last_obj;
    }

    template<typename F>
    const std::vector<double> &Simple_Opt_Problem<F>::last_hessian() const {
        return last_hassian;
    }

    template<typename F>
    bool Simple_Opt_Problem<F>::same_as_last_x(const std::vector<double> &x) {
        return x == last_x;
    }


    template<typename F>
    bool Simple_Ipopt<F>::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
                                       IndexStyleEnum &index_style) {
        n = problem.dim_x;
        m = 0;
        nnz_jac_g = 0;
        nnz_h_lag = (n + 1) * n / 2;
        index_style = TNLP::C_STYLE;
        return true;
    }

    template<typename F>
    bool
    Simple_Ipopt<F>::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
        for (Index i = 0; i < n; ++i){
            x_l[i] = problem.x_lb[i];
        }
        for (Index i = 0; i < n; ++i){
            x_u[i] = problem.x_ub[i];
        }
        return true;
    }

    template<typename F>
    Simple_Ipopt<F>::Simple_Ipopt(Simple_Opt_Problem<F> &problem,
                                  Result & result,
                                  const std::vector<double> &starting_point)
                                                  :problem(problem)
                                                  , starting_point(starting_point)
                                                  , result(result)
    {
        if (starting_point.size() == 0){
            this->starting_point = std::vector<double>(problem.dim_x);
            for (Index i = 0; i < problem.dim_x; ++i){
                this->starting_point[i] = (problem.x_lb[i] + problem.x_ub[i]) / 2;
            }
        }
    }

    template <typename F>
    bool
    Simple_Ipopt<F>::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U,
                                        Index m, bool init_lambda, Number *lambda) {
        if (init_x){
            for (Index i = 0; i < problem.dim_x; ++i){
                x[i] = starting_point[i];
            }
        }
        if (init_z){
            for (Index i = 0; i < problem.dim_x; ++i){
                z_L[i] = 0;
                z_U[i] = 0;
            }
        }
        return true;
    }

    template <typename F>
    bool Simple_Ipopt<F>::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
        std::vector<double> xx(problem.dim_x);
        for (int i = 0; i < problem.dim_x; ++i){
            xx[i] = x[i];
        }
        obj_value = problem.obj_val(xx);

        return true;
    }

    template <typename F>
    bool Simple_Ipopt<F>::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
        std::vector<double> xx(problem.dim_x);
        for (int i = 0; i < problem.dim_x; ++i){
            xx[i] = x[i];
        }
        const std::vector<double> & grad = problem.gradient(xx);
        for (int i = 0; i < problem.dim_x; ++i){
            grad_f[i] = grad[i];
        }
        return true;
    }

    template <typename F>
    bool Simple_Ipopt<F>::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
        if (m > 0){
            throw std::runtime_error("Error!");
        }
        return true;
    }

    template <typename F>
    bool Simple_Ipopt<F>::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow,
                                     Index *jCol, Number *values) {
        if (m > 0){
            throw std::runtime_error("Error!");
        }
        return true;
    }

    template <typename F>
    bool Simple_Ipopt<F>::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m,
                                 const Number *lambda, bool new_lambda, Index nele_hess, Index *iRow,
                                 Index *jCol, Number *values) {
        if (values == NULL){
            // return the structure. This is a symmetric matrix, fill the lower left
            // triangle only.
            int idx = 0;
            for (int i = 0; i < problem.dim_x; ++i){
                for (int j = 0; j <= i; ++j){
                    iRow[idx] = i;
                    jCol[idx] = j;
                    idx++;
                }
            }
            if (idx != nele_hess){
                throw std::runtime_error("Error!");
            }
        } else {
            std::vector<double> xx(problem.dim_x);
            for (int i = 0; i < problem.dim_x; ++i){
                xx[i] = x[i];
            }
            const std::vector<double> & hessian = problem.hessian(xx);

            int idx = 0;
            for (int i = 0; i < problem.dim_x; ++i){
                for (int j = 0; j <= i; ++j){
                    values[idx] = obj_factor * hessian[i * problem.dim_x + j];
                    idx++;
                }
            }
        }
        return true;
    }

    template <typename F>
    void Simple_Ipopt<F>::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L,
                                            const Number *z_U, Index m, const Number *g, const Number *lambda,
                                            Number obj_value, const IpoptData *ip_data,
                                            IpoptCalculatedQuantities *ip_cq) {
        auto & ipopt_status = result.status;
        if (status == Ipopt::SUCCESS){
            ipopt_status = success;
        } else if (status == Ipopt::CPUTIME_EXCEEDED){
            ipopt_status = max_cpu_time_reached;
        } else if (status == Ipopt::MAXITER_EXCEEDED){
            ipopt_status = max_iter_reached;
        } else if (status == Ipopt::STOP_AT_ACCEPTABLE_POINT){
            ipopt_status = acceptable;
        } else {
            ipopt_status = other_failure;
        }

        result.opt_obj = obj_value;
        result.opt_x = std::vector<double>(problem.dim_x);
        for (int i = 0; i < problem.dim_x; ++i){
            result.opt_x[i] = x[i];
        }
        result.opt_grad = problem.gradient(result.opt_x);

        result.multiplier_x_lb = std::vector<double>(problem.dim_x);
        result.multiplier_x_ub = std::vector<double>(problem.dim_x);

        for (int i = 0; i < problem.dim_x; ++i){
            result.multiplier_x_lb[i] = z_L[i];
        }
        for (int i = 0; i < problem.dim_x; ++i){
            result.multiplier_x_ub[i] = z_U[i];
        }
    }

    template <typename F>
    Simple_Ipopt<F>::~Simple_Ipopt<F>() {

    }

    template<typename F>
    Result Ipopt_Solver::optimize(Simple_Opt_Problem<F> &opt, const std::vector<double> &initial_value) {
        if (initial_value.size() != opt.dim_x && initial_value.size() != 0){
            throw std::runtime_error("wrong dimension in initial value!");
        }

        if (opt.use_hassian){
            app->Options()->SetStringValue("hessian_approximation", "exact");
        } else {
            app->Options()->SetStringValue("hessian_approximation", "limited-memory");
            app->Options()->SetIntegerValue("limited_memory_max_history", 100);
        }

        // Prepare Result
        Result result;
        result.status = other_failure;

        // Create an instance of your nlp...
        Ipopt::SmartPtr<Ipopt::TNLP> simple_ipopt = new Simple_Ipopt<F>(opt, result, initial_value);

        // Ask Ipopt to solve the problem
        auto status = app->OptimizeTNLP(simple_ipopt);

        return result;
    }

}



#endif //EXAMPLE_IPOPT_WRAPPER_H
