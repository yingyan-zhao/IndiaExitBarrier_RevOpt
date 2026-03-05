/**
 * This header implements the 'SOO' algorithm in
 *     Optimistic Optimization of a Deterministic Function without the Knowledge of its Smoothness
 *     by Remi Munos
 */


#ifndef DFS_SOO_H
#define DFS_SOO_H

#include "DFS.h"
#include <vector>
#include <array>
#include <tuple>
#include <set>
#include <map>
#include <cinttypes>

namespace DFS{
struct SOO_Option {
    int display_level = 1;  // 0 --- no display
                            // 1 --- display some progress
                            // 2 --- display details

    // 'track_all_local_optimizers' is only meaningful if you
    // want the algorithm to return local optimizers.
    //
    // if true, the algorithm will keep a record of all local optimizer,
    // but it can be very costly, when the dim_x is large.
    // if false, the algorithm will try to find some local optimizers,
    // but it might miss some local optimizer
    bool track_all_local_optimizers = false;

    // the program will stop and return as converged, if
    // there is no new local optimizer after 'n_iterations_with_same_local_optimizers' iterations.
    //
    // if you think the function has lots of local optimizer,
    // you should use a larger value for n_iterations_with_same_local_optimizers.
    unsigned n_iterations_with_same_local_optimizers = 20;
    double x_tol = 1e-13;

    unsigned max_iter = 500;
    unsigned max_cubes = std::numeric_limits<unsigned>::max();


    double* f_target = nullptr;   // stop when optimal f value is below *f_target (if minimizing)
                                  // stop when optimal f value is above *f_target (if maximizing)
};


/**
 * The following function implements the 'SOO' algorithm in
 *     Optimistic Optimization of a Deterministic Function without the Knowledge of its Smoothness
 *            by Remi Munos, (2011), NIPS'11: Proceedings of the 24th International Conference on Neural Information Processing System
 * @tparam F
 * @tparam G
 * @param problem           the optimization problem
 * @param option            option
 * @param local_optimizer   if provided and if option.track_all_local_optimizers == true, then return all the local optimizer
 *                          if provided and local_optimizer->size() > 0, then will try to find 'n' local optimizer
 *                             with n = min(local_optimizer->size(), number of local optimizer found).
 * @return  optimization result
 */
    template<typename F, typename G>
    Result SOO(Optimization_Problem<F, G> & problem, const SOO_Option & option = SOO_Option(), std::vector<std::vector<double>>* local_optimizer = nullptr);

    /********************************************************************************************************/
    /*                                      Implementation                                                  */
    /********************************************************************************************************/
    class HyperCube {
    public:
        const uint_fast16_t dim;
        HyperCube(const uint_fast16_t & dim, std::vector<double>& memory);
        HyperCube(const std::vector<double>& left_vec, const std::vector<double>& right_vec, std::vector<double>& memory);
        void          set_left(const std::vector<double>& left_update);
        void          set_right(const std::vector<double>& right_update);
        std::vector<double> center() const;
        std::vector<double> get_left() const;
        std::vector<double> get_right() const;

        void split(HyperCube& C_left, HyperCube& C_center, HyperCube& C_right, const std::vector<double> & split_weight) const;

        HyperCube(const HyperCube& cube);
        HyperCube& operator=(const HyperCube& cube);
        bool is_neighbour(const HyperCube & cube) const;
        bool lexicographical_order(const HyperCube & cube) const;

        bool contains(const std::vector<double> & point) const;
        double norm() const;

    public:
        std::vector<double>& memory;
        unsigned left;  // the position in 'memory' where cooridate of left-lower corner is stored, i.e.
        // memory[left] to memory[left + dim - 1]
        unsigned right; // the position in 'memory' where cooridate of right-upper corner is stored, i.e.
        // memory[right] to memory[right + dim -1]

        HyperCube(unsigned left_idx, const std::vector<double>& right_vec, std::vector<double>& memory);
        HyperCube(const std::vector<double>& left_vec, unsigned right_idx, std::vector<double>& memory);
    };

    template<typename F, typename G>
    Result SOO(Optimization_Problem<F, G> & problem, const SOO_Option & option, std::vector<std::vector<double>>* local_optimizer_ptr){
        // check the validity of the inputs
        unsigned dim_x = problem.dim_effective_x;
        check(dim_x < std::numeric_limits<uint_fast16_t>::max(), "dim_x should not be larger than 65535");
        if (problem.has_nonlinear_constrains){
            throw std::runtime_error("this algorithm can only deal with box constraint.");
        }
        for (unsigned i = 0; i < dim_x; ++i){
            if (problem.effective_x_lb[i] == std::numeric_limits<double>::lowest() || problem.effective_x_ub[i] == std::numeric_limits<double>::max()){
                throw std::runtime_error("SOO only works well when the feasible region is bounded.");
            }
        }
        // deal with trivial case
        if (dim_x == 0){
            Result result;
            result.optimal_x = problem.from_effective_x_to_x(std::vector<double>());
            result.optimal_f = problem.obj_for_minimization(std::vector<double>());
            if (!problem.minimize){
                result.optimal_f *= -1;
            }
            result.flag = converged;
            return result;
        }

        // prepare for the optimal x and val
        double opt_f = std::numeric_limits<double>::lowest();
        std::vector<double> opt_x(dim_x);

        // split weight
        std::vector<double> split_weight(dim_x);
        for (unsigned i = 0; i < dim_x; ++i){
            split_weight[i] = 1.0 / (problem.effective_x_ub[i] - problem.effective_x_lb[i]);
        }

        // cubes management
        std::vector<double> memory;
        std::vector<HyperCube> cubes;
        std::map<unsigned, std::set<uint64_t>> depth;
        std::vector<double> cube_value;
        std::vector<unsigned> cube_depth;
        std::vector<unsigned> cubes_to_be_deleted;

        std::vector<unsigned> cube_better_neighbor;
        std::vector<std::set<unsigned>> cube_neighbor;
        std::set<unsigned> local_optimizer;


        auto declare_new_cube = [&](){
            unsigned id = 0;
            if (cubes_to_be_deleted.size() > 0){
                id = cubes_to_be_deleted.back();
                cubes_to_be_deleted.pop_back();
            } else {
                id = cubes.size();
                cubes.emplace_back(dim_x, memory);
                cube_value.push_back(std::numeric_limits<double>::lowest());
                cube_depth.push_back(0);
                if (option.track_all_local_optimizers){
                    cube_neighbor.push_back({});
                    cube_better_neighbor.push_back(0);
                }
            }
            return id;
        };

        auto remove_cube = [&](unsigned id){
            // take care of depth
            auto it = depth.find(cube_depth[id]);
            it->second.erase(id);
            if (it->second.size() == 0){
                depth.erase(it);
            }

            // put it onto the to-be-deleted list
            cubes_to_be_deleted.push_back(id);
        };

        auto split_cube = [&](unsigned cube_id){
            unsigned idx_left = declare_new_cube();
            unsigned idx_center = declare_new_cube();
            unsigned idx_right = declare_new_cube();
            cubes[cube_id].split(cubes[idx_left], cubes[idx_center], cubes[idx_right], split_weight);

            // take care of depth
            unsigned new_depth = cube_depth[cube_id] + 1;
            auto it = depth.find(new_depth);
            if (it == depth.end()){
                depth.insert(depth.end(), {new_depth, {idx_left, idx_center, idx_right}});
            } else{
                it->second.insert(idx_left);
                it->second.insert(it->second.end(), idx_center);
                it->second.insert(it->second.end(), idx_right);
            }
            cube_depth[idx_left] = new_depth;
            cube_depth[idx_center] = new_depth;
            cube_depth[idx_right] = new_depth;

            // take care of value
            double val_left = problem.obj_for_maximization(cubes[idx_left].center());
            double val_center = cube_value[cube_id];
            double val_right = problem.obj_for_maximization(cubes[idx_right].center());
            if (val_left > opt_f){
                opt_f = val_left;
                opt_x = cubes[idx_left].center();
            }
            if (val_right > opt_f){
                opt_f = val_right;
                opt_x = cubes[idx_right].center();
            }
            cube_value[idx_left] = val_left;
            cube_value[idx_center] = val_center;
            cube_value[idx_right] = val_right;

            if (option.track_all_local_optimizers) {
                // take care of the neighbors
                cube_neighbor[idx_left].clear();
                cube_neighbor[idx_center].clear();
                cube_neighbor[idx_right].clear();
                cube_neighbor[idx_left].insert(idx_center);
                cube_neighbor[idx_right].insert(idx_center);
                cube_neighbor[idx_center].insert(idx_left);
                cube_neighbor[idx_center].insert(idx_right);

                for (const auto &id: cube_neighbor[cube_id]) {
                    cube_neighbor[id].erase(cube_id);
                    for (auto myid : {idx_left, idx_center, idx_right}){
                        if (cubes[myid].is_neighbour(cubes[id])) {
                            cube_neighbor[myid].insert(id);
                            cube_neighbor[id].insert(cube_neighbor[id].end(), myid);
                        }
                    }
                }

                // take care of the local optimizer
                if (cube_better_neighbor[cube_id] == 0){
                    local_optimizer.erase(cube_id);
                }
                for (const auto &id: cube_neighbor[cube_id]) {
                    if (cube_value[id] <= val_center){
                        --cube_better_neighbor[id];
                        if (cube_better_neighbor[id] == 0){
                            local_optimizer.insert(id);
                        }
                    }
                }

                cube_better_neighbor[idx_left] = val_left <= val_center;
                cube_better_neighbor[idx_center] = (unsigned)(val_center <= val_left) + (unsigned)(val_center <= val_right);
                cube_better_neighbor[idx_right] = val_right <= val_center;
                for (auto myid: {idx_left, idx_center, idx_right}) {
                    const double &v = cube_value[myid];
                    for (const auto &id: cube_neighbor[myid]) {
                        if (id != idx_left && id != idx_center && id != idx_right){
                            if (v >= cube_value[id]) {
                                if (cube_better_neighbor[id]++ == 0){
                                    local_optimizer.erase(id);
                                }
                            }
                            if (v <= cube_value[id]) {
                                ++cube_better_neighbor[myid];
                            }
                        }
                    }
                    if (cube_better_neighbor[myid] == 0){
                        local_optimizer.insert(myid);
                    }
                }
            }

            // take care of the old cube
            remove_cube(cube_id);
        };

        // initialize cubes
        memory.reserve(dim_x * 10000);
        cubes.emplace_back(problem.effective_x_lb, problem.effective_x_ub, memory);
        cube_value.push_back(problem.obj_for_maximization(cubes[0].center()));
        cube_depth.push_back(0);
        depth.insert({0, {0}});
        opt_f = cube_value[0];
        opt_x = cubes[0].center();
        if (option.track_all_local_optimizers){
            cube_better_neighbor.push_back(0);
            cube_neighbor.push_back({});
            local_optimizer.insert(0);
        }

        // start loop
        unsigned n_iter_with_same_local_optimizer = 0;
        Flag flag = max_iter_reached;
        const double f_target = option.f_target != nullptr ? (problem.minimize ? -(*option.f_target) : *option.f_target) : std::numeric_limits<double>::max();

        for (unsigned iter = 0; iter < option.max_iter; ++iter){
            unsigned n_cubes_old = cubes.size() - cubes_to_be_deleted.size();
            unsigned n_local_optimizer_old = local_optimizer.size();

            if (option.display_level >= 1){
                double f = problem.minimize ? -opt_f : opt_f;
                std::cout << "iter " << iter << std::endl;
                std::cout << "opt_f: " << f << " n_cubes: " << n_cubes_old << " n_local_optimizer_old: " << local_optimizer.size() << std::endl;
            }

            if (opt_f > f_target){
                flag = target_reached;
                break;
            }

            if (cubes.size() - cubes_to_be_deleted.size() > option.max_cubes){
                flag = other;
                break;
            }

            double highest_value_so_far = std::numeric_limits<double>::lowest();
            unsigned h = 0;

            while (h <= depth.rbegin()->first){
                auto it = depth.find(h);
                if (it != depth.end()){
                    // find the best value in depth 'h'
                    double high_v = std::numeric_limits<double>::lowest();
                    unsigned best_cube;
                    for (const auto & id : it->second){
                        if (cube_value[id] > high_v){
                            high_v = cube_value[id];
                            best_cube = id;
                        }
                    }

                    // compare it with the highest_value so far
                    if (high_v > highest_value_so_far){
                        // update the highest value so far
                        highest_value_so_far = high_v;
                        if (cubes[best_cube].norm() > option.x_tol){
                            // split the best cube
                            split_cube(best_cube);
                        }

                        if (opt_f > f_target){
                            flag = target_reached;
                            break;
                        }
                    }
                }
                ++h;
            }

            if (opt_f > f_target){
                flag = target_reached;
                break;
            }

            if (option.track_all_local_optimizers){
                if (n_local_optimizer_old == local_optimizer.size()){
                    if (n_iter_with_same_local_optimizer++ >= option.n_iterations_with_same_local_optimizers){
                        flag = converged;
                        break;
                    }
                } else {
                    n_iter_with_same_local_optimizer = 0;
                }
            }

            if (cubes.size() - cubes_to_be_deleted.size() == n_cubes_old){
                flag = converged;
                break;
            }
        }

        Result result;
        result.optimal_f = problem.minimize ? -opt_f : opt_f;
        result.optimal_x = opt_x;
        result.flag = flag;

        if (flag != target_reached){
            if (option.track_all_local_optimizers && local_optimizer_ptr != nullptr){
                local_optimizer_ptr->clear();
                local_optimizer_ptr->reserve(local_optimizer.size());
                for (const auto & i : local_optimizer){
                    local_optimizer_ptr->emplace_back(cubes[i].center());
                }
            } else if (local_optimizer_ptr != nullptr && local_optimizer_ptr->size() > 0){
                std::vector<unsigned> valid_id;
                valid_id.reserve(cubes.size() - cubes_to_be_deleted.size());
                bool empty_cubes_to_be_deleted = cubes_to_be_deleted.empty();
                for (unsigned i = 0; i < cubes.size(); ++i){
                    if (empty_cubes_to_be_deleted || cubes_to_be_deleted.front() != i){
                        valid_id.push_back(i);
                    }
                }
                std::sort(valid_id.begin(), valid_id.end(), [&](unsigned i, unsigned j){return cube_value[i] > cube_value[j];});

                auto distance = [&](const std::vector<double>& a, const std::vector<double> & b){
                    double d = 0;
                    for (unsigned i = 0; i < dim_x; ++i){
                        d = std::max(d, std::fabs(a[i] - b[i]) * split_weight[i]);
                    }
                    return d;
                };

                unsigned pos = 0;
                unsigned n_local_opt_found = 0;
                while (pos < valid_id.size() && n_local_opt_found < local_optimizer_ptr->size()){
                    bool local_opt = true;
                    for (unsigned i = 0; i < pos; ++i){
                        bool is_neighbor = cubes[valid_id[pos]].is_neighbour(cubes[valid_id[i]]);
                        if (is_neighbor){
                            local_opt = false;
                            break;
                        }
                    }
                    if (local_opt){
                        auto x = cubes[valid_id[pos]].center();
                        bool far = true;
                        for (unsigned j = 0; j < n_local_opt_found; ++j){
                            far = distance((*local_optimizer_ptr)[j], x) > 0.01;
                            if (!far){
                                break;
                            }
                        }
                        if (far){
                            (*local_optimizer_ptr)[n_local_opt_found].clear();
                            (*local_optimizer_ptr)[n_local_opt_found] = problem.from_effective_x_to_x(x);
                            n_local_opt_found++;
                        }
                    }
                    pos++;
                }
                local_optimizer_ptr->resize(n_local_opt_found);
            }
        } else {
            if (local_optimizer_ptr != nullptr && local_optimizer_ptr->size() > 0){
                local_optimizer_ptr->clear();
            }
        }

        return result;
    }


    inline HyperCube::HyperCube(const std::vector<double>& left_vec,
                                const std::vector<double>& right_vec,
                                std::vector<double>&       memory)
            : dim(left_vec.size())
            , memory(memory)
            , left(memory.size())
            , right(memory.size() + dim)
    {
        check(left_vec.size() == right_vec.size(), "sizes shoud be consistent");
        memory.insert(memory.end(), left_vec.begin(), left_vec.end());
        memory.insert(memory.end(), right_vec.begin(), right_vec.end());
    }

    inline HyperCube::HyperCube(unsigned left_idx, const std::vector<double>& right_vec, std::vector<double>& memory)
            : dim(right_vec.size())
            , memory(memory)
            , left(left_idx)
            , right(memory.size())
    {
        check(memory.size() >= left + dim, "memory should have enough room");
        memory.insert(memory.end(), right_vec.begin(), right_vec.end());
    }

    inline HyperCube::HyperCube(const std::vector<double>& left_vec, unsigned right_idx, std::vector<double>& memory)
            : dim(left_vec.size())
            , memory(memory)
            , left(memory.size())
            , right(right_idx)
    {
        check(memory.size() >= right + dim, "memory should have enough room");
        memory.insert(memory.end(), left_vec.begin(), left_vec.end());
    }

    inline HyperCube::HyperCube(const HyperCube& cube)
            : dim(cube.dim)
            , memory(cube.memory)
            , left(cube.left)
            , right(cube.right)
    {
    }

    inline
    HyperCube::HyperCube(const uint_fast16_t& dim, std::vector<double>& memory)
            : dim(dim)
            , memory(memory)
            , left(0)
            , right(0)
    {
    }

    inline HyperCube& HyperCube::operator=(const HyperCube& cube)
    {
        check(dim == cube.dim, "dimension should be consistent.");
        check(&memory == &cube.memory, "these cubes should point to the same memory");
        left  = cube.left;
        right = cube.right;
        return *this;
    }

    inline
    void HyperCube::set_left(const std::vector<double>& left_update)
    {
        copy(left_update.begin(), left_update.end(), memory.begin() + left);
    }

    inline
    void HyperCube::set_right(const std::vector<double>& right_update)
    {
        copy(right_update.begin(), right_update.end(), memory.begin() + right);
    }

    inline
    std::vector<double> HyperCube::get_left() const
    {
        std::vector<double> left_vec(memory.begin() + left, memory.begin() + left + dim);
        return left_vec;
    }

    inline
    std::vector<double> HyperCube::get_right() const
    {
        std::vector<double> right_vec(memory.begin() + right, memory.begin() + right + dim);
        return right_vec;
    }

    inline
    std::vector<double> HyperCube::center() const
    {
        check(dim != 0, "the cube must be initialized.");
        std::vector<double> c(memory.begin() + left, memory.begin() + left + dim);
        for (unsigned i = 0; i < dim; ++i)
        {
            c[i] = (c[i] + memory[right + i]) / 2;
        }
        return c;
    }

    inline
    void HyperCube::split(HyperCube&    C_left,
                          HyperCube&    C_center,
                          HyperCube&    C_right,
                          const std::vector<double> & split_weight) const
    {
        check(split_weight.size() == dim, "dimension of split_weight should be consistent");

        // split along the longest dimension
        double max_length = std::numeric_limits<double>::lowest();
        unsigned  sp_dim;

        for (unsigned i = 0; i < dim; ++i) {
            double temp = (memory[right + i] - memory[left + i]) * split_weight[i];
            if (temp > max_length)
            {
                max_length = temp;
                sp_dim     = i;
            }
        }

        std::vector<double> C_left_right(dim);
        std::vector<double> C_center_left(dim);
        std::vector<double> C_center_right(dim);
        std::vector<double> C_right_left(dim);

        for (unsigned i = 0; i < dim; ++i)
        {
            if (i != sp_dim)
            {
                C_left_right[i] = memory[right + i];

                C_center_left[i]  = memory[left + i];
                C_center_right[i] = memory[right + i];

                C_right_left[i] = memory[left + i];
            }
            else
            {
                C_left_right[i] = memory[left + i] + (memory[right + i] - memory[left + i]) / 3;

                C_center_left[i]  = memory[left + i] + (memory[right + i] - memory[left + i]) / 3;
                C_center_right[i] = memory[right + i] - (memory[right + i] - memory[left + i]) / 3;

                C_right_left[i] = memory[right + i] - (memory[right + i] - memory[left + i]) / 3;
            }
        }

        C_left   = HyperCube(left, C_left_right, memory);
        C_center = HyperCube(C_center_left, C_center_right, memory);
        C_right  = HyperCube(C_right_left, right, memory);
    }

    inline
    bool HyperCube::is_neighbour(const HyperCube &cube) const {
        auto this_left = [&](const unsigned & i){
            return memory[left + i];
        };
        auto this_right = [&](const unsigned & i){
            return memory[right + i];
        };
        auto other_left = [&](const unsigned & i){
            return memory[cube.left + i];
        };
        auto other_right = [&](const unsigned & i){
            return memory[cube.right + i];
        };

        bool this_small = false;
        bool other_small = false;
        for (unsigned i = 0; i < dim; ++i){
            if (this_left(i) <= other_left(i) && other_right(i) <= this_right(i)){
                other_small = true;
                break;
            } else if (this_left(i) >= other_left(i) && other_right(i) >= this_right(i)){
                this_small = true;
                break;
            }
        }

        bool neighbor = this_small || other_small;
        if (this_small){
            for (unsigned i = 0; i < dim; ++i){
                neighbor = (other_left(i) <= this_left(i)) && (this_left(i) <= other_right(i));
                neighbor = neighbor || (other_left(i) <= this_right(i)) && (this_right(i) <= other_right(i));
                if (!neighbor){
                    break;
                }
            }
        } else if (other_small){
            for (unsigned i = 0; i < dim; ++i){
                neighbor = (this_left(i) <= other_left(i)) && (other_left(i) <= this_right(i));
                neighbor = neighbor || (this_left(i) <= other_right(i)) && (other_right(i) <= this_right(i));
                if (!neighbor){
                    break;
                }
            }
        }
        return neighbor;

//        bool neighbour_1 = false;
//        unsigned neighbour_2 = 0;
//        for (unsigned i = 0; i < dim; ++i) {
//            if (memory[left + i] != memory[right + i]) {
//                if (cube.memory[cube.left + i] == memory[right + i]  or cube.memory[cube.right + i] == memory[left + i]) {
//                    neighbour_1 = true;
//                }
//                else
//                {
//                    if (cube.memory[cube.left + i] <= memory[left + i] and cube.memory[cube.right + i] >= memory[right + i])
//                    {
//                        ++neighbour_2;
//                    }
//                    else
//                    {
//                        if (cube.memory[cube.left + i] >= memory[left + i] and cube.memory[cube.right + i] <= memory[right + i])
//                        {
//                            ++neighbour_2;
//                        }
//                    }
//                }
//            }
//            else
//            {
//                ++neighbour_2;
//            }
//        }
//
//        return neighbour_1 and (neighbour_2 >= dim - 1);
    }


    inline
    bool HyperCube::lexicographical_order(const HyperCube &cube) const
    {
        bool result = true;
        for (unsigned i = 0; i < dim; ++i)
        {
            if (memory[left + i] < cube.memory[cube.left + i])
            {
                result = true;
                break;
            }
            if (memory[left + i] > cube.memory[cube.left + i])
            {
                result = false;
                break;
            }
        }
        return result;
    }

    inline
    bool HyperCube::contains(const std::vector<double> & point) const
    {
        if (point.size() != dim) throw std::runtime_error("inconsistent dim");
        bool contain = true;
        for (unsigned i = 0; i < dim; ++i)
        {
            contain = (point[i] >= memory[left + i]) and (point[i] <= memory[right + i]);
            if (not contain) break;
        }
        return contain;
    }

    inline
    double HyperCube::norm() const
    {
        double v = 0;
        for (unsigned i = 0; i < dim; ++i)
        {
            v = std::max(v, memory[right + i] - memory[left + i]);
        }
        return v;
    }





}



#endif //DFS_SOO_H
