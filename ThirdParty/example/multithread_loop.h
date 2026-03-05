#ifndef MULTITHREAD_EXAMPLES_MULTITHREAD_LOOP_H
#define MULTITHREAD_EXAMPLES_MULTITHREAD_LOOP_H

#include <thread>
#include <atomic>
#include <vector>
#include <mutex>
#include <memory>
#include <condition_variable>
#include <stdexcept>
#include <functional>
#include <queue>
#include <future>
#include <utility>
#include <chrono>

namespace MultiThreads{

    /**
     * This class handles notification that a job is done.
     */
    class Notification{
        bool done = false;
        std::mutex mutex;
        std::condition_variable condition;

    public:
        /**
         * This function is used on the receiver end. It blocks the receiver thread
         */
        inline void wait_until_done(){
            std::unique_lock<std::mutex> lock(mutex);
            // block current thread
            condition.wait(lock, [this]{return done;});
        }

        /**
         * This function is used on the sender send. It send a notification to the receiver that the job is done
         */
        inline void is_done(){
            std::unique_lock<std::mutex> lock(mutex);
            done = true;
            condition.notify_one();
        }
    };

    /**
     * This class manages a thread pool
     */
    class Threads_Management{
    public:
        /**
         * Constructor
         * @param n_threads  number of threads that are available to this object
         */
        Threads_Management(int64_t n_threads = std::thread::hardware_concurrency() - 1);

        /**
         * return number of threads availabe to this class in total
         * @return
         */
        uint16_t n_threads_in_total() const;

        /**
         * Submit a task and expect it to run immediately. If the task cannot start immediately, it will NOT be
         * added to the task queue and this function will return nullptr.
         * @tparam F        template type of functional
         * @tparam Args     template type of arguments
         * @param f         function enclosure/ functor/
         * @param args      arguments
         * @return          pointer to Notification. = nullptr if this task has not been added to the task queue.
         */
        template<class F, class... Args>
        std::shared_ptr<Notification> submit_task_that_starts_immediately(F && f, Args && ... args);

        template<class F, class... Args>
        std::shared_ptr<Notification> submit_task(F && f, Args && ... args);

        uint16_t n_idle_threads() const;

        ~Threads_Management();

        std::mutex print_mutex;
    private:
        // stop flag, if true, all worker will stop waiting for new jobs
        std::atomic_bool stop;

        // maximum threads we will use
        const uint16_t n_threads;

        // queue of tasks
        std::queue<std::function<void()>> tasks = {};

        // for synchonization
        std::mutex queue_mutex;

        // this condition variable is used to tell thread to wait or start a task
        std::condition_variable condition;

        // thread pool
        std::vector<std::thread> threads;

        // holding the busy/idle status of each thread
        std::vector<std::atomic_uint8_t> status_of_threads;


        bool at_least_one_thread_is_idle();

        bool all_threads_are_idle();

        void mark_busy(uint16_t id);

        void mark_idle(uint16_t id);

        bool is_idle(uint16_t id) const;

        void create_threads();
    };

    /**
     * as long as the order of processing f(i) does not matter, the following function is equivalent of running
     *
     *    for (uint i = 0; i < n; ++i){
     *      f(i);
     *    }
     *
     * @tparam F                    type of function
     * @param f                     function, must implement  f(i) or f(i, thread_id), where 'i' is of the same type as 'n', and thread_id is of type unsigned
     * @param n                     number of loops
     * @param threadsManagement     thread pool
     * @param maximum_threads       maximum_threads to be used
     */
    template<typename F>
    void simple_parallel_for(F & f,
                             std::int64_t n,
                             MultiThreads::Threads_Management & threadsManagement,
                             unsigned maximum_threads = std::thread::hardware_concurrency());

     /**
     * As long as the sequence of calling f(i,j) does not matter, the following code is equivalent to
     *
     *  std::vector<int8_t> result(n_outerloop, 1);
     *  for (unsigned i = 0; i < n_outerloop; ++i){
     *      for (unsigned j = 0; j < n_innerloop; ++j){
     *          bool result_ij = f(i, j);
     *          if (!result_ij){
     *              result[i] = 0;
     *              break;
     *          }
     *      }
     *  }
     *  return result
     *
     * @tparam F            type of function / functor F
     * @param f             function / functor, must implement f(outer_index, inner_index) or f(outer_index, inner_index, thread_id)
     * @param n_outerloop   number of outer loops
     * @param n_innerloop   number of inner loops
     */
    template<typename F, typename UINT_OUT, typename UINT_IN>
    std::vector<int8_t> nested_parallel_for(F & f, const UINT_OUT n_outerloop, const UINT_IN n_innerloop, Threads_Management & threadsManagement);



    /*******************************************************************************************/
    /*                  Implementation                                                         */
    /*******************************************************************************************/

    template <typename T>
    struct get_arity : get_arity<decltype(&T::operator())> {};
    template <typename R, typename... Args>
    struct get_arity<R(*)(Args...)> : std::integral_constant<unsigned, sizeof...(Args)> {};
    template <typename R, typename C, typename... Args>
    struct get_arity<R(C::*)(Args...)> : std::integral_constant<unsigned, sizeof...(Args)> {};
    template <typename R, typename C, typename... Args>
    struct get_arity<R(C::*)(Args...) const> : std::integral_constant<unsigned, sizeof...(Args)> {};
// All combinations of variadic/non-variadic, cv-qualifiers and ref-qualifiers


    template<typename F, typename UINT>
    void simple_parallel_for_templated_integer(F & f,
                                               const UINT n,
                                               MultiThreads::Threads_Management & threadsManagement,
                                               unsigned maximum_n_threads){
        constexpr size_t n_f_args = get_arity<F>{}; // get number of arguments taken by the f
        constexpr bool f_has_2_args = n_f_args == 2;

        if (n == 0){
            return;
        }
        if (maximum_n_threads == 0){
            throw std::runtime_error("number of threads should be greater or equal to 1");
        }
        if (maximum_n_threads > n){
            maximum_n_threads = n;
        }
        if (maximum_n_threads > std::thread::hardware_concurrency()){
            maximum_n_threads = std::thread::hardware_concurrency();
        }

        if (maximum_n_threads == 1){
            for (UINT i = 0; i < n; ++i){
                if constexpr(f_has_2_args){
                    f(i, 0);
                } else {
                    f(i);
                }
            }
        }
        else {
            std::atomic<UINT> progress(0);
            auto worker = [&](unsigned thread_id){
                while(true){
                    UINT i = progress.fetch_add(1, std::memory_order_relaxed);
                    if (i >= n){
                        break;
                    }
                    if constexpr(f_has_2_args){
                        f(i, thread_id);
                    } else {
                        f(i);
                    }
                }
            };

            unsigned n_threads_used = 0;
            std::vector<std::shared_ptr<Notification>> futures;

            auto manager = [&](unsigned thread_id){
                auto time_mark = std::chrono::steady_clock::now();
                while(true){
                    UINT i = progress.fetch_add(1, std::memory_order_relaxed);
                    if (i >= n){
                        break;
                    }
                    if constexpr(f_has_2_args){
                        f(i, thread_id);
                    } else {
                        f(i);
                    }

                    // get new worker
                    if (n_threads_used < maximum_n_threads){
                        auto now = std::chrono::steady_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_mark).count();
                        if (duration > 500){
                            auto future = threadsManagement.template submit_task_that_starts_immediately(worker, n_threads_used);
                            if (future != nullptr){
                                futures.template emplace_back(std::move(future));
                                ++n_threads_used;
                            }
                            time_mark = std::chrono::steady_clock::now();
                        }
                    }
                }
            };

            for (n_threads_used = 0; n_threads_used + 1 < maximum_n_threads; n_threads_used++){
                auto future = threadsManagement.template submit_task_that_starts_immediately(worker, n_threads_used);
                if (future != nullptr){
                    futures.template emplace_back(std::move(future));
                    ++n_threads_used;
                }
            }

            manager(n_threads_used++);
            for (auto & future : futures){
                future->wait_until_done();
            }
        }
    }

    template<typename F>
    void simple_parallel_for(F & f, std::int64_t n, MultiThreads::Threads_Management & threadsManagement, unsigned maximum_n_threads){
        if (n < 0){
            throw std::runtime_error("Error in simple_parallel_for. n is negative");
        }
        if (n < std::numeric_limits<std::uint8_t>::max()){
            simple_parallel_for_templated_integer(f, static_cast<std::uint8_t>(n), threadsManagement, maximum_n_threads);
        } else if (n < std::numeric_limits<std::uint16_t>::max()){
            simple_parallel_for_templated_integer(f, static_cast<std::uint16_t>(n), threadsManagement, maximum_n_threads);
        } else if (n < std::numeric_limits<std::uint32_t>::max()){
            simple_parallel_for_templated_integer(f, static_cast<std::uint32_t>(n), threadsManagement, maximum_n_threads);
        } else {
            simple_parallel_for_templated_integer(f, static_cast<std::uint64_t>(n), threadsManagement, maximum_n_threads);
        }
    }

    template<class F, class... Args>
    std::shared_ptr<Notification> Threads_Management::submit_task(F &&f, Args &&...args) {
        std::shared_ptr<Notification> future_pointer = std::make_shared<Notification>();

        std::lock_guard<std::mutex> lock(queue_mutex);

        // push task into queue
        if constexpr (sizeof...(Args) == 0){
            // if there is no argument
            tasks.template emplace([future_pointer, f = std::move(f)]{
                f();
                future_pointer->is_done();
            });
        } else {
            // if there exist at least one argument
#if (__cpp_init_captures >= 201803)
        tasks.template emplace([future_pointer, f = std::move(f), ...args = std::move(args)]{
            f(args...);
            future_pointer->is_done();
        });
#else
            auto f_function = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
            tasks.template emplace([future_pointer, f_function = std::move(f_function)]{
                f_function();
                future_pointer->is_done();
            });
#endif
        }
        condition.notify_one();
        return future_pointer;
    }

    template<class F, class... Args>
    std::shared_ptr<Notification> Threads_Management::submit_task_that_starts_immediately(F &&f, Args &&...args) {
        std::shared_ptr<Notification> future_pointer(nullptr);
        if (n_threads > 0 && at_least_one_thread_is_idle()) {
            // critical region
            std::lock_guard<std::mutex> lock(queue_mutex);
            if (n_idle_threads() > tasks.size()){
                future_pointer = std::make_shared<Notification>();
#if (__cpp_init_captures >= 201803)
                // push task into queue
                tasks.template emplace([future_pointer, f = std::move(f), ...args = std::move(args)]{
                    f(args...);
                    future_pointer->is_done();
                });
#else
                // push task into queue
                auto f_function = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
                // push task into queue
                tasks.template emplace([future_pointer, f_function = std::move(f_function)]{
                    f_function();
                    future_pointer->is_done();
                });
#endif
            }

            if (future_pointer != nullptr){
                // notify an idle thread, the existence of submit_mutex ensures that
                // no one is holding the queue_mutex
                condition.notify_one();
            }
        }
        return future_pointer;
    }


    inline
    Threads_Management::~Threads_Management() {
        // wake up all threads and tell them to stop
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            stop = true;
            condition.notify_all();
        }
        // wait for them to join the main thread
        for (auto & thread : threads){
            thread.join();
        }
    }

    inline
    uint16_t Threads_Management::n_threads_in_total() const {
        return n_threads;
    }

    inline
    uint16_t Threads_Management::n_idle_threads() const {
        uint16_t idle_threads = 0;
        for (uint16_t i = 0; i < n_threads; ++i){
            idle_threads += is_idle(i);
        }
        return idle_threads;
    }

    inline
    bool MultiThreads::Threads_Management::all_threads_are_idle() {
        bool all_idle = true;
        for (uint16_t i = 0; i < n_threads; ++i){
            all_idle = is_idle(i);
            if (!all_idle){
                break;
            }
        }
        return all_idle;
    }

    inline bool MultiThreads::Threads_Management::at_least_one_thread_is_idle() {
        bool has_idle_thread = false;
        for (uint16_t i = 0; i < n_threads; ++i){
            has_idle_thread = is_idle(i);
            if (has_idle_thread){
                break;
            }
        }
        return has_idle_thread;
    }

    inline void MultiThreads::Threads_Management::mark_busy(uint16_t id) {
        uint i = id / 8;
        uint j = id % 8;
        uint8_t mask = static_cast<uint8_t>(1) << j;
        status_of_threads[i] |= mask;
    }

    inline void MultiThreads::Threads_Management::mark_idle(uint16_t id) {
        uint i = id / 8;
        uint j = id % 8;
        uint8_t mask = 7 - (static_cast<uint8_t>(1) << j);
        status_of_threads[i] &= mask;
    }

    inline
    bool MultiThreads::Threads_Management::is_idle(uint16_t id) const {
        uint i = id / 8;
        uint j = id % 8;
        uint8_t mask = static_cast<uint8_t>(1) << j;
        return (status_of_threads[i].load() & mask) == 0;
    }

    inline MultiThreads::Threads_Management::Threads_Management(int64_t n_threads)
    : stop(false)
    , n_threads(n_threads)
    , status_of_threads((n_threads + 7) / 8)
    {
        if (n_threads > std::numeric_limits<uint16_t>::max()){
            throw std::runtime_error("n_threads is too large");
        }
        if (n_threads < 0){
            throw std::runtime_error("n_threads must be nonnegative");
        }

        for (uint i = 0; i < status_of_threads.size(); ++i){
            status_of_threads[i] = 0;
        }
        create_threads();
    }

    inline void MultiThreads::Threads_Management::create_threads() {
        threads.reserve(n_threads);

        auto worker = [&](uint16_t id) {
//            constexpr auto sleep_duration = std::chrono::seconds(1);
            while(!stop){
                bool task_assigned = false;
                std::function<void()> task;
                // critical section
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    // block current thread
                    condition.wait(lock, [this](){return stop || !tasks.empty();});
                    if (!tasks.empty()){
                        // plan to execute the first element of queue
                        task = std::move(tasks.front());
                        tasks.pop();
                        task_assigned = true;
                    }
                }
                if (task_assigned){
                    mark_busy(id);  // mark this thread as busy
                    task();         // run task
                    mark_idle(id);  // mark this thread as idle
                }
            }
        };

        for (uint16_t i = 0; i < n_threads; ++i){
            threads.emplace_back(worker, i);
        }
    }


    template<typename F, typename UINT_OUT, typename UINT_IN>
    std::vector<int8_t> nested_parallel_for(F& f, const UINT_OUT n_outerloop, const UINT_IN n_innerloop, Threads_Management & threadsManagement){
        if (n_outerloop == 0){
            return {};
        }

        constexpr size_t n_f_args = get_arity<F>{}; // get number of arguments taken by the f
        constexpr bool f_has_3_args = n_f_args == 3; // true if f takes 3 arguments

        const unsigned n_threads = std::thread::hardware_concurrency();
        // global outer index progress
        std::atomic<UINT_OUT> outer_loop_progress(0);
        // outer index that each thread is working on
        std::vector<UINT_OUT> outer_indices(n_threads, n_outerloop);

        // last stage flag(s)
        std::atomic_bool some_thread_entered_last_stage;
        some_thread_entered_last_stage = false;

        std::vector<std::atomic_bool> last_stage_flag(n_threads);
        for (unsigned i = 0; i < n_threads; ++i){
            last_stage_flag[i].store(false);
        }

        // in each thread, the range of inner loop is [0, some_thread_entered_last_stage ? n_innerloop : last_stage_inner_index_ub[thread_id])
        // the last_stage_inner_index_ub will not be changed until the last stage
        std::vector<std::atomic<UINT_IN>> last_stage_inner_index_ub(n_threads);
        for (auto & idx : last_stage_inner_index_ub){
            idx.store(n_innerloop);
        }


        // the following collects results in the last round
        std::vector<UINT_OUT> negative_results_last_round;
        std::mutex last_round_mutex; // mutex for reading / writing the above object

        std::vector<int8_t> result(n_outerloop, 1);
        auto worker = [&](unsigned thread_id){
            while(true){
                // assign outer index
                unsigned outer_idx = outer_loop_progress.fetch_add(1, std::memory_order_relaxed);

                if (outer_idx >= n_outerloop){
                    // if all outer index has been assigned, enter the last round
                    if (!last_stage_flag[thread_id]){
                        last_stage_flag[thread_id].store(true);
                    }
                    if (!some_thread_entered_last_stage.load()){
                        some_thread_entered_last_stage.store(true);
                    }

                    unsigned other_thread_id = thread_id;
                    unsigned inner_idx = n_innerloop;
                    // look for what other threads are working on and share the workload
                    for (unsigned i = 0; i < n_threads; ++i){
                        // for thread i to have workload to share, it cannot be my thread
                        bool has_workload_to_share = i != thread_id;
                        // for thread i to have workload to share, it cannot be at the last stage
                        has_workload_to_share = has_workload_to_share && !last_stage_flag[i].load();
                        // for thread i to have workload to share, it cannot have invalid inner index
                        has_workload_to_share = has_workload_to_share && last_stage_inner_index_ub[i].load() > 0;

                        if (has_workload_to_share){
                            outer_idx = outer_indices[i];
                            inner_idx = (last_stage_inner_index_ub[i]--) - 1;
                            other_thread_id = i;
                            if (inner_idx < n_innerloop && outer_idx < n_outerloop){
                                break;
                            }
                        }
                    }
                    // work on the task from the other thread
                    if (other_thread_id != thread_id && outer_idx < n_outerloop && inner_idx < n_innerloop){
                        bool test_result;
                        if constexpr(f_has_3_args){
                            test_result = f(outer_idx, inner_idx, thread_id);
                        } else {
                            test_result = f(outer_idx, inner_idx);
                        }
                        if (!test_result){
                            last_stage_flag[other_thread_id].store(true);
                            last_round_mutex.lock();
                            negative_results_last_round.push_back(outer_idx);
                            last_round_mutex.unlock();
                        }
                    } else {
                        // if there is nothing to work on, end this thread
                        break;
                    }
                } else {
                    // if this thread has been assigned a meaningful outer index, loop over all the inner indices
                    outer_indices[thread_id] = outer_idx;
                    // restart the inner index
                    UINT_IN inner_idx = 0;
                    while (true){
                        if (!some_thread_entered_last_stage.load()){
                            if (inner_idx >= n_innerloop){
                                break;
                            }
                        } else {
                            if (inner_idx >= last_stage_inner_index_ub[thread_id] || last_stage_flag[thread_id]){
                                break;
                            }
                        }

                        bool test_result;
                        if constexpr(f_has_3_args){
                            test_result = f(outer_idx, inner_idx, thread_id);
                        } else {
                            test_result = f(outer_idx, inner_idx);
                        }
                        if (!test_result){
                            result[outer_idx] = 0;
                            break;
                        }
                        ++inner_idx;
                    }
                }
            }
            return true;
        };


        unsigned i = 0;
        std::vector<std::shared_ptr<Notification>> futures;
        while (i < threadsManagement.n_threads_in_total()){
            auto future = threadsManagement.template submit_task_that_starts_immediately(worker, i);
            if (future == nullptr){
                break;
            } else{
                futures.template emplace_back(std::move(future));
                ++i;
            };
        }
        worker(i);
        for (auto & future : futures){
            future->wait_until_done();
        }

        for (auto i : negative_results_last_round){
            result[i] = 0;
        }

        return result;
    }









} // end of namespace







#endif //MULTITHREAD_EXAMPLES_MULTITHREAD_LOOP_H
