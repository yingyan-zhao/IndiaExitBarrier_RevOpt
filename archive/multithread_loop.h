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
#include <set>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace MultiThreads{

    /**
     * This class handles notification that a job is done.
     */
    class Notification{
        std::promise<void> signal;
        std::shared_future<void> fut = signal.get_future().share();
    public:
        /**
         * This function is used on the receiver end. It blocks the receiver thread
         */
        inline void wait_until_done(){
            // block current thread
            fut.wait();
        }

        /**
         * This function is used on the sender end. It send a notification to the receiver that the job is done
         */
        inline void is_done(){
            signal.set_value();
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
        std::shared_ptr<Notification> submit_task_only_when_idle(F && f, Args && ... args);

        template<class F, class... Args>
        std::shared_ptr<Notification> submit_task(F && f, Args && ... args);

        ~Threads_Management();

        std::mutex print_mutex;
    private:
        // stop flag, if true, all worker will stop waiting for new jobs
        bool stop;

        // maximum threads we will use
        const uint16_t n_threads;

        // queue of tasks
        std::queue<std::function<void()>> tasks = {};

        // for synchonization
        std::mutex queue_mutex;

        // this condition variable is used to tell thread to wait or start a task
        std::vector<std::unique_ptr<std::promise<void>>> signals;

        // thread pool
        std::vector<std::thread> threads;

        // recording idle threads
        std::vector<std::uint8_t> idle_status;
        constexpr static std::uint8_t IDLE = 1;
        constexpr static std::uint8_t BUSY = 0;

        // thread identifier
        std::set<std::thread::id> all_thread_id;

        bool at_least_one_thread_is_idle();

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
        } else {
            std::atomic<UINT> progress(0);

            auto time_stamp = std::chrono::high_resolution_clock::now();

            auto worker_part = [&](unsigned thread_id,
                                   uint64_t & finished_loop,
                                   const decltype(time_stamp) & start_time_stamp){
                unsigned n_tasks_this_loop = 1;
                if (finished_loop > 3){
                    auto elapsed = std::chrono::high_resolution_clock::now() - start_time_stamp;
                    uint64_t duration =  std::chrono::duration_cast<std::chrono::microseconds>( elapsed).count() + 1;

                    constexpr uint64_t max_continuous_loop = std::min(static_cast<uint64_t>(1000), static_cast<uint64_t>(std::numeric_limits<UINT>::max()));
                    n_tasks_this_loop = std::min(std::max((500 * finished_loop) / duration, static_cast<uint64_t>(1)),
                                                 max_continuous_loop);
                }

                UINT start_idx = progress.fetch_add(n_tasks_this_loop);
                UINT end_idx = std::min(static_cast<UINT>(start_idx + n_tasks_this_loop), n);
                for (UINT k = start_idx; k < end_idx; ++k){
                    if constexpr(f_has_2_args){
                        f(k, thread_id);
                    } else {
                        f(k);
                    }
                    finished_loop++;
                }

                bool done = end_idx >= n;
                return done;
            };

            auto worker = [&](unsigned thread_id){
                auto start_time_stamp = std::chrono::high_resolution_clock::now();
                uint64_t finished_loop = 0;
                bool done = false;
                while(!done){
                    done = worker_part(thread_id, finished_loop, start_time_stamp);
                }
            };

            std::vector<std::shared_ptr<Notification>> finish_signals;

            auto manager_part = [&](){
                // get more workers
                if (finish_signals.size() + 1 < maximum_n_threads){
                    auto now = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_stamp).count();
                    if (duration > 200){
                        auto future = threadsManagement.submit_task_only_when_idle(worker, finish_signals.size());
                        if (future != nullptr){
                            finish_signals.emplace_back(std::move(future));
                        }
                        time_stamp = std::chrono::high_resolution_clock::now();
                    }
                }
            };

            auto manager = [&](){
                auto start_time_stamp = std::chrono::high_resolution_clock::now();
                uint64_t finished_loop = 0;
                unsigned thread_id = maximum_n_threads - 1;
                bool done = false;
                while(!done){
                    done = worker_part(thread_id, finished_loop, start_time_stamp);
                    auto prog = progress.load(std::memory_order_relaxed);
                    if (!done && n > prog && static_cast<std::size_t>(n - prog) > finish_signals.size()){
                        manager_part();
                    }
                }
            };


            for (unsigned i = 0; i + 1 < maximum_n_threads; ++i){
                auto future = threadsManagement.submit_task_only_when_idle(worker, finish_signals.size());
                if (future != nullptr){
                    finish_signals.push_back(future);
                }
            }
//            worker(finish_signals.size());
            manager();
            for (auto & future : finish_signals){
                future->wait_until_done();
            }
        }
    }

    template<typename F>
    void simple_parallel_for(F & f, std::int64_t n, MultiThreads::Threads_Management & threadsManagement, unsigned maximum_n_threads){
        if (n < 0){
            throw std::runtime_error("Error in simple_parallel_for. n is negative");
        }
        if (n < std::numeric_limits<std::uint32_t>::max()){
            simple_parallel_for_templated_integer(f, static_cast<std::uint32_t>(n), threadsManagement, maximum_n_threads);
        } else {
            simple_parallel_for_templated_integer(f, static_cast<std::uint64_t>(n), threadsManagement, maximum_n_threads);
        }
    }

    template<class F, class... Args>
    std::shared_ptr<Notification> Threads_Management::submit_task(F &&f, Args &&...args) {
        std::shared_ptr<Notification> future_pointer = std::make_shared<Notification>();
        bool deadlock_possible = false;

        std::unique_lock<std::mutex> lock(queue_mutex);

        // check the possibility of deadlock
        auto this_thread_is_within_pool = all_thread_id.find(std::this_thread::get_id()) != all_thread_id.end();
        deadlock_possible = this_thread_is_within_pool && !at_least_one_thread_is_idle();

        // if deadlock is impossible, put it into the queue.
        if (!deadlock_possible){
            // push task into queue
            if constexpr (sizeof...(Args) == 0){
                // if there is no argument
                tasks.emplace([future_pointer, f = std::move(f)]{
                    f();
                    future_pointer->is_done();
                });
            } else {
                // if there exist at least one argument
#if (__cpp_init_captures >= 201803)
                tasks.emplace([future_pointer, f = std::move(f), ...args = std::move(args)]{
            f(args...);
            future_pointer->is_done();
        });
#else
                auto f_function = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
                tasks.emplace([future_pointer, f_function = std::move(f_function)]{
                    f_function();
                    future_pointer->is_done();
                });
#endif
            }

            // wake some idle thread to run it
            for (uint16_t i = 0; i < n_threads; ++i){
                if (idle_status[i] == IDLE){
                    idle_status[i] = BUSY;
                    signals[i]->set_value();
                    break;
                }
            }
        }
        lock.unlock();

        // if deadlock is possible, finish the task right now.
        if (deadlock_possible){
            f(args...);
            future_pointer->is_done();
        }

        return future_pointer;
    }

    template<class F, class... Args>
    std::shared_ptr<Notification> Threads_Management::submit_task_only_when_idle(F &&f, Args &&...args) {
        std::shared_ptr<Notification> future_pointer(nullptr);
        std::unique_lock<std::mutex> lock(queue_mutex);

        if (at_least_one_thread_is_idle()){
            future_pointer = std::make_shared<Notification>();

            // push task into queue
            if constexpr (sizeof...(Args) == 0){
                // if there is no argument
                tasks.emplace([future_pointer, f = std::move(f)]{
                    f();
                    future_pointer->is_done();
                });
            } else {
                // if there exist at least one argument
#if (__cpp_init_captures >= 201803)
                tasks.emplace([future_pointer, f = std::move(f), ...args = std::move(args)]{
                    f(args...);
                    future_pointer->is_done();
                });
#else
                auto f_function = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
                tasks.emplace([future_pointer, f_function = std::move(f_function)]{
                    f_function();
                    future_pointer->is_done();
                });
#endif
            }

            for (uint16_t i = 0; i < n_threads; ++i){
                if (idle_status[i] == IDLE){
                    idle_status[i] = BUSY;
                    signals[i]->set_value();
                    break;
                }
            }
        }
        lock.unlock();

        return future_pointer;
    }


    inline
    Threads_Management::~Threads_Management() {
        // wake up all threads and tell them to stop
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;

//        // debug
//        unsigned n_idle_threads = 0;
//        for (unsigned  i = 0; i < n_threads; ++i){
//            n_idle_threads += idle_status[i];
//        }
//        std::cout << "n_idle_thread: " << n_idle_threads << " out of " << n_threads << std::endl;
//        std::cout << "size of threads: " << threads.size() << std::endl;
//        std::cout << "stop set to true" << std::endl;
//        // debug

        for (uint16_t i = 0; i < n_threads; ++i){
            if (signals[i]) {
                try { signals[i]->set_value(); } catch (...) {}
            }
        }
        lock.unlock();

        // wait for them to join the main thread
        for (auto & thread : threads){
            thread.join();
        }
    }

    inline
    uint16_t Threads_Management::n_threads_in_total() const {
        return n_threads;
    }

    inline bool MultiThreads::Threads_Management::at_least_one_thread_is_idle() {
        bool has_idle_thread = false;
        for (uint16_t i = 0; i < n_threads; ++i){
            has_idle_thread = (idle_status[i] == IDLE);
            if (has_idle_thread){
                break;
            }
        }
        return has_idle_thread;
    }

    inline MultiThreads::Threads_Management::Threads_Management(int64_t n_threads)
    : stop(false)
    , n_threads(n_threads)
    , signals(n_threads)
    , idle_status(n_threads, IDLE)
    {
        if (n_threads > std::numeric_limits<uint16_t>::max()){
            throw std::runtime_error("n_threads is too large");
        }
        if (n_threads < 0){
            throw std::runtime_error("n_threads must be nonnegative");
        }

        create_threads();
    }

    inline void MultiThreads::Threads_Management::create_threads() {
        threads.reserve(n_threads);

        auto worker = [&](uint16_t thread_id) {
            bool keep_looping = true;

            std::unique_lock<std::mutex> lock(queue_mutex);
            all_thread_id.insert(std::this_thread::get_id());
            lock.unlock();

            while(keep_looping){
                bool task_assigned = false;
                std::function<void()> task;

                // critical section
                lock.lock();
                {
                    // block the current thread if we have nothing to do
                    if (tasks.empty() && !stop){
                        idle_status[thread_id] = IDLE;
                        signals[thread_id] = std::make_unique<std::promise<void>>();
                        lock.unlock();

                        signals[thread_id]->get_future().wait();
                        lock.lock();
                    }

                    // download the work if we have something to do
                    if (!tasks.empty()){
                        // plan to execute the first element of queue
                        task = std::move(tasks.front());
                        tasks.pop();
                        task_assigned = true;
                    }
                    keep_looping = !stop;
                }
                lock.unlock();

                if (task_assigned){
                    task();         // run task
                }
            }
        };

        for (uint16_t i = 0; i < n_threads; ++i){
            threads.emplace_back(worker, i);
        }

        // make sure all threads have been created
        bool ready = false;
        while(!ready){
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            ready = true;
            std::unique_lock<std::mutex> lock(queue_mutex);
            for (uint16_t i = 0; i < n_threads; ++i){
                ready = ready && (idle_status[i] == IDLE) && signals[i] != nullptr;
            }
            lock.unlock();
        }
    }




} // end of namespace







#endif //MULTITHREAD_EXAMPLES_MULTITHREAD_LOOP_H
