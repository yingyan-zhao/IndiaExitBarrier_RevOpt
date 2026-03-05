#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include "multithread_loop.h"

void test_simple_paralle_for(MultiThreads::Threads_Management & threadsManagement){
    std::vector<double> abc(100);
    auto f = [&](int i){
        abc[i] = std::sin(i);
    };
    simple_parallel_for(f, 100, threadsManagement);
    // check results
    for (int i = 0; i < 100; ++i){
        if (abc[i] != std::sin(i)){
            std::cout << "Error!" << std::endl;
        }
    }
}

void test_submit_task(){
    int sleep_duration = 1;
    std::vector<double> result(2);

    auto worker = [&](int i){
        std::this_thread::sleep_for(std::chrono::seconds(sleep_duration));
        result[i] = 0;
    };

    MultiThreads::Threads_Management threadsManagement;
    auto notification0 = threadsManagement.submit_task(worker, 0);
    auto notification1 = threadsManagement.submit_task(worker, 1);

    notification0->wait_until_done();
    notification1->wait_until_done();
}


void test_nested_parallel_for(MultiThreads::Threads_Management & threadsManagement){
    std::mt19937 rng(14242);
    std::uniform_real_distribution<double> uniform(0, 1);

    const unsigned n_innerloop = 1000;
    const unsigned n_outerloop = 10;
    std::vector<double> M(n_innerloop * n_outerloop);
    for (uint i = 0; i < M.size(); ++i){
        M[i] = uniform(rng);
    }

    auto f = [&](uint outer_idx, uint inner_idx){
        bool result = M[inner_idx + outer_idx * n_innerloop] * n_innerloop > 1;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        return result;
    };
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    std::vector<int8_t> result_seq(n_outerloop, 1);
    for (uint i = 0; i < n_outerloop; ++i){
        for (uint j = 0; j < n_innerloop; ++j){
            if (!f(i, j)){
                result_seq[i] = 0;
                break;
            }
        }
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    std::cout << "sequential loop finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0 << "[s]" << std::endl;

    std::vector<int8_t> result_parallel = nested_parallel_for(f, n_outerloop, n_innerloop, threadsManagement);

    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    std::cout << "sequential loop finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() / 1000.0 << "[s]" << std::endl;

    int k = 0;
    for (uint i = 0; i < n_outerloop; ++i){
        if (result_seq[i] != result_parallel[i]){
            throw std::runtime_error("Error!");
        }
        k += result_seq[i];
    }
    std::cout << "sum result = " << k << std::endl;
}

int main() {
    MultiThreads::Threads_Management threadsManagement;

    test_simple_paralle_for(threadsManagement);
    test_submit_task();
//    test_nested_parallel_for(threadsManagement);

    return 0;
}
