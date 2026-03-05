cmake_minimum_required(VERSION 3.19)
project(multithread_examples)

set(CMAKE_CXX_STANDARD 20)

find_package(Threads REQUIRED)

add_executable(example main.cpp multithread_loop.h)
target_link_libraries(example Threads::Threads)
