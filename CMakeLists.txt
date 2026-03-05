cmake_minimum_required(VERSION 3.18)
project(StructureEstimationC)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OSX_ARCHITECTURES arm64)

#find_package(Armadillo REQUIRED)

#target_include_directories(StructureEstimationC PUBLIC ${ARMADILLO_INCLUDE_DIRS})
#target_link_libraries(StructureEstimationC PUBLIC ${ARMADILLO_LIBRARIES})

# ------------------------------------------------------------------------------
# Part 1: Compiler Flags
# ------------------------------------------------------------------------------
add_compile_definitions(SOURCE_FOLDER="${PROJECT_SOURCE_DIR}")
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -mllvm -inline-threshold=20000 -mllvm -inlinehint-threshold=20000")
    else()
        set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -mllvm -inline-threshold=20000 -mllvm -inlinehint-threshold=20000")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -mllvm -inline-threshold=20000 -mllvm -inlinehint-threshold=20000")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
# using GCC
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -finline-limit=10000")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # using Visual Studio C++
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /arch")
endif()

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -fsanitize-address-use-after-scope  -fsanitize-address-use-after-return=always")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize-address-use-after-scope  -fsanitize-address-use-after-return=always")

#Find Eigen
if (DEFINED ENV{EIGEN3_PATH})
    include_directories($ENV{EIGEN3_PATH})
else()
    include_directories(${PROJECT_SOURCE_DIR}/ThirdParty/eigen-3.4.0)
endif()

#Find boost
if (DEFINED ENV{boost_PATH})
    include_directories($ENV{boost_PATH})
else()
    include_directories(${PROJECT_SOURCE_DIR}/ThirdParty/boost_1_78_0)
endif()

find_package(Threads REQUIRED)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:$ENV{HOME}/.local/lib/pkgconfig")
find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)
include_directories("$ENV{HOME}/.local/include")

if (DEFINED ENV{PCG_PATH}) # PCG
    include_directories($ENV{PCG_PATH}/include)
else()
    include_directories(${PROJECT_SOURCE_DIR}/ThirdParty/pcg-cpp-0.98/include)
endif()

add_executable(StructureEstimationC main.cpp multithread_loop.h Auxillaries.h Auxillaries.cpp Ipopt_wrapper.h
    DFS.h direct_search.h soo.h Interpolation_spline.cpp Interpolation_spline.h RandomDraw.h RandomDraw.cpp
    SimulationData.h SimulationData.cpp Basics.h Basics.cpp SolveModel.h SolveModel.cpp
    EstimationAuxiliary.cpp EstimationAuxiliary.h Estimation_IndirectInference.cpp Estimation_IndirectInference.h
    CounterfactualExercise.cpp CounterfactualExercise.h
    EstimationIndirectInference_StandardDeviation.cpp EstimationIndirectInference_StandardDeviation.h)
#
#
#

target_link_libraries(StructureEstimationC PUBLIC PkgConfig::Ipopt Threads::Threads)





