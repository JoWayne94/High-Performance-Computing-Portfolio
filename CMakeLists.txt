cmake_minimum_required(VERSION 3.1)
project (Solver)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(Solver LidDrivenCavitySolver.cpp LidDrivenCavity.cpp PoissonSolver.cpp ProgramInputs.cpp MpiImplementation.cpp)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(MPI REQUIRED)

include_directories(${LAPACK_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(Solver ${BLAS_LIBRARIES} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES}) 
