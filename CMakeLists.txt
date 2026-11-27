cmake_minimum_required(VERSION 3.24)
project(mixbench LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(LAPACK REQUIRED)
find_package(OpenSSL REQUIRED COMPONENTS Crypto)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(USE_LAPACK ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(SFLOAT OFF CACHE BOOL "" FORCE)
set(DLONG OFF CACHE BOOL "" FORCE)
FetchContent_Declare(scs URL https://codeload.github.com/cvxgrp/scs/tar.gz/refs/tags/3.2.8
  URL_HASH SHA256=22d2d785b7c7a9ee8a260d2684cf17ae4733271b8421fdbc78f281d19910ca1b)
FetchContent_Declare(spectra URL https://codeload.github.com/yixuan/spectra/tar.gz/refs/tags/v1.1.0
  URL_HASH SHA256=d29671e3d1b8036728933cadfddb05668a3cd6133331e91fc4535a9b85bedc79 SOURCE_SUBDIR unused)
FetchContent_Declare(json URL https://codeload.github.com/nlohmann/json/tar.gz/refs/tags/v3.11.3
  URL_HASH SHA256=0d8ef5af7f9794e3263480193c491549b2ba6cc74bb018906202ada498a79406 SOURCE_SUBDIR unused)
FetchContent_MakeAvailable(scs spectra json)
# In the embedded solver each call owns its workspace. Disable SCS's optional
# process-global SIGINT handler so independent OpenMP tasks can solve concurrently.
if(MSVC)
  target_compile_options(scsdir PRIVATE /UCTRLC)
else()
  target_compile_options(scsdir PRIVATE -UCTRLC)
endif()
add_library(mixing src/subgraphs.cpp src/budget.cpp src/graph.cpp src/spectrum.cpp src/directed_spectrum.cpp src/kernels.cpp src/fmmc.cpp src/analysis.cpp src/experiment.cpp)
target_include_directories(mixing PUBLIC include ${json_SOURCE_DIR}/single_include PRIVATE ${spectra_SOURCE_DIR}/include)
target_compile_definitions(mixing PUBLIC EIGEN_DONT_PARALLELIZE)
target_link_libraries(mixing PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX PRIVATE scsdir OpenSSL::Crypto)
add_executable(mixbench src/main.cpp)
target_link_libraries(mixbench PRIVATE mixing)
enable_testing()
add_executable(mix_tests tests/core.cpp)
target_link_libraries(mix_tests PRIVATE mixing)
add_test(NAME core COMMAND mix_tests)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME cli COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli.py
  $<TARGET_FILE:mixbench> ${CMAKE_CURRENT_SOURCE_DIR})
foreach(study IN ITEMS edge_ranking path_oracle nested_networks wall_budget)
  add_test(NAME studies_${study} COMMAND ${Python3_EXECUTABLE} -B
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/studies.py $<TARGET_FILE:mixbench> ${study})
endforeach()
add_test(NAME riffle COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/tests/riffle.py)
add_test(NAME bridges COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/tests/bridges.py $<TARGET_FILE:mixbench>)
add_test(NAME budget_repetitions COMMAND ${Python3_EXECUTABLE} -B
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/budget_repetitions.py $<TARGET_FILE:mixbench>)
add_test(NAME link_policy_trials COMMAND ${Python3_EXECUTABLE} -B
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/link_policy_trials.py --exe $<TARGET_FILE:mixbench>)
foreach(module IN ITEMS numpy scipy matplotlib)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c
    "from importlib.util import find_spec; import sys; sys.exit(find_spec('${module}') is None)"
    RESULT_VARIABLE missing_${module})
endforeach()
if(missing_numpy EQUAL 0 AND missing_scipy EQUAL 0)
  add_test(NAME power_grid COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/tests/power_grid.py
    --exe $<TARGET_FILE:mixbench>)
  add_test(NAME dana_parity COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_dana.py
    --exe $<TARGET_FILE:mixbench> --sizes 10 20 --jobs 2)
endif()
if(missing_numpy EQUAL 0)
  foreach(suite IN ITEMS dhn real_networks edge_policies)
    add_test(NAME ${suite} COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/tests/${suite}.py
      $<TARGET_FILE:mixbench>)
    set_tests_properties(${suite} PROPERTIES ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
  endforeach()
endif()
if(missing_matplotlib EQUAL 0)
  add_test(NAME plots COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/tests/plots.py
    $<TARGET_FILE:mixbench>)
endif()
