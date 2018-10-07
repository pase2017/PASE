# cmake project 信息
cmake_minimum_required(VERSION 2.8)
project(test_pase_matrix C)

# 编译路径
message("")
message("project source dir: ${PROJECT_SOURCE_DIR}")
message("project binary dir: ${PROJECT_BINARY_DIR}")
message("")
if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  message(FATAL_ERROR "In-place build not allowed! Please use a separate build directory.")
endif()

# 编译器和优化选项
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_C_FLAGS    "-g -Wall -O2")

# 添加 find_package(xx) 路径
if(NOT cmake_module_root)
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../../config)
else()
  set(CMAKE_MODULE_PATH ${cmake_module_root})
endif()

set(use_hypre 1)

# 设置 pase 路径
if(NOT pase_root)
  set(pase_root ${PROJECT_SOURCE_DIR}/../../pase)
endif()
find_package(pase)
include_directories(${pase_INCLUDE_DIRS})

# 设置 hypre 路径
if (NOT hypre_root)
  set(hypre_root ${PROJECT_SOURCE_DIR}/../hypre)
endif()
find_package(hypre)
include_directories(${hypre_INCLUDE_DIRS})

# 编译测试算例
file(GLOB_RECURSE test_pase_src ${PROJECT_SOURCE_DIR}/source/*.c)
set(pase_bin_name "test_pase_matrix.exe")
add_executable(${pase_bin_name} ${test_pase_src})
target_link_libraries(${pase_bin_name} ${pase_LIBRARIES} ${hypre_LIBRARIES} -lm)
