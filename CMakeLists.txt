# cmake project 信息
cmake_minimum_required(VERSION 2.8)
project(PASE C)

# 版本号
set(PASE_VERSION_MAJOR 2)
set(PASE_VERSION_MINOR 7)

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
set(CMAKE_C_FLAGS    "-g -O2 -W -Wall -Wextra -Wfloat-equal")

# 添加 find_package(xx) 路径
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# 拷贝库文件至如下路径
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/pase/lib)
# 拷贝头文件至如下路径
set(pase_include_output_path ${PROJECT_BINARY_DIR}/pase/include)

# pase_kernel 库名
set(pase_library_kernel "pase_kernel")
# 添加子路径
include_directories(${PROJECT_SOURCE_DIR}/source/kernel)
add_subdirectory   (${PROJECT_SOURCE_DIR}/source/kernel)

if(use_hypre)
  if (NOT hypre_root)
    set(hypre_root ${PROJECT_SOURCE_DIR}/../hypre)
  endif()
  find_package(hypre)
  if(hypre_FOUND)
    message("Find hypre dir: ${hypre_INCLUDE_DIRS}")
    message("Find hypre lib: ${hypre_LIBRARIES}")
    message("")
    include_directories(${hypre_INCLUDE_DIRS})
    #set(pase_extern_library "${pase_extern_library} ${hypre_LIBRARIES}")
    # pase_hypre 库名
    set(pase_library_hypre "pase_hypre")
    # 添加子路径
    include_directories(${PROJECT_SOURCE_DIR}/source/app_interface/app_hypre)
    add_subdirectory   (${PROJECT_SOURCE_DIR}/source/app_interface/app_hypre)
  endif()
endif()

# 拷贝测试文件至如下目录
message("Copy test files to ${PROJECT_BINARY_DIR}...")
file(COPY ${PROJECT_SOURCE_DIR}/test DESTINATION ${PROJECT_BINARY_DIR})
message("Copy test files to ${PROJECT_BINARY_DIR}... done.")
message("")
