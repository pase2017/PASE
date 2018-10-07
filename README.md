# Parallel Auxiliary Subspace Eigensolver

[TOC]

## 1. 编译和运行

在 `PASE/tool` 文件夹下定义了一些编译脚本，下面举例说明如何将 `PASE` 编译为静态库。

```bash
mkdir build
cd build
../tool/build-linux.sh
make
```

目前 `PASE` 需要依赖 `HYPRE` 软件，因此在编译脚本（如 `build-linux.sh` 文件）中需指定 `HYPRE` 的安装路径 `hypre_root`。程序会自动寻找该路径下的 `include` 和 `lib` 文件夹并添加至编译和链接依赖中。一个典型的编译脚本如下。

```bash
cmake -Duse_hypre=1 -Dhypre_root=~/Software/hypre-2.11.2/src/hypre \
      ../
```



由于时间关系，新版 `PASE` 尚未加入示例程序，仅在 `PASE/test` 文件夹中给出矩阵和向量操作的测试程序。类似于静态库的编译，测试程序可按如下方式编译为可执行文件。

```bash
cd test/test_pase_matrix_vector
mkdir build
cd build
../tool/build-linux.sh
make
```

对于测试程序，一个典型的编译脚本如下。

```bash
cmake -Dcmake_module_root=/media/psf/Home/Workspace/pase-2.7/config \
      -Dpase_root=/media/psf/Home/Workspace/pase-2.7/build-linux/pase \
      -Duse_hypre=1 -Dhypre_root=/home/ycg/Software/hypre-2.11.2/src/hypre \
      ../
```

> 在将 `PASE` 编译为静态库时，程序会把 `test` 文件夹打包至 `build` 文件夹。由于新版 `PASE` 尚处于早期阶段，并未在 `PASE/build/test` 文件夹下开展测试，而是直接在 `PASE/test` 中开展。以后将逐渐规范。



## 2. 程序文档

用户可在 `PASE/doc` 目录下找到 `PASE` 软件的详细介绍。



## 3. 软件架构

`PASE` 软件提供了常用软件包中矩阵和向量等数据结构的统一封装，并以 `Operation Free` 的形式，基于中国科学院计算数学与科学工程计算研究所谢和虎研究员及其研究小组提出的多重校正算法，进行代数特征值问题的高效并行求解。



在 `PASE` 软件中，核心算法的实现位于 `source/kernel` 文件夹；外部接口则位于 `source/app_interface` 文件夹，可与  `HYPRE` 等软件直接对接。



## 4. 算法实现

待补充。