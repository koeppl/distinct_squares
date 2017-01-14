Computing All Distinct Squares Efficiently
==========================================

This is a C++ implementation of the paper 
"Computing All Distinct Squares in Linear Time for Integer Alphabets"
by Hideo Bannai, Shunsuke Inenaga and Dominik Köppl, submitted to ArXiv (https://arxiv.org/abs/1610.03421)

# Goal

The goal is to report the first occurrences of all substrings of type `AA` of the input once, where `A` is an arbitrary string.

# Dependencies

- Command line tools
  - cmake
  - make
  - a C++11 compiler like gcc or clang 
  - git and svn to clone and build the external dependecies

The following dependencies will be downloaded by cmake and used to compile this project:
- glog
- gflags
- gtest
- sdsl-lite

# Running
The build system will create an executable `distinct_squares` that can process a text file.
The other files can be used for testing. `check_file` is used to check whether the output of `distinct_squares` is correct
by validating the output with a naive algorithm.
The executable `test_squares` runs different tests on random generated data.
