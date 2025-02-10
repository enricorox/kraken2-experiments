cmake_minimum_required(VERSION 3.26)
project(kraken2_experiments)

set(CMAKE_CXX_STANDARD 14)

add_executable(evaluate_calls
        evaluate_calls_marghe.cc)
