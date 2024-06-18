cmake_minimum_required(VERSION 3.4...3.18)
project(msastats)

add_subdirectory(pybind11)
pybind11_add_module(msastats 
                    src/MsaStatsCalculator.cpp 
                    src/read_seqs.cpp 
                    src/main.cpp)

target_compile_definitions(msastats
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})