cmake_minimum_required(VERSION 3.15)
project(CFD_Solver C)
set(CMAKE_C_STANDARD 99)

find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})

add_executable(cfd
  "src/cfd.c"
 "src/grid.c")

 target_include_directories(cfd PRIVATE include)
 target_link_libraries(cfd PUBLIC MKL::MKL)