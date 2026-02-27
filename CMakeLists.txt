cmake_minimum_required(VERSION 3.15)
project(CFD_Solver C)
set(CMAKE_C_STANDARD 99)

add_executable(cfd
  src/main.c
)