cmake_minimum_required(VERSION 3.15)

project(CFD_Solver C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)

if(NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "PETSC_DIR is not set. Run: export PETSC_DIR=$HOME/petsc")
endif()

if(NOT DEFINED ENV{PETSC_ARCH})
    message(FATAL_ERROR "PETSC_ARCH is not set. Run: export PETSC_ARCH=arch-c-debug")
endif()

set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

pkg_check_modules(PETSC REQUIRED PETSc)

add_executable(cfd
    src/cfd.c
    src/grid.c
    src/math_helpers.c
    src/solver.c
    src/setup.c
    include/math_helpers.h
    include/solver.h
    include/setup.h
    include/config.h
    include/grid.h
    include/cfd.h
)

target_include_directories(cfd PRIVATE
    include
    ${PETSC_INCLUDE_DIRS}
)

target_compile_options(cfd PRIVATE
    ${PETSC_CFLAGS_OTHER}
)

target_link_directories(cfd PRIVATE
    ${PETSC_LIBRARY_DIRS}
)

target_link_libraries(cfd PRIVATE
    ${PETSC_LIBRARIES}
    m
)