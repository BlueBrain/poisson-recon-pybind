cmake_minimum_required(VERSION 3.9.0)

project(_poisson_recon_pybind)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(extlib EXCLUDE_FROM_ALL)

pybind11_add_module(_poisson_recon_pybind
    ${CMAKE_CURRENT_SOURCE_DIR}/src/poisson_recon_bind.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/eigen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/poisson_surface_reconstruction.cpp
    )

find_package(Eigen3 CONFIG REQUIRED)

target_link_libraries(_poisson_recon_pybind PRIVATE Eigen3::Eigen)
target_link_libraries(_poisson_recon_pybind PRIVATE ${PROJECT_NAME}::poisson_recon)
