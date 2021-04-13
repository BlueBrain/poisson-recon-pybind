add_subdirectory(pybind11 EXCLUDE_FROM_ALL)

add_library("poisson_recon" INTERFACE)
target_include_directories("poisson_recon" SYSTEM INTERFACE "PoissonRecon")
add_library(${PROJECT_NAME}::poisson_recon ALIAS "poisson_recon")
