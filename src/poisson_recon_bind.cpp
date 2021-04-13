#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_eigen(py::module&);
void bind_create_poisson_surface(py::module&);

PYBIND11_MODULE(_poisson_recon_pybind, m)
{
    m.doc() = "Python binding for PoissonRecon";
    bind_eigen(m);
    bind_create_poisson_surface(m);
}
