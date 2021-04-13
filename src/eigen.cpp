// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2018 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------

// Python bindings for vectors of Eigen vectors: Vector3dVector and Vector3iVector
// These vectors are consumed by poisson_recon_pybind.create_poisson_surface
//
//
// Usage:
//  import numpy as np
//  import trimesh
//  from poisson_recon_pybind import Vector3dVector, Vector3iVector, create_poisson_surface
//  points = np.arrange(300).reshape((3, 100)).astype(float)
//  normal = np.zeros_like(points)
//  vertices, triangles = create_poisson_surface(
//      Vector3dVector(points), Vector3dVector(normals)
//  )
//  mesh = trimesh.Trimesh(
//    vertices=np.asarray(vertices), faces=np.asarray(triangles)
//  )
//

// The content of this file has been extracted from
// https://github.com/intel-isl/Open3D/blob/master/cpp/pybind/utility/eigen.cpp
// commit 37b1b46810c153486df840a074246e823b7665f8 (HEAD -> master, origin/master, origin/HEAD)
// Date:   Wed Mar 17 22:26:22 2021 -0400
//
// Dependencies on Open3D have been removed.


#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>  // Include first to suppress compiler warnings
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <tuple>
#include <vector>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<Eigen::Vector3d>);
PYBIND11_MAKE_OPAQUE(std::vector<Eigen::Vector3i>);

namespace pybind11 {

template <typename Vector,
          typename holder_type = std::unique_ptr<Vector>,
          typename... Args>
py::class_<Vector, holder_type> bind_vector_without_repr(
        py::module &m, std::string const &name, Args &&... args) {
    // hack function to disable __repr__ for the convenient function
    // bind_vector()
    using Class_ = py::class_<Vector, holder_type>;
    Class_ cl(m, name.c_str(), std::forward<Args>(args)...);
    cl.def(py::init<>());
    cl.def(
            "__bool__", [](const Vector &v) -> bool { return !v.empty(); },
            "Check whether the list is nonempty");
    cl.def("__len__", &Vector::size);
    return cl;
}

// - This function is used by Pybind for std::vector<SomeEigenType> constructor.
//   This optional constructor is added to avoid too many Python <-> C++ API
//   calls when the vector size is large using the default biding method.
//   Pybind matches np.float64 array to py::array_t<double> buffer.
// - Directly using templates for the py::array_t<double> and py::array_t<int>
//   and etc. doesn't work. The current solution is to explicitly implement
//   bindings for each py array types.
template <typename EigenVector>
std::vector<EigenVector> py_array_to_vectors_double(
        py::array_t<double, py::array::c_style | py::array::forcecast> array) {
    size_t eigen_vector_size = EigenVector::SizeAtCompileTime;
    if (array.ndim() != 2 || array.shape(1) != eigen_vector_size) {
        throw py::cast_error();
    }
    std::vector<EigenVector> eigen_vectors(array.shape(0));
    auto array_unchecked = array.mutable_unchecked<2>();
    for (auto i = 0; i < array_unchecked.shape(0); ++i) {
        // The EigenVector here must be a double-typed eigen vector, since only
        // open3d::Vector3dVector binds to py_array_to_vectors_double.
        // Therefore, we can use the memory map directly.
        eigen_vectors[i] = Eigen::Map<EigenVector>(&array_unchecked(i, 0));
    }
    return eigen_vectors;
}

template <typename EigenVector>
std::vector<EigenVector> py_array_to_vectors_int(
        py::array_t<int, py::array::c_style | py::array::forcecast> array) {
    size_t eigen_vector_size = EigenVector::SizeAtCompileTime;
    if (array.ndim() != 2 || array.shape(1) != eigen_vector_size) {
        throw py::cast_error();
    }
    std::vector<EigenVector> eigen_vectors(array.shape(0));
    auto array_unchecked = array.mutable_unchecked<2>();
    for (auto i = 0; i < array_unchecked.shape(0); ++i) {
        eigen_vectors[i] = Eigen::Map<EigenVector>(&array_unchecked(i, 0));
    }
    return eigen_vectors;
}
}  // namespace pybind11

template <typename EigenVector,
          typename Vector = std::vector<EigenVector>,
          typename holder_type = std::unique_ptr<Vector>,
          typename InitFunc>
py::class_<Vector, holder_type> pybind_eigen_vector_of_vector(
        py::module &m,
        const std::string &bind_name,
        const std::string &repr_name,
        InitFunc init_func) {
    typedef typename EigenVector::Scalar Scalar;
    auto vec = py::bind_vector_without_repr<std::vector<EigenVector>>(
            m, bind_name, py::buffer_protocol());
    vec.def(py::init(init_func));
    vec.def_buffer([](std::vector<EigenVector> &v) -> py::buffer_info {
        size_t rows = EigenVector::RowsAtCompileTime;
        return py::buffer_info(v.data(), sizeof(Scalar),
                               py::format_descriptor<Scalar>::format(), 2,
                               {v.size(), rows},
                               {sizeof(EigenVector), sizeof(Scalar)});
    });
    vec.def("__repr__", [repr_name](const std::vector<EigenVector> &v) {
        return repr_name + std::string(" with ") + std::to_string(v.size()) +
               std::string(" elements.\n") +
               std::string("Use numpy.asarray() to access data.");
    });
    vec.def("__copy__", [](std::vector<EigenVector> &v) {
        return std::vector<EigenVector>(v);
    });
    vec.def("__deepcopy__", [](std::vector<EigenVector> &v, py::dict &memo) {
        return std::vector<EigenVector>(v);
    });

    // py::detail must be after custom constructor
    using Class_ = py::class_<Vector, std::unique_ptr<Vector>>;
    py::detail::vector_if_copy_constructible<Vector, Class_>(vec);
    py::detail::vector_if_equal_operator<Vector, Class_>(vec);
    py::detail::vector_modifiers<Vector, Class_>(vec);
    py::detail::vector_accessor<Vector, Class_>(vec);

    return vec;
}

void bind_eigen(py::module &m) {
    auto vector3dvector = pybind_eigen_vector_of_vector<Eigen::Vector3d>(
            m, "Vector3dVector", "std::vector<Eigen::Vector3d>",
            py::py_array_to_vectors_double<Eigen::Vector3d>);

    auto vector3ivector = pybind_eigen_vector_of_vector<Eigen::Vector3i>(
            m, "Vector3iVector", "std::vector<Eigen::Vector3i>",
            py::py_array_to_vectors_int<Eigen::Vector3i>);
}
