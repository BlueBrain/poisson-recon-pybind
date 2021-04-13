# poisson-recon-pybind

 [PoissonRecon](https://github.com/mkazhdan/PoissonRecon) Python binding with [pybind11](https://pybind11.readthedocs.io)

These bindings have been extracted from https://github.com/intel-isl/Open3D and have been made independant from the Open3D library.

## Method
Reconstruction of a 3D surface mesh out of a cloud of 3D points with assigned 3D unit vectors.
The implementation is provided by the Open3D fork https://github.com/intel-isl/Open3D-PoissonRecon of M. Kazhdan's original PoissonRecon (https://github.com/mkazhdan/PoissonRecon).

## Installation
```bash
$> cd poisson-recon-pybind
$> git submodule init
$> git submodule update
$> pip install -e .
```

## Test
```bash
$> pip install tox
$> tox
```

## Requirements
* cmake > 3.0.9
* Eigen C++ library (http://eigen.tuxfamily.org/)
* C++ compiler (with C++11 support)

