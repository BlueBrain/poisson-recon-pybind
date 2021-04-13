"""
Module testing the Poisson Surface Reconstruction algorithm from M. Kazhdan.
"""
from poisson_recon_pybind import create_poisson_surface, Vector3dVector
import numpy as np
import numpy.testing as npt


def test_create_poisson_surface():
    points = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [1, 0, 0],
    ], dtype=float)
    normals = np.zeros_like(points)
    normals[:, :] = [0, 0, 1]
    vertices, triangles = create_poisson_surface(
        Vector3dVector(points), Vector3dVector(normals)
    )

    assert np.asarray(vertices).shape[1] == 3
    assert np.asarray(vertices).dtype == float
    assert np.asarray(triangles).shape[1] == 3
    assert np.asarray(triangles).dtype == np.int32
