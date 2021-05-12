import os
import sys
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

import imp


VERSION = imp.load_source(
    "poisson_recon_pybind.version", "poisson_recon_pybind/version.py"
).VERSION
REQUIRED_NUMPY_VERSION = "numpy>=1.12.0"
MIN_CPU_CORES = 2


def get_cpu_count():
    try:
        return len(os.sched_getaffinity(0))  # linux only
    except AttributeError:
        pass

    try:
        return os.cpu_count()  # python 3.4+
    except AttributeError:
        return 1  # default


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            _ = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ["--", "-j{}".format(max(MIN_CPU_CORES, get_cpu_count()))]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ['cmake', '--build', '.'] + build_args, cwd=self.build_temp
        )


install_requires = [REQUIRED_NUMPY_VERSION]


setup(
    name='poisson-recon-pybind',
    version=VERSION,
    author='BBP',
    url='https://bbpgitlab.epfl.ch/nse/poisson-recon-pybind',
    author_email='luc.guyot@epfl.ch',
    description='A Python binding for the surface reconstruction '
    'method of PoissonRecon',
    long_description='',
    packages=[
        'poisson_recon_pybind',
    ],
    ext_modules=[
        CMakeExtension('poisson_recon_pybind._poisson_recon_pybind'),
    ],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=install_requires,
)
