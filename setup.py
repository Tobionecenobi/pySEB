import os
from skbuild import setup
from skbuild.constants import CMAKE_BUILD_DIR
from skbuild.cmaker import get_cmake_version

setup(
    packages=["pyseb"],
    cmake_install_dir="pyseb",
    include_package_data=True,
    python_requires=">=3.7",
    cmake_args=[
        "-DBUILD_PYTHON=ON",
        "-DSEB_ENABLE_GINAC_BACKEND=OFF",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        f"-DCMAKE_INSTALL_RPATH=$ORIGIN",
        "-DCMAKE_BUILD_TYPE=Release",
    ],
)
