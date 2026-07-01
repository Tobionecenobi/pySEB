from skbuild import setup

setup(
    packages=["pyseb"],
    cmake_install_dir="pyseb",
    include_package_data=True,
    license_files=[],
    python_requires=">=3.9",
    cmake_args=[
        "-DBUILD_PYTHON=ON",
        "-DSEB_ENABLE_GINAC_BACKEND=OFF",
        "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
        f"-DCMAKE_INSTALL_RPATH=$ORIGIN",
        "-DCMAKE_BUILD_TYPE=Release",
    ],
)
