from skbuild import setup

setup(
    name="pyseb",
    version="0.1.0",
    description="Python bindings for the Scattering Equation Builder (SEB) library",
    author="SEB Team",
    license="MIT",
    packages=["pyseb"],
    cmake_install_dir="pyseb",
    cmake_args=["-DCMAKE_BUILD_TYPE=Release"],
    python_requires=">=3.7",
    install_requires=["numpy", "matplotlib"],
)
