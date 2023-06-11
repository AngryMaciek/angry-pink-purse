# Toolset to extract and plot genomic coverages from bam files
*Maciej Bak  
wsciekly.maciek@gmail.com*

1. Create environment: `conda env create`
2. Activate environment: `conda activate pybind11-computing`
3. [Compile Cpp into a shared library](https://pybind11.readthedocs.io/en/stable/compiling.html#building-manually):
    ```
    g++ \
    -O2 \
    -Wall \
    -shared \
    -std=c++14 \
    -undefined dynamic_lookup \
    -fPIC \
    $(python3 -m pybind11 --includes) \
    test.cpp \
    -o functions$(python3-config --extension-suffix) \
    -I/Users/maciek/miniconda3/envs/pybind11-computing/include/eigen3 \
    -I/Users/maciek/miniconda3/envs/pybind11-computing/include \
    -L/Users/maciek/miniconda3/envs/pybind11-computing/lib/ \
    -larmadillo \
    -lgsl
    ```
4. Execute Python testscript: `python test.py`
