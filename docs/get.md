---
title: Build the adapter
permalink: adapter-su2-get.html
keywords: adapter, su2, development, modules
summary: "Get SU2, get preCICE, execute adapter install script"
---

The adapter depends on [SU2](https://su2code.github.io/), [preCICE v3](https://precice.org/installation-overview.html), and the [preCICE Python bindings](https://precice.org/installation-bindings-python.html).

The script `su2AdapterInstall` replaces a few files in the SU2 source code. You then need to build SU2 from source, install it into a prefix (`SU2_RUN`) and add that to your `PATH`.

To run SU2, you can use the provided Python scripts `SU2_preCICE_CHT.py` and `SU2_preCICE_FSI.py`, which call SU2 via its Python interface.

## Get the SU2 source

Download SU2 v7.5.1 "Blackbird" source directly from the [SU2 releases on GitHub](https://github.com/su2code/SU2/releases/tag/v7.5.1). Note that both swig and mpi4py must be installed to use the SU2 Python wrapper, which needs to be enabled using the flag `-Denable-pywrapper=true` when building SU2.

## Build the adapter

In order to couple SU2 using preCICE, `python_wrapper_structure.cpp` and `CDriver.hpp` must be updated. This adapter provides the updated files. The shell script `su2AdapterInstall`, which comes with this adapter, automatically replaces the files in your SU2 directory with these updated files and provides the correct commands to re-configure and re-install SU2 with the added adjustments. This script and the build/install procedure rely on a few environment variables. Set the `SU2_HOME` to your SU2 source directory, the `SU2_RUN` to a directory where SU2 will install executable files, and add `SU2_RUN` to your `PATH`, and `PYTHONPATH` variables. For example, SU2 will show this message:

```text
Based on the input to this configuration, add these lines to your .bashrc file:

export SU2_RUN=/home/myuser/software/SU2_RUN/bin
export SU2_HOME=/home/myuser/software/SU2-7.5.1
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$PYTHONPATH:$SU2_RUN

Use './ninja -C build install' to compile and install SU2
```

which means you should set:

```shell
export SU2_RUN="/home/myuser/software/SU2_RUN"
export SU2_HOME="/home/myuser/software/SU2-7.5.1"
export PATH="${SU2_RUN}/bin:/path/to/su2-adapter/run/:$PATH"
export PYTHONPATH="${SU2_RUN}/bin:${PYTHONPATH}"
```

In particular, make sure that `SU2_RUN` points to a directory into which you have write access. You later will need to pass this to the SU2 build system (meson) with `--prefix`.

To copy the adapter files into `SU2_HOME`, run:

```shell
./su2AdapterInstall
```

The installation script will prompt you to follow commands. Check the output and follow the commands. For example:

```shell
./meson.py build -Denable-pywrapper=true --prefix=$SU2_RUN
./ninja -C build install
```

This will trigger the normal building procedure of SU2. Please refer to the installation instructions of SU2 for more details. SU2 should be built with MPI support in order to make use of parallel functionalities and must be built with pywrapper functionality enabled.

To be able to run the FSI and CHT Python scripts included in the adapter from anywhere, add to your `~/.bashrc`:

```shell
export PATH=/path/to/adapter/run:$PATH
```
