# SU2 adapter

The SU2 adapter now leverages the SU2 Python wrapper and preCICE Python bindings to couple SU2 using preCICE in a minimally invasive way. This adapter simply updates existing functions and implements new ones in the SU2 Python wrapper that allow for simple preCICE implementation with implicit coupling. In addition, Python scripts are provided that can be used to quickly run FSI or CHT, or can serve as templates for your own custom application. This adapter currently works for SU2 versions 7.5.0 and 7.5.1 "Blackbird". Both conjugate heat-transfer and fluid-structure interaction can be completed with this adapter.

This adapter has been designed to work when using the compressible solver for unsteady problems with dual-time stepping, for single-zone problems. Implicit coupling currently saves the flow solution, turbulence solution, and the mesh solution (for mesh deformation). Species transport and transition model variables at this time are not saved, but may be straightforward to implement.

**Note:** In its current state, the SU2 adapter is using the Python wrapper of SU2. The [previous implementation](https://github.com/precice/su2-adapter/tree/ab843878c1d43302a4f0c66e25dcb364b7787478) was directly editing the C++ source files of SU2.

<!-- tocstop -->
## Building the adapter

The adapter depends on [SU2](https://su2code.github.io/), [preCICE v2](https://precice.org/installation-overview.html), and the [preCICE Python bindings](https://precice.org/installation-bindings-python.html).

The script `su2AdapterInstall` replaces a few files in the SU2 source code. You then need to build SU2 from source, install it into a prefix (`SU2_RUN`) and add that to your `PATH`.

To run SU2, you can use the provided Python scripts `SU2_preCICE_CHT.py` and `SU2_preCICE_FSI.py`, which call SU2 via its Python interface.

### Get the SU2 source

Download SU2 v7.5.1 "Blackbird" source directly from the [SU2 releases on GitHub](https://github.com/su2code/SU2/releases/tag/v7.5.1). Note that both swig and mpi4py must be installed to use the SU2 Python wrapper, which needs to be enabled using the flag `-Denable-pywrapper=true` when building SU2.

### Build the adapter

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

## Running simulations

After successfully installing the adapted SU2, use the provided FSI/CHT scripts to run simulations, or build your own. Note that these scripts currently are designed for a single coupling mesh, called `interface`. However, it is extremely easy to update these scripts to handle a different BC name and/or multiple interfaces. They are provided simply for their ease of use.

### Fluid-structure interaction

#### SU2 configuration file for FSI

To set up a single-interface FSI problem for coupling with preCICE, the SU2 config file should have the following:

```text
DEFORM_MESH= YES
MARKER_DEFORM_MESH= ( interface )
```

The `interface` marker should also be set as a wall boundary, such as `MARKER_EULER` or `MARKER_ISOTHERMAL`.

#### Running an FSI simulation

By default, in the `SU2_preCICE_FSI.py` script, the following settings are automatically used for coupling with preCICE:

- preCICE Participant Name: `Fluid`
- preCICE Config File: `../precice-config.xml`
- preCICE Mesh Name: `Fluid-Mesh`
- preCICE Read Data: `Displacements`
- preCICE Write Data: `Forces`

To run with these settings:

```shell
SU2_preCICE_FSI.py -f SU2_config_file.cfg --parallel
```

The `--parallel` flag must **always** be used when SU2 is built in parallel, even if running on a single process. If you do not build SU2 with MPI, do not include it.

The read/write data variables are hardcoded, but the participant name, config file, and mesh name can be changed using flags in the call to the Python file. In general, to run an FSI case:

```shell
SU2_preCICE_FSI.py -f SU2_config_file.cfg -p participant_name -c precice_config_file -m precice_mesh_name --parallel
```

### Conjugate heat transfer

#### SU2 configuration file for CHT

To set up a single-interface CHT problem for coupling with preCICE, the SU2 config file should have the following:

```text
% For having SU2 read temperature, write heat flux:
MARKER_ISOTHERMAL= (interface, ______)
%
% For having SU2 read heat flux, write temperature (the -r flag):
MARKER_HEATFLUX= (interface, ______)
%
% And in both cases include:
MARKER_PYTHON_CUSTOM= (interface)
```

Note that the blank spots in the isothermal and heat flux markers are the initial BC values. If there is a data initialization from another solver, they will be updated and are not important.

#### Running a CHT simulation

By default in the `SU2_preCICE_CHT.py` script, the following settings are automatically used for coupling with preCICE:

- preCICE Participant Name: `Fluid`
- preCICE Config File: `../precice-config.xml`
- preCICE Mesh Name: `Fluid-Mesh`
- preCICE Read Data: `Temperature`
- preCICE Write Data: `Heat-Flux`

To run with these settings:

```shell
SU2_preCICE_CHT.py -f SU2_config_file.cfg --parallel
```

The `--parallel` flag must **always** be used when SU2 is built in parallel, even if running on a single process. If you do not build SU2 with MPI, do not include it.

The read/write data for CHT can be reversed if the preCICE config file specifies for the fluid to read heat flux and write temperature. This can easily be accomplished with the `-r` flag:

```shell
SU2_preCICE_CHT.py -f SU2_config_file.cfg -r --parallel
```

The participant name, config file, and mesh name can be changed using flags in the call to the Python file. In general, to run a CHT case:

```shell
SU2_preCICE_CHT.py -f SU2_config_file.cfg -p participant_name -c precice_config_file -m precice_mesh_name --parallel
```

### Running in parallel

The Python scripts can very easily be run in parallel by just pre-pending the Python script call like:

```shell
mpirun -n 8 python3 SU2_preCICE_CHT.py -f SU2_config_file.cfg --parallel
```

**NOTE**: As of SU2 v7.5.1: Deforming `MARKER_EULER`'s are buggy when simulations are run in parallel, leading to unexpected results. More information can be found at this discussion here: https://github.com/su2code/SU2/discussions/1931.

### Further notes

Result files (vtu) generated from SU2 might be incompatible with your ParaView version. For example, ParaView 5.11.2 on Ubuntu 22.04 is known to fail with SU2 7.5.1 result files, but ParaView 5.12 works.

The replacement files included in this repository might be long, but they only introduce minimal changes compared to the original SU2 code (mainly related to checkpointing for implicit coupling). When updating to newer SU2 versions, compare the bundled and the old unmodified files in a diff tool, and start by copying the same changes into the new source. See the `replacement_files/README.md` for more details.
