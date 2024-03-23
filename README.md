# SU2 Python Wrapper Adapter

The SU2 adapter now leverages the SU2 Python wrapper and preCICE Python bindings to couple SU2 using preCICE in a minimally invasive way. This adapter simply updates existing functions and implements new ones in the SU2 Python wrapper that allow for simple preCICE implementation with implicit coupling. In addition, Python scripts are provided that can be used to quickly run FSI or CHT, or can serve as templates for your own custom application. This adapter currently works for SU2 versions 7.5.0 and 7.5.1 "Blackbird". Both conjugate heat-transfer and fluid-structure interaction can be completed with this adapter.

This adapter has been designed to work when using the compressible solver for unsteady problems with dual-time stepping, for single-zone problems. Implicit coupling currently saves the flow solution, turbulence solution, and the mesh solution (for mesh deformation). Species transport and transition model variables at this time are not saved, but may be straightforward to implement.

## Contents
<!-- toc orderedList:0 -->

- [Contents](#contents)
- [Building the Adapter](#building-the-adapter)
  - [SU2](#su2)
  - [preCICE](#precice)
  - [Adapter](#adapter)
- [Running Simulations](#running-simulations)
  - [Fluid-Structure Interaction](#fluid-structure-interaction)
  - [Conjugate Heat Transfer](#conjugate-heat-transfer)
  - [Running in parallel](#running-in-parallel)

<!-- tocstop -->
## Building the adapter

### SU2

Download SU2 v7.5.1 "Blackbird" from directly from https://github.com/su2code/SU2/releases/tag/v7.5.1. Note that both swig and mpi4py must be installed to use the SU2 Python wrapper. After installing the adapter, the flag `-Denable-pywrapper=true` must be specified.

### preCICE

In addition to having successfully installed preCICE, the preCICE Python bindings must also be installed. For installing preCICE, please navigate to https://precice.org/installation-overview.html. After successfully installing preCICE, please follow the instructions at https://precice.org/installation-bindings-python.html to get the preCICE Python bindings. As a test, run the following command:

```shell
python3 -c "import precice"
```

If there are no errors, then preCICE and its Python bindings were successfully installed.

### Adapter

In order to couple SU2 using preCICE, *python_wrapper_structure.cpp* and *CDriver.hpp* must be updated. This adapter provides the updated files. The shell script *su2AdapterInstall*, which comes with this adapter, automatically replaces the files in your SU2 directory with these updated files and provides the correct commands to re-configure and re-install SU2 with the added adjustments. For this to work, the `SU2_HOME` variable must be set to your SU2 directory prior to running. Set also the `SU2_RUN`, `PATH`, and `PYTHONPATH` variables accordingly. For example, SU2 will show this message:

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
export SU2_HOME="/home/myuser/software/SU2-7.5.1"
export SU2_RUN="/home/myuser/software/SU2_RUN"
export PATH="${SU2_RUN}/bin:/path/to/su2-adapter/run/:$PATH"
export PYTHONPATH="${SU2_RUN}/bin:${PYTHONPATH}"
```

In particular, make sure that `SU2_RUN` points to a directory into which you have write access. You later will need to pass this to meson with `--prefix`.

SU2 will advise you to add this variable (and others) to your `~/.bashrc` (Linux) or `~/.bash_profile` (Mac) after configuring, so it may already be set if SU2 is already configured and installed on your computer. To install the adapter, run from the adapter directory:

```shell
./su2AdapterInstall
```

The script will not execute if the environment variable `SU2_HOME` is not set or empty.

If you do not want to use this script, manually copy the files to the locations given in it. The environment variable needs to be defined as stated above, nevertheless.

The installation script will prompt you to follow commands. Check the output and follow the commands. For example:

```shell
./meson.py build -Denable-pywrapper=true --prefix=$SU2_RUN
./ninja -C build install
```

After copying the adapter files to the correct locations within the SU2 package, SU2 can be configured and built just like the original version of the solver suite. Please refer to the installation instructions provided with the SU2 source code. SU2 should be built with MPI support in order to make use of parallel functionalities and must be built with pywrapper functionality enabled. The script *su2AdapterInstall* states recommended command sequences for both the configuration and the building process upon completion of the execution.

To utilize the default FSI and CHT scripts anywhere, add to your ~/.bashrc:

```shell
export PATH=/path/to/adapter/run:$PATH
```

## Running simulations

After successfully installing the adapted SU2, the default FSI/CHT scripts may be utilized. Note that these scripts currently are designed for a single coupling mesh, called *interface*. However it is extremely easy to update these scripts to handle a different BC name and/or multiple interfaces. They are provided simply for their ease of use.

### Fluid-structure interaction

#### SU2 configuration file for FSI

To set up a single-interface FSI problem for coupling with preCICE, the SU2 config file should have the following:

```text
DEFORM_MESH= YES
MARKER_DEFORM_MESH= ( interface )
```

The *interface* marker should also be set as a wall boundary, such as `MARKER_EULER` or `MARKER_ISOTHERMAL`.

#### Running an FSI simulation

By default in the *SU2_preCICE_FSI.py* script, the following settings are automatically used for coupling with preCICE:

- preCICE Participant Name: *Fluid*
- preCICE Config File: *../precice-config.xml*
- preCICE Mesh Name: *Fluid-Mesh*
- preCICE Read Data: *Displacements*
- preCICE Write Data: *Forces*

To run with these settings:

```shell
SU2_preCICE_FSI.py -f SU2_config_file.cfg --parallel
```

The `--parallel` flag must **always** be used when SU2 is built in parallel, even if running on a single process. If you do not build SU2 with MPI, do not include it.

The read/write data are hardcoded, but the participant name, config file, and mesh name can be changed using flags in the call to the Python file. In general, to run an FSI case:

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

By default in the *SU2_preCICE_CHT.py* script, the following settings are automatically used for coupling with preCICE:

- preCICE Participant Name: *Fluid*
- preCICE Config File: *../precice-config.xml*
- preCICE Mesh Name: *Fluid-Mesh*
- preCICE Read Data: *Temperature*
- preCICE Write Data: *Heat-Flux*

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

**NOTE**: Result files (vtu) might be incompatible with your ParaView version. For example, ParaView 5.11.2 on Ubuntu 22.04 is known to fail, but ParaView 5.12 works.
