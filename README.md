<a style="text-decoration: none" href="https://travis-ci.org/precice/su2-adapter" target="_blank">
    <img src="https://travis-ci.org/precice/su2-adapter.svg?branch=master" alt="Build status">
</a>

# Building and Usage of the preCICE Adapter for SU2

The adapter for the CFD-Code "Stanford University Unstructured" (SU2) was developed in the scope of the [bachelor's thesis of Alexander Rusch](https://www5.in.tum.de/pub/Rusch2016_BA.pdf).
All steps for integrating the adapter into SU2 are described in detail in the appendices of the thesis. Note that by now, the adapter has been extended with new functionalities, which are not covered in the thesis. However, the basic structure of the adapter has remained unchanged and can be studied by means of this work.
This adapter was developed for SU2 version 6.0.0 "Falcon". Other versions may be compatible as well, yet they have not been tested. Please let us know if you want to use a different version.

## Contents
<!-- toc orderedList:0 -->

- [Building and Usage of the preCICE Adapter for SU2](#building-and-usage-of-the-preCICE-adapter-for-su2)
    - [Contents](#contents)
    - [Building the Adapter](#building-the-adapter)
        - [SU2](#su2)
        - [preCICE](#precice)
        - [Adapter](#adapter)
    - [Running Simulations](#running-simulations)
        - [SU2 Configuration File](#su2-configuration-file)
        - [Running the Adapted SU2 Executable](#running-the-adapted-su2-executable)
    - [References](#references)

<!-- tocstop -->

## Building the Adapter

### SU2
Before installing the adapter SU2 itself must be installed. Download version 6.0.0 directly from https://github.com/su2code/SU2/releases/tag/v6.0.0. Using a different version is not recommended, since the adapter is only tested with this version. If necessary unpack the code and move it to your preferred location. Please do not configure and build the package before installing the adapter. In case you have already used SU2 you will need to rebuild the suite after installing the adapter.

**Troubleshooting**

* If you run `./configure --prefix=$SU2_HOME` and get the error `configure: error: cannot find python-config for /usr/bin/python`, check via `ls /usr/bin` whether there is a `python-config` and/or `python2.7-config`. If not, you can create a symbolic link via `ln /usr/bin/python3-config /usr/bin/python-config` such that `python-config` uses `python3-config`.

### preCICE
It is assumed that preCICE has been installed successfully beforehand. Concerning installation instructions for preCICE, have a look at the preCICE-wiki pages on GitHub: https://github.com/precice/precice/wiki/Building.

### Adapter
In order to run SU2 with the preCICE adapter, some SU2-native solver routines need to be changed. The altered SU2 files are provided with this adapter in the directory "replacement_files". Moreover, preCICE-specific files are contained in the directory "adapter_files". These need to be added to the source code of SU2. A simple shell script called *su2AdapterInstall* comes with the adapter, which automates this process and replaces/copies the adapted and preCICE-specific files to the correct locations within the SU2 package including the appropriately adjusted *Makefile* of SU2. For the script to work correctly, the environment variable `SU2_HOME` needs to be set to the location of SU2 (top-level directory).

It is recommended to set this variable permanently in your ~/.bashrc (Linux) or ~/.bash_profile (Mac). After setting the variable the script *su2AdapterInstall* can be run from the directory, in which it is contained:

```
./su2AdapterInstall
```

The script will not execute if the environment variable is not set or empty.

If you do not want to use this script, manually copy the files to the locations given in it. The environment variable needs to be defined as stated above, nevertheless.

After copying the adapter files to the correct locations within the SU2 package, SU2 can be configured and built just like the original version of the solver suite. Please refer to the installation instructions provided with the SU2 source code. SU2 should be built with MPI support in order to make use of parallel functionalities. The script *su2AdapterInstall* states recommended command sequences for both the configuration and the building process upon completion of the execution.

The SU2 executable is linked against the **dynamic library** of preCICE, so make sure you have built it like this.

## Running Simulations

### SU2 Configuration File
The adapter is turned on or off via the native SU2 configuration file. If it is turned off, SU2 executes in its original version. Moreover, the adapter is configured in the SU2 configuration file. The following adapter-related options are currently available (default values given in brackets):

1. `PRECICE_USAGE` (NO): Determines whether a preCICE-coupled simulation is run or not.
2. `PRECICE_VERBOSITYLEVEL_HIGH` (NO): Produces more output, mainly for debugging purposes.
3. `PRECICE_LOADRAMPING` (NO): Allows to linearly ramp up the load on the structural component at the beginning of the simulation. This may resolve stability issues due to large loads at the beginning of simulations.
4. `PRECICE_CONFIG_FILENAME` (precice-config.xml): Location and name of the preCICE configuration file.
5. `PRECICE_PARTICIPANT_NAME` (SU2): Name of the participant in the preCICE configuration file.
6. `PRECICE_MESH_NAME` (SU2-Mesh): Name of the mesh in the preCICE configuration file.
7. `PRECICE_READ_DATA_NAME` (Displacements): Name of the read data in the preCICE configuration file.
8. `PRECICE_WRITE_DATA_NAME` (Forces): Name of the write data in the preCICE configuration file.
9. `PRECICE_WETSURFACE_MARKER_NAME` (wetSurface): Name of the marker, which identifies the FSI interface in the geometry file of SU2.
10. `PRECICE_LOADRAMPING_DURATION` (1): Number of time steps, in which the load ramping is active, counting from the beginning of the simulation. The ramped load increases linearly with each time step.
11. `PRECICE_NUMBER_OF_WETSURFACES` (1): In case multiple FSI-interfaces exist, their count needs to specified here.

If multiple interfaces exist, the names of all related entries (`PRECICE_WETSURFACE_MARKER_NAME`, `PRECICE_READ_DATA_NAME`,`PRECICE_WRITE_DATA_NAME`, `PRECICE_MESH_NAME`) must be appended by consecutive numbers. Hence, the names (also in the geometry file) need to be alike differing only by the appending number, which must be successively increasing from zero. E.g. for three interfaces, the marker name could be defined as `PRECICE_WETSURFACE_MARKER_NAME= wetSurface` in the SU2 configuration file, while the markers in the geometry file would need to be named *wetSurface*, *wetSurface1* and *wetSurface2*.

Moreover, in the SU2 configuration file grid movement must be allowed: `GRID_MOVEMENT= YES` and the type of grid movement must be set correctly for preCICE-coupled simulations: `GRID_MOVEMENT_KIND= PRECICE_MOVEMENT`. Also, the boundary, which is allowed to move needs to be specified. Here the name of the FSI-interface marker including its appending identifying number as stated above needs to be used, e.g., `MARKER_MOVING= ( wetSurface0 )`. If multiple FSI-interfaces exist (as in the example above), this may look like `MARKER_MOVING= ( wetSurface, wetSurface1, wetSurface2 )`.

### Running the Adapted SU2 Executable
Since the adapter (as well as its options) is turned on or off via the SU2 configuration file, the execution procedure is just as for the original version of SU2. For execution with one process working on the fluid domain from the directory, in which both the `SU2_CFD` executable and the SU2 configuration file are located:

```
    ./SU2_CFD su2ConfigurationFile.cfg
```

The adapter is designed such that it can be executed in an intra-parallel manner meaning that the flow domain is decomposed into several parts. The execution is then as follows (again assuming that executable and configuration file are within the current directory; exemplifying a decomposition of the fluid domain with eight processes):

```
    mpirun -n 8 ./SU2_CFD su2ConfigurationFile.cfg
```

## References
[1] Alexander Rusch. Extending SU2 to fluid-structure interaction via preCICE. Bachelor's thesis, Munich School of Engineering, Technical University of Munich, 2016.
