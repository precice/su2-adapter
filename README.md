
# Building and Usage of the preCICE Adapter for SU2

The adapter for the CFD-Code "Stanford University Unstructured" (SU2) was developed in the scope of the [bachelor's thesis of Alexander Rusch](https://www5.in.tum.de/pub/Rusch2016_BA.pdf).
All steps for integrating the adapter into SU2 are described in detail in the appendices of the thesis. Note that by now, the adapter has been extended with new functionalities, which are not covered in the thesis. However, the basic structure of the adapter has remained unchanged and can be studied by means of this work.
This adapter was developed for SU2 version 4.0.2 "Cardinal". Other versions may be compatible, yet they have not been tested. Please let us know if you want to use a different version.

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
Before installing the adapter SU2 itself must be downloaded from https://github.com/su2code/SU2. If necessary unpack the code and move it to your preferred location. Please do not configure and build the package before installing the adapter. In case you have already used SU2 you will need to rebuild the suite after installing the adapter.

### preCICE
It is assumed that preCICE has been installed successfully beforehand. Concerning installation instructions for preCICE, have a look at the preCICE-wiki pages on GitHub: https://github.com/precice/precice/wiki/Building.

### Adapter
In order to run SU2 with the preCICE adapter, some SU2-native solver routines need to be changed. The altered SU2 files are provided with this adapter in the directory "replacement_files". Moreover, preCICE-specific files are contained in the directory "adapter_files". These need to be added to the source code of SU2. A simple shell script called *su2AdapterInstall* comes with the adapter, which automates this process and replaces/copies the adapted and preCICE-specific files to the correct locations within the SU2 package including the appropriately adjusted *Makefile* of SU2. For the script to work correctly, two environment variables need to be set prior to executing it:

1. `PRECICE_ROOT`: Location of preCICE (top-level directory)
2. `SU2_HOME`: Location of SU2 (top-level directory)

It is recommended to set these variables permanently in your ~/.bashrc (Linux) or ~/.bash_profile (Mac). After setting these variables the script *su2AdapterInstall* can be run from the directory, in which it is contained:

```
    ./su2AdapterInstall
```

The script will not execute if the environment variables are unset or empty.

If you do not want to use this script, manually copy the files to the locations given in it. The two environment variables need to be defined as stated above, nevertheless.

Note that the automake/autoconf tools are required to install the adapter. If they are not available, the adapter installation script aborts and points out the requested tools.

After copying the adapter files to the correct locations within the SU2 package, SU2 can be configured and built just like the original version of the solver suite. Please refer to the installation instructions provided with the SU2 source code. SU2 should be built with MPI support in order to make use of parallel functionalities. The script *su2AdapterInstall* states recommended command sequences for both the configuration and the building process upon completion of the execution. Depending on your compiler, you might need to add `CXXFLAGS="-std=c++11"` in the SU2 configuration step as preCICE uses C++-11 features.

The SU2 executable is linked against the dynamic library of preCICE. In order for the executable to find the library at runtime, the dynamic linker search path needs to be adapted manually. The easiest way is to copy the following line into ~/.bashrc (Linux).
```
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PRECICE_ROOT/build/last"
```
For Mac-users the respective line, which needs to be added to ~/.bash_profile, is:
```
    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$PRECICE_ROOT/build/last"
```
Alternatively, the correct dynamic linker search path can be defined for the SU2 executable only (and not in ~/.bashrc or ~/.bash_profile).

## Running Simulations

### SU2 Configuration File
The adapter is turned on or off via the native SU2 configuration file. If it is turned off, SU2 executes in its original version. Moreover, the adapter is configured in the SU2 configuration file. The following adapter-related options are currently available (default values given in brackets):

1. `PRECICE_USAGE` (NO): Determines whether a preCICE-coupled simulation is run or not.
2. `PRECICE_VERBOSITYLEVEL_HIGH` (NO): Produces more output, mainly for debugging purposes.
3. `PRECICE_LOADRAMPING` (NO): Allows to linearly ramp up the load on the structural component at the beginning of the simulation. This may resolve stability issues due to large loads at the beginning of simulations.
4. `PRECICE_CONFIG_FILENAME` (precice.xml): Location and name of the preCICE configuration file.
5. `PRECICE_WETSURFACE_MARKER_NAME` (wetSurface): Name of the marker, which identifies the FSI interface in the geometry file of SU2. In the geometry file the name must be appended by a zero. If multiple interfaces exist, the marker names (in the geometry file) need to be alike differing only by the appending number, which must be successively increasing from zero. E.g. for three interfaces, the marker name could be defined as `PRECICE_WETSURFACE_MARKER_NAME= wetSurface` in the SU2 configuration file, while the markers in the geometry file would need to be named *wetSurface0*, *wetSurface1* and *wetSurface2*.
6. `PRECICE_LOADRAMPING_DURATION` (1): Number of time steps, in which the load ramping is active, counting from the beginning of the simulation. The ramped load increases linearly with each time step.
7. `PRECICE_NUMBER_OF_WETSURFACES` (1): In case multiple FSI-interfaces exist, their count needs to be specified here.

Moreover, in the SU2 configuration file grid movement must be allowed: `GRID_MOVEMENT= YES` and the type of grid movement must be set correctly for preCICE-coupled simulations: `GRID_MOVEMENT_KIND= PRECICE_MOVEMENT`. Also, the boundary, which is allowed to move, needs to be specified. Here the name of the FSI-interface marker including its appending identifying number as stated above must be used, e.g., `MARKER_MOVING= ( wetSurface0 )`. If multiple FSI-interfaces exist (as in the example above), this may look like `MARKER_MOVING= ( wetSurface0, wetSurface1, wetSurface2 )`.

### Running the Adapted SU2 Executable
Since the adapter (as well as its options) is turned on or off via the SU2 configuration file, the execution procedure is just as for the original version of SU2. For execution with one process working on the fluid domain from the directory, in which both the `SU2_CFD` executable and the SU2 configuration file are located:

```
    ./SU2_CFD su2ConfigurationFile.cfg
```

The adapter is designed such that it can be executed in an intra-parallel manner meaning that the flow domain is decomposed into several parts. The execution is then as follows (again assuming that executable and configuration file are within the current directory; exemplifying a decomposition of the fluid domain with eight processes):

```
    mpirun -n 8 ./SU2_CFD su2ConfigurationFile.cfg
```

A testcase is provided with this adapter in the directory `testcase`. 

## References
[1] Alexander Rusch. Extending SU2 to fluid-structure interaction via preCICE. Bachelor's thesis, Munich School of Engineering, Technical University of Munich, 2016.
