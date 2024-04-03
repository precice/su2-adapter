---
title: Running simulations
permalink: adapter-su2-configure.html
keywords: adapter, su2, development, modules
summary: "Modify SU2 configuration file, specify interfaces by SU2 markers, run SU2 either serial or parallel"
---

After successfully installing the adapted SU2, use the provided FSI/CHT scripts to run simulations, or build your own. Note that these scripts currently are designed for a single coupling mesh, called `interface`. However, it is extremely easy to update these scripts to handle a different BC name and/or multiple interfaces. They are provided simply for their ease of use.

## Fluid-structure interaction

### SU2 configuration file for FSI

To set up a single-interface FSI problem for coupling with preCICE, the SU2 config file should have the following:

```text
DEFORM_MESH= YES
MARKER_DEFORM_MESH= ( interface )
```

The `interface` marker should also be set as a wall boundary, such as `MARKER_EULER` or `MARKER_ISOTHERMAL`.

### Running an FSI simulation

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

## Conjugate heat transfer

### SU2 configuration file for CHT

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

### Running a CHT simulation

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

## Running in parallel

The Python scripts can very easily be run in parallel by just pre-pending the Python script call like:

```shell
mpirun -n 8 python3 SU2_preCICE_CHT.py -f SU2_config_file.cfg --parallel
```

{% note %}
As of SU2 v7.5.1: Deforming `MARKER_EULER`'s are buggy when simulations are run in parallel, leading to unexpected results. More information can be found at [this SU2 discussion](https://github.com/su2code/SU2/discussions/1931).
{% endnote %}

## Important note on restarts

This code **has not been tested** for restarts using initializations *from* SU2. Any restarted simulations should have SU2 be the first participant and receive initialization data. It is possible that, if SU2 must send initialization data, that it is incorrect (it may use default values in the config file, or just be zeros if the data hasn't been computed until after/during a first iteration). Admittedly, this is from a lack of understanding of the specifics of how SU2 operates and there may not be a trivial work-around.

## Further notes

Result files (vtu) generated from SU2 might be incompatible with your ParaView version. For example, ParaView 5.11.2 on Ubuntu 22.04 is known to fail with SU2 7.5.1 result files, but ParaView 5.12 works.

The replacement files included in this repository might be long, but they only introduce minimal changes compared to the original SU2 code (mainly related to checkpointing for implicit coupling). When updating to newer SU2 versions, compare the bundled and the old unmodified files in a diff tool, and start by copying the same changes into the new source. See the `replacement_files/README.md` for more details.
