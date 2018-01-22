# Running the Fluid-Structure Interaction Testcase with SU2

1. Run the start-solution configuration file, which computes a fluid-only solution for 5s of physical time with a time step size of 0.01 s:
```
    ./SU2_CFD su2_config_Flap_startSolution.cfg
```
or in parallel:
```
    mpirun -n 4 ./SU2_CFD su2_config_Flap_startSolution.cfg
```
2. Next, run the FSI-solution configuration file, which runs the SU2 adapter for preCICE for another 5s of physical time with a time step size of 0.01s:
```
    ./SU2_CFD su2_config_Flap_FSISolution.cfg
```
or in parallel:
```
    mpirun -n 4 ./SU2_CFD su2_config_Flap_FSISolution.cfg
```

There are two differently fine meshes available (3DFlapGeoCoarse.su2 and 3DFlapGeoFine.su2). By default the coarse mesh is loaded. For using the finer mesh, uncomment the respective lines in both configuration files and comment the corresponding lines specifiying the coarse mesh. However, note that both meshes are quite coarse and are just intended to yield solutions very quickly. They are not designed for obtaining physically accurate solutions.
