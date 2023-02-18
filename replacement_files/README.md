# Implicit Coupling Setup:

Method implemented for implicit coupling was to essentially treat states as reading + writing restart files
There is no non-invasive way to write preCICE saved state restart files, and I/O may be too much.

### <ins>Main Solution Variables (member variables set within each class)</ins>

`FLOW_SOL`:
- Solution is set in CFVMFlowSolverBase::LoadRestart_impl
- Solution_time_n and Solution_time_n1 are set in CFVMFlowSolverBase::PushSolutionBackInTime


`TURB_SOL`:
- Solution is set in CTurbSolver::LoadRestart
- Solution_time_n1 and Solution_time_n1 are set in CFVMFlowSolverBase::PushSolutionBackInTime


`MESH_SOL`:
- Solution, Solution_time_n, and Solution_time_n1 is set in CMeshSolver::LoadRestart (some from function within)


### <ins>Geometry Container Variables Set</ins>

*These only matter if there is grid deformation
`[MESH_0]`:
- Coord and GridVel are set in CFVMFlowSolverBase::LoadRestart_impl
- GridVel is also set in CMeshSolver::LoadRestart (function within) -- but note that it is manually re-calculated


`[iMesh]`:
- Volume_n, Volume_nM1 are set in CFVMFlowSolverBase::PushSolutionBackInTime
- (NOTE: Coord_n and Coord_n1 also set if CConfig::GetGrid_Movement(), but not true for us!- verified)



After all variables are set, remaining communications/multigrid-interpolations/calculations were copied and pasted into appropriate functions.


Method of data saving was verified by outputting RESTART_ASCII files after saving a state and after reloading a state, and both successfully are identical when this is tested.