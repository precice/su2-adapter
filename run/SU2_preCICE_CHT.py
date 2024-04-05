#!/usr/bin/env python3

## \file SU2_preCICE_CHT.py
#  \brief Python script to launch SU2_CFD with preCICE for conjugate heat transfer simulations.
#  \author Joseph Signorelli
#  \version 7.5.1 "Blackbird"
#
# Part of the SU2-preCICE adapter: https://github.com/precice/su2-adapter
# This adapter is distributed under the GNU Lesser General Public
# License (see LICENSE file).

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys
from optparse import OptionParser # use a parser for configuration
import pysu2                      # imports the SU2 wrapped module
from math import *
import precice
import numpy
from time import sleep
# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

  # Command line options
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  # preCICE options with default settings
  parser.add_option("-p", "--precice-participant", dest="precice_name", help="Specify preCICE participant name", default="Fluid" )
  parser.add_option("-c", "--precice-config", dest="precice_config", help="Specify preCICE config file", default="../precice-config.xml")
  parser.add_option("-m", "--precice-mesh", dest="precice_mesh", help="Specify the preCICE mesh name", default="Fluid-Mesh")
  parser.add_option("-r", "--precice-reverse", action="store_true", dest="precice_reverse", help="Include flag to have SU2 write temperature, read heat flux", default=False)
  
  # Dimension
  parser.add_option("-d", "--dimension", dest="nDim", help="Dimension of fluid domain (2D/3D)", type="int", default=2)
  
  (options, args) = parser.parse_args()
  options.nZone = int(1) # Specify number of zones here (1)

  # Import mpi4py for parallel run
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
  else:
    comm = 0
    rank = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
      SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return


  # Configure preCICE:
  size = comm.Get_size()
  try:
    participant = precice.Participant(options.precice_name, options.precice_config, rank, size)
  except:
    print("There was an error configuring preCICE")
    return
  
  mesh_name = options.precice_mesh

  # Check preCICE + SU2 dimensions
  if options.nDim != participant.get_mesh_dimensions(mesh_name):
    print("SU2 and preCICE dimensions are not the same! Exiting")
    return

  CHTMarkerID = None
  CHTMarker = 'interface' # Name of CHT marker to couple

  # Get all the tags with the CHT option
  CHTMarkerList =  SU2Driver.GetAllCHTMarkersTag()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetAllBoundaryMarkers() # Returns all markers defined on this rank

  #Check if the specified marker has a CHT option and if it exists on this rank.
  if CHTMarker in CHTMarkerList and CHTMarker in allMarkerIDs.keys():
    CHTMarkerID = allMarkerIDs[CHTMarker] # So: if CHTMarkerID != None, then it exists on this rank
  
  # Number of vertices on the specified marker (per rank)
  nVertex_CHTMarker = 0         #total number of vertices (physical + halo) on this rank
  nVertex_CHTMarker_HALO = 0    #number of halo vertices
  nVertex_CHTMarker_PHYS = 0    #number of physical vertices
  iVertices_CHTMarker_PHYS = [] #indices of vertices this rank is working on
  # Note: Datatypes must be primitive as input to SU2 wrapper code, not numpy.int8, numpy.int64, etc.. So a list is used

  # If the CHT marker is defined on this rank:
  if CHTMarkerID != None:
    nVertex_CHTMarker = SU2Driver.GetNumberVertices(CHTMarkerID) #Total number of vertices on the marker
    nVertex_CHTMarker_HALO = SU2Driver.GetNumberHaloVertices(CHTMarkerID)
    nVertex_CHTMarker_PHYS = nVertex_CHTMarker - nVertex_CHTMarker_HALO

    # Obtain indices of all vertices that are being worked on on this rank
    for iVertex in range(nVertex_CHTMarker):
      if not SU2Driver.IsAHaloNode(CHTMarkerID, iVertex):
        iVertices_CHTMarker_PHYS.append(int(iVertex))

  # Get coords of vertices
  coords = numpy.zeros((nVertex_CHTMarker_PHYS, options.nDim))
  for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
    coord_passive = SU2Driver.GetInitialMeshCoord(CHTMarkerID, iVertex)
    for iDim in range(options.nDim):
      coords[i, iDim] = coord_passive[iDim]

  # Set mesh vertices in preCICE:
  try:
    vertex_ids = participant.set_mesh_vertices(mesh_name, coords)
  except:
    print("Could not set mesh vertices for preCICE. Was a (known) mesh specified in the options?")
    return

  # Get read and write data IDs
  precice_read = "Temperature"
  precice_write = "Heat-Flux"
  GetFxn = SU2Driver.GetVertexNormalHeatFlux
  SetFxn = SU2Driver.SetVertexTemperature
  GetInitialFxn = SU2Driver.GetVertexTemperature
  # Reverse coupling data read/write if -r flag included
  if options.precice_reverse:
    precice_read = "Heat-Flux"
    precice_write = "Temperature"
    GetFxn = SU2Driver.GetVertexTemperature
    SetFxn = SU2Driver.SetVertexNormalHeatFlux
    GetInitialFxn = SU2Driver.GetVertexNormalHeatFlux

  # Instantiate arrays to hold temperature + heat flux info
  read_data = numpy.zeros(nVertex_CHTMarker_PHYS)
  write_data = numpy.zeros(nVertex_CHTMarker_PHYS)

  # Retrieve some control parameters from the driver
  deltaT = SU2Driver.GetUnsteady_TimeStep()
  TimeIter = SU2Driver.GetTime_Iter()
  nTimeIter = SU2Driver.GetnTimeIter()
  time = TimeIter*deltaT

  # Set up initial data for preCICE
  if (participant.requires_initial_data()):

    for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
      write_data[i] = GetInitialFxn(CHTMarkerID, iVertex)

    participant.write_data(mesh_name, precice_write, vertex_ids, write_data)

  # Initialize preCICE
  participant.initialize()

  # Sleep briefly to allow for data initialization to be processed
  # This should only be needed on some systems and use cases
  #
  sleep(3)

  # Time loop is defined in Python so that we have access to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  precice_saved_time = 0
  precice_saved_iter = 0
  while (participant.is_coupling_ongoing()):
    # Implicit coupling
    if (participant.requires_writing_checkpoint()):
      # Save the state
      SU2Driver.SaveOldState()
      precice_saved_time = time
      precice_saved_iter = TimeIter

    # Get the maximum time step size allowed by preCICE
    precice_deltaT = participant.get_max_time_step_size()

    # Retrieve data from preCICE
    read_data = participant.read_data(mesh_name, precice_read, vertex_ids, deltaT) 

    # Set the updated values
    for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
        SetFxn(CHTMarkerID, iVertex, read_data[i])

    # Tell the SU2 drive to update the boundary conditions
    SU2Driver.BoundaryConditionsUpdate()

    if options.with_MPI == True:
      comm.Barrier()

    # Update timestep based on preCICE
    deltaT = SU2Driver.GetUnsteady_TimeStep()
    deltaT = min(precice_deltaT, deltaT)
    SU2Driver.SetUnsteady_TimeStep(deltaT)

    # Time iteration preprocessing
    SU2Driver.Preprocess(TimeIter)

    # Run one time iteration (e.g. dual-time)
    SU2Driver.Run()

    # Postprocess the solver and exit cleanly
    SU2Driver.Postprocess()

    # Update the solver for the next time iteration
    SU2Driver.Update()
    
    # Monitor the solver
    stopCalc = SU2Driver.Monitor(TimeIter)

    # Update control parameters
    TimeIter += 1
    time += deltaT

    # Loop over the vertices
    for i, iVertex in enumerate(iVertices_CHTMarker_PHYS):
      # Get heat fluxes at each vertex
      write_data[i] = GetFxn(CHTMarkerID, iVertex)
      
    # Write data to preCICE
    participant.write_data(mesh_name, precice_write, vertex_ids, write_data)

    # Advance preCICE
    participant.advance(deltaT)

    # Implicit coupling:
    if (participant.requires_reading_checkpoint()):
      # Reload old state
      SU2Driver.ReloadOldState()
      time = precice_saved_time
      TimeIter = precice_saved_iter

    if (participant.is_time_window_complete()):
      SU2Driver.Output(TimeIter)
      if (stopCalc == True):
        break

    if options.with_MPI == True:
      comm.Barrier()
      
  # Postprocess the solver and exit cleanly
  SU2Driver.Postprocessing()
  
  participant.finalize()
  
  if SU2Driver != None:
    del SU2Driver

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
