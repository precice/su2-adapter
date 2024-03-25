#!/usr/bin/env python3

## File edited to work coupled with preCICE - Joseph Signorelli

## \file flatPlate_rigidMotion.py
#  \brief Python script to launch SU2_CFD with customized unsteady boundary conditions using the Python wrapper.
#  \author David Thomas
#  \version 7.5.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys
from optparse import OptionParser	# use a parser for configuration
import pysu2			            # imports the SU2 wrapped module
from math import *
import numpy
import precice
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

    (options, args) = parser.parse_args()
    options.nDim  = int(2)
    options.nZone = int(1)

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
        participant = precice.Participant(options.precice_name, options.precice_config, rank, size)#, comm)
    except:
        print("There was an error configuring preCICE")
        return

    mesh_name = options.precice_mesh

    # Check preCICE + SU2 dimensions
    if options.nDim != participant.get_mesh_dimensions(mesh_name):
        print("SU2 and preCICE dimensions are not the same! Exiting")
        return

    MovingMarkerID = None
    MovingMarker = 'interface'       #specified by the user

    # Get all the tags with the moving option
    MovingMarkerList =  SU2Driver.GetAllDeformMeshMarkersTag()

    # Get all the markers defined on this rank and their associated indices.
    allMarkerIDs = SU2Driver.GetAllBoundaryMarkers()

    # Check if the specified marker has a moving option and if it exists on this rank.
    if MovingMarker in MovingMarkerList and MovingMarker in allMarkerIDs.keys():
        MovingMarkerID = allMarkerIDs[MovingMarker]

    # Number of vertices on the specified marker (per rank)
    nVertex_MovingMarker = 0         #total number of vertices (physical + halo)
    nVertex_MovingMarker_HALO = 0    #number of halo vertices
    nVertex_MovingMarker_PHYS = 0    #number of physical vertices
    iVertices_MovingMarker_PHYS = [] # indices of vertices this rank is working on
    # Datatypes must be primitive as input to SU2 wrapper code, not numpy.int8, numpy.int64, etc.. So a list is used
    
    if MovingMarkerID != None:
        nVertex_MovingMarker = SU2Driver.GetNumberVertices(MovingMarkerID)
        nVertex_MovingMarker_HALO = SU2Driver.GetNumberHaloVertices(MovingMarkerID)
        nVertex_MovingMarker_PHYS = nVertex_MovingMarker - nVertex_MovingMarker_HALO
        
        # Obtain indices of all vertices that are being worked on on this rank
        for iVertex in range(nVertex_MovingMarker):
            if not SU2Driver.IsAHaloNode(MovingMarkerID, iVertex):
                iVertices_MovingMarker_PHYS.append(int(iVertex))
    
    # Get coords of vertices
    coords = numpy.zeros((nVertex_MovingMarker_PHYS, options.nDim))
    for i, iVertex in enumerate(iVertices_MovingMarker_PHYS):
        coord_passive = SU2Driver.GetInitialMeshCoord(MovingMarkerID, iVertex)
        for iDim in range(options.nDim):
            coords[i, iDim] = coord_passive[iDim]

    # Set mesh vertices in preCICE:
    vertex_ids = participant.set_mesh_vertices(mesh_name, coords)

    # Set mesh vertices in preCICE:
    try:
        vertex_ids = participant.set_mesh_vertices(mesh_name, coords)
    except:
        print("Could not set mesh vertices for preCICE. Was a (known) mesh specified in the options?")
        return

    # Get read and write data IDs
    # By default:
    precice_read = "Displacement"
    precice_write = "Force"

    # Instantiate arrays to hold displacements + forces info
    displacements = numpy.zeros((nVertex_MovingMarker_PHYS,options.nDim))
    forces = numpy.zeros((nVertex_MovingMarker_PHYS,options.nDim))

    # Retrieve some control parameters from the driver
    deltaT = SU2Driver.GetUnsteady_TimeStep()
    TimeIter = SU2Driver.GetTime_Iter()
    nTimeIter = SU2Driver.GetnTimeIter()
    time = TimeIter*deltaT

    # Set up initial data for preCICE
    if (participant.requires_initial_data()):

        for i, iVertex in enumerate(iVertices_MovingMarker_PHYS):
            forces[i] = SU2Driver.GetFlowLoad(MovingMarkerID, iVertex)[:-1]

        participant.write_block_vector_data(mesh_name, precice_write, vertex_ids, forces)

    # Initialize preCICE
    participant.initialize()

    # Setup time step sizes:
    precice_deltaT = participant.get_max_time_step_size()
    deltaT = SU2Driver.GetUnsteady_TimeStep()
    deltaT = min(precice_deltaT, deltaT)

    # Sleep briefly to allow for data initialization to be processed
    # This should not be needed, but open an issue if it is:
    # https://github.com/precice/su2-adapter/issues
    #
    # sleep(3)

    # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
    if rank == 0:
        print("\n------------------------------ Begin Solver -----------------------------\n")
    sys.stdout.flush()
    if options.with_MPI == True:
        comm.Barrier()

    while (participant.is_coupling_ongoing()):#(TimeIter < nTimeIter):
        
        # Implicit coupling
        if (participant.requires_writing_checkpoint()):
            # Save the state
            SU2Driver.SaveOldState()

        # Retreive data from preCICE
        displacements = participant.read_data(mesh_name, precice_read, vertex_ids, deltaT)
        
        # Set the updated displacements
        for i, iVertex in enumerate(iVertices_MovingMarker_PHYS):
            DisplX = displacements[i][0]
            DisplY = displacements[i][1]
            DisplZ = 0 if options.nDim == 2 else displacements[i][2]

            SU2Driver.SetMeshDisplacement(MovingMarkerID, iVertex, DisplX, DisplY, DisplZ)
        
        if options.with_MPI == True:
            comm.Barrier()
            
        # Update timestep based on preCICE
        deltaT = SU2Driver.GetUnsteady_TimeStep()
        deltaT = min(precice_deltaT, deltaT)
        SU2Driver.SetUnsteady_TimeStep(deltaT)
        
        # Time iteration preprocessing (mesh is deformed here)
        SU2Driver.Preprocess(TimeIter)

        # Run one time iteration (e.g. dual-time)
        SU2Driver.Run()

        # Postprocess the solver
        SU2Driver.Postprocess()

        # Update the solver for the next time iteration
        SU2Driver.Update()

        # Monitor the solver and output solution to file if required
        stopCalc = SU2Driver.Monitor(TimeIter)


        # Loop over the vertices
        for i, iVertex in enumerate(iVertices_MovingMarker_PHYS):
            # Get forces at each vertex
            forces[i] = SU2Driver.GetFlowLoad(MovingMarkerID, iVertex)[:-1]

        # Write data to preCICE
        participant.write_data(mesh_name, precice_write, vertex_ids, forces)

        # Advance preCICE
        participant.advance(deltaT)
        precice_deltaT = participant.get_max_time_step_size()


        # Implicit coupling:
        if (participant.requires_reading_checkpoint()):
            # Reload old state
            SU2Driver.ReloadOldState()
        else: # Output and increment as usual
            SU2Driver.Output(TimeIter)
            if (stopCalc == True):
                break
            # Update control parameters
            TimeIter += 1
            time += deltaT
        
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
