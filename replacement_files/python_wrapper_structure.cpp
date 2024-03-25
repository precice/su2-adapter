/*!
 *
 * Edited to include required capabilities for coupling with preCICE
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */


#include "../include/drivers/CDriver.hpp"
#include "../include/drivers/CSinglezoneDriver.hpp"
#include "../../Common/include/toolboxes/geometry_toolbox.hpp"

void CDriver::PythonInterface_Preprocessing(CConfig **config, CGeometry ****geometry, CSolver *****solver){

  int rank = MASTER_NODE;
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  /* --- Initialize boundary conditions customization, this is achieve through the Python wrapper --- */
  for(iZone=0; iZone < nZone; iZone++){

    if (config[iZone]->GetnMarker_PyCustom() > 0){

      if (rank == MASTER_NODE) cout << endl << "----------------- Python Interface Preprocessing ( Zone "<< iZone <<" ) -----------------" << endl;

      if (rank == MASTER_NODE) cout << "Setting customized boundary conditions for zone " << iZone << endl;
      for (iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        geometry[iZone][INST_0][iMesh]->SetCustomBoundary(config[iZone]);
        // preCICE: Update input file initial custom boundary temperature/HF to be nondimensional.
        //          Was fixed in SU2 v8 here: https://github.com/su2code/SU2/pull/2078
        //------------------------------------------------------------------------------------------------
        // Basically just copied + pasted (terrible coding but as a quick/least invasive fix) CGeometry::SetCustomBoundary
        unsigned short iMarker;
        unsigned long iVertex;
        string Marker_Tag;

        for(iMarker=0; iMarker < geometry[iZone][INST_0][iMesh]->GetnMarker(); iMarker++){
          Marker_Tag = config[iZone]->GetMarker_All_TagBound(iMarker);
          if(config[iZone]->GetMarker_All_PyCustom(iMarker)){
            switch(config[iZone]->GetMarker_All_KindBC(iMarker)){
              case HEAT_FLUX:
                for(iVertex=0; iVertex < geometry[iZone][INST_0][iMesh]->GetnVertex(iMarker); iVertex++){
                  geometry[iZone][INST_0][iMesh]->SetCustomBoundaryHeatFlux(iMarker, iVertex, config[iZone]->GetWall_HeatFlux(Marker_Tag)/config[iZone]->GetHeat_Flux_Ref());
                }
                break;
              case ISOTHERMAL:
                for(iVertex=0; iVertex < geometry[iZone][INST_0][iMesh]->GetnVertex(iMarker); iVertex++){
                  geometry[iZone][INST_0][iMesh]->SetCustomBoundaryTemperature(iMarker, iVertex, config[iZone]->GetIsothermal_Temperature(Marker_Tag)/config[iZone]->GetTemperature_Ref());
                }
                break;
              case INLET_FLOW:
                // This case is handled in the solver class.
                break;
              default:
                cout << "WARNING: Marker " << Marker_Tag << " is not customizable. Using default behavior." << endl;
                break;
            }
          }
        }
          //-------------------------------------------------------------------------------------------------
      }


      geometry[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);

      if ((config[iZone]->GetKind_Solver() == MAIN_SOLVER::EULER) ||
          (config[iZone]->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES) ||
          (config[iZone]->GetKind_Solver() == MAIN_SOLVER::RANS)) {

        solver[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the global performance indices (Lift, Drag, ecc..) */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::Get_Drag() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, factor, val_Drag;

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  val_Drag = CDrag*factor;

  return SU2_TYPE::GetValue(val_Drag);
}

passivedouble CDriver::Get_Lift() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, factor, val_Lift;

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  val_Lift = CLift*factor;

  return SU2_TYPE::GetValue(val_Lift);
}

passivedouble CDriver::Get_Mx() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMx, RefLengthCoeff, factor, val_Mx;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMx = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();

  val_Mx = CMx*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mx);
}

passivedouble CDriver::Get_My() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMy, RefLengthCoeff, factor, val_My;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMy = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();

  val_My = CMy*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_My);
}

passivedouble CDriver::Get_Mz() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMz, RefLengthCoeff, factor, val_Mz;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around z-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMz = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();

  val_Mz = CMz*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mz);
}

passivedouble CDriver::Get_DragCoeff() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag;

  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  return SU2_TYPE::GetValue(CDrag);
}

passivedouble CDriver::Get_LiftCoeff() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift;

  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  return SU2_TYPE::GetValue(CLift);
}

/////////////////////////////////////////////////////////////////////////////
/* Functions to obtain information from the geometry/mesh                  */
/////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberVertices(unsigned short iMarker) const {

  return geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker];

}

unsigned long CDriver::GetNumberHaloVertices(unsigned short iMarker) const {

  unsigned long nHaloVertices, iVertex, iPoint;

  nHaloVertices = 0;
  for(iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    if(!(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint))) nHaloVertices += 1;
  }

  return nHaloVertices;

}

unsigned long CDriver::GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint, GlobalIndex;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  GlobalIndex = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetGlobalIndex(iPoint);

  return GlobalIndex;

}

bool CDriver::IsAHaloNode(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint)) return false;
  else return true;

}

vector<passivedouble> CDriver::GetInitialMeshCoord(unsigned short iMarker, unsigned long iVertex) const {

  vector<su2double> coord(3,0.0);
  vector<passivedouble> coord_passive(3, 0.0);

  auto iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  for (auto iDim = 0 ; iDim < nDim ; iDim++){
    // preCICE
   coord[iDim] = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetCoord(iPoint,iDim);
                  //solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord(iPoint,iDim);
    // CSolver object only instantiates coordinates if DEFORM_MESH= YES. This above works regardless, which is handy for CHT
  }

  coord_passive[0] = SU2_TYPE::GetValue(coord[0]);
  coord_passive[1] = SU2_TYPE::GetValue(coord[1]);
  coord_passive[2] = SU2_TYPE::GetValue(coord[2]);

  return coord_passive;
}

vector<passivedouble> CDriver::GetVertexNormal(unsigned short iMarker, unsigned long iVertex, bool unitNormal) const {

  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);
  vector<passivedouble> ret_Normal_passive(3, 0.0);

  Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

  if (!unitNormal) {

    ret_Normal_passive[0] = SU2_TYPE::GetValue(Normal[0]);
    ret_Normal_passive[1] = SU2_TYPE::GetValue(Normal[1]);
    if(nDim>2) ret_Normal_passive[2] = SU2_TYPE::GetValue(Normal[2]);

    return ret_Normal_passive;
  }

  Area = GeometryToolbox::Norm(nDim, Normal);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  ret_Normal_passive[0] = SU2_TYPE::GetValue(ret_Normal[0]);
  ret_Normal_passive[1] = SU2_TYPE::GetValue(ret_Normal[1]);
  ret_Normal_passive[2] = SU2_TYPE::GetValue(ret_Normal[2]);

  return ret_Normal_passive;
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, ecc...) */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetnTimeIter() const {

  return config_container[ZONE_0]->GetnTime_Iter();
}

unsigned long CDriver::GetTime_Iter() const{

  return TimeIter;
}

passivedouble CDriver::GetUnsteady_TimeStep() const {

  return SU2_TYPE::GetValue(config_container[ZONE_0]->GetDelta_UnstTime());
  // preCICE: Changed to GetDelta_UnstTime(), as this is not the initial time step but the ACTUAL time step that is used
}

string CDriver::GetSurfaceFileName() const {

  return config_container[ZONE_0]->GetSurfCoeff_FileName();
}
//////////////////////////////////////////////////////////////////////////////////
/* Functions specifically created for use with preCICE */
//////////////////////////////////////////////////////////////////////////////////

// preCICE:
void CDriver::SetUnsteady_TimeStep(passivedouble val_delta_unsttime) {
    config_container[ZONE_0]->SetDelta_UnstTimeND(val_delta_unsttime / config_container[ZONE_0]->GetTime_Ref());
}

// preCICE:
void CDriver::ReloadOldState() {

  // Get the number of solution variables, points, and dimension
  const unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  const unsigned long nPoint_Local = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain();
  const unsigned short nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  
  // Get if RANS
  const bool rans = config_container[ZONE_0]->GetKind_Turb_Model() != TURB_MODEL::NONE;
  const unsigned short TURB_nVar = (rans) ? solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar() : 0;

  // Get if this is dynamic grid (for unsteady FSI problems)
  const bool dynamic_grid = config_container[ZONE_0]->GetDynamic_Grid();
  const unsigned short MESH_nVar = (dynamic_grid) ? solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetnVar() : 0;

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {

  // Loop through everything and set all necessary variables to current state
  for (unsigned long iPoint_Local = 0; iPoint_Local < nPoint_Local; iPoint_Local++) {

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->SetSolution(iPoint_Local, iVar, preCICE_Solution(iPoint_Local, iVar));
        solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->Set_Solution_time_n(iPoint_Local, iVar, preCICE_Solution_time_n(iPoint_Local, iVar));
        solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->Set_Solution_time_n1(iPoint_Local, iVar, preCICE_Solution_time_n1(iPoint_Local, iVar));
    }
    if (rans) {
      for (unsigned short TURB_iVar = 0; TURB_iVar < TURB_nVar; TURB_iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetNodes()->SetSolution(iPoint_Local, TURB_iVar, preCICE_TURB_Solution(iPoint_Local, TURB_iVar));
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetNodes()->Set_Solution_time_n(iPoint_Local, TURB_iVar, preCICE_TURB_Solution_time_n(iPoint_Local, TURB_iVar));
        solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetNodes()->Set_Solution_time_n1(iPoint_Local, TURB_iVar, preCICE_TURB_Solution_time_n1(iPoint_Local, TURB_iVar));
      }
    }    
    if (dynamic_grid) {

      for (unsigned short MESH_iVar = 0; MESH_iVar < MESH_nVar; MESH_iVar++) {
        solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetSolution(iPoint_Local, MESH_iVar, preCICE_MESH_Solution(iPoint_Local, MESH_iVar));
        solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n(iPoint_Local, MESH_iVar, preCICE_MESH_Solution_time_n(iPoint_Local, MESH_iVar));
        solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n1(iPoint_Local, MESH_iVar, preCICE_MESH_Solution_time_n1(iPoint_Local, MESH_iVar));
      }

      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetCoord(iPoint_Local,iDim, preCICE_Coord(iPoint_Local,iDim));
        geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetGridVel(iPoint_Local, iDim, preCICE_GridVel(iPoint_Local, iDim));
      }
      
      
      //Temporarily must set volume and then set appropriate n, n1, then reset Volume
      // Order may seem awkward, but look at CPoint::SetVolume_____ functions to understand why
      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetVolume(iPoint_Local, preCICE_Volume_nM1(iPoint_Local));
      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetVolume_n();
      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetVolume_nM1();

      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetVolume(iPoint_Local, preCICE_Volume_n(iPoint_Local));
      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetVolume_n();

      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetVolume(iPoint_Local, preCICE_Volume(iPoint_Local));
    }

  }

  }  // end safe global access, pre and postprocessing are thread-safe.
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  FinalizeFLOW_SOL();
  if (rans) FinalizeTURB_SOL();
  if (dynamic_grid) FinalizeMESH_SOL();
}

//preCICE: Finalize FLOW reloads
void CDriver::FinalizeFLOW_SOL() {

  // Get if RANS
  const bool rans = config_container[ZONE_0]->GetKind_Turb_Model() != TURB_MODEL::NONE;

  // Get if this is dynamic grid (for unsteady FSI problems)
  bool dynamic_grid = config_container[ZONE_0]->GetDynamic_Grid();


  /*--- Update the geometry for flows on deforming meshes. ---*/
  if (dynamic_grid) {
    CGeometry::UpdateGeometry(geometry_container[ZONE_0][INST_0], config_container[ZONE_0]);

    for (auto iMesh = 0u; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {

      /*--- Compute the grid velocities on the coarser levels. ---*/
      if (iMesh) geometry_container[ZONE_0][INST_0][iMesh]->SetRestricted_GridVelocity(geometry_container[ZONE_0][INST_0][iMesh - 1]);
      else {
        geometry_container[ZONE_0][INST_0][MESH_0]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], GRID_VELOCITY);
        geometry_container[ZONE_0][INST_0][MESH_0]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], GRID_VELOCITY);
      }
    }
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We also call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION);
  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION);

  /*--- For turbulent/species simulations the flow preprocessing is done by the turbulence/species solver
   *    after it loads its variables (they are needed to compute flow primitives). In case turbulence and species, the
   *    species solver does all the Pre-/Postprocessing. ---*/
  if (!rans &&
      config_container[ZONE_0]->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_0][INST_0][MESH_0], solver_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

    /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
    CSolver::MultigridRestriction(*geometry_container[ZONE_0][INST_0][iMesh - 1], solver_container[ZONE_0][INST_0][iMesh - 1][FLOW_SOL]->GetNodes()->GetSolution(),
                         *geometry_container[ZONE_0][INST_0][iMesh], solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->GetNodes()->GetSolution());
    solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], SOLUTION);
    solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], SOLUTION);

    if (!rans &&
        config_container[ZONE_0]->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->Preprocessing(geometry_container[ZONE_0][INST_0][iMesh], solver_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
  }
}

// preCICE: Finalize TURB reloads
void CDriver::FinalizeTURB_SOL() {
  /*--- MPI solution and compute the eddy viscosity ---*/
  solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION);
  solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION);

  /*--- For turbulent+species simulations the solver Pre-/Postprocessing is done by the species/transition solver. ---*/
  if (config_container[ZONE_0]->GetKind_Species_Model() == SPECIES_MODEL::NONE && config_container[ZONE_0]->GetKind_Trans_Model() == TURB_TRANS_MODEL::NONE) {
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_0][INST_0][MESH_0], solver_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], MESH_0, NO_RK_ITER,
                                            RUNTIME_FLOW_SYS, false);
    solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->Postprocessing(geometry_container[ZONE_0][INST_0][MESH_0], solver_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], MESH_0);
  } else {
    SU2_MPI::Error("Invalid configuration for using implicit coupling! Species and transition models not implemented.", CURRENT_FUNCTION);
    return;
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
    CSolver::MultigridRestriction(*geometry_container[ZONE_0][INST_0][iMesh - 1], solver_container[ZONE_0][INST_0][iMesh - 1][TURB_SOL]->GetNodes()->GetSolution(),
                        *geometry_container[ZONE_0][INST_0][iMesh], solver_container[ZONE_0][INST_0][iMesh][TURB_SOL]->GetNodes()->GetSolution());
    solver_container[ZONE_0][INST_0][iMesh][TURB_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], SOLUTION);
    solver_container[ZONE_0][INST_0][iMesh][TURB_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], SOLUTION);

    if (config_container[ZONE_0]->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->Preprocessing(geometry_container[ZONE_0][INST_0][iMesh], solver_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                            false);
      solver_container[ZONE_0][INST_0][iMesh][TURB_SOL]->Postprocessing(geometry_container[ZONE_0][INST_0][iMesh], solver_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0], iMesh);
    }
  }
} 

// preCICE: Finalize MESH reloads
void CDriver::FinalizeMESH_SOL() {

  // Get the number of solution points and dimension
  const unsigned long nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
  const unsigned short nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  
  /*--- Communicate the loaded displacements. ---*/
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION);

  /*--- Init the linear system solution. ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->LinSysSol(iPoint, iDim) = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution(iPoint, iDim);
    }
  }

  /*--- For time-domain problems, we need to compute the grid velocities. ---*/
  
  /*--- Update the old geometry (coordinates n and n-1) ---*/
  //Only relevant functions from RestartOldGeometry pasted below
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION_TIME_N);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION_TIME_N);
  
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION_TIME_N1);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], SOLUTION_TIME_N1);


  /*--- Once Displacement_n and Displacement_n1 are filled we can compute the Grid Velocity ---*/
  /*--- Compute the velocity of each node. ---*/

  const bool firstOrder = config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST;
  const bool secondOrder = config_container[ZONE_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND;
  const su2double invTimeStep = 1.0 / config_container[ZONE_0]->GetDelta_UnstTimeND();

  //SU2_OMP_FOR_STAT(omp_chunk_size) - omp_chunk_size part of CMeshSolver
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Coordinates of the current point at n+1, n, & n-1 time levels. ---*/

    const su2double* Disp_nM1 = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution_time_n1(iPoint);
    const su2double* Disp_n   = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution_time_n(iPoint);
    const su2double* Disp_nP1 = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution(iPoint);

    /*--- Compute mesh velocity for this point with 1st or 2nd-order approximation. ---*/

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {

      su2double GridVel = 0.0;
      if (firstOrder)
        GridVel = (Disp_nP1[iDim] - Disp_n[iDim]) * invTimeStep;
      else if (secondOrder)
        GridVel = (1.5*Disp_nP1[iDim] - 2.0*Disp_n[iDim] + 0.5*Disp_nM1[iDim]) * invTimeStep;

      geometry_container[ZONE_0][INST_0][MESH_0]->nodes->SetGridVel(iPoint, iDim, GridVel);
    }
  }
  //END_SU2_OMP_FOR

  for (auto iMGlevel = 1u; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++)
    geometry_container[ZONE_0][INST_0][iMGlevel]->SetRestricted_GridVelocity(geometry_container[ZONE_0][INST_0][iMGlevel-1]);



  /*--- Store the boundary displacements at the Bound_Disp variable. ---*/
  for (unsigned short iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {

    if ((config_container[ZONE_0]->GetMarker_All_Deform_Mesh(iMarker) == YES) ||
        (config_container[ZONE_0]->GetMarker_All_Moving(iMarker) == YES)) {

      for (unsigned long iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++) {

        /*--- Get node index. ---*/
        auto iNode = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

        /*--- Set boundary solution. ---*/
        solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iNode, solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution(iNode));
      }
    }
  }

}

// preCICE:
void CDriver::SaveOldState() {

  // Get the number of solution variables, points, and dimension
  // Problem: am looping through global number of points and indexing as such. Not local.
  const unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  const unsigned long nPoint_Local = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain();
  const unsigned short nDim = geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
  
  // Get if RANS
  const bool rans = config_container[ZONE_0]->GetKind_Turb_Model() != TURB_MODEL::NONE;
  const unsigned short TURB_nVar = (rans) ? solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar() : 0;

  // Get if this is dynamic grid (for unsteady FSI problems)
  const bool dynamic_grid = config_container[ZONE_0]->GetDynamic_Grid();
  const unsigned short MESH_nVar = (dynamic_grid) ? solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetnVar() : 0;
  
  // Instantiate all required member variables if they aren't already
  if (preCICE_Solution.empty()) preCICE_Solution.resize(nPoint_Local, nVar) = su2double(0.0);
  if (preCICE_Solution_time_n.empty()) preCICE_Solution_time_n.resize(nPoint_Local,nVar) = su2double(0.0);
  if (preCICE_Solution_time_n1.empty()) preCICE_Solution_time_n1.resize(nPoint_Local,nVar) = su2double(0.0);

  if (rans) {
    if (preCICE_TURB_Solution.empty()) preCICE_TURB_Solution.resize(nPoint_Local,TURB_nVar) = su2double(0.0);
    if (preCICE_TURB_Solution_time_n.empty()) preCICE_TURB_Solution_time_n.resize(nPoint_Local,TURB_nVar) = su2double(0.0);
    if (preCICE_TURB_Solution_time_n1.empty()) preCICE_TURB_Solution_time_n1.resize(nPoint_Local,TURB_nVar) = su2double(0.0);
  }

  if (dynamic_grid) {
    if (preCICE_MESH_Solution.empty()) preCICE_MESH_Solution.resize(nPoint_Local,MESH_nVar) = su2double(0.0);
    if (preCICE_MESH_Solution_time_n.empty()) preCICE_MESH_Solution_time_n.resize(nPoint_Local,MESH_nVar) = su2double(0.0);
    if (preCICE_MESH_Solution_time_n1.empty()) preCICE_MESH_Solution_time_n1.resize(nPoint_Local,MESH_nVar) = su2double(0.0);
  

    if (preCICE_Coord.empty()) preCICE_Coord.resize(nPoint_Local, nDim) = su2double(0.0);
    if (preCICE_GridVel.empty()) preCICE_GridVel.resize(nPoint_Local,nDim) = su2double(0.0);
    if (preCICE_Volume.empty()) preCICE_Volume.resize(nPoint_Local) = su2double(0.0);
    if (preCICE_Volume_n.empty()) preCICE_Volume_n.resize(nPoint_Local) = su2double(0.0);
    if (preCICE_Volume_nM1.empty()) preCICE_Volume_nM1.resize(nPoint_Local) = su2double(0.0);
  }


  // Loop through everything and save all necessary variables to reload state
  for (unsigned long iPoint_Local = 0; iPoint_Local < nPoint_Local; iPoint_Local++) {
    
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      preCICE_Solution(iPoint_Local, iVar) = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSolution(iPoint_Local, iVar);
      preCICE_Solution_time_n(iPoint_Local, iVar) = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint_Local, iVar);
      preCICE_Solution_time_n1(iPoint_Local, iVar) = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint_Local, iVar);
    }

    if (rans) {
      for (unsigned short TURB_iVar = 0; TURB_iVar < TURB_nVar; TURB_iVar++) {
        preCICE_TURB_Solution(iPoint_Local, TURB_iVar) = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint_Local, TURB_iVar);
        preCICE_TURB_Solution_time_n(iPoint_Local, TURB_iVar) = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetNodes()->GetSolution_time_n(iPoint_Local, TURB_iVar);
        preCICE_TURB_Solution_time_n1(iPoint_Local, TURB_iVar) = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetNodes()->GetSolution_time_n1(iPoint_Local, TURB_iVar);
      }
    }

    if (dynamic_grid) {
      for (unsigned short MESH_iVar = 0; MESH_iVar < MESH_nVar; MESH_iVar++) {
        preCICE_MESH_Solution(iPoint_Local, MESH_iVar) = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution(iPoint_Local, MESH_iVar);
        preCICE_MESH_Solution_time_n(iPoint_Local, MESH_iVar) = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution_time_n(iPoint_Local, MESH_iVar);
        preCICE_MESH_Solution_time_n1(iPoint_Local, MESH_iVar) = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetSolution_time_n1(iPoint_Local, MESH_iVar);
      }


      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        preCICE_Coord(iPoint_Local,iDim) = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetCoord(iPoint_Local,iDim);
        preCICE_GridVel(iPoint_Local, iDim) = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetGridVel(iPoint_Local)[iDim];
      }
      
      preCICE_Volume(iPoint_Local) = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetVolume(iPoint_Local);
      preCICE_Volume_n(iPoint_Local) = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetVolume_n(iPoint_Local);
      preCICE_Volume_nM1(iPoint_Local) = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetVolume_nM1(iPoint_Local);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
/* Functions related to CHT solver                                           */
///////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetVertexTemperature(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  su2double vertexWallTemp(0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint) && compressible){
    vertexWallTemp = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
  }

  //preCICE: re-dimensionalize before returning
  return SU2_TYPE::GetValue(vertexWallTemp * config_container[ZONE_0]->GetTemperature_Ref());

}

void CDriver::SetVertexTemperature(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallTemp){

  // preCICE: non-dimensionalize before setting
  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryTemperature(iMarker, iVertex, val_WallTemp / config_container[ZONE_0]->GetTemperature_Ref());
}

vector<passivedouble> CDriver::GetVertexHeatFluxes(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  unsigned short iDim;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;
  vector<su2double> GradT (3,0.0);
  vector<su2double> HeatFlux (3,0.0);
  vector<passivedouble> HeatFluxPassive (3,0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(compressible){
    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
      HeatFlux[iDim] = -thermal_conductivity*GradT[iDim];
    }
  }

  //preCICE: re-dimensionalize before returning
  HeatFluxPassive[0] = SU2_TYPE::GetValue(HeatFlux[0] * config_container[ZONE_0]->GetHeat_Flux_Ref());
  HeatFluxPassive[1] = SU2_TYPE::GetValue(HeatFlux[1] * config_container[ZONE_0]->GetHeat_Flux_Ref());
  HeatFluxPassive[2] = SU2_TYPE::GetValue(HeatFlux[2] * config_container[ZONE_0]->GetHeat_Flux_Ref());

  return HeatFluxPassive;
}

passivedouble CDriver::GetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex) const{

  unsigned long iPoint;
  unsigned short iDim;
  su2double vertexWallHeatFlux;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Area;
  su2double laminar_viscosity, thermal_conductivity, dTdn;
  su2double *Normal, GradT[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  vertexWallHeatFlux = 0.0;
  dTdn = 0.0;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint) && compressible){
    Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

    Area = GeometryToolbox::Norm(nDim, Normal);

    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    /*Compute wall heat flux (normal to the wall) based on computed temperature gradient*/
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
      dTdn += GradT[iDim]*UnitNormal[iDim];
    }

    vertexWallHeatFlux = -thermal_conductivity*dTdn;
  }

   //preCICE: re-dimensionalize before returning
  return SU2_TYPE::GetValue(vertexWallHeatFlux * config_container[ZONE_0]->GetHeat_Flux_Ref());
}

void CDriver::SetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallHeatFlux){

  // preCICE: non-dimensionalize before setting
  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryHeatFlux(iMarker, iVertex, val_WallHeatFlux / config_container[ZONE_0]->GetHeat_Flux_Ref());
}

passivedouble CDriver::GetThermalConductivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
  thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);

  return SU2_TYPE::GetValue(thermal_conductivity);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers                             */
////////////////////////////////////////////////////////////////////////////////

vector<string> CDriver::GetAllBoundaryMarkersTag() const {

  vector<string> boundariesTagList;
  unsigned short iMarker,nBoundariesMarkers;
  string Marker_Tag;

  nBoundariesMarkers = config_container[ZONE_0]->GetnMarker_All();
  boundariesTagList.resize(nBoundariesMarkers);

  for(iMarker=0; iMarker < nBoundariesMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    boundariesTagList[iMarker] = Marker_Tag;
  }

  return boundariesTagList;
}

vector<string> CDriver::GetAllDeformMeshMarkersTag() const {

  vector<string> interfaceBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Deform_Mesh();
  interfaceBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Deform_Mesh_TagBound(iMarker);
    interfaceBoundariesTagList[iMarker] = Marker_Tag;
  }

  return interfaceBoundariesTagList;
}

vector<string> CDriver::GetAllCHTMarkersTag() const {

  vector<string> CHTBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();

  //The CHT markers can be identified as the markers that are customizable with a BC type HEAT_FLUX or ISOTHERMAL.
  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    if((config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) && config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker)){
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      CHTBoundariesTagList.push_back(Marker_Tag);
    }
  }

  return CHTBoundariesTagList;
}

vector<string> CDriver::GetAllInletMarkersTag() const {

  vector<string> BoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    bool isCustomizable = config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker);
    bool isInlet = (config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
    if(isCustomizable && isInlet) {
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      BoundariesTagList.push_back(Marker_Tag);
    }
  }

  return BoundariesTagList;
}

map<string, int> CDriver::GetAllBoundaryMarkers() const {

  map<string, int>  allBoundariesMap;
  unsigned short iMarker, nBoundaryMarkers;
  string Marker_Tag;

  nBoundaryMarkers = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker < nBoundaryMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    allBoundariesMap[Marker_Tag] = iMarker;
  }

  return allBoundariesMap;
}

map<string, string> CDriver::GetAllBoundaryMarkersType() const {

  map<string, string> allBoundariesTypeMap;
  unsigned short iMarker, KindBC;
  string Marker_Tag, Marker_Type;

  for(iMarker=0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    KindBC = config_container[ZONE_0]->GetMarker_All_KindBC(iMarker);
    switch(KindBC){
      case EULER_WALL:
        Marker_Type = "EULER_WALL";
        break;
      case FAR_FIELD:
        Marker_Type = "FARFIELD";
        break;
      case ISOTHERMAL:
        Marker_Type = "ISOTHERMAL";
        break;
      case HEAT_FLUX:
        Marker_Type = "HEATFLUX";
        break;
      case INLET_FLOW:
        Marker_Type = "INLET_FLOW";
        break;
      case OUTLET_FLOW:
        Marker_Type = "OUTLET_FLOW";
        break;
      case SYMMETRY_PLANE:
        Marker_Type = "SYMMETRY";
        break;
      case SEND_RECEIVE:
        Marker_Type = "SEND_RECEIVE";
        break;
      default:
        Marker_Type = "UNKNOWN_TYPE";
    }
    allBoundariesTypeMap[Marker_Tag] = Marker_Type;
  }

  return allBoundariesTypeMap;
}

void CDriver::SetHeatSource_Position(passivedouble alpha, passivedouble pos_x, passivedouble pos_y, passivedouble pos_z){

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][RAD_SOL];

  config_container[ZONE_0]->SetHeatSource_Rot_Z(alpha);
  config_container[ZONE_0]->SetHeatSource_Center(pos_x, pos_y, pos_z);

  solver->SetVolumetricHeatSource(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0]);

}

void CDriver::SetInlet_Angle(unsigned short iMarker, passivedouble alpha){

  su2double alpha_rad = alpha * PI_NUMBER/180.0;

  unsigned long iVertex;

  for (iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, 0, cos(alpha_rad));
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, 1, sin(alpha_rad));
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, ecc...) */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
      integration_container[iZone][INST_0][FLOW_SOL]->SetConvergence(false);
      if (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) integration_container[iZone][INST_0][TURB_SOL]->SetConvergence(false);
      if(config_container[iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) integration_container[iZone][INST_0][TRANS_SOL]->SetConvergence(false);
      break;

    case MAIN_SOLVER::FEM_ELASTICITY:
      integration_container[iZone][INST_0][FEA_SOL]->SetConvergence(false);
      break;

    case MAIN_SOLVER::ADJ_EULER: case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      integration_container[iZone][INST_0][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) || (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) )
        integration_container[iZone][INST_0][ADJTURB_SOL]->SetConvergence(false);
      break;

    default:
      break;
    }
  }

}

void CSinglezoneDriver::SetInitialMesh() {

  DynamicMeshUpdate(0);

  SU2_OMP_PARALLEL {
    // Overwrite fictious velocities
    for (iMesh = 0u; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      SU2_OMP_FOR_STAT(roundUpDiv(geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(),omp_get_max_threads()))
      for (unsigned long iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(); iPoint++) {

        /*--- Overwrite fictitious velocities ---*/
        su2double Grid_Vel[3] = {0.0, 0.0, 0.0};

        /*--- Set the grid velocity for this coarse node. ---*/
        geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetGridVel(iPoint, Grid_Vel);
      }
      END_SU2_OMP_FOR
      /*--- Push back the volume. ---*/
      geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetVolume_n();
      geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetVolume_nM1();
    }
    /*--- Push back the solution so that there is no fictious velocity at the next step. ---*/
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n();
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n1();
  }
  END_SU2_OMP_PARALLEL
}

void CDriver::BoundaryConditionsUpdate(){

  int rank = MASTER_NODE;
  unsigned short iZone;

  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  if(rank == MASTER_NODE) cout << "Updating boundary conditions." << endl;
  for(iZone = 0; iZone < nZone; iZone++){
    geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
  }
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to finite elements                                       */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetFEA_Loads(unsigned short iMarker, unsigned long iVertex, passivedouble LoadX,
                       passivedouble LoadY, passivedouble LoadZ) {

  unsigned long iPoint;
  su2double NodalForce[3] = {0.0,0.0,0.0};
  NodalForce[0] = LoadX;
  NodalForce[1] = LoadY;
  NodalForce[2] = LoadZ;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->GetNodes()->Set_FlowTraction(iPoint,NodalForce);

}

vector<passivedouble> CDriver::GetFEA_Displacements(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Displacements(3, 0.0);
  vector<passivedouble> Displacements_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  Displacements[0] = solver->GetNodes()->GetSolution(iPoint, 0);
  Displacements[1] = solver->GetNodes()->GetSolution(iPoint, 1);
  if (geometry->GetnDim() == 3)
    Displacements[2] = solver->GetNodes()->GetSolution(iPoint, 2);
  else
    Displacements[2] = 0.0;

  Displacements_passive[0] = SU2_TYPE::GetValue(Displacements[0]);
  Displacements_passive[1] = SU2_TYPE::GetValue(Displacements[1]);
  Displacements_passive[2] = SU2_TYPE::GetValue(Displacements[2]);

  return Displacements_passive;
}


vector<passivedouble> CDriver::GetFEA_Velocity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Velocity(3, 0.0);
  vector<passivedouble> Velocity_passive(3,0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC){
    Velocity[0] = solver->GetNodes()->GetSolution_Vel(iPoint, 0);
    Velocity[1] = solver->GetNodes()->GetSolution_Vel(iPoint, 1);
    if (geometry->GetnDim() == 3)
      Velocity[2] = solver->GetNodes()->GetSolution_Vel(iPoint, 2);
    else
      Velocity[2] = 0.0;
  }

  Velocity_passive[0] = SU2_TYPE::GetValue(Velocity[0]);
  Velocity_passive[1] = SU2_TYPE::GetValue(Velocity[1]);
  Velocity_passive[2] = SU2_TYPE::GetValue(Velocity[2]);

  return Velocity_passive;
}

vector<passivedouble> CDriver::GetFEA_Velocity_n(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Velocity_n(3, 0.0);
  vector<passivedouble> Velocity_n_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC){
    Velocity_n[0] = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, 0);
    Velocity_n[1] = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, 1);
    if (geometry->GetnDim() == 3)
      Velocity_n[2] = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, 2);
    else
      Velocity_n[2] = 0.0;
  }

  Velocity_n_passive[0] = SU2_TYPE::GetValue(Velocity_n[0]);
  Velocity_n_passive[1] = SU2_TYPE::GetValue(Velocity_n[1]);
  Velocity_n_passive[2] = SU2_TYPE::GetValue(Velocity_n[2]);

  return Velocity_n_passive;

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to adjoint simulations                                   */
////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetMeshDisp_Sensitivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Disp_Sens(3, 0.0);
  vector<passivedouble> Disp_Sens_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver =  solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  Disp_Sens[0] = solver->GetNodes()->GetBoundDisp_Sens(iPoint, 0);
  Disp_Sens[1] = solver->GetNodes()->GetBoundDisp_Sens(iPoint, 1);
  if (geometry->GetnDim() == 3)
    Disp_Sens[2] = solver->GetNodes()->GetBoundDisp_Sens(iPoint, 2);
  else
    Disp_Sens[2] = 0.0;

  Disp_Sens_passive[0] = SU2_TYPE::GetValue(Disp_Sens[0]);
  Disp_Sens_passive[1] = SU2_TYPE::GetValue(Disp_Sens[1]);
  Disp_Sens_passive[2] = SU2_TYPE::GetValue(Disp_Sens[2]);

  return Disp_Sens_passive;

}

vector<passivedouble> CDriver::GetFlowLoad_Sensitivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> FlowLoad_Sens(3, 0.0);
  vector<passivedouble> FlowLoad_Sens_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  FlowLoad_Sens[0] = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, 0);
  FlowLoad_Sens[1] = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, 1);
  if (geometry->GetnDim() == 3)
    FlowLoad_Sens[2] = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, 2);
  else
    FlowLoad_Sens[2] = 0.0;

  FlowLoad_Sens_passive[0] = SU2_TYPE::GetValue(FlowLoad_Sens[0]);
  FlowLoad_Sens_passive[1] = SU2_TYPE::GetValue(FlowLoad_Sens[1]);
  FlowLoad_Sens_passive[2] = SU2_TYPE::GetValue(FlowLoad_Sens[2]);

  return FlowLoad_Sens_passive;

}

void CDriver::SetFlowLoad_Adjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                  passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 0, val_AdjointX);
  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 2, val_AdjointZ);

}

void CDriver::SetSourceTerm_DispAdjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                        passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  unsigned long iPoint;

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];
  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, 0, val_AdjointX);
  solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, 2, val_AdjointZ);

}

void CDriver::SetSourceTerm_VelAdjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                        passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];
  const auto iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, 0, val_AdjointX);
  solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, 2, val_AdjointZ);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to mesh deformation */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetMeshDisplacement(unsigned short iMarker, unsigned long iVertex, passivedouble DispX, passivedouble DispY, passivedouble DispZ) {

  unsigned long iPoint;
  su2double MeshDispl[3] =  {0.0,0.0,0.0};

  MeshDispl[0] = DispX;
  MeshDispl[1] = DispY;
  MeshDispl[2] = DispZ;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iPoint,MeshDispl);

}

void CDriver::CommunicateMeshDisplacement(void) {

  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0],
                                                                    config_container[ZONE_0], MESH_DISPLACEMENTS);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0],
                                                                    config_container[ZONE_0], MESH_DISPLACEMENTS);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to flow loads */
////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetFlowLoad(unsigned short iMarker, unsigned long iVertex) const {

  vector<su2double> FlowLoad(3, 0.0);
  vector<passivedouble> FlowLoad_passive(3, 0.0);

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetSolid_Wall(iMarker)) {
    FlowLoad[0] = solver->GetVertexTractions(iMarker, iVertex, 0);
    FlowLoad[1] = solver->GetVertexTractions(iMarker, iVertex, 1);
    if (geometry->GetnDim() == 3)
      FlowLoad[2] = solver->GetVertexTractions(iMarker, iVertex, 2);
    else
      FlowLoad[2] = 0.0;
  }

  FlowLoad_passive[0] = SU2_TYPE::GetValue(FlowLoad[0]);
  FlowLoad_passive[1] = SU2_TYPE::GetValue(FlowLoad[1]);
  FlowLoad_passive[2] = SU2_TYPE::GetValue(FlowLoad[2]);

  return FlowLoad_passive;

}
