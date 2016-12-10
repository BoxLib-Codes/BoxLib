#!/usr/bin/env bash

BOXLIB_HEADERS="\
Amr.H \
AmrCore.H \
AmrData.H \
AmrLevel.H \
AmrParGDB.H \
AmrParticles.H \
AmrvisConstants.H \
Arena.H \
Array.H \
ArrayLim.H \
AuxBoundaryData.H \
BArena.H \
BCRec.H \
BC_TYPES.H \
BLBackTrace.H \
BLFort.H \
BLPgas.H \
BLProfiler.H \
BLassert.H \
BaseFab.H \
BaseFab_f.H \
BndryData.H \
BndryRegister.H \
BoundCond.H \
Box.H \
BoxArray.H \
BoxDomain.H \
BoxLib.H \
BoxList.H \
CArena.H \
CONSTANTS.H \
COORDSYS_F.H \
Cluster.H \
CoordSys.H \
DataServices.H \
Derive.H \
DistributionMapping.H \
ErrorList.H \
Extrapolater.H \
FArrayBox.H \
FLUSH_F.H \
FLUXREG_F.H \
FPC.H \
FabArray.H \
FabConv.H \
FabSet.H \
FillPatchUtil.H \
FluxRegister.H \
Geometry.H \
IArrayBox.H \
INTERPBNDRYDATA_F.H \
INTERP_F.H \
IndexType.H \
IntVect.H \
InterpBndryData.H \
Interpolater.H \
LO_BCTYPES.H \
Lazy.H \
LevelBld.H \
Looping.H \
MAKESLICE_F.H \
MacBndry.H \
Mask.H \
MemPool.H \
MemProfiler.H \
MultiFab.H \
MultiFabUtil.H \
MultiFabUtil_F.H \
MultiMask.H \
NFiles.H \
Orientation.H \
PList.H \
PROB_AMR_F.H \
ParGDB.H \
ParallelDescriptor.H \
ParmParse.H \
ParticleInit.H \
Particles.H \
Particles_F.H \
Periodicity.H \
PhysBCFunct.H \
PlotFileUtil.H \
Pointers.H \
REAL.H \
RealBox.H \
SLABSTAT_F.H \
SPACE.H \
SPACE_F.H \
SlabStat.H \
StateData.H \
StateDescriptor.H \
StationData.H \
TagBox.H \
TinyProfiler.H \
TracerParticles.H \
Tuple.H \
UseCount.H \
Utility.H \
VisMF.H \
ccse-mpi.H \
iMultiFab.H \
winstd.H"

NUM_HEADERS=102
I_HEADER=0

for header in ${BOXLIB_HEADERS}; do
  find . -path .git -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)'"${header}"'\(\s*\(>\|\"\)\)/\1AMReX_'"${header}"'\3/g' {} +
  I_HEADER=$((I_HEADER+1))
  percent=$(awk "BEGIN { pc=100*${I_HEADER}/${NUM_HEADERS}; i=int(pc); print (pc-i<0.5)?i:i+1 }")
  echo -ne "Progress: ${percent}%"\\r
done
