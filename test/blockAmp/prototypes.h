#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "amps.h"
#include <string>

//File: aux
void conModel(Model_S *Model);
void conGrid(Grid_S *Grid);
void conDimension( Real *Xbd, int Nxp, Real *dx, int *inds, Real *xc, Real *xi);
void initialConds(RealP4 State, Grid_S Grid, Model_S Model);
Real CalcDT(RealP4 State, Grid_S Grid, Model_S Model);
void EnforceBCs(RealP4 State, Grid_S Grid, Model_S Model);


//File: amp_array.cpp
RealP4 Create4Array(int N1, int N2, int N3, int N4);
RealP4 Map4Array(Real *Start, int N1, int N2, int N3, int N4);
void Kill4Array(RealP4 ToDie);
void Wipe4Array(RealP4 Ar, int N1, int N2, int N3, int N4);
void Copy4Array(RealP4 A, RealP4 B, int N1, int N2, int N3, int N4);
MICTYPE void WipeBlockCC(BlockCC A, Block_S Block);
MICTYPE void WipeBlockIC(BlockIC A, Block_S Block);
MICTYPE void PrintBlockCC(BlockCC A, Block_S Block);
MICTYPE void PrintBlockIC(BlockIC A, Block_S Block);
MICTYPE void PrintBlockR(BlockR A);

//File: amp_output.cpp
void toConsole(RealP4 State, Grid_S Grid, Model_S Model);
void toVTK(RealP4 State, Grid_S Grid, Model_S Model);
void bswap(void *vdat, int len, int cnt); //Borrowed from athena
void WriteVar(std::string varname, int varnum, Real ****State, Grid_S Grid, FILE *VTKFile, float *data);

//File: phintegrate.cpp
MICTYPE void AdvanceFluid(BlockCC State, Block_S Block, Model_S Model, Real dt);
MICTYPE void FluxUpdate(BlockCC Prim, BlockIC Fx, BlockIC Fy, BlockIC Fz, Real dt, Block_S Grid);

//File: blockdecomp.cpp
void BlockAdvance(RealP4 State, Grid_S Grid, Model_S Model, Real dt);
void InitializeIntegrator(Grid_S Grid, Model_S Model);
void DestroyIntegrator(Grid_S Grid, Model_S Model);
void CopyinBlock (RealP4 Q, BlockCC Qblk, Grid_S Grid, Block_S Block);
void CopyoutBlock(RealP4 Q, BlockCC Qblk, Grid_S Grid, Block_S Block);
void InitBlock(Block_S *Block, Grid_S Grid, int iBlk, int jBlk, int kBlk);
Block_S*** MapBlocks(Grid_S Grid, int Bx, int By, int Bz);

//File: flux.cpp
MICTYPE void Flux_PCM(BlockCC State,BlockIC Flux_x,BlockIC Flux_y,BlockIC Flux_z,Block_S Block,Model_S Model);
MICTYPE void LRs2Flux(BlockCC lW,BlockCC rW, BlockIC Flx, int d,  Block_S Grid,Model_S Model);

//File: rsolve.cpp
MICTYPE void RiemannFluxHLLE(BlockR LeftW,BlockR RightW,BlockR FluxLR,Real Gam);
MICTYPE void Roes_Vec(BlockR LeftW,BlockR RightW,BlockR RoeLR,BlockR evals,Real Gam);
MICTYPE inline Real inlineEnthalpy(Real d, Real Vx, Real Vy, Real Vz, Real P,Real Gam);
MICTYPE void CalcLR_Fluxes(BlockR LeftW, BlockR RightW, BlockR Fl, BlockR Fr, Real bm[VECBUFF], Real bp[VECBUFF], Real Gam );

#endif //PROTOTYPES_H