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

//File: amp_output.cpp
void toConsole(RealP4 State, Grid_S Grid, Model_S Model);
void toVTK(RealP4 State, Grid_S Grid, Model_S Model);
void bswap(void *vdat, int len, int cnt); //Borrowed from athena
void WriteVar(std::string varname, int varnum, Real ****State, Grid_S Grid, FILE *VTKFile, float *data);

//File integrator.cpp
void AdvanceFluid(RealP4 State, Grid_S Grid, Model_S Model, Real dt);
void BlockAdvance(RealP4 State, Grid_S Grid, Model_S Model, Real dt);
void InitializeIntegrator(Grid_S Grid, Model_S Model);
void DestroyIntegrator(Grid_S Grid, Model_S Model);
void FluxUpdate(RealP4 Prim, RealP4 Fx, RealP4 Fy, RealP4 Fz, Real dt, Grid_S Grid);
void CopyinBlock (RealP4 Q, Real Qblk[NVAR][NZBLK][NYBLK][NXBLK], Grid_S Grid, int iBlk, int jBlk, int kBlk);
void CopyoutBlock(RealP4 Q, Real Qblk[NVAR][NZBLK][NYBLK][NXBLK], Grid_S Grid, int iBlk, int jBlk, int kBlk);



#endif //PROTOTYPES_H