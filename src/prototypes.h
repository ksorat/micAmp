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
void InitializeIntegrator(Grid_S Grid, Model_S Model);
void DestroyIntegrator(Grid_S Grid, Model_S Model);
void FluxUpdate(RealP4 Prim, RealP4 Fx, RealP4 Fy, RealP4 Fz, Real dt, Grid_S Grid);


//File recon.cpp
void Flux_PCM(RealP4 W, RealP4 Fx, RealP4 Fy, RealP4 Fz, Grid_S Grid, Model_S Model);
void Flux_PLM(RealP4 W, RealP4 Fx, RealP4 Fy, RealP4 Fz, Grid_S Grid, Model_S Model);
inline void Recon_PLM(Real Wmp[2], Real Wl, Real Wc, Real Wr);
inline Real SlopeLimit(Real dWl, Real dWr,Real dWc);
void LRs2Flux(RealP4 lW, RealP4 rW, RealP4 Flx, int d,  Grid_S Grid);

//File flux.cpp
void RiemannFluxHLLE(Real LeftW[NVAR][VECBUFF], Real RightW[NVAR][VECBUFF], Real FluxLR[NVAR][VECBUFF]);
void Roes_Vec(Real LeftW[NVAR][VECBUFF],Real RightW[NVAR][VECBUFF],Real RoeLR[NVAR][VECBUFF],Real evals[NVAR][VECBUFF]);
void CalcLR_Fluxes(Real LeftW[NVAR][VECBUFF],Real RightW[NVAR][VECBUFF], Real Fl[NVAR][VECBUFF], Real Fr[NVAR][VECBUFF], Real bm[VECBUFF], Real bp[VECBUFF] );
inline Real inlineEnthalpy(Real d, Real Vx, Real Vy, Real Vz, Real P,Real Gam);

//File viscosity.cpp
void Flux_Viscous(RealP4 W, RealP4 Fx, RealP4 Fy, RealP4 Fz, Grid_S Grid, Model_S Model);
void iJacV(RealP4 W, int ib0, int j0, int k0, int iLim, Real JacV_L[NDIM][NDIM][VECBUFFV], Real JacV_R[NDIM][NDIM][VECBUFFV], Real JacV_I[NDIM][NDIM][VECBUFFV], int intDir, Real dx, Real dy, Real dz);
Real Sig(Real Wl, Real Wc, Real Wr, Real dx);
Real ShearViscosity(Real x, Real y, Real z, Real P);
Real BulkViscosity(Real x, Real y, Real z, Real P);

#endif //PROTOTYPES_H