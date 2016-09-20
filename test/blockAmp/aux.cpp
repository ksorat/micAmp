//Handles all the routines for the shell of the driver
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amps.h"
#include "prototypes.h"


void conModel(Model_S *Model) {
	//For now just put values in

	Model->C0 = 0.4;
	for (int i=0;i<3;i++) {
		Model->gVec[i] = 0.0;
	}
	Model->TsOut = TSOUT;

	Model->Gam = 5.0/3.0;
}

//Initialize problem


void initialConds(RealP4 State, Grid_S Grid, Model_S Model) {
	int i, j, k;

	//Initialize a spherical blast
	Real P0, Vx0, Vy0, Vz0, d0; // Simple default values
	Real rad, Rcrit = 0.25;
	Real DelP = 4.0;

	Vx0 = 0.0; Vy0 = 0.0; Vz0 = 0.0;
	P0 = 1.0; d0 = 1.0;
	for (k=Grid.ksd;k<=Grid.ked;k++) {
		for (j=Grid.jsd;j<=Grid.jed;j++) {
			for (i=Grid.isd;i<=Grid.ied;i++) {
				State[PRESSURE][k][j][i] = P0;
				State[VELX][k][j][i] = Vx0;
				State[VELY][k][j][i] = Vy0;
				State[VELZ][k][j][i] = Vz0;
				State[DEN][k][j][i] = d0;
				rad = sqrt( Grid.xc[i]*Grid.xc[i] + Grid.yc[j]*Grid.yc[j] + Grid.zc[k]*Grid.zc[k] );
				if (rad <= Rcrit) {
					State[PRESSURE][k][j][i] *= DelP;
				}
			}
		}
	}

}

Real CalcDT(RealP4 State, Grid_S Grid, Model_S Model) {
	int i,j,k;

	Real SMALL_DT = 1.0e-10;
	Real dx_min, dx_max, dx_mean, Vmax, Vx, Vy, Vz, Vsq, Vabs, Cs;
	Real dt, dt_V;
	//Real Re, ReMin;

	dx_min = std::min(Grid.dx, std::min(Grid.dy, Grid.dz));
	dx_max = std::max(Grid.dx, std::max(Grid.dy, Grid.dz));
	//dx_mean = pow(3*0.25*Grid.dx*Grid.dy*Grid.dz/PI,1.0/3);

	//Calculate timestep
	// dt = CFL * min{ dx/v+c, rho*dx2/mu, rho*cp*dx2/lam } : ACM
	// dt = CFL * min{ dx/v+c } : MUSCL


	//Find largest velocity on grid
	Vmax = 0.0; //ReMin = HUGE; Re = 0.0;
	for (k=Grid.ks;k<=Grid.ke;k++) {
		for (j=Grid.js;j<=Grid.je;j++) {
			for (i=Grid.is;i<=Grid.ie;i++) {
				//Not optimal, but simple for ACM/MUSCL
				Cs = sqrt( Model.Gam*State[PRESSURE][k][j][i]/State[DEN][k][j][i] );
				Vx = State[VELX][k][j][i];
				Vy = State[VELY][k][j][i];
				Vz = State[VELZ][k][j][i];

				Vsq = sqrt( Vx*Vx + Vy*Vy + Vz*Vz );
				Vabs = Vsq + Cs;
				Vmax = fmax(Vmax, Vabs);
				//Re = State[DEN][k][j][i]*Vabs*dx_mean/SHEARVISC;
				//ReMin = fmin(Re,ReMin);
			}
		}
	}


	dt_V = dx_min/Vmax;

	dt = Model.C0*dt_V;
	//dt = VISCSIG*dt/(1+ (2/ReMin));
	return dt;
}

void EnforceBCs(RealP4 State, Grid_S Grid, Model_S Model) {

	int i,j,k,v;
	
	//x boundaries (ibx/obx)
	for (v=0;v<NVAR;v++) {
		for (k=Grid.ksd;k<=Grid.ked;k++) {
			for (j=Grid.jsd;j<=Grid.jed;j++) {
				for (i=1;i<=Grid.Ng;i++) {
					State[v][k][j][Grid.is-i] = State[v][k][j][Grid.ie-i+1];
					State[v][k][j][Grid.ie+i] = State[v][k][j][Grid.is+i-1];
				}
			}
		}	
	}

	//y boundaries (iby/oby)
	for (v=0;v<NVAR;v++) { 
		for (k=Grid.ksd;k<=Grid.ked;k++) {
			for (j=1;j<=Grid.Ng;j++) {
				for (i=Grid.isd;i<=Grid.ied;i++) {
					State[v][k][Grid.js-j][i] = State[v][k][Grid.je-j+1][i];
					State[v][k][Grid.je+j][i] = State[v][k][Grid.js+j-1][i];
				}
			}
		}	
	}

	//z boundaries (ibz/obz)
	for (v=0;v<NVAR;v++) { 
		for (k=1;k<=Grid.Ng;k++) {
			for (j=Grid.jsd;j<=Grid.jed;j++) {
				for (i=Grid.isd;i<=Grid.ied;i++) {
					State[v][Grid.ks-k][j][i] = State[v][Grid.ke-k+1][j][i];
					State[v][Grid.ke+k][j][i] = State[v][Grid.ks+k-1][j][i];
				}
			}
		}
	}

}
void conGrid(Grid_S *Grid) {

	int inds[4];
	Real dxi;

	//Grab some values from defs (should be changed to input deck)
	Grid->Ng = NUMGHOST;
	Grid->Nv = NVAR;
	//Physical cells
	Grid->Nxp = NXP;
	Grid->Nyp = NYP;
	Grid->Nzp = NZP;

	Grid->Nx = Grid->Nxp + 2*NUMGHOST;
	Grid->Ny = Grid->Nyp + 2*NUMGHOST;
	Grid->Nz = Grid->Nzp + 2*NUMGHOST;

	Grid->t = 0.0;
	Grid->dt = 0.0;
	Grid->Ts = 0;

	//Grab domain bounds from defs (use same for all dims)
	Grid->Xbd[0] = XMIN; Grid->Xbd[1] = XMAX;
	Grid->Ybd[0] = XMIN; Grid->Ybd[1] = XMAX;
	Grid->Zbd[0] = XMIN; Grid->Zbd[1] = XMAX;

	Grid->Tfin = TFIN;

	//Create X-dimension
	conDimension(Grid->Xbd, Grid->Nxp, &dxi, inds, Grid->xc,Grid->xi);
	Grid->is = inds[0]; Grid->ie = inds[1]; Grid->isd = inds[2]; Grid->ied = inds[3];
	Grid->dx = dxi;

	//Create Y-dimension
	conDimension(Grid->Ybd, Grid->Nyp, &dxi, inds, Grid->yc,Grid->yi);
	Grid->js = inds[0]; Grid->je = inds[1]; Grid->jsd = inds[2]; Grid->jed = inds[3];
	Grid->dy = dxi;

	//Create Z-dimension
	conDimension(Grid->Zbd, Grid->Nzp, &dxi, inds, Grid->zc,Grid->zi);
	Grid->ks = inds[0]; Grid->ke = inds[1]; Grid->ksd = inds[2]; Grid->ked = inds[3];
	Grid->dz = dxi;

	//Store vector info
	Grid->N[0] = Grid->Nv; Grid->N[1] = Grid->Nz; Grid->N[2] = Grid->Ny; Grid->N[3] = Grid->Nx;
	for (int i=1;i<=3;i++) {
		Grid->Ni[i] = Grid->N[i]+1;
	}
}

//Create an arbitrary dimension, notation is for x
void conDimension( Real *Xbd, int Nxp, Real *dx, int *inds, Real *xc, Real *xi) {
	Real xB, xE, x0i, dxc;
	int Nx, Ng, i, is, ie, isd, ied;
	

	Ng = NUMGHOST; //Could change this from defined constant

	xB = Xbd[0]; xE = Xbd[1];
	dxc = (xE-xB)/Nxp;

	Nx = Nxp+2*Ng;

	//Create inds=[is,ie,isd,ied]
	isd = 0; ied = Nx-1;
	is = isd+Ng;
	ie = ied-Ng;
	inds[0] = is; inds[1] = ie; inds[2] = isd; inds[3] = ied;

	//Create xc (cell centers) and interfaces
	x0i = xB-Ng*dxc;
	for (i=isd;i<=ied+1;i++) {
		xi[i] = x0i + i*dxc;
		
	}
	for (i=isd;i<=ied;i++) {
		xc[i] = 0.5*( xi[i] + xi[i+1] );
	}
	*dx = dxc;
	
}