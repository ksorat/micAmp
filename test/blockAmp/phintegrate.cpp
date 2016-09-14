#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

//RealP4 Flux_x, Flux_y, Flux_z;
RealP4 advState;

//Break up input 4D State into blocks and send each to advance routine
void BlockAdvance(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	RealP4 swpState;
	Real Qblk[NVAR][NZBLK][NYBLK][NXBLK];
	Grid_S SubGrid;

	int iblk,jblk,kblk,nblk;
	int i,j,k,nv,kP,jP,iP;
	

	int Ng = Grid.Ng;
	printf("Breaking grid into %d (%d,%d,%d) blocks\n",BX*BY*BZ,NXPBLK,NYPBLK,NZPBLK);

	nblk = 0;
	//Loop over blocks
	for (kblk=0;kblk<BZ;kblk++) {
		for (jblk=0;jblk<BY;jblk++) {
			for (iblk=0;iblk<BX;iblk++) {
				printf("Block %d (%d,%d,%d)\n",nblk,iblk,jblk,kblk);

				//Copy from State->Block
				CopyinBlock(State,Qblk,Grid,iblk,jblk,kblk);

				//Advance sub-block

				//Copy advanced sub-block back into advState holder
				//Avoid ghosts
				CopyoutBlock(advState,Qblk,Grid,iblk,jblk,kblk);

				nblk++;
			}
		}
	} //Loop over blocks

	//Now swap advState and State pointers so that State is updated
	Copy4Array(State,advState,Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);
}


//MUSCL update
void AdvanceFluid(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	//1-Step
	//FluxUpdate(State,Flux_x,Flux_y,Flux_z,dt,Grid);

}

void FluxUpdate(RealP4 Prim, RealP4 Fx, RealP4 Fy, RealP4 Fz, Real dt, Grid_S Grid) {
	__assume_aligned(Prim, ALIGN);
	__assume_aligned(Fx, ALIGN);
	__assume_aligned(Fy, ALIGN);
	__assume_aligned(Fz, ALIGN);

	int i,j,k,n;
	Real dFx, dFy, dFz, dtox, dtoy, dtoz;
	Real rho, E, Mx,My,Mz, P;

	const Real Gam = Model.Gam;
	//Use dt from argument instead of Grid for multi-step methods
	dtox = dt/Grid.dx; dtoy = dt/Grid.dy; dtoz = dt/Grid.dz;

 	#pragma omp parallel for collapse(2) \
 		default(shared) private(rho,E,Mx,My,Mz,P,dFx,dFy,dFz)
	for (k=Grid.ksd;k<=Grid.ked;k++) {
		for (j=Grid.jsd;j<=Grid.jed;j++) {
			#pragma omp simd
			for (i=Grid.isd;i<=Grid.ied;i++) {
				//In - Out, downward located fluxes

				//primitive -> conserved
				rho = fmax(Prim[DEN][k][j][i],TINY);

				//E = P/(Gam-1) + KinE
				E = ( fmax(Prim[PRESSURE][k][j][i],TINY) / (Gam-1) ) + 0.5*rho*( SQR(Prim[VELX][k][j][i]) + SQR(Prim[VELY][k][j][i]) + SQR(Prim[VELZ][k][j][i]) );
				Mx = rho*Prim[VELX][k][j][i];
				My = rho*Prim[VELY][k][j][i];
				Mz = rho*Prim[VELZ][k][j][i];

				//Update conserved
				rho += dtox*( Fx[DEN][k][j][i] - Fx[DEN][k][j][i+1] )
					+  dtoy*( Fy[DEN][k][j][i] - Fy[DEN][k][j+1][i] )
					+  dtoz*( Fz[DEN][k][j][i] - Fz[DEN][k+1][j][i] );

				Mx  += dtox*( Fx[MOMX][k][j][i] - Fx[MOMX][k][j][i+1] )
					+  dtoy*( Fy[MOMX][k][j][i] - Fy[MOMX][k][j+1][i] )
					+  dtoz*( Fz[MOMX][k][j][i] - Fz[MOMX][k+1][j][i] );

				My  += dtox*( Fx[MOMY][k][j][i] - Fx[MOMY][k][j][i+1] )
					+  dtoy*( Fy[MOMY][k][j][i] - Fy[MOMY][k][j+1][i] )
					+  dtoz*( Fz[MOMY][k][j][i] - Fz[MOMY][k+1][j][i] );

				E   += dtox*( Fx[TOTE][k][j][i] - Fx[TOTE][k][j][i+1] )
					+  dtoy*( Fy[TOTE][k][j][i] - Fy[TOTE][k][j+1][i] )
					+  dtoz*( Fz[TOTE][k][j][i] - Fz[TOTE][k+1][j][i] );

				//Store back into prim
				rho = fmax(rho,TINY);
				Prim[DEN][k][j][i] = rho;
				Prim[VELX][k][j][i] = Mx/rho;
				Prim[VELY][k][j][i] = My/rho;
				Prim[VELZ][k][j][i] = Mz/rho;
				P = (Gam-1)*( E - 0.5*rho*( SQR(Prim[VELX][k][j][i]) + SQR(Prim[VELY][k][j][i]) + SQR(Prim[VELZ][k][j][i]) ) );
				Prim[PRESSURE][k][j][i] =  fmax(P,TINY);

			}
		}
	}

}


//Initialize big data storage objects for fluxes/delta-state etc
//Can eventually make this a function pointer
void InitializeIntegrator(Grid_S Grid, Model_S Model) {

	//Holds advancing state
	advState = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);

}

//Clean up after yourself
void DestroyIntegrator(Grid_S Grid, Model_S Model) {

	Kill4Array(advState);
}

void CopyinBlock(RealP4 Q, Real Qblk[NVAR][NZBLK][NYBLK][NXBLK], Grid_S Grid, int iBlk, int jBlk, int kBlk) {

	int isB,ieB,jsB,jeB,ksB,keB;
	int nv,i,j,k;
	int iP,jP,kP;

	int Ng=Grid.Ng;

	//Copy into sub-block (include ghosts)
	isB = Grid.is + iBlk*NXPBLK;
	ieB = isB + NXPBLK;
	jsB = Grid.js + jBlk*NYPBLK;
	jeB = jsB + NYPBLK;
	ksB = Grid.ks + kBlk*NZPBLK;
	keB = ksB + NZPBLK;

	printf("\t(is,ie) = (%d,%d)\n",isB,ieB);
	printf("\t(js,je) = (%d,%d)\n",jsB,jeB);
	printf("\t(ks,ke) = (%d,%d)\n",ksB,keB);	

	for (nv=0;nv<NVAR;nv++) {
		for (k=0;k<NZBLK;k++) {
			for (j=0;j<NYBLK;j++) {
				for (i=0;i<NXBLK;i++) {

					kP = (ksB-Ng)+k;
					jP = (jsB-Ng)+j;
					iP = (isB-Ng)+i;
					Qblk[nv][k][j][i] = Q[nv][kP][jP][iP];
					Qblk[nv][k][j][i] = kBlk*BX*BY + jBlk*BX + iBlk;
				}
			}
		}
	} //4D block loop

}

void CopyoutBlock(RealP4 Q, Real Qblk[NVAR][NZBLK][NYBLK][NXBLK], Grid_S Grid, int iBlk, int jBlk, int kBlk) {
	int isB,ieB,jsB,jeB,ksB,keB;
	int nv,i,j,k;

	int Ng=Grid.Ng;

	//Copy out sub-block into main grid (don't include ghosts)
	isB = Grid.is + iBlk*NXPBLK;
	ieB = isB + NXPBLK;
	jsB = Grid.js + jBlk*NYPBLK;
	jeB = jsB + NYPBLK;
	ksB = Grid.ks + kBlk*NZPBLK;
	keB = ksB + NZPBLK;

	for (nv=0;nv<NVAR;nv++) {
		for (k=0;k<NZPBLK;k++) {
			for (j=0;j<NYPBLK;j++) {
				for (i=0;i<NXPBLK;i++) {
					Q[nv][ksB+k][jsB+j][isB+i] = Qblk[nv][k+Ng][j+Ng][i+Ng];
				}
			}
		}
	} //4D block loop

}

