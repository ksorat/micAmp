#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

extern RealP4 lW, rW;

//Given input state W, calculates fluxes based on PCM LR reconstruction
void Flux_PCM(RealP4 W, RealP4 Fx, RealP4 Fy, RealP4 Fz, Grid_S Grid, Model_S Model) {

	//Calculate fluxes, PCM means lW=rW=W
	LRs2Flux(W,W,Fx,DIR_X,Grid);
	LRs2Flux(W,W,Fy,DIR_Y,Grid);
	LRs2Flux(W,W,Fz,DIR_Z,Grid);

}

//Given input state W, calculates fluxes based on PLM LR reconstruction
//Uses MC limiter, but simple enough to try other options (maybe superbee?)
//Assumes global pointers lW,rW in recon.h have been allocated in integrator initialization

void Flux_PLM(RealP4 W, RealP4 Fx, RealP4 Fy, RealP4 Fz, Grid_S Grid, Model_S Model) {
	__assume_aligned(W, ALIGN);
	__assume_aligned(Fx, ALIGN);
	__assume_aligned(Fy, ALIGN);
	__assume_aligned(Fz, ALIGN);

	int i,j,k,n;
	int d,di,dj,dk;
	Real Wmp[2]; //Holds mW,pW (minus/positive interface values of cell)
	RealP4 Fi; //Pointer to flux for current dimension

	for (d=0;d<NDIM;d++) {
		switch(d) {
			//Setup offsets for directionality
			case DIR_X :
				//X-Dir
				di = 1; dj = 0; dk = 0;
				Fi = Fx;
				break;
			case DIR_Y :
				//Y-Dir
				di = 0; dj = 1; dk = 0;
				Fi = Fy;
				break;
			case DIR_Z :
				//Z-Dir
				di = 0; dj = 0; dk = 1;
				Fi = Fz;
				break;		
		}
		//Big loop, collapse nk
		#pragma omp parallel for collapse(3) private(Wmp)
		for (n=0;n<Grid.Nv;n++) {
			for (k=Grid.ksd+1;k<=Grid.ked-1;k++) {
				for (j=Grid.jsd+1;j<=Grid.jed-1;j++) {
					#pragma omp simd
					for (i=Grid.isd+1;i<=Grid.ied-1;i++) {
						//Reconstruct profile to get *CELL* LR states
						Recon_PLM(Wmp, W[n][k-dk][j-dj][i-di],W[n][k][j][i],W[n][k+dk][j+dj][i+di]);
						lW[n][k][j][i] = Wmp[0];
						rW[n][k][j][i] = Wmp[1];
					}
				}
			}
		}
		
		//Done with cell LR states for given direction, calculate flux
		LRs2Flux(lW,rW,Fi,d,Grid);
	}

}

//Takes *CELL* L/R values and a direction, returns flux
void LRs2Flux(RealP4 lW, RealP4 rW, RealP4 Flx, int d,  Grid_S Grid) {
	__assume_aligned(lW, ALIGN);
	__assume_aligned(rW, ALIGN);
	__assume_aligned(Flx, ALIGN);
	
	int i,j,k,n,iG,iLim;
	int di,dj,dk;
	int Vn, Vt1, Vt2;


	int iblk;
	//These hold *INTERFACE* L/R values
	Real LeftW[NVAR][VECBUFF], RightW[NVAR][VECBUFF], FluxLR[NVAR][VECBUFF] __attribute__((aligned(ALIGN)));

	switch(d) {
		//Setup offsets for directionality
		//Setup rotated triad for twist/untwist
		case DIR_X :
			//X-Dir
			di = 1; dj = 0; dk = 0;
			Vn = VELX; Vt1 = VELY; Vt2 = VELZ; 
			break;
		case DIR_Y :
			//Y-Dir
			di = 0; dj = 1; dk = 0;
			Vn = VELY; Vt1 = VELZ; Vt2 = VELX;	
			break;
		case DIR_Z :
			//Z-Dir
			di = 0; dj = 0; dk = 1;
			Vn = VELZ; Vt1 = VELX; Vt2 = VELY;
			break;		
	}
	//Wipe flux array
	Wipe4Array(Flx,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);

	//Asymmetric bounds b/c of flux centering
 	#pragma omp parallel for collapse(3) \
 		default(shared) private(i,iLim,iG,iblk,LeftW,RightW,FluxLR)
	for (k=Grid.ksd+1;k<=Grid.ked;k++) {
		for (j=Grid.jsd+1;j<=Grid.jed;j++) {
			for (iblk=Grid.isd+1;iblk<=Grid.ied;iblk+=VECBUFF) {
				//Inner vector loops

				//Limit inner loops in case of bad divisibility
				//Last i-bound is ied
				iLim = IMIN(VECBUFF,Grid.ied-iblk+1);

				//Loads data into aligned LR buffer, hands off to vector Riemann solver
				//LR buffer is interface LR not cell LR
				//"Twist" coordinate triad so that x is in direction of interface normal
				#pragma omp simd 
				for (i=0;i<iLim;i++) {
					iG = iblk+i;
					LeftW [DEN][i]      = rW[DEN][k-dk][j-dj][iG-di];
					RightW[DEN][i]      = lW[DEN][k][j][iG];
					//Twist, Vx<->Vnormal / Vy,Vz=Vt1,Vt2
					LeftW [VELX][i]     = rW[Vn][k-dk][j-dj][iG-di];
					RightW[VELX][i]     = lW[Vn][k][j][iG];
					LeftW [VELY][i]     = rW[Vt1][k-dk][j-dj][iG-di];
					RightW[VELY][i]     = lW[Vt1][k][j][iG];
					LeftW [VELZ][i]     = rW[Vt2][k-dk][j-dj][iG-di];
					RightW[VELZ][i]     = lW[Vt2][k][j][iG];
					
					LeftW [PRESSURE][i] = rW[PRESSURE][k-dk][j-dj][iG-di];
					RightW[PRESSURE][i] = lW[PRESSURE][k][j][iG];

				}
				//Call Riemann solver
				RiemannFluxHLLE(LeftW,RightW,FluxLR);

				//Unpack into fluxes
				//Untwist back to original coordinate system
				//Untwist, Vn<->Vx etc
				#pragma omp simd 
				for (i=0;i<iLim;i++) {
					iG = iblk+i;
					Flx[DEN][k][j][iG] = FluxLR[DEN][i];
					Flx[Vn][k][j][iG] = FluxLR[VELX][i];
					Flx[Vt1][k][j][iG] = FluxLR[VELY][i];
					Flx[Vt2][k][j][iG] = FluxLR[VELZ][i];
					Flx[PRESSURE][k][j][iG] = FluxLR[PRESSURE][i];
				}		

			}
		}
	}

}


inline void Recon_PLM(Real Wmp[2], Real Wl, Real Wc, Real Wr) {
	Real DelWl, DelWc, DelWr, DelWm, sgn;

	//Calculate different local slopes
	DelWl = Wc-Wl;
	DelWr = Wr-Wc;
	DelWc = 0.5*( Wr-Wl );

	//Limit slope
	DelWm = SlopeLimit(DelWl,DelWr,DelWc);

	//Calculate directed states using limited slope
	Wmp[0] = Wc - 0.5*DelWm;
	Wmp[1] = Wc + 0.5*DelWm;

}

inline Real SlopeLimit(Real dWl, Real dWr,Real dWc) {
	Real dWm, sgn;

	if (dWl*dWr <= 0) {
		dWm = 0.0;
	} else {
		sgn = (dWl > 0) - (dWl < 0); //Cheap sgn function w/o branch
		dWm = sgn*fmin( fabs(dWc),fmin(2.0*fabs(dWl),2.0*fabs(dWr)) );
	}
	return dWm;
}