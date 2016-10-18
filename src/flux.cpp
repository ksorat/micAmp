#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

//Given input state W (prim), calculates fluxes based on PCM LR reconstruction
void Flux_PCM(BlockCC W,BlockIC Fx,BlockIC Fy,BlockIC Fz,Block_S Block,Model_S Model) {
	ISALIGNED(W);
	ISALIGNED(Fx);
	ISALIGNED(Fy);
	ISALIGNED(Fz);

    LRs2Flux(W,W,Fx,DIR_X,Block,Model);
    LRs2Flux(W,W,Fy,DIR_Y,Block,Model);
    LRs2Flux(W,W,Fz,DIR_Z,Block,Model);
	
}

//Given input state W (prim), calculates fluxes based on PLM LR reconstruction
//Using MC limiter
void Flux_PLM(BlockCC W,BlockIC Fx,BlockIC Fy,BlockIC Fz,Block_S Block,Model_S Model) {
	ISALIGNED(W);
	ISALIGNED(Fx);
	ISALIGNED(Fy);
	ISALIGNED(Fz);

	BlockCC lW, rW DECALIGN;

	Real Wmp[2]; //Holds mW,pW (minus/positive interface values of cell)

	int d,n,j,k,i;
	int di,dj,dk;
	for (d=0;d<NDIM;d++) {
		di = (d == DIR_X) ? 1 : 0;
		dj = (d == DIR_Y) ? 1 : 0;
		dk = (d == DIR_Z) ? 1 : 0;

		//Loop over grid
		for (n=0;n<Block.Nv;n++) {
			for (k=Block.ksd+1;k<=Block.ked-1;k++) {
				for (j=Block.jsd+1;j<=Block.jed-1;j++) {
					for (i=Block.isd+1;i<=Block.ied-1;i++) {

						//Reconstruct profile to get *CELL* LR states
						Recon_PLM(Wmp, W[n][k-dk][j-dj][i-di],W[n][k][j][i],W[n][k+dk][j+dj][i+di]);

						lW[n][k][j][i] = Wmp[0];
						rW[n][k][j][i] = Wmp[1];
					}
				}
			}
		} //n loop

		//Done with *CELL* LR states, calculate fluxes
		//Using switch for now, figure out better pointer method
		switch(d) {
			case DIR_X :
				LRs2Flux(lW,rW,Fx,d,Block,Model);	
				break;
			case DIR_Y :
				LRs2Flux(lW,rW,Fy,d,Block,Model);	
				break;
			case DIR_Z :
				LRs2Flux(lW,rW,Fz,d,Block,Model);	
				break;
		}
		
	} //Direction loop
}

//Takes *CELL* L/R values and a direction, returns flux
void LRs2Flux(BlockCC lW,BlockCC rW, BlockIC Flx, int d,  Block_S Grid, Model_S Model) {
	ISALIGNED(lW);
	ISALIGNED(rW);
	ISALIGNED(Flx);

	int i,j,k,n,iG,iLim;
	int di,dj,dk;
	int Vn, Vt1, Vt2;


	int iblk;
	Real Gam = Model.Gam;
	
	//These hold *INTERFACE* L/R values
	BlockR LeftW, RightW, FluxLR DECALIGN;

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
	WipeBlockIC(Flx,Grid);
	//Asymmetric bounds b/c of flux centering
 	#pragma omp parallel for collapse(2) \
 		num_threads(TpSB) \
 		default(shared) private(i,j,k,iLim,iG,iblk,LeftW,RightW,FluxLR)
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
				RiemannFluxHLLE(LeftW,RightW,FluxLR,Gam);
				//RiemannFluxHLLC(LeftW,RightW,FluxLR,Gam);

				//Unpack into fluxes
				//Untwist back to original coordinate system
				//Untwist, Vn<->Vx etc
				#pragma omp simd 
				for (i=0;i<iLim;i++) {
					iG = iblk+i;
					Flx[DEN][k][j][iG]      = FluxLR[DEN][i];
					Flx[Vn][k][j][iG]       = FluxLR[VELX][i];
					Flx[Vt1][k][j][iG]      = FluxLR[VELY][i];
					Flx[Vt2][k][j][iG]      = FluxLR[VELZ][i];
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
