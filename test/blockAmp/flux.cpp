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

	// WipeBlockIC(Fx,Block);
	// WipeBlockIC(Fy,Block);
	// WipeBlockIC(Fz,Block);

	printf("Flux Gam = %f\n",Model.Gam);
    LRs2Flux(W,W,Fx,DIR_X,Block);
    //LRs2Flux(W,W,Fy,DIR_Y,Block);
    //LRs2Flux(W,W,Fz,DIR_Z,Block);
	
}

//Takes *CELL* L/R values and a direction, returns flux
void LRs2Flux(BlockCC lW,BlockCC rW, BlockIC Flx, int d,  Block_S Grid) {
	ISALIGNED(lW);
	ISALIGNED(rW);
	ISALIGNED(Flx);

	int i,j,k,n,iG,iLim;
	int di,dj,dk;
	int Vn, Vt1, Vt2;


	int iblk;
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
				RiemannFluxHLLE(LeftW,RightW,FluxLR);
				//PrintBlockR(FluxLR);

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