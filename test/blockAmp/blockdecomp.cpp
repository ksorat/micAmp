#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

RealP4 advState;

//Break up input 4D State into blocks and send each to advance routine
void BlockAdvance(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	RealP4 swpState;
	Real Qblk[NVAR][NZBLK][NYBLK][NXBLK];

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

//Initializes Block data structure for a given block index
void InitBlock(Block_S *Block, Grid_S Grid, int iBlk, int jBlk, int kBlk) {

	int Ng=Grid.Ng; //Use variable for generality, ie varying # of ghosts
	Block->Nv = NVAR;
	Block->Ng = Ng;

	Block->Nxp = NXPBLK;
	Block->Nyp = NYPBLK;
	Block->Nzp = NZPBLK;

	Block->Nx = NXBLK + 2*Ng;
	Block->Ny = NYBLK + 2*Ng;
	Block->Nz = NZBLK + 2*Ng;

	Block->t = Grid.t;
	Block->dt = Grid.dt;
	Block->Ts = Grid.Ts;

	//Set ?Bds, x/y/z-c/i
	Block->isd = 0;               Block->ied = NXBLK-1;
	Block->is  = Block->isd + Ng; Block->ie  = Block->ied - Ng; 

	Block->jsd = 0;               Block->jed = NYBLK-1;
	Block->js  = Block->jsd + Ng; Block->je  = Block->jed - Ng; 

	Block->ksd = 0;               Block->ked = NZBLK-1;
	Block->ks  = Block->ksd + Ng; Block->ke  = Block->ked - Ng; 

}