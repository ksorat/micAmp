#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

RealP4 advState;

//Break up input 4D State into blocks and send each to advance routine
void BlockAdvance(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	RealP4 swpState;
	BlockCC Qblk DECALIGN;
	Block_S Block;

	int iblk,jblk,kblk;
	

	//printf("Inside block advance\n");
	int Ng = Grid.Ng;
	//printf("Breaking grid into %d (%d,%d,%d) blocks\n",BX*BY*BZ,NXPBLK,NYPBLK,NZPBLK);

	//Loop over blocks
	for (kblk=0;kblk<BZ;kblk++) {
		for (jblk=0;jblk<BY;jblk++) {
			for (iblk=0;iblk<BX;iblk++) {
				
				InitBlock(&Block,Grid,iblk,jblk,kblk);
				
				//Copy from State->Block
				CopyinBlock(State,Qblk,Grid,Block);

				//Advance sub-block
				AdvanceFluid(Qblk,Block,Model,Grid.dt);

				//Copy advanced sub-block back into advState holder
				//Avoid ghosts
				CopyoutBlock(advState,Qblk,Grid,Block);
			}
		}
	} //Loop over blocks
	//printf("Finished block update, x-fer back\n");
	//Now swap advState and State pointers so that State is updated
	Copy4Array(State,advState,Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);
	//printf("Finished x-fer\n");
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


//Copy from Grid into sub-block (include ghosts)
void CopyinBlock(RealP4 Q, BlockCC Qblk, Grid_S Grid, Block_S Block) {
	int nv,i,j,k;
	int DelI,DelJ,DelK;

	DelK = Block.ked-Block.ksd;
	DelJ = Block.jed-Block.jsd;
	DelI = Block.ied-Block.isd;

	//printf("Del-(ijk) = %d, %d, %d\n",DelI,DelJ,DelK);
	for (nv=0;nv<NVAR;nv++) {
		for (k=0;k<=DelK;k++) {
			for (j=0;j<=DelJ;j++) {
				for (i=0;i<=DelI;i++) {	
					//printf("%d %d %d %d\n",nv,k,j,i);
					Qblk[nv][Block.ksd+k][Block.jsd+j][Block.isd+i] = Q[nv][Block.ksdG+k][Block.jsdG+j][Block.isdG+i];
				}
			}
		}
	}
	//printf("Done copy\n");
}

//Copy out sub-block into main grid (don't include ghosts)
void CopyoutBlock(RealP4 Q, BlockCC Qblk, Grid_S Grid, Block_S Block) {

	int nv,i,j,k;
	int DelI,DelJ,DelK;

	
	DelK = Block.ke-Block.ks;
	DelJ = Block.je-Block.js;
	DelI = Block.ie-Block.is;

	//printf("Copy out\n");
	for (nv=0;nv<NVAR;nv++) {
		for (k=0;k<=DelK;k++) {
			for (j=0;j<=DelJ;j++) {
				for (i=0;i<=DelI;i++) {	
					 Q[nv][Block.ksG+k][Block.jsG+j][Block.isG+i] = Qblk[nv][Block.ks+k][Block.js+j][Block.is+i];
				}
			}
		}
	}

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
	Block->dx = Grid.dx; Block->dy = Grid.dy; Block->dz = Grid.dz;

	Block->Ts = Grid.Ts;

	//Set ?Bds, x/y/z-c/i
	Block->isd = 0;               Block->ied = NXBLK-1;
	Block->is  = Block->isd + Ng; Block->ie  = Block->ied - Ng; 

	Block->jsd = 0;               Block->jed = NYBLK-1;
	Block->js  = Block->jsd + Ng; Block->je  = Block->jed - Ng; 

	Block->ksd = 0;               Block->ked = NZBLK-1;
	Block->ks  = Block->ksd + Ng; Block->ke  = Block->ked - Ng; 

	//Set global bounds for reverse map
	Block->isG  = Grid.is + iBlk*NXPBLK; Block->ieG  = Block->isG + NXPBLK-1;
	Block->isdG = Block->isG - Ng;       Block->iedG = Block->ieG + Ng;

	Block->jsG  = Grid.js + jBlk*NYPBLK; Block->jeG  = Block->jsG + NYPBLK-1;
	Block->jsdG = Block->jsG - Ng;       Block->jedG = Block->jeG + Ng;

	Block->ksG  = Grid.ks + kBlk*NZPBLK; Block->keG  = Block->ksG + NZPBLK-1;
	Block->ksdG = Block->ksG - Ng;       Block->kedG = Block->keG + Ng;


	// printf("Block (%d,%d,%d)\n",iBlk,jBlk,kBlk);
	// printf("\tLocal <-> Global\n");
	// printf("\tX: [%d,%d] <-> [%d,%d]\n",Block->is,Block->ie,Block->isG,Block->ieG);
	// printf("\tY: [%d,%d] <-> [%d,%d]\n",Block->js,Block->je,Block->jsG,Block->jeG);
	// printf("\tZ: [%d,%d] <-> [%d,%d]\n",Block->ks,Block->ke,Block->ksG,Block->keG);
	// printf("\t\tXd: [%d,%d] <-> [%d,%d]\n",Block->isd,Block->ied,Block->isdG,Block->iedG);
	// printf("\t\tYd: [%d,%d] <-> [%d,%d]\n",Block->jsd,Block->jed,Block->jsdG,Block->jedG);
	// printf("\t\tZd: [%d,%d] <-> [%d,%d]\n",Block->ksd,Block->ked,Block->ksdG,Block->kedG);

}