#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"
#include <omp.h>

#ifdef DOPHI
	#include <offload.h>
#endif

//Host variables to manage decomposition
RealP4 advState;
Block_S ***SubBlocks;

int NumDevs; //Number of devices
int TpC; //Threads per MIC core
int SBpDev; //Sub-blocks per device

int NumSBs; //Number of simultaneous sub-blocks
int TpSB; //Number of threads per sub-block

//Break up input 4D State into blocks and send each to advance routine
void BlockAdvance(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	RealP4 swpState;
	BlockCC Qblk DECALIGN;
	Block_S *myBlock;

	int iblk,jblk,kblk;
	int tID,devID; //Thread id/offload device #

	//printf("Inside block advance\n");
	int Ng = Grid.Ng;

	//Loop over blocks
	#pragma omp parallel for collapse(3) \
		num_threads(NumSBs) default(shared) \
		private(iblk,jblk,kblk,tID,devID,Qblk,myBlock)
	for (kblk=0;kblk<BZ;kblk++) {
		for (jblk=0;jblk<BY;jblk++) {
			for (iblk=0;iblk<BX;iblk++) {
				#pragma omp critical
				{
				
				tID = omp_get_thread_num();
				devID = (NumDevs>0) ? (tID % NumDevs) : 0;
				
				//printf("tID/devID = %d %d\n", tID,devID);

				myBlock = &(SubBlocks[kblk][jblk][iblk]);

				//Copy from State->Block
				CopyinBlock(State,Qblk,Grid,*myBlock);

				#ifdef DOPHI
				#pragma offload target(mic:devID) \
					inout(Qblk), in(myBlock : length(1)) \
					in(Model,Grid.dt)
				#endif
				{
					//Advance sub-block
					AdvanceFluid(Qblk,*myBlock,Model,Grid.dt);
				}
				//Copy advanced sub-block back into advState holder
				//Avoid ghosts
				CopyoutBlock(advState,Qblk,Grid,*myBlock);

				}//Critical
			}
		}
	} //Loop over blocks

	//printf("Finished block update, x-fer back\n");
	//Now swap advState and State pointers so that State is updated
	Copy4Array(State,advState,Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);
	
}

//Initialize big data storage objects for fluxes/delta-state etc
//Can eventually make this a function pointer
void InitializeIntegrator(Grid_S Grid, Model_S Model) {

	printf("Initializing Integrator ...\n");
	//Holds advancing state
	advState = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);

	//Note, BZ-BX can be variables, they don't need to be compile-time defined
	SubBlocks = MapBlocks(Grid,BX,BY,BZ);
	omp_set_nested(1); //Allow nested
	omp_set_max_active_levels(3); //Limit to 2-deep

#ifdef DOPHI
	NumDevs = _Offload_number_of_devices();
	//These are defs now, but will be changed to variables
	TpC = TPC;
	SBpDev = SBPDEV;
	NumSBs = NumDevs*SBpDev;

	TpSB = (CPMIC/SBpDev)*TpC; //Threads per sub-block
	printf("\tFound %d devices\n", NumDevs);
	printf("\tRunning %d Threads/Core and %d Sub-blocks/Device\n", TpC, SBpDev);
	printf("\tRunning %d Threads/Device\n",TpC*CPMIC);
#else
	NumDevs = 0;
	TpC = 0;
	SBpDev = 0;
	NumSBs = NUMSBS;
	TpSB = TPSB;
	printf("\tNo devices found, running on host\n");
	
#endif
	printf("\n");
	printf("\tSimultaneous # of blocks = %d\n",NumSBs);
	printf("\tNumber of Threads/Sub-Block = %d\n",TpSB);
}

//Clean up after yourself
void DestroyIntegrator(Grid_S Grid, Model_S Model) {

	Kill4Array(advState);

	//Kill SubBlocks array
	delete[] SubBlocks[0][0];
	delete[] SubBlocks[0];
	delete[] SubBlocks;
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
					//printf("Qblk[%d][%d][%d][%d] = %f\n", nv,k,j,i,Qblk[nv][Block.ks+k][Block.js+j][Block.is+i]);
					Q[nv][Block.ksG+k][Block.jsG+j][Block.isG+i] = Qblk[nv][Block.ks+k][Block.js+j][Block.is+i];
				}
			}
		}
	}

}

//Creates 3D block array and initializes it
//Note reverse order, Map[k][j][i]
Block_S*** MapBlocks(Grid_S Grid, int Bx, int By, int Bz) {
	int i,j,k;
	int ind1,ind2;

	printf("\tBreaking grid into %d [%d x %d x %d] blocks of size (%d,%d,%d)\n",Bx*By*Bz,Bx,By,Bz,NXPBLK,NYPBLK,NZPBLK);

	Block_S*** Map3D = new Block_S** [Bz];
	Map3D[0]         = new Block_S*  [Bz*By];
	Map3D[0][0]      = new Block_S   [BZ*BY*BX];

	for (k=0;k<Bz;k++) {
		ind1 = k*By;
		Map3D[k] = &Map3D[0][ind1];
		for (j=0;j<By;j++) {
			ind2 = k*(By*Bx) + j*Bx;
			Map3D[k][j] = &Map3D[0][0][ind2];
		}
	}

	for (k=0;k<Bz;k++) {
		for (j=0;j<By;j++) {
			for (i=0;i<Bx;i++) {
				InitBlock( &(Map3D[k][j][i]), Grid, i,j,k);
			}
		}
	}
	return Map3D;
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