//Test of memory allocation/transfer to MIC

//Allocate contiguous 2D array, x-fer to MIC, do simple calculation and bring back
//For SuperMIC
//
#define DIMX 4
#define DIMY 3
#define DIMZ 5
#define DIMV 2

#define ALIGN 64

#define MIC0 mic:0
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define FREE  alloc_if(0) free_if(1)
#define RETAIN alloc_if(1) free_if(0)

#define MICTYPE __declspec(target(mic))

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double Real;
typedef Real **** RealP4;

MICTYPE RealP4 Create4Array(int N1, int N2, int N3, int N4);
MICTYPE RealP4 Map4Array(Real *Start, int N1, int N2, int N3, int N4);
MICTYPE void Kill4Array(RealP4 ToDie);

MICTYPE RealP4 Data;
MICTYPE RealP4 DataPhi;
MICTYPE Real *Start;
MICTYPE Real *StartPhi;

int main() {
	int i,j,k,n, nCum, Ntot;
	Real s, c, cumsum;
	//Setup

	Ntot = DIMX*DIMY*DIMZ*DIMV;

	printf("Alloc on host\n");
	Data = Create4Array(DIMV,DIMZ,DIMY,DIMX);
	Start = &(Data[0][0][0][0]);
	printf("Finish alloc on host\n");

	printf("Start alloc on MIC\n");
	#pragma offload target(MIC0) nocopy(DataPhi : REUSE) nocopy(StartPhi : length(Ntot) RETAIN)
	{
		DataPhi = Map4Array(StartPhi,DIMV,DIMZ,DIMY,DIMX);
		DataPhi[0][0][0][0] = -1.0;

	}
	printf("Finish alloc on MIC\n");


	printf("Initialize data on host\n");
	nCum = 0;
	for (n=0;n<DIMV;n++) {
		for (k=0;k<DIMZ;k++) {
			for (j=0;j<DIMY;j++) {
				for (i=0;i<DIMX;i++) {
					//Data[n][k][j][i] = 1.0*n*k*j*i;
					Data[n][k][j][i] = 1.0*nCum;
					nCum++;
				}
			}
		}
	}
	printf("Finish initialization on host\n");

	printf("Begin transfer to device\n");
	#pragma offload_transfer target(MIC0) in( Start : into(StartPhi) length(Ntot) REUSE )
	printf("End transfer to device\n");

	printf("Begin computation on device\n");
	#pragma offload target(MIC0) nocopy(DataPhi : REUSE) nocopy(StartPhi : REUSE)
	{
		for (n=0;n<DIMV;n++) {
			for (k=0;k<DIMZ;k++) {
				for (j=0;j<DIMY;j++) {
					#pragma omp simd
					for (i=0;i<DIMX;i++) {
						s = sin(DataPhi[n][k][j][i]);
						c = cos(DataPhi[n][k][j][i]);
						DataPhi[n][k][j][i] = s*s + c*c;
					}
				}
			}	
		}
	}
	printf("End computation on device\n");
	#pragma offload_transfer target(MIC0) out( StartPhi : into(Start) length(Ntot) FREE)


	// //Check sum
	cumsum = 0;
	for (n=0;n<DIMV;n++) {
		for (k=0;k<DIMZ;k++) {
			for (j=0;j<DIMY;j++) {
				for (i=0;i<DIMX;i++) {
					cumsum += Data[n][k][j][i];
				}
			}
		}
	}

	printf("Checksum = %e\n", cumsum/(DIMX*DIMY*DIMZ*DIMV) - 1.0);
	Kill4Array(Data);
	return 0;
}



//Input: Start, points to contiguous allocated memory
//Function creates 4D structure w/ mapped data
RealP4 Map4Array(Real *Start, int N1, int N2, int N3, int N4) {
	int n,i,j;
	int ind1,ind2,ind3;

	RealP4 Data4D   = (Real ****) _mm_malloc(sizeof(Real****)*N1,       ALIGN);
	Data4D[0]       = (Real *** ) _mm_malloc(sizeof(Real***)*N1*N2,     ALIGN);
	Data4D[0][0]    = (Real **  ) _mm_malloc(sizeof(Real**)*N1*N2*N3,   ALIGN);
	Data4D[0][0][0] = (Real *   ) Start;

	for(n=0;n<N1;n++) {
		ind1 = n*N2;
		Data4D[n] = &Data4D[0][ind1];
		for(i=0;i<N2;i++) {
			ind2 = n*(N2*N3) + i*N3;
			Data4D[n][i] = &Data4D[0][0][ind2];
			for(j=0;j<N3;j++) {
				ind3 = n*(N2*N3*N4) + i*(N3*N4) + j*N4;
				Data4D[n][i][j] = &Data4D[0][0][0][ind3];

			}
		}
	}
	return Data4D;
}

//Creates array of dimension Data4D[N1][N2][N3][N4]
RealP4 Create4Array(int N1, int N2, int N3, int N4) {
	Real *Start;
	RealP4 Data4D;

	//Allocate memory
	Start = (Real *) _mm_malloc(sizeof(Real*)*N1*N2*N3*N4, ALIGN);
	Data4D = Map4Array(Start,N1,N2,N3,N4);

	return Data4D;
}

void Kill4Array(RealP4 ToDie) {
	_mm_free (***ToDie);
	_mm_free (**ToDie);
	_mm_free (*ToDie);
	_mm_free (ToDie);
}
