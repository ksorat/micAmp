//Various routines to create/kill contiguous dynamic arrays

#include <stdio.h>
#include <stdlib.h>

#include "amps.h"
#include "prototypes.h"

//Creates array of dimension Data4D[N1][N2][N3][N4]
RealP4 Create4Array(int N1, int N2, int N3, int N4) {

	int n,i,j;
	int ind1, ind2, ind3;

	RealP4 Data4D = (Real ****)_mm_malloc(sizeof(Real****)*N1, ALIGN);
    Data4D[0] = (Real ***)_mm_malloc(sizeof(Real***)*N1*N2, ALIGN);
    Data4D[0][0] = (Real **)_mm_malloc(sizeof(Real**)*N1*N2*N3, ALIGN);
    Data4D[0][0][0] = (Real *)_mm_malloc(sizeof(Real*)*N1*N2*N3*N4, ALIGN);	

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

void Kill4Array(RealP4 ToDie) {
    _mm_free(ToDie[0][0][0]);
    _mm_free(ToDie[0][0]);
    _mm_free(ToDie[0]);
    _mm_free(ToDie);	
}

void Wipe4Array(RealP4 Ar, int N1, int N2, int N3, int N4) {
	__assume_aligned(Ar, ALIGN);
	int i,j,k,l;
	
	#pragma omp parallel for collapse(3)
	for (i=0;i<N1;i++) {
		for (j=0;j<N2;j++) {
			for (k=0;k<N3;k++) {
				#pragma omp simd
				for (l=0;l<N4;l++) {
					Ar[i][j][k][l] = 0.0;
				}
			}
		}
	}

}

//A = B
//Consider using memcpy, although compiler should do the transformation automagically
void Copy4Array(RealP4 A, RealP4 B, int N1, int N2, int N3, int N4) {
	__assume_aligned(A, ALIGN);
	__assume_aligned(B, ALIGN);	
	int i,j,k,l;

	#pragma omp parallel for collapse(3)
	for (i=0;i<N1;i++) {
		for (j=0;j<N2;j++) {
			for (k=0;k<N3;k++) {
				#pragma omp simd
				for (l=0;l<N4;l++) {
					A[i][j][k][l] = B[i][j][k][l];
				}
			}
		}
	}
}


//Zeros out block
void WipeBlockCC(BlockCC A, Block_S Block) {
	int n,i,j,k;
	#pragma omp parallel for collapse(3) num_threads(TpSB)
	for (n=0;n<NVAR;n++) {
		for (k=Block.ksd;k<=Block.ked;k++) {
			for (j=Block.jsd;j<=Block.jed;j++) {
				#pragma omp simd
				for (i=Block.isd;i<=Block.ied;i++) {
					A[n][k][j][i] = 0.0;
				}
			}
		}
	} //Block loop
}

//Zeros out block
void WipeBlockIC(BlockIC A, Block_S Block) {
	int n,i,j,k;
	#pragma omp parallel for collapse(3) num_threads(TpSB)
	for (n=0;n<NVAR;n++) {
		for (k=Block.ksd;k<=Block.ked+1;k++) {
			for (j=Block.jsd;j<=Block.jed+1;j++) {
				#pragma omp simd
				for (i=Block.isd;i<=Block.ied+1;i++) {
					A[n][k][j][i] = 0.0;
				}
			}
		}
	} //Block loop
}

//Prints out block
void PrintBlockCC(BlockCC A, Block_S Block) {
	int n,i,j,k;
	for (n=0;n<NVAR;n++) {
		for (k=Block.ksd;k<=Block.ked;k++) {
			for (j=Block.jsd;j<=Block.jed;j++) {
				for (i=Block.isd;i<=Block.ied;i++) {
					printf("A[%d][[%d][%d][%d] = %f\n", n,k,j,i,A[n][k][j][i]);
				}
			}
		}
	} //Block loop
}

