//Test of memory allocation/transfer to MIC

//Allocate contiguous 2D array, x-fer to MIC, do simple calculation and bring back

#define DIMX 1000
#define DIMY 1000

#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define FREE  alloc_if(0) free_if(1)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double Real;

Real **Create2Array(int N1, int N2);
void Kill2Array(Real **ToDie);

int main() {
	Real **Data;
	int i,j;
	Real s, c, cumsum;
	//Setup
	Data = Create2Array(DIMX,DIMY);
	for (i=0;i<DIMX;i++) {
		for (j=0;j<DIMY;j++) {
			Data[i][j] = 1.0*i*j;
		}
	}

	//Xfer and calculate
	#pragma offload target(mic:0) inout(Data[0:DIMX][0:DIMY])
	{
		#pragma omp parallel for private(s,c)
		for (i=0;i<DIMX;i++) {
			#pragma omp simd
			for (j=0;j<DIMY;j++) {
				s = sin(Data[i][j]);
				c = cos(Data[i][j]);
	
				Data[i][j] = s*s + c*c;
			}
		}
	} //Data region

	//Check sum
	cumsum = 0;
	for (i=0;i<DIMX;i++) {
		for (j=0;j<DIMY;j++) {
			cumsum += Data[i][j];
		}
	}

	printf("Checksum = %f\n", cumsum/(DIMX*DIMY));
	Kill2Array(Data);
	return 0;
}

Real **Create2Array(int N1, int N2) {
	int i;
	int ind1;

	Real **Data2D = new Real*[N1];

	Data2D[0] = new Real [N1*N2];
	for (i=0;i<N1;i++) {
		ind1 = i*N2;
		Data2D[i] = &Data2D[0][ind1];
	}
	return Data2D;
}

void Kill2Array(Real **ToDie) {
	delete[] ToDie[0];
	delete[] ToDie;	
}