#ifndef AMPS_H
#define AMPS_H
/* Defines overall data structures for the entire module

*/


//Define precision
typedef double Real;
typedef Real **** RealP4;


//#define VISCOSITY
//#define DEBUG
#define SHEARVISC 0.05
#define VISCSIG 0.8 //Safety factor for viscous timestep

#define PI 3.14159
//Define handles for Gas[NVAR][NZ][NY][NX]
#define NVAR 5
#define NUMGHOST 4
#define VECBUFF 16 //Buffer size for vector functions (MUSCL)
#define ALIGN 64
#define TINY 1.0e-6

//For now, simply define resolution/domain
#define NXP 128
#define NYP 64
#define NZP 64
#define NX (NXP+2*NUMGHOST)
#define NY (NYP+2*NUMGHOST)
#define NZ (NZP+2*NUMGHOST)

//Block decomposition information
//In reality, define N?PBLK @ compile-time, rest is derived at run-time

#define BX (2)
#define BY (4)
#define BZ (4)

#define NXPBLK (NXP/BX)
#define NYPBLK (NYP/BY)
#define NZPBLK (NZP/BZ)
#define NXBLK (NXPBLK+2*NUMGHOST)
#define NYBLK (NYPBLK+2*NUMGHOST)
#define NZBLK (NZPBLK+2*NUMGHOST)

//Space and time domain information
#define XMIN -1.0
#define XMAX 1.0
#define TFIN 10.0
#define TSOUT 10


typedef struct {

	int Nx,  Ny,  Nz, Ng, Nv;
	int Nxp, Nyp, Nzp;

	int is, ie, isd, ied;
	int js, je, jsd, jed;
	int ks, ke, ksd, ked;
	//Loop i=is;i<=ie is physical domain
	//Loop i=isd;i<=ied is entire domain (ie, include ghosts)

	int Ts;
	Real t, dt, dx, dy, dz;
	Real Xbd[2], Ybd[2], Zbd[2];
	Real Tfin;

	//These should be switched to allocated memory
	Real xc[NX];
	Real yc[NY];
	Real zc[NZ];

	Real xi[NX+1];
	Real yi[NY+1];
	Real zi[NZ+1];


	int N[4], Ni[4]; //4-vectors for centered/interface quantities
} Grid_S;

//Block data structure
//Similar to grid, but relies on known-size memory for offload simplicity
typedef struct {
	int Nx,  Ny,  Nz, Ng, Nv;
	int Nxp, Nyp, Nzp;

	//Loop i=is;i<=ie is physical domain
	//Loop i=isd;i<=ied is entire domain (ie, include ghosts)
	//This is local info
	int is, ie, isd, ied;
	int js, je, jsd, jed;
	int ks, ke, ksd, ked;

	//How to map from local->global
	int isG, ieG, isdG, iedG;
	int jsG, jeG, jsdG, jedG;
	int ksG, keG, ksdG, kedG;

	int Ts;
	Real t, dt, dx, dy, dz;
	Real Xbd[2], Ybd[2], Zbd[2];
	Real Tfin;

	Real xc[NXBLK];
	Real yc[NYBLK];
	Real zc[NZBLK];

	Real xi[NXBLK+1];
	Real yi[NYBLK+1];
	Real zi[NZBLK+1];

} Block_S;


typedef struct {
	Real C0; //Courant Number
	Real gVec[3]; //Gravity vectory
	int TsOut; //Output cadence for command line
	Real Gam; //For MUSCL equation of state
} Model_S;

//Cell-centered and interface-centered state blocks
typedef Real BlockCC[NVAR][NZBLK][NYBLK][NXBLK];
typedef Real BlockIC[NVAR][NZBLK+1][NYBLK+1][NXBLK+1];
//Block type for Riemann solver
typedef Real BlockR[NVAR][VECBUFF];

extern Model_S Model;

#define SQR(x) ((x)*(x))
#define IMIN(a,b) ((a) < (b) ? a : b)


//Prim vars (Main variable)
#define DEN 0
#define VELX 1
#define VELY 2
#define VELZ 3
#define PRESSURE 4

//Conserved vars (Secondary)
#define MOMX 1
#define MOMY 2
#define MOMZ 3
#define TOTE 4

#define DIR_X 0
#define DIR_Y 1
#define DIR_Z 2

#define NDIM 3

//Declare local array as aligned
#define DECALIGN __attribute__((aligned(ALIGN)))
//Direct compiler to assume alignment for passed array
#define ISALIGNED(x) __assume_aligned(x, ALIGN)

//Various defs for Intel Phi
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define FREE  alloc_if(0) free_if(1)
#define RETAIN alloc_if(1) free_if(0)
#define TPC 3 //Threads/core on device (1-4)
#define SBPDEV 4 //Sub-blocks per device

#ifdef DOPHI
	//Decore variable as offload-able
	#define MICTYPE __declspec(target(mic))
	#define CPMIC 60
#else
	//Define MICTYPE to empty string
	#define MICTYPE
#endif


#endif //AMPS_H