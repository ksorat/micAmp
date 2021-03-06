#ifndef AMPS_H
#define AMPS_H
/* Defines overall data structures for the entire module

*/


//Define precision
typedef double Real;
typedef Real **** RealP4;
#define MICTYPE __declspec(target(mic))

//#define VISCOSITY
//#define DEBUG
#define SHEARVISC 0.05
#define VISCSIG 0.8 //Safety factor for viscous timestep

#define PI 3.14159
//Define handles for Gas[NVAR][NZ][NY][NX]
#define NVAR 5
#define NUMGHOST 4
#define VECBUFF 16 //Buffer size for vector functions (MUSCL)
#define VECBUFFV 16 //Buffer size for vector functions (Viscous)
#define ALIGN 32
#define TINY 1.0e-6

//For now, simply define resolution/domain
#define NXP 64
#define NYP 64
#define NZP 32
#define NX NXP+2*NUMGHOST
#define NY NYP+2*NUMGHOST
#define NZ NZP+2*NUMGHOST

#define XMIN -1.0
#define XMAX 1.0
#define TFIN 10.0
#define TSOUT 1


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

	Real xc[NX+2*NUMGHOST];
	Real yc[NY+2*NUMGHOST];
	Real zc[NZ+2*NUMGHOST];

	Real xi[NX+2*NUMGHOST+1];
	Real yi[NY+2*NUMGHOST+1];
	Real zi[NZ+2*NUMGHOST+1];


	int N[4], Ni[4]; //4-vectors for centered/interface quantities
} Grid_S;

typedef struct {
	Real C0; //Courant Number
	Real gVec[3]; //Gravity vectory
	int TsOut; //Output cadence for command line
	Real Gam; //For MUSCL equation of state
} Model_S;

extern MICTYPE Model_S Model;

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

//Various defs for Intel Phi
#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define FREE  alloc_if(0) free_if(1)
#define RETAIN alloc_if(1) free_if(0)

#endif //AMPS_H