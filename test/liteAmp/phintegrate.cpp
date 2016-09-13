#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

MICTYPE RealP4 Flux_x, Flux_y, Flux_z;

//MUSCL update
void AdvanceFluid(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	//1-Step
	FluxUpdate(State,Flux_x,Flux_y,Flux_z,dt,Grid);

}

void FluxUpdate(RealP4 Prim, RealP4 Fx, RealP4 Fy, RealP4 Fz, Real dt, Grid_S Grid) {
	__assume_aligned(Prim, ALIGN);
	__assume_aligned(Fx, ALIGN);
	__assume_aligned(Fy, ALIGN);
	__assume_aligned(Fz, ALIGN);

	int i,j,k,n;
	Real dFx, dFy, dFz, dtox, dtoy, dtoz;
	Real rho, E, Mx,My,Mz, P;

	const Real Gam = Model.Gam;
	//Use dt from argument instead of Grid for multi-step methods
	dtox = dt/Grid.dx; dtoy = dt/Grid.dy; dtoz = dt/Grid.dz;

 	#pragma omp parallel for collapse(2) \
 		default(shared) private(rho,E,Mx,My,Mz,P,dFx,dFy,dFz)
	for (k=Grid.ksd;k<=Grid.ked;k++) {
		for (j=Grid.jsd;j<=Grid.jed;j++) {
			#pragma omp simd
			for (i=Grid.isd;i<=Grid.ied;i++) {
				//In - Out, downward located fluxes

				//primitive -> conserved
				rho = fmax(Prim[DEN][k][j][i],TINY);

				//E = P/(Gam-1) + KinE
				E = ( fmax(Prim[PRESSURE][k][j][i],TINY) / (Gam-1) ) + 0.5*rho*( SQR(Prim[VELX][k][j][i]) + SQR(Prim[VELY][k][j][i]) + SQR(Prim[VELZ][k][j][i]) );
				Mx = rho*Prim[VELX][k][j][i];
				My = rho*Prim[VELY][k][j][i];
				Mz = rho*Prim[VELZ][k][j][i];

				//Update conserved
				rho += dtox*( Fx[DEN][k][j][i] - Fx[DEN][k][j][i+1] )
					+  dtoy*( Fy[DEN][k][j][i] - Fy[DEN][k][j+1][i] )
					+  dtoz*( Fz[DEN][k][j][i] - Fz[DEN][k+1][j][i] );

				Mx  += dtox*( Fx[MOMX][k][j][i] - Fx[MOMX][k][j][i+1] )
					+  dtoy*( Fy[MOMX][k][j][i] - Fy[MOMX][k][j+1][i] )
					+  dtoz*( Fz[MOMX][k][j][i] - Fz[MOMX][k+1][j][i] );

				My  += dtox*( Fx[MOMY][k][j][i] - Fx[MOMY][k][j][i+1] )
					+  dtoy*( Fy[MOMY][k][j][i] - Fy[MOMY][k][j+1][i] )
					+  dtoz*( Fz[MOMY][k][j][i] - Fz[MOMY][k+1][j][i] );

				E   += dtox*( Fx[TOTE][k][j][i] - Fx[TOTE][k][j][i+1] )
					+  dtoy*( Fy[TOTE][k][j][i] - Fy[TOTE][k][j+1][i] )
					+  dtoz*( Fz[TOTE][k][j][i] - Fz[TOTE][k+1][j][i] );

				//Store back into prim
				rho = fmax(rho,TINY);
				Prim[DEN][k][j][i] = rho;
				Prim[VELX][k][j][i] = Mx/rho;
				Prim[VELY][k][j][i] = My/rho;
				Prim[VELZ][k][j][i] = Mz/rho;
				P = (Gam-1)*( E - 0.5*rho*( SQR(Prim[VELX][k][j][i]) + SQR(Prim[VELY][k][j][i]) + SQR(Prim[VELZ][k][j][i]) ) );
				Prim[PRESSURE][k][j][i] =  fmax(P,TINY);

			}
		}
	}

}


//Initialize big data storage objects for fluxes/delta-state etc
//Can eventually make this a function pointer
void InitializeIntegrator(Grid_S Grid, Model_S Model) {

	int Ntot 
	int m=0;
	
	RealP4 FxH;
	RealP4 FyH;
	RealP4 FzH;
	Real *FxH0, *FyH0, *FzH0;
	Real *Fx0, *Fy0, *Fz0;

	//Do all allocating on card, but can't use uninitialized pointers
	//Stupid hack, allocate 4d on host, mirror on device, then kill host
	Ntot = Grid.Nv*(Grid.Nz+1)*(Grid.Ny+1)*(Grid.Nx+1);
	FxH = Create4Array(Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
	FyH = Create4Array(Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
	FzH = Create4Array(Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);

	Fx0 = &(FxH[0][0][0][0]);
	Fy0 = &(FyH[0][0][0][0]);
	Fz0 = &(FzH[0][0][0][0]);

	#pragma offload target(mic:m) \
		nocopy(Flux_x:REUSE) nocopy(Fx0:length(Ntot) ALLOC) \
		nocopy(Flux_y:REUSE) nocopy(Fy0:length(Ntot) ALLOC) \
		nocopy(Flux_z:REUSE) nocopy(Fz0:length(Ntot) ALLOC) 
	{
		Flux_x = Map4Array(Fx0,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
		Flux_y = Map4Array(Fy0,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
		Flux_z = Map4Array(Fz0,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);

		//Initialize
		Wipe4Array(Flux_x,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
		Wipe4Array(Flux_y,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
		Wipe4Array(Flux_z,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
	}

}

//Clean up after yourself
void DestroyIntegrator(Grid_S Grid, Model_S Model) {
	int m=0;
	
	//Clean up arrays on card
	#pragma offload target(mic:m) \
		nocopy(Flux_x,Flux_y,Flux_z)
	{
		Kill4Array(Flux_x);
		Kill4Array(Flux_y);
		Kill4Array(Flux_z);
	}
}


