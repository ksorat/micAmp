#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

//Globals internal to integrator
//Should be allocated on the phi

RealP4 Wmid, lW,rW; //Half timestep, CELL L/R
RealP4 Flux_x, Flux_y, Flux_z;

//MUSCL update
void AdvanceFluid(RealP4 State, Grid_S Grid, Model_S Model, Real dt) {

	//Update State->Wmid (0.5*dt) using PCM LR states

	//Get PCM-fluxes
	
	Flux_PCM(State, Flux_x, Flux_y, Flux_z, Grid, Model);
	Copy4Array(Wmid,State,Grid.N[0],Grid.N[1],Grid.N[2],Grid.N[3]);
	FluxUpdate(Wmid, Flux_x, Flux_y, Flux_z, 0.5*dt, Grid);

	//Get PLM-fluxes
	Flux_PLM(Wmid, Flux_x, Flux_y, Flux_z, Grid, Model);
	//Finish update, State->State (dt)
	FluxUpdate(State, Flux_x, Flux_y, Flux_z, dt, Grid);

#ifdef VISCOSITY
	//Apply viscosity in an operator-split manner for now
	
	Flux_Viscous(State, Flux_x, Flux_y, Flux_z, Grid, Model);
	FluxUpdate(State, Flux_x, Flux_y, Flux_z, dt, Grid);	
#endif

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

	//Wmid is the half timestep update to the input state
	Wmid = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);

	//Create holders for directed states
	lW = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);
	rW = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);

	Flux_x = Create4Array(Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
	Flux_y = Create4Array(Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
	Flux_z = Create4Array(Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);

}
//Clean up after yourself
void DestroyIntegrator(Grid_S Grid, Model_S Model) {

	Kill4Array(Wmid);
	Kill4Array(lW);
	Kill4Array(rW);

	Kill4Array(Flux_x);
	Kill4Array(Flux_y);
	Kill4Array(Flux_z);

}


