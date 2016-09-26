#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

//MUSCL update
void AdvanceFluid(BlockCC State, Block_S Block, Model_S Model, Real dt) {
	BlockIC Flux_x, Flux_y, Flux_z DECALIGN;

	printf("dt = %f\n",dt);
	//Get PCM fluxes
	DEBUG_MSG("Calculating fluxes\n");

	Flux_PCM(State,Flux_x,Flux_y,Flux_z,Block,Model);
	//PrintBlockIC(Flux_x,Block);

	DEBUG_MSG("Applying fluxes\n");
	FluxUpdate(State,Flux_x,Flux_y,Flux_z,dt,Block,Model);
	
	DEBUG_MSG("Fluid advance complete\n");
	
}	

void FluxUpdate(BlockCC Prim, BlockIC Fx, BlockIC Fy, BlockIC Fz, Real dt, Block_S Grid, Model_S Model) {
	ISALIGNED(Prim);
	ISALIGNED(Fx);
	ISALIGNED(Fy);
	ISALIGNED(Fz);

	int i,j,k,n;
	Real dFx, dFy, dFz, dtox, dtoy, dtoz;
	Real rho, E, Mx,My,Mz, P;

	const Real Gam = Model.Gam;

	//Use dt from argument instead of Grid for multi-step methods
	dtox = dt/Grid.dx; dtoy = dt/Grid.dy; dtoz = dt/Grid.dz;

 	#pragma omp parallel for collapse(2) \
 		num_threads(TpSB) \
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

				Mz  += dtox*( Fx[MOMZ][k][j][i] - Fx[MOMZ][k][j][i+1] )
					+  dtoy*( Fy[MOMZ][k][j][i] - Fy[MOMZ][k][j+1][i] )
					+  dtoz*( Fz[MOMZ][k][j][i] - Fz[MOMZ][k+1][j][i] );

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



