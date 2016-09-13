//Main driver file for MIC-amp
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amps.h"
#include "prototypes.h"

MICTYPE Model_S Model;

int main(int argc, char *argv[]) {

	
	Grid_S Grid;
	RealP4 State;
	RealP4 StatePhi;
	Real *StatePhi0, *State0; 

	int m=0; //Index of mic card
	int Ntot;

	//Setup
	// - Fill in Model data structure
	// - Fill in Grid data structure
	printf("Begin host initialization\n");
	conModel(&Model);
	conGrid(&Grid);

	Ntot = Grid.Nv*Grid.Nz*Grid.Ny*Grid.Nx;
	State = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);

	//Initialize system on host
	initialConds(State,Grid,Model);
	State0 = &(State[0][0][0][0]);

	printf("Begin device initialization\n");
	//Create data container on Phi
	#pragma offload_transfer target(mic:m) nocopy(StatePhi0 : length(Ntot) ALLOC)

	/*#pragma offload target(mic:m) in(Ntot,Grid) nocopy(StatePhi:REUSE) nocopy(StatePhi0:length(Ntot) ALLOC)
	{
		StatePhi = Map4Array(StatePhi0,Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);
	}*/
	printf("Initializing integrator\n");
	//InitializeIntegrator(Grid,Model);

	//Enforce BC's
	EnforceBCs(State, Grid, Model);

	//Calculate initial timestep
	Grid.dt = CalcDT(State,Grid,Model);

	//Do initial outputs
	if ( (Grid.Ts) % Model.TsOut == 0 ) {
		toConsole(State,Grid,Model);
		toVTK(State,Grid,Model); 
	}	

	while (Grid.t < Grid.Tfin) {
		//Evolve system
		#pragma offload target(mic:m) \
			in (State0: into(StatePhi0) length(Ntot) REUSE) \
			out(StatePhi0: into(State0) length(Ntot) REUSE) \
			nocopy(StatePhi)
		{
			AdvanceFluid(StatePhi, Grid, Model, Grid.dt);
		}

		//Enforce BCs
		EnforceBCs(State, Grid, Model);

		//Update time
		Grid.Ts++;
		Grid.t = Grid.t + Grid.dt;

		//Calculate new timestep
		Grid.dt = CalcDT(State,Grid,Model);

		//Output if necessary
		if ( (Grid.Ts) % Model.TsOut == 0 ) {
			toConsole(State,Grid,Model);
			toVTK(State,Grid,Model); 
		}

	}


	Kill4Array(State);
	#pragma offload target(mic:m) \
		nocopy(StatePhi) nocopy(StatePhi0: length(Ntot) FREE)
	{
		Kill4Array(StatePhi);
	}

	DestroyIntegrator(Grid,Model);
}