//Main driver file for MIC-amp
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amps.h"
#include "prototypes.h"

Model_S Model;

int main(int argc, char *argv[]) {

	
	Grid_S Grid;
	RealP4 State;


	//Setup
	// - Fill in Model data structure
	// - Fill in Grid data structure
	conModel(&Model);
	conGrid(&Grid);

	State = Create4Array(Grid.Nv,Grid.Nz,Grid.Ny,Grid.Nx);

	//Initialize system
	initialConds(State,Grid,Model);
	InitializeIntegrator(Grid,Model);

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
		tic();
		BlockAdvance(State,Grid,Model,Grid.dt);
		toc();

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
	DestroyIntegrator(Grid,Model);
}