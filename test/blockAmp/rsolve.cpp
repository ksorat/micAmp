#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

//Standard flux functions
//Input: LeftW[NVAR][VECBUFF], RightW[][] the Left/Right *INTERFACE* states
//Output: Flux @ the interface
//General HLL-style flux calculation
//Calculate Roe LR states/eigenvalues
//Calculate wave-speeds, use these to calculate weights/scales
//Use outer-most speeds of fan & LR states to get L-R fluxes (two sets of fluxes based on lW/rW)
//Use scaling numbers to combine Fl & Fr => FluxLR

//Takes LR states (VECBUFF of them) and returns Riemann fluxes
void RiemannFluxHLLE(BlockR LeftW,BlockR RightW,BlockR FluxLR) {
	ISALIGNED(LeftW);
	ISALIGNED(RightW);
	ISALIGNED(FluxLR);

	int i;
	Real Gam;
	Real cfl, cfr, ar, al;
	Real Scl[VECBUFF], bp[VECBUFF], bm[VECBUFF] DECALIGN;
	BlockR RoeLR, evals, Fl, Fr DECALIGN;

	Gam = Model.Gam;

	//Calculate Roe averages & eigenvalues
	Roes_Vec(LeftW,RightW,RoeLR,evals);
	
	
	//Calculate wave speeds/scale factor
	//#pragma omp simd private(cfl,cfr,al,ar)
	//#pragma omp simd
	for (i=0;i<VECBUFF;i++) {
		cfl = sqrt( Gam*LeftW [PRESSURE][i]/LeftW [DEN][i] );
		cfr = sqrt( Gam*RightW[PRESSURE][i]/RightW[DEN][i] );

		al = fmin( evals[0]     [i], LeftW [VELX][i] - cfl );
		ar = fmax( evals[NVAR-1][i], RightW[VELX][i] + cfr );

		bp[i] = fmax(ar,0.0);
		bm[i] = fmin(al,0.0);
		Scl[i] = 0.5*(bp[i] + bm[i])/(bp[i] - bm[i]);
		printf("bp/bm/Scl = %f %f %f", bp[i],bm[i],Scl[i]);
	}

	CalcLR_Fluxes(LeftW,RightW,Fl,Fr,bm,bp);

	//Use scale factor to combine L-R fluxes into interface flux
	//Assuming NVAR ordering: DEN/Vel-xyz/TotalE
	//#pragma omp simd
	for (i=0;i<VECBUFF;i++) {
		FluxLR[0][i] = 0.5*( Fl[0][i] + Fr[0][i] ) + Scl[i]*( Fl[0][i] - Fr[0][i] );
		FluxLR[1][i] = 0.5*( Fl[1][i] + Fr[1][i] ) + Scl[i]*( Fl[1][i] - Fr[1][i] );
		FluxLR[2][i] = 0.5*( Fl[2][i] + Fr[2][i] ) + Scl[i]*( Fl[2][i] - Fr[2][i] );
		FluxLR[3][i] = 0.5*( Fl[3][i] + Fr[3][i] ) + Scl[i]*( Fl[3][i] - Fr[3][i] );
		FluxLR[4][i] = 0.5*( Fl[4][i] + Fr[4][i] ) + Scl[i]*( Fl[4][i] - Fr[4][i] );
	}

}

void Roes_Vec(BlockR LeftW,BlockR RightW,BlockR RoeLR,BlockR evals) {
	ISALIGNED(LeftW );
	ISALIGNED(RightW);
	ISALIGNED(RoeLR );
	ISALIGNED(evals );

	int i;
	Real invD,hL,hR, vsq, asq, a;
	const Real Gam = Model.Gam;

	//#pragma omp simd private(invD,hL,hR,vsq,asq,a)
	//#pragma omp simd
	for (i=0;i<VECBUFF;i++) {
		RoeLR[DEN][i] = sqrt( LeftW[DEN][i]*RightW[DEN][i] );
		invD = 1.0/( sqrt(LeftW[DEN][i]) + sqrt(RightW[DEN][i]) );

		RoeLR[VELX][i] = invD*( sqrt(LeftW[DEN][i])*LeftW[VELX][i] + sqrt(RightW[DEN][i])*RightW[VELX][i] );
		RoeLR[VELY][i] = invD*( sqrt(LeftW[DEN][i])*LeftW[VELY][i] + sqrt(RightW[DEN][i])*RightW[VELY][i] );
		RoeLR[VELZ][i] = invD*( sqrt(LeftW[DEN][i])*LeftW[VELZ][i] + sqrt(RightW[DEN][i])*RightW[VELZ][i] );

		//Calculate L/R enthalpy
		hL = inlineEnthalpy(LeftW [DEN][i], LeftW [VELX][i], LeftW [VELY][i], LeftW [VELZ][i], LeftW [PRESSURE][i], Gam);
		hR = inlineEnthalpy(RightW[DEN][i], RightW[VELX][i], RightW[VELY][i], RightW[VELZ][i], RightW[PRESSURE][i], Gam);

		//Enthalpy really
		RoeLR[PRESSURE][i] = invD*( sqrt(LeftW[DEN][i])*hL + sqrt(RightW[DEN][i])*hR );

		//Calculate eigenvalues
		vsq = SQR(RoeLR[VELX][i]) + SQR(RoeLR[VELY][i]) + SQR(RoeLR[VELZ][i]);
		asq = (Gam-1)*fmax(RoeLR[PRESSURE][i] - 0.5*vsq, TINY);
		a = sqrt(asq);
		evals[0][i] = RoeLR[VELX][i] - a;
		evals[1][i] = RoeLR[VELX][i];
		evals[2][i] = RoeLR[VELX][i];
		evals[3][i] = RoeLR[VELX][i];
		evals[4][i] = RoeLR[VELX][i] + a;
	}
}

//Use outer-most wave speeds of fan to calculate fluxes of lW/rW (outer states)
//Will later use weights to turn these fluxes into the fluxes in the intermediate regions of the fan
void CalcLR_Fluxes(BlockR LeftW, BlockR RightW, BlockR Fl, BlockR Fr, Real bm[VECBUFF], Real bp[VECBUFF] ) {
	ISALIGNED(LeftW );
	ISALIGNED(RightW);
	ISALIGNED(Fl    );
	ISALIGNED(Fr    );

	Real vL,vR, El, Er;
	int i;
	const Real Gam = Model.Gam;
	

	//#pragma omp simd private(vL,vR,El,Er)
	//#pragma omp simd
	for (i=0;i<VECBUFF;i++) {
		vL = LeftW [VELX][i] - bm[i];
		vR = RightW[VELX][i] - bp[i];

		//Mass flux (Mx-bm*d)
		Fl[DEN][i] = LeftW [DEN][i]*vL;
		Fr[DEN][i] = RightW[DEN][i]*vR;

		//Mx flux ( Mx*(Vx-bm) )
		Fl[MOMX][i] = LeftW [DEN][i]*LeftW [VELX][i]*vL;
		Fr[MOMX][i] = RightW[DEN][i]*RightW[VELX][i]*vR;

		//My flux
		Fl[MOMY][i] = LeftW [DEN][i]*LeftW [VELY][i]*vL;
		Fr[MOMY][i] = RightW[DEN][i]*RightW[VELY][i]*vR;

		//Mz flux
		Fl[MOMZ][i] = LeftW [DEN][i]*LeftW [VELZ][i]*vL;
		Fr[MOMZ][i] = RightW[DEN][i]*RightW[VELZ][i]*vR;

		//Add pressure correction in x (normal to interface) direction
		Fl[MOMX][i] = Fl[MOMX][i] + LeftW[PRESSURE][i];
		Fr[MOMX][i] = Fr[MOMX][i] + RightW[PRESSURE][i];

		//Energy flux [ E*(Vx-bm) + P*Vx ]
		El = ( LeftW [PRESSURE][i]/(Gam-1) )+ 0.5*LeftW [DEN][i]*( SQR(LeftW [VELX][i]) + SQR(LeftW [VELY][i]) + SQR(LeftW [VELZ][i]) );
		Er = ( RightW[PRESSURE][i]/(Gam-1) )+ 0.5*RightW[DEN][i]*( SQR(RightW[VELX][i]) + SQR(RightW[VELY][i]) + SQR(RightW[VELZ][i]) );

		Fl[TOTE][i] = El*vL + LeftW [PRESSURE][i]*LeftW [VELX][i];
		Fr[TOTE][i] = Er*vR + RightW[PRESSURE][i]*RightW[VELX][i];

	}	
}

inline Real inlineEnthalpy(Real d, Real Vx, Real Vy, Real Vz, Real P,Real Gam) {
	Real K, e, H;

	e = P/(Gam-1);
	K = 0.5*d*( SQR(Vx) + SQR(Vy) + SQR(Vz) );
	H = (P + e + K)/d;
	return H;

}
