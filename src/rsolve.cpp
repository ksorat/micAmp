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
void RiemannFluxHLLE(BlockR LeftW,BlockR RightW,BlockR FluxLR,Real Gam) {
	ISALIGNED(LeftW);
	ISALIGNED(RightW);
	ISALIGNED(FluxLR);

	int i;
	
	Real cfl, cfr, ar, al;
	Real Scl[VECBUFF], bp[VECBUFF], bm[VECBUFF] DECALIGN;
	BlockR RoeLR, evals, Fl, Fr DECALIGN;

	//Calculate Roe averages & eigenvalues
	Roes(LeftW,RightW,RoeLR,evals,Gam);
	
	//Calculate wave speeds/scale factor
	#pragma omp simd private(cfl,cfr,al,ar)
	for (i=0;i<VECBUFF;i++) {
		cfl = sqrt( Gam*LeftW [PRESSURE][i]/LeftW [DEN][i] );
		cfr = sqrt( Gam*RightW[PRESSURE][i]/RightW[DEN][i] );

		al = fmin( evals[0]     [i], LeftW [VELX][i] - cfl );
		ar = fmax( evals[NVAR-1][i], RightW[VELX][i] + cfr );

		bp[i] = fmax(ar,0.0);
		bm[i] = fmin(al,0.0);
		Scl[i] = 0.5*(bp[i] + bm[i])/(bp[i] - bm[i]);
	}

	CalcLR_Fluxes(LeftW,RightW,Fl,Fr,bm,bp,Gam);

	//Use scale factor to combine L-R fluxes into interface flux
	//Assuming NVAR ordering: DEN/Vel-xyz/TotalE
	#pragma omp simd
	for (i=0;i<VECBUFF;i++) {
		FluxLR[DEN ][i] = 0.5*( Fl[DEN ][i] + Fr[DEN ][i] ) + Scl[i]*( Fl[DEN ][i] - Fr[DEN ][i] );
		FluxLR[MOMX][i] = 0.5*( Fl[MOMX][i] + Fr[MOMX][i] ) + Scl[i]*( Fl[MOMX][i] - Fr[MOMX][i] );
		FluxLR[MOMY][i] = 0.5*( Fl[MOMY][i] + Fr[MOMY][i] ) + Scl[i]*( Fl[MOMY][i] - Fr[MOMY][i] );
		FluxLR[MOMZ][i] = 0.5*( Fl[MOMZ][i] + Fr[MOMZ][i] ) + Scl[i]*( Fl[MOMZ][i] - Fr[MOMZ][i] );
		FluxLR[TOTE][i] = 0.5*( Fl[TOTE][i] + Fr[TOTE][i] ) + Scl[i]*( Fl[TOTE][i] - Fr[TOTE][i] );
	}

}

//HLLC Riemann solver
void RiemannFluxHLLC(BlockR LeftW,BlockR RightW,BlockR FluxLR,Real Gam) {
	ISALIGNED(LeftW);
	ISALIGNED(RightW);
	ISALIGNED(FluxLR);

	int i;

	Real cfl, cfr, ar, al;
	Real dvl, dvr, tl, tr, dl, dr;
	Real bp[VECBUFF], bm[VECBUFF], am[VECBUFF], cp[VECBUFF] DECALIGN;
	Real SclL[VECBUFF], SclR[VECBUFF], SclM[VECBUFF] DECALIGN;
	BlockR RoeLR, evals, Fl, Fr DECALIGN;

	//Calculate Roe averages & eigenvalues
	Roes(LeftW,RightW,RoeLR,evals,Gam);

	//Calculate wave speeds/scale factor
	#pragma omp simd private(cfl,cfr,al,ar,dvl,dvr,tl,tr,dl,dr)
	for (i=0;i<VECBUFF;i++) {
		cfl = sqrt( Gam*LeftW [PRESSURE][i]/LeftW [DEN][i] );
		cfr = sqrt( Gam*RightW[PRESSURE][i]/RightW[DEN][i] );

		al = fmin( evals[0]     [i], LeftW [VELX][i] - cfl );
		ar = fmax( evals[NVAR-1][i], RightW[VELX][i] + cfr );

		//Force al/ar to be directional
		bp[i] = fmax(ar,0.0);
		bm[i] = fmin(al,0.0);

		//Calculate extra wave speeds for middle region of fan
		//Residual velocities
		dvl = LeftW [VELX][i] - al;
		dvr = RightW[VELX][i] - ar;

		//Pressure in contact wave
		tl = LeftW [PRESSURE][i] + LeftW [DEN][i]* LeftW [VELX][i]*dvl;
		tr = RightW[PRESSURE][i] + RightW[DEN][i]* RightW[VELX][i]*dvr;

		//rho*(v-a)
		dl = LeftW [DEN][i]*dvl;
		dr = -1.0*RightW[DEN][i]*dvr;

		//Contact wave speed and pressure @ moving surface
		am[i] = (tl-tr) / (dl + dr);
		cp[i] = fmax( (dl*tr + dr*tl)/(dl+dr), 0.0 );

		if ( am[i] >= 0) { //Right-moving contact
			SclL[i] = am[i] / ( am[i] - bm[i] );
			SclR[i] = 0.0;
			SclM[i] = -bm[i]/( am[i] - bm[i] );
		} else {
			SclL[i] = 0.0;
			SclR[i] = -am[i]/( bp[i] - am[i] );
			SclM[i] = bp[i]/( bp[i] - am[i] );
		}

	}

	CalcLR_Fluxes(LeftW,RightW,Fl,Fr,bm,bp,Gam);

	//Use scale factor to combine L-R fluxes into interface flux
	#pragma omp simd
	for (i=0;i<VECBUFF;i++) {
		FluxLR[DEN ][i] = SclL[i]*Fl[DEN ][i] + SclR[i]*Fr[DEN ][i];
		FluxLR[MOMX][i] = SclL[i]*Fl[MOMX][i] + SclR[i]*Fr[MOMX][i];
		FluxLR[MOMY][i] = SclL[i]*Fl[MOMY][i] + SclR[i]*Fr[MOMY][i];
		FluxLR[MOMZ][i] = SclL[i]*Fl[MOMZ][i] + SclR[i]*Fr[MOMZ][i];
		FluxLR[TOTE][i] = SclL[i]*Fl[TOTE][i] + SclR[i]*Fr[TOTE][i];

		//Add fixes
		//Correct Mx w/ contact pressure
		FluxLR[MOMX][i] += SclM[i]*cp[i];
		FluxLR[TOTE][i] += SclM[i]*cp[i]*am[i]; 

	}

}
void Roes(BlockR LeftW,BlockR RightW,BlockR RoeLR,BlockR evals,Real Gam) {
	ISALIGNED(LeftW );
	ISALIGNED(RightW);
	ISALIGNED(RoeLR );
	ISALIGNED(evals );

	int i;
	Real invD,hL,hR, vsq, asq, a;

	#pragma omp simd private(invD,hL,hR,vsq,asq,a)
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
void CalcLR_Fluxes(BlockR LeftW, BlockR RightW, BlockR Fl, BlockR Fr, Real bm[VECBUFF], Real bp[VECBUFF], Real Gam ) {
	ISALIGNED(LeftW );
	ISALIGNED(RightW);
	ISALIGNED(Fl    );
	ISALIGNED(Fr    );

	Real vL,vR, El, Er;
	int i;
	
	

	#pragma omp simd private(vL,vR,El,Er)
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
