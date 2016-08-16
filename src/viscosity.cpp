#include <stdio.h>
#include <math.h>

#include "amps.h"
#include "prototypes.h"

void Flux_Viscous(RealP4 W, RealP4 Fx, RealP4 Fy, RealP4 Fz, Grid_S Grid, Model_S Model) {
	__assume_aligned(W, ALIGN);
	__assume_aligned(Fx, ALIGN);
	__assume_aligned(Fy, ALIGN);
	__assume_aligned(Fz, ALIGN);

	int i,iblk,iG,j,k,d;
	int di,dj,dk, iLim;

	Real ds, dx, dy, dz;
	Real DivV[VECBUFFV], Eta[VECBUFFV], Zeta[VECBUFFV] __attribute__((aligned(ALIGN)));
	Real Vi[NDIM][VECBUFFV], ViscPi[NDIM][VECBUFFV] __attribute__((aligned(ALIGN)));
	Real JacV_L[NDIM][NDIM][VECBUFFV], JacV_R[NDIM][NDIM][VECBUFFV], JacV_I[NDIM][NDIM][VECBUFFV] __attribute__((aligned(ALIGN)));

	RealP4 Flx; //Worth switching to this type?

	Real x,y,z,P; //Spatial position and pressure at interface for viscosity function

	dx = Grid.dx; dy = Grid.dy; dz = Grid.dz;
	//Outer loop over dimensions to set each component of flux vector
	for (d=0;d<NDIM;d++) {
		switch(d) {
			//Setup offsets for directionality
			case DIR_X :
				//X-Dir
				di = 1; dj = 0; dk = 0;
				ds = dx;
				Flx = Fx;
				break;
			case DIR_Y :
				//Y-Dir
				di = 0; dj = 1; dk = 0;
				ds = dy;
				Flx = Fy;
				break;
			case DIR_Z :
				//Z-Dir
				di = 0; dj = 0; dk = 1;
				ds = dz;
				Flx = Fz;
				break;		
		}
		Wipe4Array(Flx,Grid.Nv,Grid.Nz+1,Grid.Ny+1,Grid.Nx+1);
		//Now calculate individual viscous flux comp.
		//Asymmetric bounds b/c of flux centering
		//Fx @ i,j,k is flux between cell i-1,j,k and i,j,k
		//Flx @ i,j,k is flux between cell i-di,j-dj,k-dk and i,j,k

		#pragma omp parallel for default(shared) collapse(2)  \
			private(iblk,i,j,k,iLim,iG) \
			private(x,y,z,P,Eta,Zeta,DivV) \
			private(JacV_L, JacV_R, JacV_I, Vi, ViscPi)
		for (k=Grid.ks;k<=Grid.ke+1;k++) {
			for (j=Grid.js;j<=Grid.je+1;j++) {
				for (iblk=Grid.is;iblk<=Grid.ie+1;iblk+=VECBUFFV) {
					//Inner vector loops
					
					//Limit inner loops in case of bad divisibility
					//Last i-bound is ied
					iLim = IMIN(VECBUFFV,Grid.ie+1-iblk+1);

					iJacV(W, iblk, j, k, iLim, JacV_L, JacV_R, JacV_I, d,dx,dy,dz);
					//Calculate interface velocities, average slope-limited L/R estimates
					#pragma omp simd 
					for (i=0;i<iLim;i++) {
						iG = iblk+i;
						Vi[DIR_X][i] = 0.5*( W[VELX][k-dk][j-dj][iG-di] + 0.5*ds*JacV_L[DIR_X][d][i] 
											+ W[VELX][k][j][iG] - 0.5*ds*JacV_R[DIR_X][d][i] );
	
						Vi[DIR_Y][i] = 0.5*( W[VELY][k-dk][j-dj][iG-di] + 0.5*ds*JacV_L[DIR_Y][d][i] 
											+ W[VELY][k][j][iG] - 0.5*ds*JacV_R[DIR_Y][d][i] );
	
						Vi[DIR_Z][i] = 0.5*( W[VELZ][k-dk][j-dj][iG-di] + 0.5*ds*JacV_L[DIR_Z][d][i] 
										+ W[VELZ][k][j][iG] - 0.5*ds*JacV_R[DIR_Z][d][i] );

					}

					//Calculate interface viscosities
					#pragma omp simd 
					for (i=0;i<iLim;i++) {
						iG = iblk+i;
						//Calculate interface pressure (just average for now, only used in temp-dep. viscosities)
						P = 0.5*( W[PRESSURE][k][j][iG] + W[PRESSURE][k-dk][j-dj][iG-di] );
						x = 0; y = 0; z = 0; //Set these to actual values if visc dep. on space
						Eta[i] = ShearViscosity(x,y,z,P);
						Zeta[i] = BulkViscosity(x,y,z,P);
					}

					//Calculate viscous flux tensor (only need 1 row, \Pi_{d,:})
					#pragma omp simd 
					for (i=0;i<iLim;i++) {
						iG = iblk+i;

						DivV[i] = JacV_I[DIR_X][DIR_X][i] + JacV_I[DIR_Y][DIR_Y][i] + JacV_I[DIR_Z][DIR_Z][i];

						//Calculate non-Identity part of Pi, Eta*[ \grad V \grad V^{T} ]
						ViscPi[DIR_X][i] = Eta[i]*( JacV_I[d][DIR_X][i] + JacV_I[DIR_X][d][i] );
						ViscPi[DIR_Y][i] = Eta[i]*( JacV_I[d][DIR_Y][i] + JacV_I[DIR_Y][d][i] );
						ViscPi[DIR_Z][i] = Eta[i]*( JacV_I[d][DIR_Z][i] + JacV_I[DIR_Z][d][i] );
						//Now add in identity part to \Pi_{d,d}
						ViscPi[d][i] += Zeta[i]*DivV[i] - (2.0/3.0)*Eta[i]*DivV[i];

					}
					//Calculate flux in direction d
					//F_{d} = [0,\Pi_{d,x},\Pi_{d,y},\Pi_{d,z},V_{i}\Pi_{d,i}]
					//Add negative sign b/c math
					#pragma omp simd 
					for (i=0;i<iLim;i++) {
						iG = iblk+i;
						Flx[DEN][k][j][iG] = 0.0;
						Flx[VELX][k][j][iG] = -1*ViscPi[DIR_X][i];
						Flx[VELY][k][j][iG] = -1*ViscPi[DIR_Y][i];
						Flx[VELZ][k][j][iG] = -1*ViscPi[DIR_Z][i];
						Flx[PRESSURE][k][j][iG] = -1*( Vi[DIR_X][i]*ViscPi[DIR_X][i] + Vi[DIR_Y][i]*ViscPi[DIR_Y][i] + Vi[DIR_Z][i]*ViscPi[DIR_Z][i] );

					}


				} //iblock loop
			}
		}
	} //Direction loop

}

//Calculate slope-limited Jacobian of velocity @ interface
//Using convention that J[i,j] = \partial_{x_j} V_{x_i}
//Vectorized version, start at ib0/j0/k0 and do until ib0+iLim
//Find interface Jacobians in direction intDir

void iJacV(RealP4 W, int ib0, int j0, int k0, int iLim, Real JacV_L[NDIM][NDIM][VECBUFFV], Real JacV_R[NDIM][NDIM][VECBUFFV], Real JacV_I[NDIM][NDIM][VECBUFFV], int intDir, Real dx, Real dy, Real dz) {
	__assume_aligned(JacV_L, ALIGN);
	__assume_aligned(JacV_R, ALIGN);
	__assume_aligned(JacV_I, ALIGN);

	const int vOff = VELX; // Offset for velocity, ie d+vOff = VEL-D
	int vDIR, iL, iR, jL, jR, kL, kR, diL;
	int d2,d1,i, di, dj, dk;
	Real ds, jvL,jvR,jvI;

	jR = j0; kR = k0;
	jL = jR - (int)(intDir == DIR_Y); //Cast bool to get 0/1
	kL = kR - (int)(intDir == DIR_Z);
	diL = (int)(intDir == DIR_X);

	for (d2=0;d2<NDIM;d2++) { //Velocity component
		vDIR = d2 + vOff;

		for (d1=0;d1<NDIM;d1++) { //Derivative wrt
			dk = (int) (d1 == DIR_Z);
			dj = (int) (d1 == DIR_Y);
			di = (int) (d1 == DIR_X);
			ds = dx*di + dy*dj + dz*dk; //Pick out relevant grid spacing
			#pragma omp simd 
			for (i=0;i<iLim;i++) {
				iR = ib0+i;
				iL = iR - diL;
				//Calculate JacV[d2][d1][i]
				jvL = Sig( W[vDIR][kL-dk][jL-dj][iL-di], W[vDIR][kL][jL][iL], W[vDIR][kL+dk][jL+dj][iL+di], ds );
				jvR = Sig( W[vDIR][kR-dk][jR-dj][iR-di], W[vDIR][kR][jR][iR], W[vDIR][kR+dk][jR+dj][iR+di], ds );
				//Construct interface Jacobian componenet
				//For now simply averaging, could use upwinding or weighted-correction
				//See Loh/Puigt
				jvI = 0.5*( jvL + jvR );
				JacV_L[d2][d1][i] = jvL;
				JacV_R[d2][d1][i] = jvR;
				JacV_I[d2][d1][i] = jvI;
			}

		}
	}
}


//Calculate limited slope given oriented L/C/R values and grid spacing
#pragma omp declare simd uniform(dx)
Real Sig(Real Wl, Real Wc, Real Wr, Real dx) {
	Real dL, dC, dR, sgn, mSig;

	mSig = 0.0;
	dL = (Wc-Wl);
	dR = (Wr-Wc);
	dC = 0.5*(Wr-Wl);

	//Use MC limiter but could be anything (and indep. of slope lim in recon)
	if (dL*dR <= 0) {
		mSig = 0.0;
	} else {
		sgn = (dL > 0) - (dL < 0); //Cheap sgn function w/o branch
		mSig = sgn*fmin( fabs(dC), 2.0*fmin(fabs(dL),fabs(dR)) );
	}
	mSig = mSig/dx; //For correct units
	return mSig;

}

#pragma omp declare simd
Real ShearViscosity(Real x, Real y, Real z, Real P) {
	Real Eta = SHEARVISC;

	return Eta;
}

#pragma omp declare simd
inline Real BulkViscosity(Real x, Real y, Real z, Real P) {
	Real Zeta = 0.0;

	return Zeta;
}

