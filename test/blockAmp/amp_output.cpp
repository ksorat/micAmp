//Various routines to output data
#include <stdio.h>
#include <cstdlib>
#include <string>

#include "amps.h"
#include "prototypes.h"

void toConsole(RealP4 State, Grid_S Grid, Model_S Model) {
	printf("\n");
	printf("Simulation: t=%4.4f -/- Ts = %d -/- dt = %4.4e\n", Grid.t,Grid.Ts,Grid.dt);

}

void toVTK(RealP4 State, Grid_S Grid, Model_S Model) {


	//Output to data file
	char fname[20]="Output.vtk";
	float *data; //Always output in single precision
	float scldata; //Hold scalar data
	FILE *VTKFile;
	int i,j,k, tsdat;
	static int VTKOutNum = 0;

	sprintf(fname, "Output.%04d.vtk", VTKOutNum);
	printf("\tWriting VTK (%d) ...\n", VTKOutNum);

	//Allocate buffer for writing
	if ( (data = (float *)malloc(Grid.Nxp*sizeof(float))) == NULL) {
		printf("Error allocating VTK write-out buffer\n");
		exit(1);
	}

	//Open file for editing
	if ( (VTKFile = fopen(fname,"w")) == NULL ) {
		printf("Error opening VTK output file\n");
		exit(1);
	}

	fprintf(VTKFile,"# vtk DataFile Version 2.0\n");
	fprintf(VTKFile,"VTK file output\n");
	fprintf(VTKFile,"BINARY\n");

	fprintf(VTKFile, "DATASET STRUCTURED_POINTS\n");

	//Add some scalar data
	//Current time
	scldata = (float)Grid.t; bswap(&scldata,sizeof(float),1);
	fprintf(VTKFile, "FIELD FieldData 2\n");
	fprintf(VTKFile, "TIME 1 1 float\n");
	fwrite(&scldata,sizeof(float),1,VTKFile);

	//Current timestep
	tsdat = (int)Grid.Ts; bswap(&tsdat,sizeof(int),1);
	fprintf(VTKFile, "CYCLE 1 1 int\n");
	fwrite(&tsdat,sizeof(int),1,VTKFile);

	//Finish describing grid	
	fprintf(VTKFile, "DIMENSIONS %d %d %d \n", Grid.Nxp+1,Grid.Nyp+1,Grid.Nzp+1);
	fprintf(VTKFile, "ORIGIN %e %e %e \n", Grid.Xbd[0],Grid.Ybd[0],Grid.Zbd[0]);
	fprintf(VTKFile, "SPACING %e %e %e \n", Grid.dx,Grid.dy,Grid.dz);

	fprintf(VTKFile, "CELL_DATA %d \n", Grid.Nxp*Grid.Nyp*Grid.Nzp);


	//Start outputting variables
	WriteVar("Density",DEN, State, Grid, VTKFile, data);
	WriteVar("Pressure",PRESSURE, State, Grid, VTKFile, data);
	WriteVar("Vx",VELX, State, Grid, VTKFile, data);
	WriteVar("Vy",VELY, State, Grid, VTKFile, data);
	WriteVar("Vz",VELZ, State, Grid, VTKFile, data);
	
	fclose(VTKFile);
	free(data);

	VTKOutNum++;
}

void WriteVar(std::string varname, int varnum, Real ****State, Grid_S Grid, FILE *VTKFile, float *data) {

	int i, j, k;

	fprintf(VTKFile, "SCALARS %s float\n", varname.c_str());
	fprintf(VTKFile, "LOOKUP_TABLE default\n");

	for (k=Grid.ks;k<=Grid.ke;k++) {
		for (j=Grid.js;j<=Grid.je;j++) {
			for (i=Grid.is;i<=Grid.ie;i++) {
				data[i-Grid.is] = (float)State[varnum][k][j][i];
			}
			bswap(data,sizeof(float),Grid.Nxp); //Flip for endian-ness
			fwrite(data,sizeof(float),Grid.Nxp,VTKFile);
		}
	}

}
void bswap(void *vdat, int len, int cnt)
{
  char tmp, *dat = (char *) vdat;
  int k;
 
  if (len==1)
    return;
  else if (len==2)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
      dat += 2;
    }
  else if (len==4)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
      tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
      dat += 4;
    }
  else if (len==8)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
      tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
      tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
      tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
      dat += 8;
    }
  else {  /* the general SLOOOOOOOOOW case */
    for(k=0; k<len/2; k++) {
      tmp = dat[k];
      dat[k] = dat[len-1-k];
      dat[len-1-k] = tmp;
    }
  }
}
