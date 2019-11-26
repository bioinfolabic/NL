/*****************************************************************
defines.h: defines, includes and globals
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#ifndef TRUE
#define TRUE    1
#define FALSE   0
#endif

#define FAIL    0
#define NDIM 3     										//2 --> two-dimensional space definition
														//3 --> three-dimensional space definition


#define MASK 2147483647									//para o gerador de numeros aleatorios--> trocar pelo Mersenne Twister!
#define SCALE 0.4656612873e-9
#define IADD 453806245
#define IMUL 314159269


//select a numerical method to integrate Newton's equations of motion: 		 
#define VELOCITYVERLET									
//#define LEAPFROG
//#define VERLET

#define RESCALEVELOCITIES 								//rescale velocities enable

#define SITESMOL	4									//number of interaction sites in each molecule


#define Max(x1, x2) (((x1) > (x2)) ? (x1) : (x2))
#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * (x) * (x))


//Defines Linear Equations solver
#define SLE_NMAX 200  				//max number of variables

//Protein max length
#define SEQ_MAX_LENGTH 200

#define TAMPATHWAY 1001


typedef unsigned int BOOL;
typedef double real; //synonymous for "float" or "double"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											/** Para salvar os apenas dados no final **/
typedef struct 										//struct of one conformation
{
	double x[SEQ_MAX_LENGTH];
	double y[SEQ_MAX_LENGTH];
	double z[SEQ_MAX_LENGTH];
	double potEnergy;
	int stepCount;
	double rgAll;
	double rgH;
	double rgP;
}conformation;

conformation conformations[TAMPATHWAY];


/***********************************************************************************************************************************/
typedef struct 										//struct of best conformation
{
	real potEnergy;
	int step;
}Bestconformation; 

Bestconformation best;


//vector definition:
#if NDIM == 2
	typedef struct{ real x,y;}VectorR; 					//struct definition to represent a two-dimensional vector
#else 
	typedef struct{ real x,y,z;}VectorR;				//struct definition to represet a three-dimensional vector
#endif

//vector definition for lattice coordinates:
#if NDIM == 2
	typedef struct{int x,y;}VectorI; 					//struct definition to represent a two-dimensional vector for square lattice
#else 
	typedef struct{int x,y,z;}VectorI;					//struct definition to represet a three-dimensional vector for cubic lattice
#endif


//particle's variables definition (molecule) 
typedef struct {
	VectorR r, rv, ra;
}Particle; 			//r, v, a correspond to the coordinate, velocity and acceleration vectors of the particle


//Definitions for rigid molecules:

typedef struct {real u1, u2, u3, u4;}Quat;	//Hamilton quaternions components				


//molecule's variables definition (Mol) 
typedef struct {
	VectorR r, rv, ra, ra1, ra2, ro, rvo;
	Quat q, qv, qa, qa1, qa2, qo, qvo;
	VectorR torq;
}Mol; 			
//r, v, a correspond to the coordinate, velocity and acceleration vectors of the particle

typedef struct {VectorR f, r;} IntSite; 	//interaction sites structure
typedef struct {real u[9];} RotMat; 		//rotation matrix structure 3x3
typedef struct {VectorR r; int typeF;} MSite;//interaction site coordinates. 
											//typeF is used to distinguish individual sites for the force evaluation.
											//for example, typeF=0 -->Oxygen, typeF=1 --> Hydrogen 


typedef struct {VectorR v; real BondLen2; real dist2; int ik; int jk;} constr;


//typedef struct {real val, sum, sum2;} Prop; 			//some measurement properties


//Globals:
IntSite *site;
MSite *mSite;
VectorR mInert;
real bCon; 
constr *constraint;
int *mMat;
real *lMat;
real *consV;

real consPrec; //tolerance used in establishing convergence used in the relaxation method (shake)
int maxCycle;  // number of iterations used in the relaxation method (shake)
int nCycleR, nCycleV;


int moreCycles, ID, nMol, stepAvg, stepCount, displayInterval, stepLimit, step2rescVels, step2shake, step2report;
char savepathways;
int pathwaysstep;
int printsummaryInterval;
char printsummary2file;
real LV; //linear size of cubical ou square volume
//int LV2; // 

time_t ti, tp;							//initial time and processing time

Mol *mol;
Mol *beststructure; 					//best structure coordinates, velocities, among others.

VectorR region, vSum;
real dt, rCut, temperature, temperature_steps, velMag, current_temperature, mass;
int BondLen, ProtLen, nC;				//bond length, protein size, number of constraints
real uSum;								//total interaction energy
real uLJ;								//total Lennard-Jones potential
real uTorsion; 							//total Torsion potential
real uChainAngles;						//total Chain angles potential
char report_file, printsummary;
VectorR center_mass;

//properties
real kineticEnergy;
real potEnergy;
real pressure;
real virial;
real density;
real TotalEnergy;							//internal energy = Ekinetic + Epotencial

real rgH, rgP, rgAll, bondavg;

 
char sequence[SEQ_MAX_LENGTH]; //protein sequence (max length = SEQ_MAX_LENGTH)

//variables used to draw the 3D picture of the system:
#define xWindowSize  640				  // window size in screen pixels
#define yWindowSize  640 

#define PI 4 * atan(1.0)
#define RADIUS 0.3	                      // radius of molecule
#define ANGLE  5                          // rotation angle in degrees

int main_w;						  		  //main window e sub-windows (w1, w2, ...)


double minExtent[3], maxExtent[3];        // extent of system volume



int zoomfactor;						  	  //zoomfactor
int phi, theta;                           // to rotate system using arrow keys


//thermostat
real cT; 								//Berendsen thermostat parameter
