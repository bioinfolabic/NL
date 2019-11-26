#include "defines_CPU.h"


void Initialize_MD();
void AllocArrays();
void freeArrays();
void SetParams();
void Coords_Init();
void Vel_Init();
void Accel_Init();
void Initialize_MD();
void EvalProps();
void PrintSummary (FILE *fp);
void AccumProps (int icode);
int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer);
void GetParameters(char **argv);
void PrintParameters();
void show_coordinates(int part_id); 
void create_PDB(int iter);

//Force and Accelarations, quaternions
void velocityVerlet();
void computeAccelerations();	
void computeTorsionForces();  	//compute chain torsion forces 
void computeAngleForces(); //compute bond angle forces
void computeLennard_Jones();	//compute Lennard Jones potential
void verify_boundary(); 
void rescaleVelocities(); 	//velocity rescaling thermostat
void berendsen_thermostat(); //berendsen thermostat
void calc_currentTemp();
void calc_AngVelocities(int n, VectorR *w);
void calc_AccelQuat();
void calc_Torqs();
void BuildRotMatrix(RotMat *rMat, Quat *q, int transpose);
void GenSiteCoords();
void PredictorStep();
void PredictorStepQ();
void CorrectorStep();
void CorrectorStepQ();
void thermostat();
void InitAngCoords();
void EulerToQuat(Quat *qe, real *eAng);
void InitAngVels();
void InitAngAccels();
void VRand(VectorR *p);
void Calc_SiteForces();
void DefineMol();
void AdjustQuat();
void ApplyBoundaryCond();
real RandR(); 				//trocar por Mersenne Twister
void calc_centermass();
void makebackbone();		//draw the backbone
void shake_relaxation();
real RgH();
real RgP();
real RgAll();
double randdouble(double max); 


void evaluate();
void calc_kinetic_energy();
void calc_currentTemp();
void calc_pressure();
void calc_density();

//functions for geometrically constrained molecules:
void buildConstMatrix();
void calc_constraints();
void solve_linearEq(real *a, real *x, int n);
real BondAverage();


void print_summary();
real calc_Total_energy();
void report();
void time_ex();
void initial_message();
void saveBestcoordinates();
void savePathways(); 			//this function saves protein structure coordinates in several steps. 
								//these files can be used to show the folding pathways.
void savePathwaysVetor(int pos);
void savePathwaysVetorArquivo();



int Step(); //main function



