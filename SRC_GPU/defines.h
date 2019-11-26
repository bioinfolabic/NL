#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define FAIL 0
#define SUCCESS 1
#define N_DIM 3
#define SEQ_MAX_LENGTH 80000	// Max Sequence
#define RAND_MAX_MT 4294967295  // Mersenne Twister maximum 2^32-1
#define BLOCKSIZE 256 			// Block size of the functions in gpu
#define BLOCKSIZESUM 1024		// Block size specifically to the function of vector sum
#define sqr(x) ((x) * (x))
#define min(x, y) (((x) < (y))? (x) : (y))

//typedef double double; // Problem with atomic operations

/***   Structures   ***/
typedef struct {
	double x;
	double y;
	double z;
} VectorR; // 3D vectors


typedef struct {
	VectorR v_r;
	VectorR v_v;
	VectorR v_a;
} Particle;	// Particle stores the position, velocity and acceleration of a amino-acid



typedef struct {
	VectorR v;
	int ik;
	int jk;
} Constraint; // Stores a position and the position of the particle in the vector

Particle *particles;		// Amino-acids vector
Particle *best_structure;	// Best Structure
Constraint * constraint;	// Constraint vector (For computation of geometric constraints)

/***   From file   ***/
char sequence[SEQ_MAX_LENGTH];	// AB Sequence
double mass;				// mass of each particle
int n_mol;					// Number of Molecules
int bond_len;				// Bond Length
int prot_len;				// Protein Length
int n_c;					// Number of Constraints
int n_blocks;				// Number of blocks that will be used in the GPUs computations
int n_blocks_sum;			// Number of blocks that will be used in the sum of vector in GPU
double LV;					// Dimension of the Box
double dt;					// Time Step
double c_T;					
double temperature;
double temperature_steps;

int step2neighbour;			// Steps between each update of neighbour list
int display_interval;
int step_limit;
int step2resc_vels;
char report_file;
int step2report;
char print_summary;
int print_summary_interval;
char print_summary2file;
char save_pathways;
int pathways_step;

double shake_cons_prec;		// Shake constratint precision
int shake_max_cycle;	
int shake_step2shake;		// Stepes between each update of the geometrics constraints

/***   Simulation control   ***/
time_t t_i, delta_t;
struct timeval t_now, last_t, result_t;
int i_step;
double vel_mag;				// Velocity magnetude(Temperature dependent)
int en_update;				// Enable update


/***   Energies   ***/
double *auxSumLJ;			// Lennard-Jones Energy Vector
double *auxSumB;			// Bond Angle Energy Vector
double *auxSumT;			// Torsion Force Energy Vector
double uBond;				// Bond Angle energy of the structure
double uTorsion;			// Torsion force energy of the structure
double uLJ;					// Lennard-Jones energy of the structure
double uSum;				// Sum of the three partial energies
double kinetic_energy;
double total_energy;

/***   SHAKE   ***/
int nCycleR;
int nCycleV;

/***   Other   ***/
double current_temperature;		// Actual temperature of the structure
double density;					// Density of the structure
double bond_avg;				// Bond distance average
VectorR center_mass;			// Center of mass
double rGH, rGP, rG;			// Gyration radius

/***   Best   ***/
double best_potencial_energy;	// Best potential energy obtained
int best_step;					// Best step

/***   Device variables   ***/
unsigned int blockSize;
unsigned int blockSizeSum;
Particle *d_particles;
double *d_LV;
double *d_uB;
double *d_uT;
double *d_uLJ;
double *d_uLJVector;
double *d_uBVector;
double *d_uTVector;
double *d_auxSumLJ;
double *d_auxSumB;
double *d_auxSumT;
int device;
char *d_sequence;			// AB sequence

Constraint * d_constraint;


