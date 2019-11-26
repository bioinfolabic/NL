#include "defines.h"


/***   Timer   ***/
void initTimer();
void finTimer();

/***   Files   ***/
int getParameter(char *field_name, char *token, char *format, void *variable, FILE *fi);
int loadFile(char **argv);
void putParameters();

/***   Utilities   ***/
double randdouble(double max);
int isUnique(Particle *p, int last);
void verifyBoundary1(VectorR *v);
void verifyBoundary2(VectorR *v);
void setVelMag();
void buildConstMatrix();
void calcStatus();
void calcRG();
void evaluate();
void calcCenterMass();
void printSummary(char **argv);
void report(char **argv);
void savePathways(char **argv);
void saveBestCoordinates(char **argv);


/***   Initialize   ***/
void allocArrays();
void initCoords();
void unitaryVectors();
void initVels();
void initAccs();
void get_structure();
void initMD();

/***   MD   ***/
void updatePos();
void computeForces();
void updateVelocities();
void step();

/***   Neighbourhood List   ***/
void *updateNeighbourList();

/***   SHAKE   ***/
void shakeRelaxation();


/***   Thermostat   ***/
void berendsenThermostat();

/***   End   ***/
void freeArrays();
void finishSim(char **argv);
void freeDevice();
