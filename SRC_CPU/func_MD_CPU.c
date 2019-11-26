#include "func_MD_CPU.h"


/* Dynamic array allocation */
/* mol and sites */
void AllocArrays()	
{
	mol = (Mol *) malloc ((nMol) * sizeof (Mol));
	beststructure = (Mol *) malloc ((nMol) * sizeof (Mol));				//best structure features

	site = (IntSite *)	 malloc (nMol * SITESMOL * sizeof(IntSite));  //water molecules, TIP4P model
	mSite = (MSite *)	 malloc (SITESMOL * sizeof(MSite));

	constraint = (constr *) malloc ((nC) * sizeof (constr));		//constraints (for geometrically constrained molecules)

	mMat = (int *) malloc (ProtLen * nC * sizeof (int));

	lMat = (real *) malloc (nC * nC * sizeof (real));

	consV = (real *) malloc (nC * sizeof (real));

}


/* Free arrays */
/* free mol and sites */
void freeArrays()
{
	free(mol);
	free(beststructure);
	free(site);
	free(mSite);
	free(constraint);
	free(mMat);
	free(lMat);
	free(consV);
}

/* Parameters setting */
/* velMag =magnitude of initial velocities */
/* velMag needs to be updated  when the temperature changes 
(for constant-temperature simulations). See Functions  mainw_menu() and rescaleVelocities()*/
void SetParams()
{
	//rCut = pow (2., 1./6.);

	velMag = sqrt(NDIM * (1. - 1. / ProtLen) * temperature);  //velMag = sqrt ((NDIM * (1. - 1. / ProtLen) - (real) nC / ProtLen) * temperature);
}


/* Coordinates initialization */
/* Molecular Coordinates are initialized in a cubic (or square) lattice*/
/* The Face-centered cubic (FCC) lattice will be implemented in another version*/
void Coords_Init()
{
	int i, n;
	double alpha, theta;
	int conf_OK; //flag
	real dist;
	//random protein conformation:

	mol[0].r.x = LV/2;	//first amino acid
	mol[0].r.y = LV/2;	
	mol[0].r.z = LV/2;

	for (i = 1; i < ProtLen; i++)
	{
		
		do
		{
			alpha = randdouble(360)*PI/180;
			theta = randdouble(180)*PI/180;		
					
			mol[i].r.x = mol[i - 1].r.x + cos (theta) * sin (alpha);
			mol[i].r.y = mol[i - 1].r.y + sin (theta) * sin (alpha);
			mol[i].r.z = mol[i - 1].r.z + cos (alpha); //mol[i - 1].r.z + cos (theta);

	
			//verification
			conf_OK = 1; //OK
			for (n=0; n < i; n++)
			{
				dist = sqrt((mol[i].r.x - mol[n].r.x) * (mol[i].r.x - mol[n].r.x) + (mol[i].r.y - mol[n].r.y) * (mol[i].r.y - mol[n].r.y) + (mol[i].r.z - mol[n].r.z) * (mol[i].r.z - mol[n].r.z));
				if (dist < 1)
					conf_OK = 0; //NOK
			}
		}while(conf_OK != 1); //NOK?
	}

}

/* Initial velocities */
/* The initial velocities are set to velocity magnitude that depends on the
temperature (see function SetParams() ). Velocity directions are also randomly assigned.
rand() can be replaced by the MErsenne-Twister number generator */
void Vel_Init()
{
	int n;

	VectorR sumrv;
	
	sumrv.x = sumrv.y = sumrv.z = 0;

	for (n = 0; n < ProtLen; n ++)
	{
		VRand (&mol[n].rv);			//it generates a randomly oriented vector of unit length
		mol[n].rv.x *=  velMag;
		mol[n].rv.y *=  velMag;
		mol[n].rv.z *=  velMag; 
		sumrv.x = sumrv.x + mol[n].rv.x;
		sumrv.y = sumrv.y + mol[n].rv.y;
		sumrv.z = sumrv.z + mol[n].rv.z;
	
	}

	for (n = 0; n < ProtLen; n ++)							//uniform distribution --> Scale to zero total linear momentum 
	{	
		mol[n].rv.x = mol[n].rv.x - sumrv.x / ProtLen;
		mol[n].rv.y = mol[n].rv.y - sumrv.y / ProtLen;
		mol[n].rv.z = mol[n].rv.z - sumrv.z / ProtLen;  

	}
}

/* Initial Accelerations*/
/* The initial accelerations are initialized to zero */
void Accel_Init()
{
	int n;

	for (n = 0; n < nMol; n ++)
	{
		mol[n].ra.x  =  0;
		mol[n].ra1.x =  0;
		mol[n].ra2.x =  0;

		mol[n].ra.y  =  0;
		mol[n].ra1.y =  0;
		mol[n].ra2.y =  0;

		mol[n].ra.z  =  0;
		mol[n].ra1.z =  0;
		mol[n].ra2.z =  0;
	}

}


/* Initialize the MD simulation */
void Initialize_MD()
{
	stepCount = 0; 				//step counter initialization
	printf("Alocando arrays\n");
	AllocArrays(); 				//arrays allocation
	printf("Coordes init\n");
	Coords_Init(); 		
	printf("Vel init\n");		
	Vel_Init();				//velocities initialization
	printf("aceel init");
	Accel_Init();
	printf("matriz");
	buildConstMatrix();			//Constraint Matrix




}


/*Input file reading*/
void GetParameters(char **argv)
{
	FILE *file = fopen( argv[1], "r" );   // precisamos indicar a pasta na hora de compilar
    

	if (file == 0)
	{
    	printf( "Could not open file!\n" );
	}
    else 
    {	
		ffscanf("sequence", file, "%s", &sequence);
		ffscanf("ProtLen", file, "%d", &ProtLen);
		ffscanf("LV", file, "%lf", &LV);
		ffscanf("stepLimit", file, "%d", &stepLimit);
		ffscanf("temperature", file, "%lf", &temperature);
		ffscanf("savepathways", file, "%c", &savepathways);
		ffscanf("pathwaysstep", file, "%d", &pathwaysstep);	
		
		mass = 1.0; 
		nMol = ProtLen;
		BondLen = 1.0;
		nC = (ProtLen-1);
		dt = 0.0001;
		rCut = 2.38;
		cT = 0.01;
		displayInterval = 100;
		step2rescVels = 1;
		temperature_steps = 0.10;
		report_file = 'y';
		step2report = 'n';
		printsummary = 160;
		printsummaryInterval = 'y';
		printsummary2file = 'y';
		consPrec = 1.0e-06;
		maxCycle = 1; 
		step2shake = 10;

    }
    fclose(file);

}


int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer)
{
	char buffer[600];
	int len;
    int commentflag = 0;
	char *pch;
	char *pch2;
	
	do
	{
		if(fgets(buffer, 599, fp) == 0) return FAIL;
		buffer[599] = '\0';
		len = strlen(buffer);
	    if (buffer[len - 1] == '\n') buffer[len - 1] = '\0';

	    switch (commentflag) 
		{
		    case 0:
				if (strstr(buffer, "/*") != 0) 
				{
			    	commentflag = 1;
			    	if (strstr(buffer, "*/") != 0)
						commentflag = 2;
				}
				break;

		    case 1:
				if (strstr(buffer, "*/") != 0)
				    commentflag = 2;
				break;

			    case 2:
				if (strstr(buffer, "/*") != 0) 
				{
				    commentflag = 1;
				    if (strstr(buffer, "*/") != 0)
					commentflag = 2;
				}
				else
				    commentflag = 0;
				break;
			    }	
	}while(commentflag != 0);	

	//separate field name: token = "="
	if (strstr (buffer, fieldname) != 0)
	{
		pch = strtok (buffer,"=");

		while (pch != NULL)
		{
			pch2 = pch;
    		pch = strtok (NULL, "= ");
		}
		sscanf(pch2, format, inbuffer);
		return 1;//ok
	}
	else return 0; 
}


void PrintParameters()
{
	printf("\n--> DM Parameters from input file:\n");
	printf("ID = %d\n", ID);
	printf("sequence = %s\n", sequence);
	printf("nMol = %d\n", nMol);
	printf("BondLen = %d\n", BondLen);
	printf("ProtLen = %d\n", ProtLen);
	printf("nC = %d\n", nC);
	printf("LV = %lf\n", LV);
	printf("dt = %lf\n", dt);
	printf("cT = %lf\n", cT); 
	printf("displayInterval = %d\n", displayInterval);
	printf("stepLimit = %d\n", stepLimit);
	printf("step2rescVels = %d\n", step2rescVels);
	printf("temperature = %lf\n", temperature);
	printf("temperature_steps = %lf\n", temperature_steps);
	printf("report_file = %c\n", report_file);
	printf("step2report=%d\n", step2report);
	printf("printsummary = %c\nprintsummaryInterval=%d\nprintsummary2file=%c\n", printsummary, printsummaryInterval, printsummary2file);
	printf("Shake relaxation's Parameters:\n");
	printf("\tmaxCycle = %d\n\tconsPrec = %lf\n\tstep2shake = %d\n", maxCycle, consPrec, step2shake);
}


void initial_message()
{
	system("clear");
	printf("*************************************************************************\n");
	printf(" UNIVERSIDADE TECNOLOGICA FEDERAL DO PARANA (UTFPR)                  \n");
	printf(" LABORATORIO DE BIOINFORMATICA									   \n");
	printf(" CESAR M. VARGAS BENITEZ       									   \n");
	printf("*************************************************************************\n");

}


void show_coordinates(int part_id) 
{
	//This function show the coordinate of all particles or of a specific particle selected through the parameter part_id
	//if part_id < 0 --> show the coordinate of all particle
	int n;

	if (part_id < 0)
	{
		for (n=0; n < nMol; n++)
		{
			printf("Particule %d -- coordinates: \n\tx=%lf\n\ty=%lf\n\tz=%lf\n\t\n", n,  mol[n].r.x, mol[n].r.y, mol[n].r.z);
		}
	}
	else
	{
			printf("Particule %d -- coordinates: \n\tx=%lf\n\ty=%lf\n\tz=%lf\n\t\n", part_id,  mol[part_id].r.x, mol[part_id].r.y, mol[part_id].r.z);
	}
}


/*Creates de PDB file for a specific iteration */
void create_PDB(int iter)
{
	int n;
	FILE *file;	
	char *file_name;					
	char x_str[8], y_str[8], z_str[8];

	file_name = (char *)malloc(500*sizeof(char));
   	sprintf(file_name,"argon_iter%d.pdb", iter);
	file = fopen(file_name, "w");
	free(file_name);  			
	
	fprintf(file,"REMARK ARGON GAS COORDINATES - ITERATION %d\n", iter);

	for (n=0; n< nMol; n++)
	{
		if (mol[n].r.x<0)
		{
			if(mol[n].r.x>-10)	
				sprintf(x_str, "-0%2.3f", mol[n].r.x*(-1));
			else
				sprintf(x_str, "-%3.3f", mol[n].r.x*(-1));
		}
		else
		{
			if (mol[n].r.x<10)
				sprintf(x_str, " 0%2.3f", mol[n].r.x);
			else
				sprintf(x_str, " %3.3f", mol[n].r.x);	
		}	
		if (mol[n].r.y<0)
		{
			if(mol[n].r.y>-10)	
				sprintf(y_str, "-0%2.3f", mol[n].r.y*(-1));
			else
				sprintf(y_str, "-%3.3f", mol[n].r.y*(-1));
		}		
		else
		{
			if (mol[n].r.y<10)		
				sprintf(y_str, " 0%2.3f", mol[n].r.y);
			else
				sprintf(y_str, " %3.3f", mol[n].r.y);
		}

		if (mol[n].r.z<0)
		{
			if(mol[n].r.y>-10)	
				sprintf(z_str, "-0%2.3f", mol[n].r.z*(-1));
			else
				sprintf(z_str, "-%3.3f", mol[n].r.z*(-1));
		}		
		else
		{
			if (mol[n].r.z<10)		
				sprintf(z_str, " 0%2.3f", mol[n].r.z);
			else
				sprintf(z_str, " %3.3f", mol[n].r.z);
		}

		fprintf(file, "ATOM    %3d   AR ARG  %3d     %s  %s  %s  1.00  0.00           SOL\n", n+1, n+1, x_str, y_str, z_str); 	

	}
	fprintf(file, "END");

	fclose(file);

}



/* Step - Iteration process*/
/* function that handles the processing for a single step (or iteration) */
/* report and makeMolecules are also executed here */
int Step() 
{
	int finish;

	if (stepCount <= stepLimit)			
	{
		velocityVerlet();

		nCycleR	= 0;
		nCycleV = 0;
		
		if (stepCount % step2shake == 0){
			shake_relaxation();
		}

		berendsen_thermostat();	//Berendsen thermostat
		
		evaluate(); 

		if((printsummary == 'y' || printsummary == 'Y') && stepCount % printsummaryInterval == 0){
			print_summary();
		}
		if(   (report_file == 'y' || report_file == 'Y')   &&   (stepCount % step2report == 0 || stepCount == 0)   ){
			report();
		}
		if (stepCount % pathwaysstep == 0 || stepCount == 0){
			savePathwaysVetor(stepCount/3000);
		}
		stepCount++;
		finish = 0;
	}
	else	
	{
		tp = time(NULL) - ti; 	

		report();

		savePathwaysVetorArquivo();

		finish = 1;
	}
	return finish;
}

/*print summary*/
/*it shows the current temperature, potential energy and center of mass */
void print_summary()
{
	char *nome_arquivo;						
	FILE *results;

//	real rgH, rgP, rgAll, bondavg;
	calc_currentTemp();
	calc_centermass();
	//uSum = calc_Total_energy();

	uSum = uLJ + uTorsion + uChainAngles;
	
	rgH = RgH();
	rgP = RgP();
	rgAll = RgAll();
	bondavg = BondAverage();
	

	system("clear");
	printf("Step %d:\n\tTemp=%f\n\tTotal Lennard-Jones potential = %lf\n\tTotal Torsion potential = %lf\n\tTotal Chain angle potential = %lf\n\tTotal Energy - U: uSum = %lf\n", stepCount, current_temperature, uLJ, uTorsion, uChainAngles , uSum);
	printf("\tCenter of mass: = (%.2lf, %.2lf, %.2lf)\n", center_mass.x, center_mass.y, center_mass.z);
	printf("\tRgH = %lf\n\tRgP = %lf\n\tRg_All = %lf\n", rgH, rgP, rgAll);

	printf("\tBond Lenght Average = %lf\n", bondavg);


	if (printsummary2file == 'y')
	{
    	nome_arquivo = (char *)malloc(500*sizeof(char));
	    sprintf(nome_arquivo,"MD_ProteinAB_%s.txt", sequence); 
	    results = fopen(nome_arquivo, "a+");
	    free(nome_arquivo); 
	
		if (stepCount == 0) fprintf(results, "Step\tE_LJ\tE_torsion\tE_chainAngles\tE_total\tRgH\tRgP\tRgAll\n");
		fprintf(results, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", stepCount, uLJ, uTorsion, uChainAngles, uSum, rgH, rgP, rgAll);
	
		fclose(results);
	}
}

/* Report file generator*/
void report()
{
	FILE *file;	
	static char *file_name;					

	file_name = (char *)malloc(50*sizeof(char));
	sprintf(file_name,"report_%daminoacids_%d.txt", nMol, ID); 
	file = fopen(file_name, "a+");
	free(file_name);  	

	if (stepCount == 0) 
		fprintf(file, "Step\tt\tCurrent temperature\tKineticEnergy\tPotential Energy\tInternal Energy\tPressure\tDensity\tRgH\tRgP\tRgAll\tBond length Average\n");


	if (bondavg <= 1) 
			fprintf(file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", stepCount, stepCount * dt, current_temperature, kineticEnergy, potEnergy,  TotalEnergy, pressure, density, rgH, rgP, rgAll, bondavg);

	fclose(file);
}




void time_ex()
{	
	FILE *arqv;	
	arqv = fopen("time.txt", "a");				
												
	
	fprintf(arqv, "%d\t%d\t%d\n",nMol, ID, (int)(tp));

	fclose(arqv);
}




/* Center of mass calculation */
void calc_centermass()
{
	int n;
	VectorR r;

	r.x = 0;
	r.y = 0;
	r.z = 0;
    for (n = 0; n < nMol; n++)
	{
	    r.x += mol[n].r.x; 
		r.y +=  mol[n].r.y; 
	    r.z += mol[n].r.z;
	}

	center_mass.x = r.x / nMol; //r.x * mass/(nMol * mass);
	center_mass.y = r.y / nMol;
	center_mass.z = r.z / nMol;
	
}

real calc_Total_energy()
{
    VectorR r_ij;               
	int n, j;
	real r, U_LJ, U_LJ_pair;


	U_LJ = 0; //Lennard-Jones potential

    for (n = 0; n < ProtLen-2; n++)        
        for (j = n+2; j < ProtLen; j++) 
		{
            r_ij.x = mol[n].r.x - mol[j].r.x;
			r_ij.y = mol[n].r.y - mol[j].r.y;
			r_ij.z = mol[n].r.z - mol[j].r.z;
////////////
			if (r_ij.x >= 0.5 * LV) r_ij.x -= LV;
			else if (r_ij.x < -0.5 * LV) r_ij.x += LV;

			if (r_ij.y >= 0.5 * LV) r_ij.y -= LV;
			else if (r_ij.y < -0.5 * LV) r_ij.y += LV;
			if (r_ij.z >= 0.5 * LV) r_ij.z -= LV;
			else if (r_ij.z < -0.5 * LV) r_ij.z += LV;
////////////
	        r  = r_ij.x * r_ij.x + r_ij.y * r_ij.y;
	        r += r_ij.z * r_ij.z;

			r = sqrt(r);			
			            
			//Lennard-Jones potential:
			U_LJ_pair = 4 * (pow(r, -12) - pow(r, -6));	//from Lennard-Jones equation

			//if (sequence[n] == 'A' && sequence[j] == 'A')	//AA or HH
			//	U_LJ_pair *= 1;
			//else			
			if ((sequence[n] == 'A' && sequence[j] == 'B') || (sequence[n] == 'B' && sequence[j] == 'A')) //AB or HP
				U_LJ_pair *= 0.5; 
			else if (sequence[n] == 'B' && sequence[j] == 'B') //BB or PP
				U_LJ_pair *= 0.5;
	
			U_LJ += U_LJ_pair;
        }

	return U_LJ;
}

/* Velocity verlet*/
/* This function is not used for DM of water molecules */
void velocityVerlet() 
{
	int n;


    for (n= 0; n < nMol; n++)
	{
        mol[n].r.x += mol[n].rv.x * dt + 0.5 * mol[n].ra.x * dt * dt;
        mol[n].r.y += mol[n].rv.y * dt + 0.5 * mol[n].ra.y * dt * dt;
        mol[n].r.z += mol[n].rv.z * dt + 0.5 * mol[n].ra.z * dt * dt;

		//verify boundary conditions
		if (mol[n].r.x < 0)
		{
            mol[n].r.x += LV;//= 0; 
		}
        if (mol[n].r.x >= LV)
            mol[n].r.x  -= LV; //LV
		if (mol[n].r.y < 0)
            mol[n].r.y += LV; //0
        if (mol[n].r.y >= LV)
            mol[n].r.y -= LV;//LV
		if (mol[n].r.z < 0)
            mol[n].r.z += LV; //0
        if (mol[n].r.z >= LV)
            mol[n].r.z -= LV; //LV
        mol[n].rv.x += 0.5 * mol[n].ra.x * dt;
        mol[n].rv.y += 0.5 * mol[n].ra.y * dt;
        mol[n].rv.z += 0.5 * mol[n].ra.z * dt;
	}


    for (n = 0; n < nMol; n++)     	//initialize accelerations     
	{
    	mol[n].ra.x = 0;
		mol[n].ra.y = 0;
	   	mol[n].ra.z = 0;
	}


	computeAngleForces();	//chain angles forces

	computeTorsionForces();  //chain torsion forces

	computeLennard_Jones();  //Lennard-Jones potential


    for (n= 0; n < nMol; n++)					//update velocities
	{
        mol[n].rv.x += 0.5 * mol[n].ra.x * dt;
        mol[n].rv.y += 0.5 * mol[n].ra.y * dt;
        mol[n].rv.z += 0.5 * mol[n].ra.z * dt;
	}

}

/* Boundary verification */
/* Periodic boundary function*/
void verify_boundary() 
{
	int n;
  	for (n = 0; n < nMol; n++) 
 	{
		if (mol[n].r.x < LV) //
		{
            mol[n].r.x += LV;
		}
        if (mol[n].r.x > LV) //
            mol[n].r.x -= LV;
		if (mol[n].r.y < LV) //
            mol[n].r.y += LV;
        if (mol[n].r.y > LV) //
            mol[n].r.y -= LV;
		if (mol[n].r.z < LV) //
            mol[n].r.z += LV;
        if (mol[n].r.z > LV) //
            mol[n].r.z -= LV;
    }

}

/* Accelerations calc */
/* from AB 3D protein model Energy function */
/* The Energy function considerates bond, torsion and Lennard-Jones potential */
/* See the paper "Local interactions and protein folding:
A three-dimensional off-lattice approach" - (Irback et al, 1997) */
void computeAccelerations() 
{
    VectorR r_ij;               
	int n, j;
	real r, forceLJ;

    for (n = 0; n < ProtLen; n++)          
	{
    	mol[n].ra.x = 0;
		mol[n].ra.y = 0;
	   	mol[n].ra.z = 0;
	}


    for (n = 0; n < ProtLen-2; n++)        
        for (j = n+2; j < ProtLen; j++) 
		{
            r_ij.x = mol[n].r.x - mol[j].r.x;
			r_ij.y = mol[n].r.y - mol[j].r.y;
			r_ij.z = mol[n].r.z - mol[j].r.z;
////////////
			if (r_ij.x >= 0.5 * LV) r_ij.x -= LV;
			else if (r_ij.x < -0.5 * LV) r_ij.x += LV;

			if (r_ij.y >= 0.5 * LV) r_ij.y -= LV;
			else if (r_ij.y < -0.5 * LV) r_ij.y += LV;
			if (r_ij.z >= 0.5 * LV) r_ij.z -= LV;
			else if (r_ij.z < -0.5 * LV) r_ij.z += LV;
////////////
	        r  = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
			r = sqrt(r);			
			            
			//Lennard-Jones potential:
			forceLJ = 24 * (2 * pow(r, -13) - pow(r, -7));	//from Lennard-Jones equation

			if (sequence[n] == 'A' && sequence[j] == 'A')	//AA or HH
				forceLJ *= 1;
			else if ((sequence[n] == 'A' && sequence[j] == 'B') || (sequence[n] == 'B' && sequence[j] == 'A')) //AB or HP
				forceLJ *= 1/2;
			else if (sequence[n] == 'B' && sequence[j] == 'B') //BB or PP
				forceLJ *= 1/2;

            mol[n].ra.x += r_ij.x * forceLJ;
            mol[n].ra.y += r_ij.y * forceLJ;
            mol[n].ra.z += r_ij.z * forceLJ;
            mol[j].ra.x -= r_ij.x * forceLJ;
            mol[j].ra.y -= r_ij.y * forceLJ;
            mol[j].ra.z -= r_ij.z * forceLJ;
        }
}


/* Compute Lennard Jones Potential */
void computeLennard_Jones() 
{
    VectorR r_ij;               
	int n, j;
	real r2, forceLJ, U_LJ_pair;

	uLJ = 0;

    for (n = 0; n < ProtLen-2; n++)        
        for (j = n+2; j < ProtLen; j++) 
		{
            r_ij.x = mol[n].r.x - mol[j].r.x;
			r_ij.y = mol[n].r.y - mol[j].r.y;
			r_ij.z = mol[n].r.z - mol[j].r.z;

			if (r_ij.x >= 0.5 * LV) r_ij.x -= LV;
			else if (r_ij.x < -0.5 * LV) r_ij.x += LV;

			if (r_ij.y >= 0.5 * LV) r_ij.y -= LV;
			else if (r_ij.y < -0.5 * LV) r_ij.y += LV;
			if (r_ij.z >= 0.5 * LV) r_ij.z -= LV;
			else if (r_ij.z < -0.5 * LV) r_ij.z += LV;

	        r2  = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
			//r = sqrt(r);			
			            
			//Lennard-Jones potential:
			forceLJ = 24 * (2 * pow(r2, -7) - pow(r2, -4));	//from Lennard-Jones equation

			//Lennard-Jones potential:
			U_LJ_pair = 4 * (pow(r2, -6) - pow(r2, -3));	//from Lennard-Jones equation

			if ((sequence[n] == 'A' && sequence[j] == 'B') || (sequence[n] == 'B' && sequence[j] == 'A')) //AB or HP
			{
				forceLJ *= 0.5; //-0.5
				U_LJ_pair *= 0.5; 
			}
			else if (sequence[n] == 'B' && sequence[j] == 'B') //BB or PP
			{
				forceLJ *= 0.5;
				U_LJ_pair *= 0.5;
			}

            mol[n].ra.x += r_ij.x * forceLJ;
            mol[n].ra.y += r_ij.y * forceLJ;
            mol[n].ra.z += r_ij.z * forceLJ;
            mol[j].ra.x -= r_ij.x * forceLJ;
            mol[j].ra.y -= r_ij.y * forceLJ;
            mol[j].ra.z -= r_ij.z * forceLJ;


			uLJ += U_LJ_pair;
			//
        }


}


/* Compute chain torsion forces */
void computeTorsionForces() 
{
	VectorR dr1, dr2, dr3;
	
	real c11, c12, c13, c22, c23, c33, pi, qia, qib;
	real cr1, cr2;
	real f;
	real t1, t2, t3, t4, t5, t6;
	VectorR fi1, fi2;
	int i;
	
	uTorsion = 0;

	for (i = 0; i < ProtLen -3; i++)
	{
		dr1.x = mol[i + 1].r.x - mol[i].r.x;
		dr1.y = mol[i + 1].r.y - mol[i].r.y;
		dr1.z = mol[i + 1].r.z - mol[i].r.z;
		dr2.x = mol[i + 2].r.x - mol[i + 1].r.x;
		dr2.y = mol[i + 2].r.y - mol[i + 1].r.y;
		dr2.z = mol[i + 2].r.z - mol[i + 1].r.z;
		dr3.x = mol[i + 3].r.x - mol[i + 2].r.x;
		dr3.y = mol[i + 3].r.y - mol[i + 2].r.y;
		dr3.z = mol[i + 3].r.z - mol[i + 2].r.z;


		//boundary conditions
		if (dr1.x > 0.5 * LV)  dr1.x -= LV;
		else if (dr1.x < -0.5 * LV)  	dr1.x += LV;
		if (dr1.y > 0.5 * LV)  dr1.y -= LV;
		else if (dr1.y < -0.5 * LV)  	dr1.y += LV;
		if (dr1.z > 0.5 * LV)  dr1.z -= LV;
		else if (dr1.z < -0.5 * LV)  	dr1.z += LV;


		if (dr2.x > 0.5 * LV)  dr2.x -= LV;
		else if (dr2.x < -0.5 * LV)  	dr2.x += LV;
		if (dr2.y > 0.5 * LV)  dr2.y -= LV;
		else if (dr2.y < -0.5 * LV)  	dr2.y += LV;
		if (dr2.z > 0.5 * LV)  dr2.z -= LV;
		else if (dr2.z < -0.5 * LV)  	dr2.z += LV;
	
		if (dr3.x > 0.5 * LV)  dr3.x -= LV;
		else if (dr3.x < -0.5 * LV)  	dr3.x += LV;
		if (dr3.y > 0.5 * LV)  dr3.y -= LV;
		else if (dr3.y < -0.5 * LV)  	dr3.y += LV;
		if (dr3.z > 0.5 * LV)  dr3.z -= LV;
		else if (dr3.z < -0.5 * LV)  	dr3.z += LV;


		c11 = dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z;
		c12 = dr1.x * dr2.x + dr1.y * dr2.y + dr1.z * dr2.z;	
		c22 = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;
		c13 = dr1.x * dr3.x + dr1.y * dr3.y + dr1.z * dr3.z;	
		c23 = dr2.x * dr3.x + dr2.y * dr3.y + dr2.z * dr3.z;
		c33 = dr3.x * dr3.x + dr3.y * dr3.y + dr3.z * dr3.z;		
		pi = c13 * c22 - c12 * c23;
		qia = c11 * c22 - c12 * c12;
		qib = c22 * c33 - c23 * c23;
		cr1 = c12 / c22;
		cr2 = c23 / c22;

		t1 = pi; 
		t2 = c11 * c23 - c12 * c13;		
		t3 = -qia;
		t4 = qib;
		t5 = c13 * c23 - c12 * c33;
		t6 = -t1;


		f = 0.5; 

		fi1.x = f * c22 * (t1 * dr1.x  + t2 * dr2.x + t3 * dr3.x) / ((sqrt(qia * qib)) * qia);  
		fi1.y = f * c22 * (t1 * dr1.y  + t2 * dr2.y + t3 * dr3.y) / ((sqrt(qia * qib)) * qia);
		fi1.z = f * c22 * (t1 * dr1.z  + t2 * dr2.z + t3 * dr3.z) / ((sqrt(qia * qib)) * qia);

		mol[i].ra.x = mol[i].ra.x + fi1.x;														
		mol[i].ra.y = mol[i].ra.y + fi1.y;
		mol[i].ra.z = mol[i].ra.z + fi1.z;

		fi2.x = f * c22 * (t4 * dr1.x  + t5 * dr2.x + t6 * dr3.x) / ((sqrt(qia * qib)) * qib); //f(i+1)
		fi2.y = f * c22 * (t4 * dr1.y  + t5 * dr2.y + t6 * dr3.y) / ((sqrt(qia * qib)) * qib);
		fi2.z = f * c22 * (t4 * dr1.z  + t5 * dr2.z + t6 * dr3.z) / ((sqrt(qia * qib)) * qib);

		mol[i + 1].ra.x = mol[i + 1].ra.x + (- (1. + cr1) * fi1.x) + cr2 * fi2.x;				//mol[i+1] :  (i-1)-atom
		mol[i + 1].ra.y = mol[i + 1].ra.y + (- (1. + cr1) * fi1.y) + cr2 * fi2.y;
		mol[i + 1].ra.z = mol[i + 1].ra.z + (- (1. + cr1) * fi1.z) + cr2 * fi2.z;

		mol[i + 2].ra.x = mol[i + 2].ra.x + cr1 * fi1.x + (- (1. + cr2)) * fi2.x;				//mol[i+2] : (i)-atom
		mol[i + 2].ra.y = mol[i + 2].ra.y + cr1 * fi1.y + (- (1. + cr2)) * fi2.y;
		mol[i + 2].ra.z = mol[i + 2].ra.z + cr1 * fi1.z + (- (1. + cr2)) * fi2.z;

		mol[i + 3].ra.x = mol[i + 3].ra.x + fi2.x;												//mol[i+3] : (i+1)-atom
		mol[i + 3].ra.y = mol[i + 3].ra.y + fi2.y;
		mol[i + 3].ra.z = mol[i + 3].ra.z + fi2.z;


		//Torsion potential
		uTorsion +=  (-0.5) * c13; 		//-k2 * bi . bi+2		

	}

}


/* Compute chain angle forces */
void computeAngleForces() 
{
	VectorR dr1, dr2;
	real c, cd, c11, c22, c12, f;
	int i;
	VectorR fi1, fi2;

	uChainAngles = 0;

	for (i = 0; i < ProtLen - 2; i++)
	{
		dr1.x = mol[i + 1].r.x - mol[i].r.x;
		dr1.y = mol[i + 1].r.y - mol[i].r.y;
		dr1.z = mol[i + 1].r.z - mol[i].r.z;
		dr2.x = mol[i + 2].r.x - mol[i + 1].r.x;
		dr2.y = mol[i + 2].r.y - mol[i + 1].r.y;
		dr2.z = mol[i + 2].r.z - mol[i + 1].r.z;


		//boundary conditions
		if (dr1.x > 0.5 * LV)  dr1.x -= LV;
		else if (dr1.x < -0.5 * LV)  	dr1.x += LV;
		if (dr1.y > 0.5 * LV)  dr1.y -= LV;
		else if (dr1.y < -0.5 * LV)  	dr1.y += LV;
		if (dr1.z > 0.5 * LV)  dr1.z -= LV;
		else if (dr1.z < -0.5 * LV)  	dr1.z += LV;


		if (dr2.x > 0.5 * LV)  dr2.x -= LV;
		else if (dr2.x < -0.5 * LV)  	dr2.x += LV;
		if (dr2.y > 0.5 * LV)  dr2.y -= LV;
		else if (dr2.y < -0.5 * LV)  	dr2.y += LV;
		if (dr2.z > 0.5 * LV)  dr2.z -= LV;
		else if (dr2.z < -0.5 * LV)  	dr2.z += LV;
	
		c11 = dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z;
		c12 = dr1.x * dr2.x + dr1.y * dr2.y + dr1.z * dr2.z;	
		c22 = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;

		cd = sqrt (c11 * c22);
		c  = c12 / cd;

		f = -1 ; //* sqrt(c11) * sqrt(c22); // - (- k1 * abs(bi) * abs(b_(i+1)) * 1) ;  component
		
		fi1.x = (f / cd ) * ((c12 / c11) * dr1.x  - dr2.x);
		fi1.y = (f / cd ) * ((c12 / c11) * dr1.y  - dr2.y);
		fi1.z = (f / cd ) * ((c12 / c11) * dr1.z  - dr2.z);

		fi2.x = (f / cd) * (dr1.x - (c12 / c22) * dr2.x); 
		fi2.y = (f / cd) * (dr1.y - (c12 / c22) * dr2.y); 
		fi2.z = (f / cd) * (dr1.z - (c12 / c22) * dr2.z); 

		mol[i].ra.x = mol[i].ra.x + fi1.x;
		mol[i].ra.y = mol[i].ra.y + fi1.y;
		mol[i].ra.z = mol[i].ra.z + fi1.z;

		mol[i + 1].ra.x = mol[i + 1].ra.x - fi1.x - fi2.x;
		mol[i + 1].ra.y = mol[i + 1].ra.y - fi1.y - fi2.y;
		mol[i + 1].ra.z = mol[i + 1].ra.z - fi1.z - fi2.z;		

		mol[i + 2].ra.x = mol[i + 2].ra.x + fi2.x;
		mol[i + 2].ra.y = mol[i + 2].ra.y + fi2.y;
		mol[i + 2].ra.z = mol[i + 2].ra.z + fi2.z;

		//Chain Angle potential
		uChainAngles += c12; // -k1 *  (bi . bi+1) ; k1 = -1

		//
	}
}

//Berendsen Thermostat
void berendsen_thermostat()
{	
    real sum;
	real lambda;
	int i;
	real temp;

	sum = 0;
    for (i = 0; i < nMol; i++)
	{
        sum += mol[i].rv.x * mol[i].rv.x + mol[i].rv.y * mol[i].rv.y + mol[i].rv.z * mol[i].rv.z;
	}

	temp = sum / (3 * (nMol-1));
	lambda = sqrt(1 + (dt/ cT)  * (temperature / temp   - 1)); 

    for (i = 0; i < nMol; i++)
	{
        mol[i].rv.x *= lambda;
        mol[i].rv.y *= lambda;
        mol[i].rv.z *= lambda;
	}

	/*printf("lambda=%lf, temperature=%lf, temp=%lf\n", lambda, temperature, temp);*/

}

/* rescaleVelocities from a given temperature */
//velocity rescaling thermostat
void rescaleVelocities() 
{
    real sum = 0;
	real lambda;
	int i;
	
	VectorR w;    


	velMag = sqrt(NDIM * (1. - 1. / nMol) * temperature);  //update velMag 
														   	

    for (i = 0; i < nMol; i++)
	{
		sum += mol[i].rv.x * mol[i].rv.x + mol[i].rv.y * mol[i].rv.y + mol[i].rv.z * mol[i].rv.z;
	}

	lambda = velMag / sqrt(sum / nMol);

    for (i = 0; i < nMol; i++)
	{
        mol[i].rv.x *= lambda;
        mol[i].rv.y *= lambda;
        mol[i].rv.z *= lambda;
	}


	sum = 0;
    for (i = 0; i < nMol; i++)
	{
		calc_AngVelocities(i, &w);		
        sum += mInert.x * w.x * w.x +  mInert.y * w.y * w.y + mInert.z * w.z * w.z;
	}

	lambda = velMag / sqrt(sum / nMol);
    
    for (i = 0; i < nMol; i++)
	{
        mol[i].qv.u1 *= lambda;
        mol[i].qv.u2 *= lambda;
        mol[i].qv.u3 *= lambda;
        mol[i].qv.u4 *= lambda;
	}
}


/* Function calc_AngVelocities */
/*this function calculates the angular velocities (w) for each molecule*/
void calc_AngVelocities(int n, VectorR *w)
{
	Quat qt, qvt;

	qvt = mol[n].qv;
	qvt.u4 *= -1;

	qt.u1 = qvt.u4 * mol[n].q.u1 - qvt.u3 * mol[n].q.u2 + qvt.u2 * mol[n].q.u3 + qvt.u1 * mol[n].q.u4; 
	qt.u2 = qvt.u3 * mol[n].q.u1 + qvt.u4 * mol[n].q.u2 - qvt.u1 * mol[n].q.u3 + qvt.u2 * mol[n].q.u4; 
	qt.u3 = -qvt.u2 * mol[n].q.u1 + qvt.u1 * mol[n].q.u2 + qvt.u4 * mol[n].q.u3 + qvt.u3 * mol[n].q.u4; 
	
	//scale
	qt.u1 = 2 * qt.u1;
	qt.u2 = 2 * qt.u2;
	qt.u3 = 2 * qt.u3;

	(*w).x = qt.u1;
	(*w).y = qt.u2;
	(*w).z = qt.u3;

}

/*Function calc_AccelQuat */
/*this function calculates the quaternions' accelerations (angular)*/
void calc_AccelQuat()
{
	Quat qs;
	VectorR w;
	int n;

    for (n = 0; n < nMol; n++)
	{
		calc_AngVelocities(n, &w);			

		qs.u1 = (mol[n].torq.x + (mInert.y - mInert.z) * w.y * w.z) / mInert.x;
		qs.u2 = (mol[n].torq.y + (mInert.z - mInert.x) * w.z * w.x) / mInert.y;
		qs.u3 = (mol[n].torq.z + (mInert.x - mInert.y) * w.x * w.y) / mInert.z;
		qs.u4 = -2. * (mol[n].qv.u1 * mol[n].qv.u1 + mol[n].qv.u2 * mol[n].qv.u2 + mol[n].qv.u3 * mol[n].qv.u3 + mol[n].qv.u4 * mol[n].qv.u4);

		mol[n].qa.u1 =   mol[n].q.u4 * qs.u1 - mol[n].q.u3 * qs.u2 + mol[n].q.u2 * qs.u3 + mol[n].q.u1 * qs.u4;
		mol[n].qa.u2 = 	 mol[n].q.u3 * qs.u1 + mol[n].q.u4 * qs.u2 - mol[n].q.u1 * qs.u3 + mol[n].q.u2 * qs.u4;
		mol[n].qa.u3 = - mol[n].q.u2 * qs.u1 + mol[n].q.u1 * qs.u2 + mol[n].q.u4 * qs.u3 + mol[n].q.u3 * qs.u4;
		mol[n].qa.u4 = - mol[n].q.u1 * qs.u1 - mol[n].q.u2 * qs.u2 - mol[n].q.u3 * qs.u3 + mol[n].q.u4 * qs.u4;

		mol[n].qa.u1 = mol[n].qa.u1 * 0.5;
		mol[n].qa.u2 = mol[n].qa.u2 * 0.5;
		mol[n].qa.u3 = mol[n].qa.u3 * 0.5;
		mol[n].qa.u4 = mol[n].qa.u4 * 0.5;

	}	

}

/*Function: calc_Torqs */
/*This function calculates the torques and translational accelerations*/
void calc_Torqs()
{
	RotMat rMat; 	
	VectorR dr, t, torqS;
	int i, n;

    for (n = 0; n < nMol; n++)
	{
		mol[n].ra.x = 0;
		mol[n].ra.y = 0;
		mol[n].ra.z = 0;

		torqS.x = 0;
		torqS.y = 0;
		torqS.z = 0;

		for (i = 0; i < SITESMOL; i++) 
		{
			mol[n].ra.x = mol[n].ra.x + site[n * SITESMOL + i].f.x;  //m=1
			mol[n].ra.y = mol[n].ra.y + site[n * SITESMOL + i].f.y;
			mol[n].ra.z = mol[n].ra.z + site[n * SITESMOL + i].f.z;

			dr.x =  site[n * SITESMOL + i].r.x - mol[n].r.x;
			dr.y =  site[n * SITESMOL + i].r.y - mol[n].r.y;
			dr.z =  site[n * SITESMOL + i].r.z - mol[n].r.z;

			t.x  =  dr.y * site[n *SITESMOL + i].f.z - dr.z * site[n *SITESMOL + i].f.y;
			t.y  =  dr.z * site[n *SITESMOL + i].f.x - dr.x * site[n *SITESMOL + i].f.z;
			t.z  =  dr.x * site[n *SITESMOL + i].f.y - dr.y * site[n *SITESMOL + i].f.x;


			torqS.x = torqS.x + t.x;
			torqS.y = torqS.y + t.y;
			torqS.z = torqS.z + t.z;
		}

		BuildRotMatrix (&rMat, &mol[n].q, 0);

		mol[n].torq.x = rMat.u[0] * torqS.x + rMat.u[3] * torqS.y + rMat.u[6] * torqS.z;
		mol[n].torq.y = rMat.u[1] * torqS.x + rMat.u[4] * torqS.y + rMat.u[7] * torqS.z;
		mol[n].torq.z = rMat.u[2] * torqS.x + rMat.u[5] * torqS.y + rMat.u[8] * torqS.z;

	}
}


/*BuildRotMatrix: generates the transpose matrix for each molecule from the quaternion components*/
void BuildRotMatrix(RotMat *rMat, Quat *q, int transpose)
{
	real p[10], tq[4], s;
	int k, k1, k2;

	tq[0] = q->u1;
	tq[1] = q->u2;
	tq[2] = q->u3;
	tq[3] = q->u4;

	for(k=0, k2=0; k2<4; k2++)
	{
		for(k1 = k2; k1<4; k1++, k++)	
		{	
			p[k] = 2 * tq[k1] * tq[k2];
		}
	}
	
	rMat->u[0] = p[0] + p[9] - 1;
	rMat->u[4] = p[4] + p[9] - 1;
	rMat->u[8] = p[7] + p[9] - 1;
	s = transpose ? 1 : -1;
	rMat->u[1] = p[1] + s * p[8];
	rMat->u[3] = p[1] - s * p[8];
	rMat->u[2] = p[2] - s * p[6];
	rMat->u[6] = p[2] + s * p[6];
	rMat->u[5] = p[5] + s * p[3];
	rMat->u[7] = p[5] - s * p[3];

}

/* Sites coordinates generations */
void GenSiteCoords()
{
	RotMat rMat;
	VectorR t;
	int i, n;

	for (n = 0; n < nMol; n++)
	{
		BuildRotMatrix (&rMat, &mol[n].q, 1);

		for (i=0; i<SITESMOL; i++)
		{
			t.x = rMat.u[0] * mSite[i].r.x + rMat.u[3] * mSite[i].r.y + rMat.u[6] * mSite[i].r.z;
			t.y = rMat.u[1] * mSite[i].r.x + rMat.u[4] * mSite[i].r.y + rMat.u[7] * mSite[i].r.z;
			t.z = rMat.u[2] * mSite[i].r.x + rMat.u[5] * mSite[i].r.y + rMat.u[8] * mSite[i].r.z;

			site[SITESMOL * n + i].r.x =  mol[n].r.x + t.x;
			site[SITESMOL * n + i].r.y =  mol[n].r.y + t.y;
			site[SITESMOL * n + i].r.z =  mol[n].r.z + t.z;
		}
	}

}

/* PC method : predictor-corrector functions */
/* Predictor step */
/* Adams–Bashforth formula */
/* Second-order equation */
void PredictorStep()
{
	real cr[] = {19., -10., 3.};
	real cv[] = {27., -22., 7.};
	real div = 24.;
	real wr, wv;
	int n; 

	wr = (dt * dt) / div;
	wv = dt / div;

	
	for (n = 0; n < nMol; n++)
	{
		mol[n].ro = mol[n].r; 
		mol[n].rvo = mol[n].rv;

		mol[n].r.x = mol[n].r.x + dt * mol[n].rv.x + wr * (cr[0] * mol[n].ra.x + cr[1] * mol[n].ra1.x + cr[2] * mol[n].ra2.x);
		mol[n].rv.x = (mol[n].r.x - mol[n].ro.x) / dt + wv * (cv[0] * mol[n].ra.x + cv[1] * mol[n].ra1.x + cv[2] * mol[n].ra2.x);
		mol[n].r.y = mol[n].r.y + dt * mol[n].rv.y + wr * (cr[0] * mol[n].ra.y + cr[1] * mol[n].ra1.y + cr[2] * mol[n].ra2.y);
		mol[n].rv.y = (mol[n].r.y - mol[n].ro.y) / dt + wv * (cv[0] * mol[n].ra.y + cv[1] * mol[n].ra1.y + cv[2] * mol[n].ra2.y);
		mol[n].r.z = mol[n].r.z + dt * mol[n].rv.z + wr * (cr[0] * mol[n].ra.z + cr[1] * mol[n].ra1.z + cr[2] * mol[n].ra2.z);
		mol[n].rv.z = (mol[n].r.z - mol[n].ro.z) / dt + wv * (cv[0] * mol[n].ra.z + cv[1] * mol[n].ra1.z + cv[2] * mol[n].ra2.z);
		mol[n].ra2 = mol[n].ra1;
		mol[n].ra1 = mol[n].ra;

	} 

}


/* Predictor step for quaternions*/
/* Adams–Bashforth formula */
/* Second-order equation */
void PredictorStepQ()
{
	real cr[] = {19., -10., 3.};
	real cv[] = {27., -22., 7.}; 
	real div = 24.;
	real wr, wv;
	int n;

	wr = (dt * dt) / div;
	wv = dt / div;

	for (n = 0; n < nMol; n++)
	{
		mol[n].qo = mol[n].q;
		mol[n].qvo = mol[n].qv;

		mol[n].q.u1 = mol[n].q.u1 + dt * mol[n].qv.u1 + wr * (cr[0] * mol[n].qa.u1 + cr[1] * mol[n].qa1.u1 + cr[2] * mol[n].qa2.u1);
		mol[n].q.u2 = mol[n].q.u2 + dt * mol[n].qv.u2 + wr * (cr[0] * mol[n].qa.u2 + cr[1] * mol[n].qa1.u2 + cr[2] * mol[n].qa2.u2);
		mol[n].q.u3 = mol[n].q.u3 + dt * mol[n].qv.u3 + wr * (cr[0] * mol[n].qa.u3 + cr[1] * mol[n].qa1.u3 + cr[2] * mol[n].qa2.u3);
		mol[n].q.u4 = mol[n].q.u4 + dt * mol[n].qv.u4 + wr * (cr[0] * mol[n].qa.u4 + cr[1] * mol[n].qa1.u4 + cr[2] * mol[n].qa2.u4);

		mol[n].qv.u1 = (mol[n].q.u1 - mol[n].qo.u1) / dt + wv * (cv[0] * mol[n].qa.u1 + cv[1] * mol[n].qa1.u1 + cv[2] * mol[n].qa2.u1);
		mol[n].qv.u2 = (mol[n].q.u2 - mol[n].qo.u2) / dt + wv * (cv[0] * mol[n].qa.u2 + cv[1] * mol[n].qa1.u2 + cv[2] * mol[n].qa2.u2);
		mol[n].qv.u3 = (mol[n].q.u3 - mol[n].qo.u3) / dt + wv * (cv[0] * mol[n].qa.u3 + cv[1] * mol[n].qa1.u3 + cv[2] * mol[n].qa2.u3);
		mol[n].qv.u4 = (mol[n].q.u4 - mol[n].qo.u4) / dt + wv * (cv[0] * mol[n].qa.u4 + cv[1] * mol[n].qa1.u4 + cv[2] * mol[n].qa2.u4);
		
		mol[n].qa2 = mol[n].qa1;
		mol[n].qa1 = mol[n].qa;

	}

}

/* Corrector step */
/* Adams–Bashforth formula */
/* Second-order equation */
void CorrectorStep()
{
	real cr[] = {3., 10., -1.};
	real cv[] = {7., 6., -1.}; 
	real div = 24.; 
	real wr, wv;
	int n;

	wr = (dt*dt) / div;
	wv = dt / div; 

	for (n = 0; n < nMol; n++)
	{
		mol[n].r.x = mol[n].ro.x + dt * mol[n].rvo.x + wr * (cr[0] * mol[n].ra.x + cr[1] * mol[n].ra1.x + cr[2] * mol[n].ra2.x);
		mol[n].rv.x = (mol[n].r.x -  mol[n].ro.x) / dt + wv * (cv[0] * mol[n].ra.x + cv[1] * mol[n].ra1.x + cv[2] * mol[n].ra2.x);
		mol[n].r.y = mol[n].ro.y + dt * mol[n].rvo.y + wr * (cr[0] * mol[n].ra.y + cr[1] * mol[n].ra1.y + cr[2] * mol[n].ra2.y);
		mol[n].rv.y = (mol[n].r.y -  mol[n].ro.y) / dt + wv * (cv[0] * mol[n].ra.y + cv[1] * mol[n].ra1.y + cv[2] * mol[n].ra2.y);
		mol[n].r.z = mol[n].ro.z + dt * mol[n].rvo.z + wr * (cr[0] * mol[n].ra.z + cr[1] * mol[n].ra1.z + cr[2] * mol[n].ra2.z);
		mol[n].rv.z = (mol[n].r.z -  mol[n].ro.z) / dt + wv * (cv[0] * mol[n].ra.z + cv[1] * mol[n].ra1.z + cv[2] * mol[n].ra2.z);
	}
}

/* Corrector step for quaternions */
/* Adams–Bashforth formula */
/* Second-order equation */
void CorrectorStepQ()
{
	real cr[] = {3., 10., -1.};
	real cv[] = {7., 6., -1.};
	real div = 24.;
	real wr, wv;
	int n;

	wr = (dt * dt) / div;
	wv = dt / div;

	for (n = 0; n < nMol; n++)
	{
		mol[n].q.u1 = mol[n].qo.u1 + dt * mol[n].qvo.u1 + wr * (cr[0] * mol[n].qa.u1 + cr[1] * mol[n].qa1.u1 + cr[2] * mol[n].qa2.u1);
		mol[n].q.u2 = mol[n].qo.u2 + dt * mol[n].qvo.u2 + wr * (cr[0] * mol[n].qa.u2 + cr[1] * mol[n].qa1.u2 + cr[2] * mol[n].qa2.u2);
		mol[n].q.u3 = mol[n].qo.u3 + dt * mol[n].qvo.u3 + wr * (cr[0] * mol[n].qa.u3 + cr[1] * mol[n].qa1.u3 + cr[2] * mol[n].qa2.u3);
		mol[n].q.u4 = mol[n].qo.u4 + dt * mol[n].qvo.u4 + wr * (cr[0] * mol[n].qa.u4 + cr[1] * mol[n].qa1.u4 + cr[2] * mol[n].qa2.u4);

		mol[n].qv.u1 = (mol[n].q.u1 - mol[n].qo.u1) / dt + wv * (cv[0] * mol[n].qa.u1 + cv[1] * mol[n].qa1.u1 + cv[2] * mol[n].qa2.u1);
		mol[n].qv.u2 = (mol[n].q.u2 - mol[n].qo.u2) / dt + wv * (cv[0] * mol[n].qa.u2 + cv[1] * mol[n].qa1.u2 + cv[2] * mol[n].qa2.u2);
		mol[n].qv.u3 = (mol[n].q.u3 - mol[n].qo.u3) / dt + wv * (cv[0] * mol[n].qa.u3 + cv[1] * mol[n].qa1.u3 + cv[2] * mol[n].qa2.u3);
		mol[n].qv.u4 = (mol[n].q.u4 - mol[n].qo.u4) / dt + wv * (cv[0] * mol[n].qa.u4 + cv[1] * mol[n].qa1.u4 + cv[2] * mol[n].qa2.u4);

	}
}



//AdjustQuat: normalization of the quaternions to prevent numerical error
void AdjustQuat()
{
	real qi;
	int n;

	for (n = 0; n < nMol; n++)
	{
		qi = 1. / sqrt (mol[n].q.u1 * mol[n].q.u1 + mol[n].q.u2 * mol[n].q.u2 + mol[n].q.u3 * mol[n].q.u3 + mol[n].q.u4 * mol[n].q.u4);
		mol[n].q.u1 = mol[n].q.u1 * qi;
		mol[n].q.u2 = mol[n].q.u2 * qi;
		mol[n].q.u3 = mol[n].q.u3 * qi;
		mol[n].q.u4 = mol[n].q.u4 * qi;
	}
}


/* Mathematical thermostat */
void thermostat()
{
	real s1, s2, vFac;
	VectorR w;
	int n; 

	s1 = 0.;
	s2 = 0.;

	for (n = 0; n < nMol; n++)
	{
		s1 += mol[n].rv.x * mol[n].ra.x + mol[n].rv.y * mol[n].ra.y + mol[n].rv.z * mol[n].ra.z;
		s2 += mol[n].rv.x * mol[n].rv.x + mol[n].rv.y * mol[n].rv.y + mol[n].rv.z * mol[n].rv.z;
	}

	for (n = 0; n < nMol; n++)
	{
		calc_AngVelocities(n, &w);
		s1 += w.x * mol[n].torq.x + w.y * mol[n].torq.y + w.z * mol[n].torq.z;
		s2 += mInert.x * w.x * w.x + mInert.y * w.y * w.y + mInert.z * w.z * w.z;
	}
	
	vFac = - s1 / s2;

	for (n = 0; n < nMol; n++)
	{
		mol[n].ra.x = mol[n].ra.x + vFac * mol[n].rv.x; 
		mol[n].ra.y = mol[n].ra.y + vFac * mol[n].rv.y;
		mol[n].ra.z = mol[n].ra.z + vFac * mol[n].rv.z;

		mol[n].qa.u1 = mol[n].qa.u1 + vFac * mol[n].qv.u1; 
		mol[n].qa.u2 = mol[n].qa.u2 + vFac * mol[n].qv.u2;
		mol[n].qa.u3 = mol[n].qa.u3 + vFac * mol[n].qv.u3;
		mol[n].qa.u4 = mol[n].qa.u4 + vFac * mol[n].qv.u4;
	}

}


/* Angular Coordinates initialization */
void InitAngCoords()
{
	VectorR e;
	real eAng[3];
	int n; 
	
	for (n = 0; n < nMol; n++)
	{
		VRand (&e);
		eAng[0] = atan2 (e.x, e.y);
		eAng[1] = acos (e.z); 
		eAng[2] = 2. * PI * RandR();
		EulerToQuat (&mol[n].q, eAng);
	}

}


//Random Number generator --> TODO: Mersenne Twister

real RandR()
{
	int randSeedP = rand() % 20; //17;
	
	randSeedP = (randSeedP * IMUL + IADD) & MASK; 
	return (randSeedP * SCALE);
}


/* Euler to quaternions */
/* Angular coordinates are coverted to quaternions form */
void EulerToQuat(Quat *qe, real *eAng)
{
	real a1, a2, a3;
	a1 = 0.5 * eAng[1]; 
	a2 = 0.5 * (eAng[0] - eAng[2]);
	a3 = 0.5 * (eAng[0] + eAng[2]);

	(*qe).u1 = sin (a1) * cos (a2);
	(*qe).u2 = sin (a1) * sin (a2);
	(*qe).u3 = cos (a1) * sin (a3);
	(*qe).u4 = cos (a1) * cos (a3);
}

/* Angular velocities initialization based on temperature (or velMag) */
void InitAngVels()
{
	Quat qe;
	VectorR e; 
	real f;
	int n;

	for (n = 0; n < nMol; n++)
	{
		VRand (&e); 
		qe.u1 = e.x;
		qe.u2 = e.y;
		qe.u3 = e.z;
		qe.u4 = 0;

		mol[n].qv.u1 =   mol[n].q.u4 * qe.u1 - mol[n].q.u3 * qe.u2 + mol[n].q.u2 * qe.u3 + mol[n].q.u1 * qe.u4; 
		mol[n].qv.u2 = 	 mol[n].q.u3 * qe.u1 + mol[n].q.u4 * qe.u2 - mol[n].q.u1 * qe.u3 + mol[n].q.u2 * qe.u4; 
		mol[n].qv.u3 = - mol[n].q.u2 * qe.u1 + mol[n].q.u1 * qe.u2 + mol[n].q.u4 * qe.u3 + mol[n].q.u3 * qe.u4;
		mol[n].qv.u4 = - mol[n].q.u1 * qe.u1 - mol[n].q.u2 * qe.u2 - mol[n].q.u3 * qe.u3 + mol[n].q.u4 * qe.u4;

		f = 0.5 * velMag / sqrt(mInert.x * e.x * e.x + mInert.y * e.y * e.y + mInert.z * e.z * e.z);

		mol[n].qv.u1 = mol[n].qv.u1 * f;
		mol[n].qv.u2 = mol[n].qv.u2 * f;
		mol[n].qv.u3 = mol[n].qv.u3 * f;
		mol[n].qv.u4 = mol[n].qv.u4 * f;
	} 
}

/* Angular accelerations initialization  (to zero) */
void InitAngAccels()
{
	int n; 
	for (n = 0; n < nMol; n++)
	{
		mol[n].qa.u1 = 0;
		mol[n].qa.u2 = 0;
		mol[n].qa.u3 = 0;
		mol[n].qa.u4 = 0;
		mol[n].qa1.u1 = 0;
		mol[n].qa1.u2 = 0;
		mol[n].qa1.u3 = 0;
		mol[n].qa1.u4 = 0;
		mol[n].qa2.u1 = 0;
		mol[n].qa2.u2 = 0;
		mol[n].qa2.u3 = 0;
		mol[n].qa2.u4 = 0;
	}
}

/* Random number generator with uniform distribution */
/* This function generates random unit vectors */
void VRand (VectorR *p)
{
	real s, x, y;
	s = 2.; 
	while (s > 1.)
	{
		x = 2. * RandR() - 1.;
		y = 2. * RandR() - 1.;
		s = x*x + y*y;
	} 
	p->z = 1. - 2. * s;
	s = 2. * sqrt (1. - s);
	p->x = s * x;
	p->y = s * y;
}

/* Sites forces force and energy calculation */
void Calc_SiteForces()
{
	VectorR dr, shift;
	real fcVal, rr, rrCut, rri, rri3, uVal;
	int j1, j2, m1, m2, ms1, ms2, n, type;

	rrCut = rCut * rCut;
	
	for (n = 0; n < nMol * SITESMOL; n++)
	{
		site[n].f.x = 0;
		site[n].f.y = 0;
		site[n].f.z = 0;
	}
	
	uSum = 0;

	for (m1 = 0; m1 < nMol - 1; m1++)		
	{
		for (m2 = m1 + 1; m2 < nMol; m2++)
		{
			dr.x = mol[m1].r.x - mol[m2].r.x;
			dr.y = mol[m1].r.y - mol[m2].r.y;
			dr.z = mol[m1].r.z - mol[m2].r.z;

			shift.x = 0;
			shift.y = 0;
			shift.z = 0;

			if (dr.x >= 0.5 * LV) shift.x -= LV; 		//shift
			else if (dr.x < -0.5 * LV) shift.x += LV;

			if (dr.y >= 0.5 * LV) shift.y -= LV; 		//shift
			else if (dr.y < -0.5 * LV) shift.y += LV;

			if (dr.z >= 0.5 * LV) shift.z -= LV; 		//shift
			else if (dr.z < -0.5 * LV) shift.z += LV;

	
			dr.x = dr.x + shift.x;
			dr.y = dr.y + shift.y;
			dr.z = dr.z + shift.z;

			rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
				
			if (rr < rrCut)
			{
				ms1 = m1 * SITESMOL;
				ms2 = m2 * SITESMOL;
			
				for (j1 = 0; j1 < SITESMOL; j1++)
				{
					for (j2 = 0; j2 < SITESMOL; j2++)
					{
						type = mSite[j1].typeF + mSite[j2].typeF;

						if (mSite[j1].typeF == mSite[j2].typeF || type == 5)
						{
							dr.x = site[ms1 + j1].r.x - site[ms2 + j2].r.x;
							dr.y = site[ms1 + j1].r.y - site[ms2 + j2].r.y;
							dr.z = site[ms1 + j1].r.z - site[ms2 + j2].r.z;

							dr.x = dr.x + shift.x;
							dr.y = dr.y + shift.y;
							dr.z = dr.z + shift.z;

							rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
							rri = 1/rr;

							switch(type)
							{
								case 2: //OO interaction
									rri3  = rri * rri * rri;
									uVal  =  4 * rri3 * (rri3 - 1);
									fcVal = 48 * rri3 * (rri3 - 0.5) * rri;
									break;
								
								case 4: //MM interaction
									uVal  = 4 * bCon * sqrt(rri);
									fcVal = uVal * rri;
									break;
								
								case 5: //MH interaction
									uVal = -2 * bCon * sqrt(rri);
									fcVal = uVal * rri;
									break; 

								case 6: //HH interaction
									uVal = bCon * sqrt(rri);
									fcVal = uVal * rri;
									break; 
							}

							site[ms1 + j1].f.x = site[ms1 + j1].f.x + fcVal * dr.x;  // + fcVal * r(x, y, z)
							site[ms1 + j1].f.y = site[ms1 + j1].f.y + fcVal * dr.y;
							site[ms1 + j1].f.z = site[ms1 + j1].f.z + fcVal * dr.z;

							site[ms2 + j2].f.x = site[ms2 + j1].f.x - fcVal * dr.x;
							site[ms2 + j2].f.y = site[ms2 + j1].f.y - fcVal * dr.y;
							site[ms2 + j2].f.z = site[ms2 + j1].f.z - fcVal * dr.z;

							uSum += uVal;

						}
					}
				}
			}
		}
	}		

}


/********************************
DefineMol
This functiones defines 
the details of the molecule
*********************************/
void DefineMol()
{
 	mSite[0].r.x = 0;	
 	mSite[1].r.x = 0;	
 	mSite[2].r.x = 0;	
 	mSite[3].r.x = 0;	
	mSite[0].r.z = -0.0206;
	mSite[1].r.z = 0.0274;
	mSite[2].r.y = 0.240;
	mSite[2].r.z = 0.165;
	mSite[3].r.y = - mSite[2].r.y; 
	mSite[3].r.z = mSite[2].r.z;
	mInert.x = 0.00980;
	mInert.y = 0.00340;
	mInert.z = 0.00640;
	bCon = 183.5;
	mSite[0].typeF = 1;   //Oxygen
	mSite[1].typeF = 2;   //M
	mSite[2].typeF = 3;   //hydrogen
	mSite[3].typeF = 3;   //hydrogen
}


/* Functions for geometrically constrained molecules */
/*for example, proteins using the AB off-lattice model */

/* Constraint matrix building */
void buildConstMatrix()
{
	int i;
	
	for (i = 0; i < ProtLen * nC; i++) mMat[i] = 0;	//ProtLen = protein size (number of amino acids); 
													//nC = number of constraints (nC = ProtLen - 1) 	
	
	for (i = 0; i < ProtLen; i++)
	{
		if ( i-1 > 0) mMat[(i-1) * ProtLen + i] = 2;  	//mMat[(i-1) * ProtLen + i] = 2; //else mMat[(i-1) * ProtLen + i] = 0; 

		if (i < nC) mMat[i * ProtLen + i] = -2; 		//mMat[i * ProtLen + i] = -2; //else mMat[(i-1) * ProtLen + i] = 0;
	}

	for (i = 0; i < nC ; i++)
	{
		constraint[i].dist2 = BondLen * BondLen;  		//BondLen = BondLength 
														//BondLen = 1 for the AB off-lattice model  
		constraint[i].ik = i;							//ik and jk atoms (neighboors )
		constraint[i].jk = i + 1;
	}	
}


/* Geometrical constraints calculation */

void calc_constraints()
{
	VectorR dv, da;
	int i, j, k;
	int d;

	for (i = 0; i < nC; i++)
	{
		constraint[i].v.x =  mol[constraint[i].ik].r.x - mol[constraint[i].jk].r.x;
		constraint[i].v.y =  mol[constraint[i].ik].r.y - mol[constraint[i].jk].r.y;
		constraint[i].v.z =  mol[constraint[i].ik].r.z - mol[constraint[i].jk].r.z;
		//
		if (constraint[i].v.x >= 0.5 * LV) constraint[i].v.x -= LV;
		else if (constraint[i].v.x < -0.5 * LV) constraint[i].v.x += LV;

		if (constraint[i].v.y >= 0.5 * LV) constraint[i].v.y -= LV;
		else if (constraint[i].v.y < -0.5 * LV) constraint[i].v.y += LV;

		if (constraint[i].v.z >= 0.5 * LV) constraint[i].v.z -= LV;
		else if (constraint[i].v.z < -0.5 * LV) constraint[i].v.z += LV;
		//
	}

	k = 0;

	for(i = 0; i < nC; i++)
	{
		for(j = 0; j < nC; j++)
		{
			lMat[k] = 0;
			d = mMat[i * ProtLen +  constraint[j].ik]  - mMat[i * ProtLen +  constraint[j].jk];
			if(d != 0)
				lMat[k] = (real) (d * (constraint[i].v.x  * constraint[j].v.x  + constraint[i].v.y  * constraint[j].v.y + constraint[i].v.z  * constraint[j].v.z));

			k++;	
	
		}
	}	
	
	for (i = 0; i < nC; i++)
	{
		dv.x = mol[constraint[i].ik].rv.x - mol[constraint[i].jk].rv.x;
		dv.y = mol[constraint[i].ik].rv.y - mol[constraint[i].jk].rv.y;
		dv.z = mol[constraint[i].ik].rv.z - mol[constraint[i].jk].rv.z;

		da.x = mol[constraint[i].ik].ra.x - mol[constraint[i].jk].ra.x;
		da.y = mol[constraint[i].ik].ra.y - mol[constraint[i].jk].ra.y;
		da.z = mol[constraint[i].ik].ra.z - mol[constraint[i].jk].ra.z;		

		consV[i] = - (da.x * constraint[i].v.x + da.y * constraint[i].v.y + da.z * constraint[i].v.z) - (dv.x * dv.x + dv.y * dv.y + dv.z * dv.z);
	}

	solve_linearEq(lMat, consV, nC);
		
	for (i = 0; i < nC; i++)
	{
		for (j = 0; j < ProtLen; j++)
		{
			if (mMat[i * ProtLen + j] != 0) 
			{
				mol[j].ra.x = mol[j].ra.x + mMat[i * ProtLen + j] * consV[i] * constraint[i].v.x;
				mol[j].ra.y = mol[j].ra.y + mMat[i * ProtLen + j] * consV[i] * constraint[i].v.y;
				mol[j].ra.z = mol[j].ra.z + mMat[i * ProtLen + j] * consV[i] * constraint[i].v.z;
			}
		}
	}	
}



/* Relaxation method for maintaining rigid bonds: 
this function restores the coordinates and
velocities to their constrained values
-- adapted version from Rappaport, 2008
*/
void shake_relaxation()
{
	VectorR dr, dv;	
	real cDev, cDevR, cDevV, g, ga;			
	int changed, m, m1, m2, maxCycle, n;
	
	maxCycle = 200;// --> to define ?
	cDevR = cDevV = 0;

	for(n = 0; n < ProtLen; n++)	
	{
		nCycleR = 0;
		changed = 1;

		while (nCycleR < maxCycle && changed)
		{
			nCycleR++;
			changed = 0;
			cDev = 0;
			for (m=0; m <nC; m++)		
			{
				m1 = constraint[m].ik; //n * ProtLen +
				m2 = constraint[m].jk;
				
				dr.x = mol[m1].r.x - mol[m2].r.x;
				dr.y = mol[m1].r.y - mol[m2].r.y;
				dr.z = mol[m1].r.z - mol[m2].r.z;

				if (dr.x >= 0.5 * LV) dr.x -= LV;
				else if (dr.x < -0.5 * LV) dr.x += LV;

				if (dr.y >= 0.5 * LV) dr.y -= LV;
				else if (dr.y < -0.5 * LV) dr.y += LV;

				if (dr.z >= 0.5 * LV) dr.z -= LV;
				else if (dr.z < -0.5 * LV) dr.z += LV;				

				g 		= ( (dr.x * dr.x + dr.y * dr.y + dr.z * dr.z) - constraint[m].dist2) / (4 * constraint[m].dist2);
				ga 		= fabs(g); 
				cDev 	= (cDev > ga) ? cDev : ga;

				if (ga > consPrec) 
				{
					changed = 1;
					mol[m1].r.x = mol[m1].r.x - g * dr.x;
					mol[m1].r.y = mol[m1].r.y - g * dr.y;
					mol[m1].r.z = mol[m1].r.z - g * dr.z;

					mol[m2].r.x = mol[m2].r.x + g * dr.x;
					mol[m2].r.y = mol[m2].r.y + g * dr.y;
					mol[m2].r.z = mol[m2].r.z + g * dr.z;


				}			
			}
		}
		cDevR = (cDevR > cDev) ? cDevR:cDev;
		nCycleV = 0;
		changed = 1;
		
		while (nCycleV < maxCycle && changed)
		{
			nCycleV++;
			changed = 0;
			cDev = 0;
			for (m = 0; m < nC; m++)
			{
				m1 = constraint[m].ik; 					//n * ProtLen
				m2 = constraint[m].jk;

				dr.x = mol[m1].r.x - mol[m2].r.x;
				dr.y = mol[m1].r.y - mol[m2].r.y;
				dr.z = mol[m1].r.z - mol[m2].r.z;

				if (dr.x >= 0.5 * LV) dr.x -= LV;
				else if (dr.x < -0.5 * LV) dr.x += LV;

				if (dr.y >= 0.5 * LV) dr.y -= LV;
				else if (dr.y < -0.5 * LV) dr.y += LV;

				if (dr.z >= 0.5 * LV) dr.z -= LV;
				else if (dr.z < -0.5 * LV) dr.z += LV;	

				dv.x = mol[m1].rv.x - mol[m2].rv.x;
				dv.y = mol[m1].rv.y - mol[m2].rv.y;
				dv.z = mol[m1].rv.z - mol[m2].rv.z;

				g = (dv.x * dr.x + dv.y * dr.y + dv.z * dr.z) / (2 * constraint[m].dist2);
				ga = fabs (g);
				cDev = (cDev > ga) ? cDev:ga;

				if (ga > consPrec)
				{
					changed = 1;
					mol[m1].rv.x = mol[m1].rv.x - g * dr.x;
					mol[m1].rv.y = mol[m1].rv.y - g * dr.y;
					mol[m1].rv.z = mol[m1].rv.z - g * dr.z;

					mol[m2].rv.x = mol[m2].rv.x + g * dr.x;
					mol[m2].rv.y = mol[m2].rv.y + g * dr.y;
					mol[m2].rv.z = mol[m2].rv.z + g * dr.z;
				}
			}
		}
		cDevV = (cDevV > cDev) ? cDevV:cDev;

	}
}

/* Linear equations Solver  using the standard LU decomposition method */
/* LU decomposition is a matrix decomposition which writes a matrix as 
the product of a lower triangular matrix and an upper triangular matrix */
/* based on the Crout version of the LU method */

void solve_linearEq(real *a, real *x, int n)
{	
	real vMax[SLE_NMAX], v, max;
	int ptr[SLE_NMAX], i, j, k, m;

	for (i = 0; i < n; i++)
	{
		max = 0;

		for(j = 0; j < n; j++)
		{
			v = fabs (a[i + n * j]);
			if (v > max) max = v;
		}

		vMax[i] = 1 / max;
	}

	for (m = 0; m < n; m++)
	{
		max = 0;

		for (i = m; i < n; i++)
		{
			for (k = 0; k < m; k++) 
			{
				a[i + n * m] = a[i + n * m] - (a[i + n * k] * a[k + n * m]);
			}	

			v = fabs(a[i + n * m]) * vMax[i];	

			if (v > max)
			{
				max = v;
				ptr[m] = i;
			}	
		}

		if(m != ptr[m])
		{
			for (k = 0; k < n; k++)
			{
				v = a[i + n * m];
				a[i + n * m] = a[ptr[m] + n * k];
				a[ptr[m] + n * k] = v;
			}
			
			vMax[ptr[m]] = vMax[m];
		}

		for (j = m + 1; j < n; j++)
		{
			for (k = 0; k < m; k++)	
			{
				a[m + n * j]  = a[m + n * j] - (a[m + n * k] * a[k + n * j]);
			}
/*
printf("m=%d, n=%d, j=%d\n", m, n, j);
printf("a[m + n * j] = a[%d] = %lf\n", m + n * j, a[m + n * j]);
printf("a[m + n * m] = a[%d] = %lf\n", m + n * m, a[m + n * m]);
printf("---\n");
*/

			a[m + n * j] = a[m + n * j] / a[m + n * m];
		}
	}	
	for (i = 0; i < n; i++)
	{
		v = x[ptr[i]];
		x[ptr[i]] = x[i];
		x[i] = v;
	
		for (j = 0; j < i; j++)			
		{
			x[i] = x[i] - (a[i + n * j] * x[j]);
		}

		x[i] = x[i] / a[i + n * i];

	}
			
	for(i = n - 2; i >= 0; i--)
	{
		for(j = i + 1; j < n; j++)
		{
			x[i] = x[i] - (a[i + n * j] * x[j]);
		}
	}
		
}



/* Bond Lenght (between 2 molecules - amino acids) average  */
real BondAverage()
{
	real sum;
	int i;

	sum = 0;
	for(i = 0; i < ProtLen-1; i++)
	{
		sum = sum + sqrt((mol[i+1].r.x - mol[i].r.x) * (mol[i+1].r.x - mol[i].r.x) + (mol[i+1].r.y - mol[i].r.y) * (mol[i+1].r.y - mol[i].r.y) + (mol[i+1].r.z - mol[i].r.z) * (mol[i+1].r.z - mol[i].r.z));
	}
	return sum/(ProtLen-1);
}


/* Gyration Radius */
//Hydrophobic Gyration Radius:
real RgH()
{
	int i, n;

	real RgH;
	real x_avg, y_avg, z_avg;

	
	x_avg = 0;
	y_avg = 0;
	z_avg = 0;
	RgH = 0;
	
	n = 0;

	for (i = 0; i < ProtLen; i++)
	{
		if (sequence[i] == 'A')
		{
			x_avg = x_avg + mol[i].r.x;
			y_avg = y_avg + mol[i].r.y;
			z_avg = z_avg + mol[i].r.z;

			n++;
		}
	}
		x_avg = x_avg / n;
		y_avg = y_avg / n;
		z_avg = z_avg / n;

	for (i = 0; i < ProtLen; i++)
	{
		if (sequence[i] == 'A')
		{
			RgH =  RgH + (mol[i].r.x - x_avg) * (mol[i].r.x - x_avg) + (mol[i].r.y - y_avg) * (mol[i].r.y - y_avg) + (mol[i].r.z - z_avg) * (mol[i].r.z - z_avg);
		}
	}
	
	RgH = sqrt(RgH/n);
	
	return RgH;

}

//Hydrophilic gyration radius
real RgP()
{
	int i, n;

	real RgP;
	real x_avg, y_avg, z_avg;
	
	x_avg = 0;
	y_avg = 0;
	z_avg = 0;
	RgP = 0;
	
	n = 0;

	for (i = 0; i < ProtLen; i++)
	{
		if (sequence[i] == 'B')
		{
			x_avg = x_avg + mol[i].r.x;
			y_avg = y_avg + mol[i].r.y;
			z_avg = z_avg + mol[i].r.z;

			n++;
		}
	}
		x_avg = x_avg / n;
		y_avg = y_avg / n;
		z_avg = z_avg / n;

	for (i = 0; i < ProtLen; i++)
	{
		if (sequence[i] == 'B')
		{
			RgP = RgP + (mol[i].r.x - x_avg) * (mol[i].r.x - x_avg) + (mol[i].r.y - y_avg) * (mol[i].r.y - y_avg) + (mol[i].r.z - z_avg) * (mol[i].r.z - z_avg);
		}
	}
	
	RgP = sqrt(RgP/n);
	
	return RgP;

}

//all amino acids gyration radius
real RgAll()
{
	int i;

	real Rg;
	real x_avg, y_avg, z_avg;
	
	x_avg = 0;
	y_avg = 0;
	z_avg = 0;
	Rg = 0;
	
	for (i = 0; i < ProtLen; i++)
	{
		x_avg = x_avg + mol[i].r.x;
		y_avg = y_avg + mol[i].r.y;
		z_avg = z_avg + mol[i].r.z;	
	}
	x_avg = x_avg / ProtLen;
	y_avg = y_avg / ProtLen;
	z_avg = z_avg / ProtLen;

	for (i = 0; i < ProtLen; i++)
	{
		Rg = Rg + (mol[i].r.x - x_avg) * (mol[i].r.x - x_avg) + (mol[i].r.y - y_avg) * (mol[i].r.y - y_avg) + (mol[i].r.z - z_avg) * (mol[i].r.z - z_avg);
	}
	
	Rg = sqrt(Rg / ProtLen);
	
	return Rg;
}






/*Generates double random number between 0 and max */
double randdouble(double max) 
{
	return ((rand()/((double)(RAND_MAX)+1)) * max);
}


void evaluate()							
{
	int n;

	calc_kinetic_energy();
	potEnergy = uTorsion + uChainAngles + uLJ;
	TotalEnergy = kineticEnergy + potEnergy;
	calc_currentTemp();
	calc_pressure(); 
	calc_density();

	bondavg = BondAverage();

	//salva o melhor
	if (bondavg <= 1) 
	{
		if (stepCount == 0)
		{
			best.potEnergy = potEnergy;
			best.step = stepCount;
		}
		else
		{
			if (potEnergy < best.potEnergy)
			{
				best.potEnergy = potEnergy;
				best.step = stepCount;

				for (n= 0; n < nMol; n++)
				{
		       	 	beststructure[n].r.x = mol[n].r.x; 
		        	beststructure[n].r.y = mol[n].r.y; 
		        	beststructure[n].r.z = mol[n].r.z; 
				}
			}
		}


		rgH = RgH();
		rgP = RgP();
		rgAll = RgAll();
	}


}



void calc_kinetic_energy() 
{
	int n;

	kineticEnergy = 0;

	//Translational
	for (n = 0; n < nMol; n++)
	{	
		kineticEnergy += mol[n].rv.x * mol[n].rv.x + mol[n].rv.y * mol[n].rv.y + mol[n].rv.z * mol[n].rv.z;
	}
	
	kineticEnergy *=  0.5; 
}


/* Current temperature calculation using Boltzmann equipartition*/
void calc_currentTemp() 
{
	int n;
    real sum = 0;
    for (n = 0; n < nMol; n++)
	{
        sum += mol[n].rv.x * mol[n].rv.x + mol[n].rv.y * mol[n].rv.y + mol[n].rv.z * mol[n].rv.z;
	}

    current_temperature = sum / (3 * (nMol - 1));
}


void calc_pressure() 
{
	pressure = kineticEnergy * 2 + virial;
	pressure = pressure / (3 * Cube(LV)); //	pressure = pressure / (3 * Cube(LV));
}


void calc_density()
{
	density = nMol  / Cube(LV);  
}



/*TODO: adapt createsRASMOL function for AB model */
void saveBestcoordinates()
{
	int i;
	FILE *coordFile;
	char *nome_arquivo;						
   	nome_arquivo = (char *)malloc(500*sizeof(char));
   	sprintf(nome_arquivo,"coordFile_%daminoacids_%d.txt", ProtLen, ID);
	coordFile = fopen(nome_arquivo, "w");					
	free(nome_arquivo);  										
	
	
	
	//fprintf(coordFile, "Coordinates:\nAminoacid\tx\ty\tz\n");
	for (i=0;i<ProtLen;i++)          
	{
		fprintf(coordFile, "%d\t%lf\t%lf\t%lf\n", i, beststructure[i].r.x, beststructure[i].r.y, beststructure[i].r.z);
	}

	fprintf(coordFile, "BestEnergy = %lf\nStep = %d\n", best.potEnergy, best.step);
	fprintf(coordFile, "RgAll = %lf\nrgH:%lf\nrgP = %lf\n", rgAll, rgH, rgP);
	fprintf(coordFile, "tp = %lfs\n", (double) tp);	//processing time
	fclose(coordFile);
}

void savePathwaysVetor(int pos)
{
	int i;
	

	for (i=0;i<ProtLen;i++) 
	{ 
		(conformations[pos]).x[i] = mol[i].r.x;
		(conformations[pos]).y[i] = mol[i].r.y;
		(conformations[pos]).z[i] = mol[i].r.z;
	}
	(conformations[pos]).potEnergy = potEnergy;
	(conformations[pos]).stepCount = stepCount;
	(conformations[pos]).rgAll = rgAll;
	(conformations[pos]).rgH = rgH;
	(conformations[pos]).rgP = rgP;
}

void savePathwaysVetorArquivo() 
{
	int i, k;
	FILE *coordFile;

   	char *file_name;					
	file_name = (char *)malloc(500*sizeof(char));
   	sprintf(file_name,"pathways%d_%d.txt", nMol, ID);
	
	coordFile = fopen(file_name, "a");
   	 free(file_name);
	
	for(k = 0 ; k < TAMPATHWAY ; k++)
	{	
		fprintf(coordFile, "N\tx\ty\tz\n");
		for (i=0;i<ProtLen;i++) 
		{
			fprintf(coordFile, "%d\t%lf\t%lf\t%lf\n", i, (conformations[k]).x[i], (conformations[k]).y[i], (conformations[k]).z[i]);
		}
		fprintf(coordFile, "\n\n");
		fprintf(coordFile, "Potential Energy = %lf\nStep = %d\n", (conformations[k]).potEnergy, (conformations[k]).stepCount);
		fprintf(coordFile, "rGAll = %lf\nrGH = %lf\nrGP = %lf\n", (conformations[k]).rgAll, (conformations[k]).rgH, (conformations[k]).rgP);
		
		fprintf(coordFile, "\n\n");
	}
	fclose(coordFile);
}
