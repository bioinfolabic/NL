#include "functions.h"


int main(int argc, char **argv) {

	initTimer();
	int i_seed = (int)*argv[3] - 48;
	printf("i_seed = %d\n", i_seed);
	srand(time(NULL));
	seedMT();

	device = (int)*argv[4] - 48;
	cudaSetDevice(device);
	printf("GPU: %d\n", device);

    if(loadFile(argv) == SUCCESS) {
        putParameters();

   		/***   Initiallize atoms' status and LJ potential tables   ***/
    	/***   Set parameters controlling the simulation   ***/
    	initMD();

		/***   Iterations   ***/
		for(i_step = 0; i_step < step_limit; i_step++) {
			step();
			if(i_step % shake_step2shake == 0)
			{

				/***   Shake   ***/
				shakeRelaxation();
						    
			}
			/***   Thermostat   ***/
			berendsenThermostat();

			evaluate();

			if((save_pathways == 'y' || save_pathways == 'Y') && (i_step % pathways_step == 0) && (i_step != 0) || (i_step == 1) )
				savePathways(argv);

		}

		finishSim(argv);

    } else {
        printf("Error: could not load file.");
    }
	return 0;
}
