#include "func_MD_CPU.h"


int main (int argc, char **argv)
{
	sscanf(argv[2],"%d",&ID);
   
	int finish, k = 0;
	srand(time(0));

	ti = time(NULL);				


	initial_message();
	GetParameters(argv);			
	SetParams();					
	PrintParameters();				
	Initialize_MD();				

	finish = 0;
	while (finish == 0)
	{
		finish = Step();		
		k++;
	}	

	freeArrays();					
	
	return 1;
}

