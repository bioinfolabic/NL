# PathMolD-AB(Neighbourhood List Version)

PathMolD-AB(Neighbourhood List Version) or PathMolD-AB(NL Version) is a software package for generating biomolecular folding pathway datasets using massively parallel MD with neighbourhood list methods in the canonical ensemble. The software describes proteins by a coarse-grained AB model, and it is designed for extensive processing in GPU architectures, therefore being ready for data production in modern High-Performance Computing (HPC) facilities. Also, PathMolD-AB(NL Version) incorporates modules for the analyses of its output folding trajectories.

## Directories Structure

```
PATHMOLD-AB/

||INPUT
||OUTPUT
||-- img
||SRC_CPU
||SRC_GPU
||SRC_GPU_NL
||DOC
```

## Directories Description

INPUT - directory with input files, such as the AB sequence, parameters file for the folding simulation, and FASTA file.

OUTPUT - directory with output files, such as the protein folding data, folding simulation images, and the video of the protein folding simulation.

SRC_CPU - directory with the simulation of Molecular Dynamics software in the sequential/CPU version (main.c, define.h, function.c, function.h, mt.h).

SRC_GPU - directory with the parallel Molecular Dynamics software (main.c, define.h, function.cu, function.h, mt.h).

SRC_GPU_NL - directory with the parallel Molecular Dynamics with Neighbourhood List software (main.c, define.h, function.cu, function.h, mt.h), software that converts amino acid sequence to the AB sequence (AB_sequence .py), and the software for generating the folding simulation video (AB_sequence.py).

DOC - directory with a step-by-step manual for the new users.

*the root directory contains makefile to run the programs. 

	$ make ab               # convert the AB sequence from the amino acid sequence
	
	$ make cpu              # compile the sequential molecular dynamic version

	$ make gpu	 	        # compile the parallel molecular dynamic version

	$ make gpu_nl 	        # compile the parallel molecular dynamic with neighbourhood list version

	$ make run_cpu          # execute the sequential molecular dynamic version

	$ make run_gpu          # execute the parallel molecular dynamic version

	$ make run_gpu_nl       # execute the parallel molecular dynamic with neighbourhood list version

	$ make vizualize        # generating the folding simulation video


## Computing Capability of the PathMolD-AB(NL Version)
The programs were tested in the following configurations.

PathMolD-AB(NL Version) capability (ubuntu 18 LTS)
Computing Capability (CC)

| GPU model| CC         | CUDA7 | CUDA8 | CUDA9 |

| GTX660   | 3          |   NO  |   NO  |   NO  |

| K40      | 3.5        |   NO  |   NO  |   NO  |

| GTX750   | 5          |   NO  |   NO  |   NO  |

| Titan X  | 5.2        |   NO  |  YES  |  YES  |

| GTX 1080 | 6.1        |   NO  |  YES  |  YES  |

| Titan Xp | 6.1        |   NO  |  YES  |  YES  |

| GCC/G                 |  4.8  |  5.3  |  6.5  |

| Python                |        2.7/3.6        |

for any doubt send a email to ​leandrotakeshihattori@gmail.com

Laboratory of Bioinformatics and Computational Intelligence
Federal Technological University of Paraná – UTFPR
http://labic.utfpr.edu.br/
