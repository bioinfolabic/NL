# Copyright (c) 2019, Universidade Tecnológica Federal do Paraná - Bioinformatics and Computacional Intelligence Laboratory
#All rights reserved.
#
#The 3-clause BSD license is applied to this software.
ab:
	cd SRC_GPU_NL && python AB_sequence.py 2gb1 

cpu:
	cd SRC_CPU && make -f Makefile_CPU 

gpu:
	cd SRC_GPU && make -f Makefile_PATHMOLD-AB

gpu_nl:
	cd SRC_GPU_NL && make -f Makefile_PATHMOLD-AB_NL

run_cpu:
	cd SRC_CPU && ./a.out ../INPUT/2GB1_56.in 0 && mv pathways56_0.txt ../OUTPUT/pathways56.txt

run_gpu:
	cd SRC_GPU && ./a.out ../INPUT/2GB1_56.in ../OUTPUT/pathways56 0 0

run_gpu_nl:
	cd SRC_GPU_NL && ./a.out ../INPUT/2GB1_56.in ../OUTPUT/pathways56 0 0

visualize:
	cd SRC_GPU_NL && python pathway_print_multi-subplot.py && mv folding.mp4 ../OUTPUT/folding.mp4 && mv img ../OUTPUT/


