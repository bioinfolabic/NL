# -*- coding: utf-8 -*-

import csv
import pickle

import re
import imageio
import glob
import os


import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

version = 2
try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO
	version = 3
#######################################################################################################################################################################
################################################################### CHANGE HERE #######################################################################################
#######################################################################################################################################################################
path_pathways = '../OUTPUT/' 	                   # folder in which the datasets are
filename = 'pathway_2GB1_0.txt'                    # name of the dataset from which the images will be created
filesequencia = '../INPUT/2GB1.txt'                # file containing the AB sequence of the protein
path_save = 'img'                            	   # folder where the images will be saved
########################################################################################################################################################################
########################################################################################################################################################################

def create_folding_img_directiory():
	
	# Create target Directory if don't exist
	if not os.path.exists(path_save):
		os.mkdir(path_save)
		print("Directory " , path_save ,  " Created ")
	else:    
		print("Directory " , path_save ,  " already exists")
        	
def open_file(filename):
	f = open(filename, 'r')
	return f 

def close_file(f):
	f.close()

def read_file(f):
	count = 0
	data    = [] 
	features = []
	first = 0
	f.readline()
	
	for lines_read, line in enumerate(f): # lines_read - the number of the line
		                                  # line       - the content in the current lines_read variable
		test  = re.findall(r'[a-zA-Z]+', line) #Check if there are letters  
		test2 = re.findall(r'[0-9]+', line)    #Check if there are numbers   
		# Begin Test
		#print("Lineread: ",lines_read)
		#print("Test letters content: ", test)
		#print("Test letters size ", len(test))
		#print("Test number contant: ", test2) 
		#print("Test number size ", len(test2)) 
        # End Test  		
		if(len(test) == 0 and len(test2) > 0): # copy the line that has the coord information   
			if(version < 3):
    				aux = (np.genfromtxt(StringIO.StringIO(line), delimiter="\t"))
			else:
    				aux = (np.genfromtxt(StringIO(line), delimiter="\t"))
			
			data = np.append(aux,data)	
		if(len(test) > 0 and len(test2) > 0): # copy the line that has the coord information 
			if(version < 3):
				aux = (np.genfromtxt(StringIO.StringIO(line), dtype=[('mystring','S5'),('myfloat','f8')], delimiter=" = "))
			else:
				aux = (np.genfromtxt(StringIO(line), dtype=[('mystring','S5'),('myfloat','f8')], delimiter=" = "))
				#print(aux)
				#print(type(aux))
			features.append(aux)
	return np.array(data), np.array(features)

def get_shape(dataset):
    sequence_size = int(dataset[0]) + 1
    dimension     = 4 # amino acid sequence; x; y; z     
    total_element = dataset.shape[0]
    folding       = total_element / (sequence_size * dimension) 
    return [int(folding), int(sequence_size), int(dimension)]

def format_dataset(dataset,shape):		     # dataset = raw dataset, shape = (numere of sample, sequence, coord x-y-z)
	print(np.array(dataset).shape)
	#shape = int(shape[0]), int(shape[1]), int(shape[2])
	dataset = dataset.reshape((shape))       # make the reshape 
	
	
	dataset = dataset[:,::-1,:]              # invert to first aminio acid to the last
	dataset = dataset[::-1,:,:]              # invert to the first fold to the last fold
	dataset = dataset[:,:,1:4]               # remove the amino acid position 
	print("dataset format shape: ", dataset.shape)
	return dataset

def get_feature(feature):
	#print("feature.shape[0]: ", feature.shape[0])
	number_feature = 5
	if(feature.shape[0]/number_feature != 1001):
		number_feature = 8
	new_shape = int(feature.shape[0]*1.0/number_feature)

	print("get_feature shape: ", new_shape)
	feature_Pot = np.zeros(new_shape)
	feature_Step = np.zeros(new_shape)
	feature_rGAll = np.zeros(new_shape)
	feature_rGH = np.zeros(new_shape)
	feature_rGP = np.zeros(new_shape)

	i_Pot = 0
	i_Step = 0
	i_rGAll = 0
	i_rGH = 0
	i_rGP = 0
	print("feature shape ", feature.shape)

	for i in range(0,feature.shape[0]):
		#print(feature[i][0])
		#print(type(print(feature[i][0])))
		feat = str(feature[i][0].decode("utf-8")) 
		#print(feat)

		#print(feature[i][1])
		if("Poten" == feat):
			feature_Pot[i_Pot] = feature[i][1]
			i_Pot = i_Pot +1	
		elif("Step" == feat):
			feature_Step[i_Step] = feature[i][1]
			i_Step = i_Step +1
		elif("rGAll" == feat):
			feature_rGAll[i_rGAll] = feature[i][1]
			i_rGAll = i_rGAll +1
		elif("rGH" == feat):
			feature_rGH[i_rGH] = feature[i][1]
			i_rGH = i_rGH +1
		elif("rGP" == feat):
			feature_rGP[i_rGP] = feature[i][1]
			i_rGP = i_rGP +1

	return feature_Pot, feature_Step, feature_rGAll, feature_rGH, feature_rGP 


def generate_protein_folding_video(frames):
	os.system('ffmpeg -r 3 -f image2 -s 720x480 -start_number 0 -i ' + path_save + '/%d.png -vframes '+ str(frames) +' -vcodec libx264 -crf 25 -pix_fmt yuv420p folding.mp4')

def main():

	print("reading...")
	print("READING: " + path_pathways + filename)
	print("READING: " + path_pathways + filesequencia)
	
	f = open_file(path_pathways + filename)
	fs = open_file(path_pathways + filesequencia)
	dataset, feat = read_file(f)
	print("Dataset shape: " ,dataset)
	print("feat shape: ", feat)
	close_file(f)
	
	print("pre-processing...")
	shape = get_shape(dataset)
	print("get a shape: ", shape)
	dataset = format_dataset(dataset, shape)
	
	#print "new shape: ", dataset.shape
	
	coords_max_min = np.zeros(shape=(2,3))
	for i in range(0,3):
		coords_max_min[0,i] = dataset[:,:,i].min()
		coords_max_min[1,i] = dataset[:,:,i].max()
	
	ab_seq = []
	sequencia = fs.read()
	sequencia = sequencia.replace('\n','')
	for  i in sequencia:
		ab_seq.append(i)
	ab_seq = np.asarray(ab_seq)
	close_file(fs)
	
	hydrophobic_pos = np.where(ab_seq == 'A')[0]
	polar_pos       = np.where(ab_seq == 'B')[0]

	feature_Pot, feature_Step, feature_rGAll, feature_rGH, feature_rGP  = get_feature(feat)
	
	create_folding_img_directiory()
	print("ploting...")
	#print("DATASET: ", dataset.shape)
	#print("POT: ", feature_Pot.shape)
	'''		
	for i in range(0,dataset.shape[0]):
		print("iteration: ", i)

		fig = plt.figure(figsize=(7, 8))
		gs = gridspec.GridSpec(2, 2)
		
		# First subplot
		# ===========================================================================

		ax1 = fig.add_subplot(gs[0, 0])

		ax1.plot(np.linspace(0, feature_Pot.shape[0], dataset.shape[0]), feature_Pot, color='cornflowerblue')
		ax1.plot(np.array([i,i]), np.array([feature_Pot.min(), feature_Pot.max()]), color='orange', dashes=[6, 2], label='current energy')
		
		ax1.grid(True)
		ax1.set_ylabel('Potencial Energy')
		ax1.set_xlabel('Iteration')
		ax1.legend(loc='upper right')

		# Second subplot
		#===========================================================================
		ax2 = fig.add_subplot(gs[0, 1])

		ax2.plot(np.linspace(0, feature_rGH.shape[0], dataset.shape[0]), feature_rGH, color='salmon', label='RgH')
		ax2.plot(np.linspace(0, feature_rGP.shape[0], dataset.shape[0]), feature_rGP, color='cornflowerblue', label='RgP')

		ax2.plot(np.array([i,i]), np.array([np.asarray([feature_rGH.min(), feature_rGP.min()]).min(), np.asarray([feature_rGH.max(), feature_rGP.max()]).max()]) , color='orange', dashes=[6, 2], label='current Rg')


		ax2.grid(True)
		ax2.set_ylabel('Rg')
		ax2.set_xlabel('Iteration')
		ax2.legend(loc='upper right')

		# Third subplot
		#===========================================================================
		
		
		ax3 = fig.add_subplot(gs[1, :], projection='3d')
		ax3.set_xlim(coords_max_min[0,0],coords_max_min[1,0])
		ax3.set_ylim(coords_max_min[0,1],coords_max_min[1,1])
		ax3.set_zlim(coords_max_min[0,2],coords_max_min[1,2])
		m = 'o'
		hydrophobic_color = 'r'
		polar_color       = 'b'

		coordenates = dataset[i]
		hydrophobic_aa = np.delete(coordenates, polar_pos, 0)
		polar_aa       = np.delete(coordenates, hydrophobic_pos, 0)

		lines  = ax3.plot(coordenates[:,0], coordenates[:,1], coordenates[:,2], label='backbone', color='black')
		points = ax3.scatter(hydrophobic_aa[:,0], hydrophobic_aa[:,1], hydrophobic_aa[:,2], c=hydrophobic_color, marker=m, label='hydrophobic')
		pontts2 = ax3.scatter(polar_aa[:,0], polar_aa[:,1], polar_aa[:,2], c=polar_color, marker=m, label='polar')

		ax3.legend(bbox_to_anchor=(1.1,0.0))
		plot = plt.savefig( path_save + '/' + str(i) + '.png', pad_inches=1.0, dpi=80)
		plt.close('all')
	'''
	generate_protein_folding_video(dataset.shape[0])
	
if __name__ == "__main__":
    main()
