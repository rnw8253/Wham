#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl

zeros = np.zeros

meta_file1 = sys.argv[1]
wham_file1 = sys.argv[2]
system1 = sys.argv[3]

meta_file2 = sys.argv[4]
wham_file2 = sys.argv[5]
system2 = sys.argv[6]

meta_file3 = sys.argv[7]
wham_file3 = sys.argv[8]
system3 = sys.argv[9]

nBins = 100
k = 0.001987 # Kcal K^-1 mol^-1
T = 298. # K
kT = k*T
boltz = 2*kT
four_pi = 4*np.pi

# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------

# READ META DATA FILE TO SUBSEQUENTLY ANALYZE THE CORRECT DATASET
file_list1 = []
data_list1 = []
with open(meta_file1,'r') as f:
	for line in f:
		temp1 = line.split()
		if temp1[0] == '#':
			continue
		else:
			file_list1.append(temp1[0])
			data_list1.append([float(temp1[1]),float(temp1[2])])

data_list1 = np.array(data_list1)
nProds1 = len(file_list1)

file_list2 = []
data_list2 = []
with open(meta_file2,'r') as f:
	for line in f:
		temp2 = line.split()
		if temp2[0] == '#':
			continue
		else:
			file_list2.append(temp2[0])
			data_list2.append([float(temp2[1]),float(temp2[2])])

data_list2 = np.array(data_list2)
nProds2 = len(file_list2)

file_list3 = []
data_list3 = []
with open(meta_file3,'r') as f:
	for line in f:
		temp3 = line.split()
		if temp3[0] == '#':
			continue
		else:
			file_list3.append(temp3[0])
			data_list3.append([float(temp3[1]),float(temp3[2])])

data_list3 = np.array(data_list3)
nProds3 = len(file_list3)

# ----------------------------------------
# Read in WHAM results and plot PMF;
bin_centers1 = []
free_energy1 = []
fe_err1 = []
C_values1 = []
with open(wham_file1,'r') as f:
	for line in f:
		temp1 = line.split()
		if temp1[0] == '#Coor' or temp1[0] == '#Window':
			continue
		elif temp1[1] == 'inf' or temp1[2] == 'nan':
			continue
		elif temp1[0][0] != '#':
			bin_centers1.append(float(temp1[0]))
			free_energy1.append(float(temp1[1]))
			fe_err1.append(float(temp1[2]))
		# GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT IS USED TO ALIGN THE WINDOWS WITH EACH OTHER;
		else:
			C_values1.append(float(temp1[1]))

bin_centers1 = np.array(bin_centers1)
free_energy1 = np.array(free_energy1)
fe_err1 = np.array(fe_err1)
C_values1 = np.array(C_values1)

if len(C_values1) != nProds1:
	print 'number of windows (in meta data file) do not match the number of values at the bottom of the wham file...'
	sys.exit()

wham_bins1 = len(bin_centers1)
for i in range(wham_bins1):
	free_energy1[i] += kT*np.log(four_pi*bin_centers1[i]**2.0)	# SHODDY VOLUME CORRECTION, BUT ALL I CAN DO
free_energy1 -= np.ndarray.min(free_energy1)

#######

bin_centers2 = []
free_energy2 = []
fe_err2 = []
C_values2 = []
with open(wham_file2,'r') as f:
	for line in f:
		temp2 = line.split()
		if temp2[0] == '#Coor' or temp2[0] == '#Window':
			continue
		elif temp2[1] == 'inf' or temp2[2] == 'nan':
			continue
		elif temp2[0][0] != '#':
			bin_centers2.append(float(temp2[0]))
			free_energy2.append(float(temp2[1]))
			fe_err2.append(float(temp2[2]))
		# GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT IS USED TO ALIGN THE WINDOWS WITH EACH OTHER;
		else:
			C_values2.append(float(temp2[1]))

bin_centers2 = np.array(bin_centers2)
free_energy2 = np.array(free_energy2)
fe_err2 = np.array(fe_err2)
C_values2 = np.array(C_values2)

if len(C_values2) != nProds2:
	print 'number of windows (in meta data file) do not match the number of values at the bottom of the wham file...'
	sys.exit()

wham_bins2 = len(bin_centers2)
for i in range(wham_bins2):
	free_energy2[i] += kT*np.log(four_pi*bin_centers2[i]**2.0)	# SHODDY VOLUME CORRECTION, BUT ALL I CAN DO
free_energy2 -= np.ndarray.min(free_energy2)

#######

bin_centers3 = []
free_energy3 = []
fe_err3 = []
C_values3 = []
with open(wham_file3,'r') as f:
	for line in f:
		temp3 = line.split()
		if temp3[0] == '#Coor' or temp3[0] == '#Window':
			continue
		elif temp3[1] == 'inf' or temp3[2] == 'nan':
			continue
		elif temp3[0][0] != '#':
			bin_centers3.append(float(temp3[0]))
			free_energy3.append(float(temp3[1]))
			fe_err3.append(float(temp3[2]))
		# GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT IS USED TO ALIGN THE WINDOWS WITH EACH OTHER;
		else:
			C_values3.append(float(temp3[1]))

bin_centers3 = np.array(bin_centers3)
free_energy3 = np.array(free_energy3)
fe_err3 = np.array(fe_err3)
C_values3 = np.array(C_values3)

if len(C_values3) != nProds3:
	print 'number of windows (in meta data file) do not match the number of values at the bottom of the wham file...'
	sys.exit()

wham_bins3 = len(bin_centers3)
for i in range(wham_bins3):
	free_energy3[i] += kT*np.log(four_pi*bin_centers3[i]**2.0)	# SHODDY VOLUME CORRECTION, BUT ALL I CAN DO
free_energy3 -= np.ndarray.min(free_energy3)

plt.errorbar(bin_centers1[:],free_energy1[:],yerr=fe_err1[:], label='%s' %(system1))
plt.errorbar(bin_centers2[:],free_energy2[:],yerr=fe_err2[:], label='%s' %(system2))
plt.errorbar(bin_centers3[:],free_energy3[:],yerr=fe_err3[:], label='%s' %(system3))
plt.legend(loc=4, borderaxespad=0.)
#plt.legend(handles=[blue_line])
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylim((-0.5,14.0))
plt.ylabel(r'Relative Free Energy (kCal mol$^{-1}$)',size=14)
plt.xlabel(r'Distance ($\AA$)',size=14)
plt.savefig('comparison.Free_energy.png',dpi=300)
plt.close()


## ----------------------------------------
## LOOP THROUGH ALL DATA FILES, COLLECT DATA, HISTOGRAM DATA INTO FREQ, PROB DENSITY, AND FREE ENERGY COUNTERS
#for i in range(nProds1):
#	with open('%s' %(file_list1[i]),'r') as f:	# loading file into a numpy array
#		temp1 = np.loadtxt(f,dtype=np.float)
#	
#	# collecting data to be used for creating the histograms
#	x_min1 = np.ndarray.min(temp1[:,1])
#	x_max1 = np.ndarray.max(temp1[:,1])
#	delta_x1 = (x_max1 - x_min1)/nBins1
#	nValues1 = len(temp1)
#	prob_density_divisor1 = nValues1*delta_x1
#	fe_divisor1 = prob_density_divisor1*four_pi
#
#	half_bins1 = zeros(nBins)
#	for j in range(nBins):
#		half_bins1[j] = x_min1 + delta_x1*(j+0.5)
#
#	counts1 = zeros(nBins)		# binning data with no weighting to use later as a counter
#	prob_density1 = zeros(nBins)	# binning data with prob density weighting to observe shape of distribution and overlap between windows
#	fe_counts1 = zeros(nBins)	# binning data with boltzmann weighting to subsequently calc free energy within a window
#	for j in range(nValues1):
#		exponent1 = data_list1[i][1]*(temp1[j][1] - data_list1[i][0])**2/boltz
#		index1 = int((temp1[j][1]-x_min1)/delta_x1)
#		if index1 == nBins:
#			counts1[-1] += 1
#			prob_density1[-1] += 1/prob_density_divisor1
#			fe_counts1[-1] += 1/(fe_divisor1*temp1[j][1]**2*np.exp(-exponent1))
#		else:
#			counts1[index1] += 1
#			prob_density1[index1] += 1/prob_density_divisor1
#			fe_counts[index1] += 1/(fe_divisor1*temp1[j][1]**2*np.exp(-exponent1))
#
#	#c = i/float(nProds)
#	# HISTOGRAM PROB DENSITY OF ALL DATA FILES ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION OF OVERLAP BETWEEN WINDOWS AND VARIATIONS IN PRODUCTION RUN DISTRIBUTIONS WITHIN EACH WINDOW
#	plt.figure(1)
##	plt.plot(half_bins[:],prob_density[:],color=(c,0,0,1))
#	plt.plot(half_bins1[:],prob_density1[:])
#
#	for j in range(nBins):
#		fe_counts1[j] = -kT*np.log(fe_counts1[j])		# taking negative log of the boltzmann weighted fe counter;
#		fe_counts1[j] += C_values1[i] 	# subtracting out the constant used in WHAM to align the windows with each other; this will align the unstitched free energy surfaces, allowing for comparison of overlap of windows
#
#	# PLOT THE FREE ENERGY SURFACE FOR EACH WINDOW ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION/COMPARISON OF FREE ENERGY SURFACES BEFORE BEING STITCHED TOGETHER BY WHAM; IF GAPS OR LARGE VALUE DISPARITIES BETWEEN WINDOWS ARE PRESENT, THIS INDICATES THAT WINDOWS ARE SAMPLING DIFFERENT DISTRIBUTIONS (AKA NOT ERGODIC)
#	plt.figure(2)
##	plt.plot(half_bins[counts > 10],fe_counts[counts > 10],color=(c,0,0,1))
#	plt.plot(half_bins1[counts1 > 10],fe_counts1[counts1 > 10])
#
#plt.figure(1)
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylabel('Probability Density')
#plt.xlabel(r'Distance ($\AA$)',size=14)
#plt.savefig('%s.data_histogram.png' %(system1),dpi=300)
#plt.close()
#
#plt.figure(2)
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylabel(r'Relative Free Energy (kCal mol$^{-1}$)',size=14)
#plt.xlabel(r'Distance ($\AA$)',size=14)
#plt.savefig('%s.unstitched_fe.png' %(system1),dpi=300)
#plt.close()
#
#
###########
#
## LOOP THROUGH ALL DATA FILES, COLLECT DATA, HISTOGRAM DATA INTO FREQ, PROB DENSITY, AND FREE ENERGY COUNTERS
#for i in range(nProds2):
#	with open('%s' %(file_list2[i]),'r') as f:	# loading file into a numpy array
#		temp2 = np.loadtxt(f,dtype=np.float)
#	
#	# collecting data to be used for creating the histograms
#	x_min2 = np.ndarray.min(temp2[:,1])
#	x_max2 = np.ndarray.max(temp2[:,1])
#	delta_x2 = (x_max2 - x_min2)/nBins
#	nValues2 = len(temp2)
#	prob_density_divisor2 = nValues2*delta_x2
#	fe_divisor2 = prob_density_divisor2*four_pi
#
#	half_bins2 = zeros(nBins)
#	for j in range(nBins):
#		half_bins2[j] = x_min2 + delta_x2*(j+0.5)
#
#	counts2 = zeros(nBins)		# binning data with no weighting to use later as a counter
#	prob_density2 = zeros(nBins)	# binning data with prob density weighting to observe shape of distribution and overlap between windows
#	fe_counts2 = zeros(nBins)	# binning data with boltzmann weighting to subsequently calc free energy within a window
#	for j in range(nValues):
#		exponent = data_list[i][1]*(temp[j][1] - data_list[i][0])**2/boltz
#		index = int((temp[j][1]-x_min)/delta_x)
#		if index == nBins:
#			counts[-1] += 1
#			prob_density[-1] += 1/prob_density_divisor
#			fe_counts[-1] += 1/(fe_divisor*temp[j][1]**2*np.exp(-exponent))
#		else:
#			counts[index] += 1
#			prob_density[index] += 1/prob_density_divisor
#			fe_counts[index] += 1/(fe_divisor*temp[j][1]**2*np.exp(-exponent))
#
#	#c = i/float(nProds)
#	# HISTOGRAM PROB DENSITY OF ALL DATA FILES ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION OF OVERLAP BETWEEN WINDOWS AND VARIATIONS IN PRODUCTION RUN DISTRIBUTIONS WITHIN EACH WINDOW
#	plt.figure(1)
##	plt.plot(half_bins[:],prob_density[:],color=(c,0,0,1))
#	plt.plot(half_bins[:],prob_density[:])
#
#	for j in range(nBins):
#		fe_counts[j] = -kT*np.log(fe_counts[j])		# taking negative log of the boltzmann weighted fe counter;
#		fe_counts[j] += C_values[i] 	# subtracting out the constant used in WHAM to align the windows with each other; this will align the unstitched free energy surfaces, allowing for comparison of overlap of windows
#
#	# PLOT THE FREE ENERGY SURFACE FOR EACH WINDOW ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION/COMPARISON OF FREE ENERGY SURFACES BEFORE BEING STITCHED TOGETHER BY WHAM; IF GAPS OR LARGE VALUE DISPARITIES BETWEEN WINDOWS ARE PRESENT, THIS INDICATES THAT WINDOWS ARE SAMPLING DIFFERENT DISTRIBUTIONS (AKA NOT ERGODIC)
#	plt.figure(2)
##	plt.plot(half_bins[counts > 10],fe_counts[counts > 10],color=(c,0,0,1))
#	plt.plot(half_bins[counts > 10],fe_counts[counts > 10])
#
#plt.figure(1)
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylabel('Probability Density')
#plt.xlabel(r'Distance ($\AA$)',size=14)
#plt.savefig('%s.data_histogram.png' %(system),dpi=300)
#plt.close()
#
#plt.figure(2)
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylabel(r'Relative Free Energy (kCal mol$^{-1}$)',size=14)
#plt.xlabel(r'Distance ($\AA$)',size=14)
#plt.savefig('%s.unstitched_fe.png' %(system),dpi=300)
#plt.close()
#
################
## LOOP THROUGH ALL DATA FILES, COLLECT DATA, HISTOGRAM DATA INTO FREQ, PROB DENSITY, AND FREE ENERGY COUNTERS
#for i in range(nProds):
#	with open('%s' %(file_list[i]),'r') as f:	# loading file into a numpy array
#		temp = np.loadtxt(f,dtype=np.float)
#	
#	# collecting data to be used for creating the histograms
#	x_min = np.ndarray.min(temp[:,1])
#	x_max = np.ndarray.max(temp[:,1])
#	delta_x = (x_max - x_min)/nBins
#	nValues = len(temp)
#	prob_density_divisor = nValues*delta_x
#	fe_divisor = prob_density_divisor*four_pi
#
#	half_bins = zeros(nBins)
#	for j in range(nBins):
#		half_bins[j] = x_min + delta_x*(j+0.5)
#
#	counts = zeros(nBins)		# binning data with no weighting to use later as a counter
#	prob_density = zeros(nBins)	# binning data with prob density weighting to observe shape of distribution and overlap between windows
#	fe_counts = zeros(nBins)	# binning data with boltzmann weighting to subsequently calc free energy within a window
#	for j in range(nValues):
#		exponent = data_list[i][1]*(temp[j][1] - data_list[i][0])**2/boltz
#		index = int((temp[j][1]-x_min)/delta_x)
#		if index == nBins:
#			counts[-1] += 1
#			prob_density[-1] += 1/prob_density_divisor
#			fe_counts[-1] += 1/(fe_divisor*temp[j][1]**2*np.exp(-exponent))
#		else:
#			counts[index] += 1
#			prob_density[index] += 1/prob_density_divisor
#			fe_counts[index] += 1/(fe_divisor*temp[j][1]**2*np.exp(-exponent))
#
#	#c = i/float(nProds)
#	# HISTOGRAM PROB DENSITY OF ALL DATA FILES ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION OF OVERLAP BETWEEN WINDOWS AND VARIATIONS IN PRODUCTION RUN DISTRIBUTIONS WITHIN EACH WINDOW
#	plt.figure(1)
##	plt.plot(half_bins[:],prob_density[:],color=(c,0,0,1))
#	plt.plot(half_bins[:],prob_density[:])
#
#	for j in range(nBins):
#		fe_counts[j] = -kT*np.log(fe_counts[j])		# taking negative log of the boltzmann weighted fe counter;
#		fe_counts[j] += C_values[i] 	# subtracting out the constant used in WHAM to align the windows with each other; this will align the unstitched free energy surfaces, allowing for comparison of overlap of windows
#
#	# PLOT THE FREE ENERGY SURFACE FOR EACH WINDOW ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION/COMPARISON OF FREE ENERGY SURFACES BEFORE BEING STITCHED TOGETHER BY WHAM; IF GAPS OR LARGE VALUE DISPARITIES BETWEEN WINDOWS ARE PRESENT, THIS INDICATES THAT WINDOWS ARE SAMPLING DIFFERENT DISTRIBUTIONS (AKA NOT ERGODIC)
#	plt.figure(2)
##	plt.plot(half_bins[counts > 10],fe_counts[counts > 10],color=(c,0,0,1))
#	plt.plot(half_bins[counts > 10],fe_counts[counts > 10])
#
#plt.figure(1)
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylabel('Probability Density')
#plt.xlabel(r'Distance ($\AA$)',size=14)
#plt.savefig('%s.data_histogram.png' %(system),dpi=300)
#plt.close()
#
#plt.figure(2)
#plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylabel(r'Relative Free Energy (kCal mol$^{-1}$)',size=14)
#plt.xlabel(r'Distance ($\AA$)',size=14)
#plt.savefig('%s.unstitched_fe.png' %(system),dpi=300)
#plt.close()
#
