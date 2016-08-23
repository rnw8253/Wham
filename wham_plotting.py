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

meta_file = sys.argv[1]
wham_file = sys.argv[2]
system = sys.argv[3]

nBins = 100
k = 0.001987 # Kcal K^-1 mol^-1
T = 300. # K
kT = k*T
boltz = 2*kT
four_pi = 4*np.pi

# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------

# READ META DATA FILE TO SUBSEQUENTLY ANALYZE THE CORRECT DATASET
file_list = []
data_list = []
with open(meta_file,'r') as f:
	for line in f:
		temp = line.split()
		if temp[0] == '#':
			continue
		else:
			file_list.append(temp[0])
			data_list.append([float(temp[1]),float(temp[2])])

data_list = np.array(data_list)
nProds = len(file_list)

# ----------------------------------------
# Read in WHAM results and plot PMF;
bin_centers = []
free_energy = []
fe_err = []
C_values = []
with open(wham_file,'r') as f:
	for line in f:
		temp = line.split()
		if temp[0] == '#Coor' or temp[0] == '#Window':
			continue
		elif temp[1] == 'inf' or temp[2] == 'nan':
			continue
		elif temp[0][0] != '#':
			bin_centers.append(float(temp[0]))
			free_energy.append(float(temp[1]))
			fe_err.append(float(temp[2]))
		# GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT IS USED TO ALIGN THE WINDOWS WITH EACH OTHER;
		else:
			C_values.append(float(temp[1]))

bin_centers = np.array(bin_centers)
free_energy = np.array(free_energy)
fe_err = np.array(fe_err)
C_values = np.array(C_values)

if len(C_values) != nProds:
	print 'number of windows (in meta data file) do not match the number of values at the bottom of the wham file...'
	sys.exit()

wham_bins = len(bin_centers)
for i in range(wham_bins):
	free_energy[i] += kT*np.log(four_pi*bin_centers[i]**2.0)	# SHODDY VOLUME CORRECTION, BUT ALL I CAN DO
free_energy -= np.ndarray.min(free_energy)

plt.errorbar(bin_centers[:],free_energy[:],yerr=fe_err[:])
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#plt.ylim((-0.5,14.0))
plt.ylabel(r'Relative Free Energy (kCal mol$^{-1}$)',size=14)
plt.xlabel(r'Distance ($\AA$)',size=14)
plt.savefig('%s.Free_energy.png' %(system),dpi=300)
plt.close()

# ----------------------------------------
# LOOP THROUGH ALL DATA FILES, COLLECT DATA, HISTOGRAM DATA INTO FREQ, PROB DENSITY, AND FREE ENERGY COUNTERS
for i in range(nProds):
	with open('%s' %(file_list[i]),'r') as f:	# loading file into a numpy array
		temp = np.loadtxt(f,dtype=np.float)
	
	# collecting data to be used for creating the histograms
	x_min = np.ndarray.min(temp[:,1])
	x_max = np.ndarray.max(temp[:,1])
	delta_x = (x_max - x_min)/nBins
	nValues = len(temp)
	prob_density_divisor = nValues*delta_x
	fe_divisor = prob_density_divisor*four_pi

	half_bins = zeros(nBins)
	for j in range(nBins):
		half_bins[j] = x_min + delta_x*(j+0.5)

	counts = zeros(nBins)		# binning data with no weighting to use later as a counter
	prob_density = zeros(nBins)	# binning data with prob density weighting to observe shape of distribution and overlap between windows
	fe_counts = zeros(nBins)	# binning data with boltzmann weighting to subsequently calc free energy within a window
	for j in range(nValues):
		exponent = data_list[i][1]*(temp[j][1] - data_list[i][0])**2/boltz
		index = int((temp[j][1]-x_min)/delta_x)
		if index == nBins:
			counts[-1] += 1
			prob_density[-1] += 1/prob_density_divisor
			fe_counts[-1] += 1/(fe_divisor*temp[j][1]**2*np.exp(-exponent))
		else:
			counts[index] += 1
			prob_density[index] += 1/prob_density_divisor
			fe_counts[index] += 1/(fe_divisor*temp[j][1]**2*np.exp(-exponent))

	#c = i/float(nProds)
	# HISTOGRAM PROB DENSITY OF ALL DATA FILES ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION OF OVERLAP BETWEEN WINDOWS AND VARIATIONS IN PRODUCTION RUN DISTRIBUTIONS WITHIN EACH WINDOW
	plt.figure(1)
#	plt.plot(half_bins[:],prob_density[:],color=(c,0,0,1))
	plt.plot(half_bins[:],prob_density[:])

	for j in range(nBins):
		fe_counts[j] = -kT*np.log(fe_counts[j])		# taking negative log of the boltzmann weighted fe counter;
		fe_counts[j] += C_values[i] 	# subtracting out the constant used in WHAM to align the windows with each other; this will align the unstitched free energy surfaces, allowing for comparison of overlap of windows

	# PLOT THE FREE ENERGY SURFACE FOR EACH WINDOW ONTO THE SAME PLOT; ALLOWS FOR VISUALIZATION/COMPARISON OF FREE ENERGY SURFACES BEFORE BEING STITCHED TOGETHER BY WHAM; IF GAPS OR LARGE VALUE DISPARITIES BETWEEN WINDOWS ARE PRESENT, THIS INDICATES THAT WINDOWS ARE SAMPLING DIFFERENT DISTRIBUTIONS (AKA NOT ERGODIC)
	plt.figure(2)
#	plt.plot(half_bins[counts > 10],fe_counts[counts > 10],color=(c,0,0,1))
	plt.plot(half_bins[counts > 10],fe_counts[counts > 10])

plt.figure(1)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.ylabel('Probability Density')
plt.xlabel(r'Distance ($\AA$)',size=14)
plt.savefig('%s.data_histogram.png' %(system),dpi=300)
plt.close()

plt.figure(2)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.ylabel(r'Relative Free Energy (kCal mol$^{-1}$)',size=14)
plt.xlabel(r'Distance ($\AA$)',size=14)
plt.savefig('%s.unstitched_fe.png' %(system),dpi=300)
plt.close()

