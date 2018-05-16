"""
2/24/18
Small script that determines the index of the beams derived from a calibration grid observation. It reads in aggregated_weights*.mat
file, loops through the XEL/EL offset pairs, and finds the Antenna pointing that is the shortest angular distance from the offsets. 
Users Inputs:
Project - GBO project name (e.g., AGBT17B_360_01)

Usage:
ipython genSensitivityMap.py Project
Example:
ipython genSensitivityMap.py Project 
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""

## imports
from astropy.io import fits
import numpy as np
import scipy.io
import matplotlib.pyplot as pyplot
import matplotlib 
import glob
import os 
import sys
import scipy.io

## get projectID 
projectID  = sys.argv[1]

if projectID == 'AGBT17B_360_01':
	onScans = np.array([19, 20, 21, 22, 23,25, 26, 27, 28, 29, 31,32, 33, 34, 35, 37, 38, 39, 40, 41,43, 44, 45, 46, 47,49, 50, 51, 52, 53, 55, 56, 57, 58])
	wtFile = 'w_AGBT17B_360_01_grid_A.FITS'
	matFileName = 'AGBT17B_360_01_aggregated_weights_X_grid.mat'
wtFile = 'w_' + projectID + '_A_FullGrid.FITS'
matFileName = projectID + '_aggregated_weights_X_grid.mat'

## read in weight FITS file to get beam offsets
#wtHdu = fits.open('/home/scratch/npingel/FLAG/data/AGBT16B_400/' + projectID + '/weight_files/fits_files/' + wtFile)
#bmXelOff = wtHdu[1].data['BeamOff_XEL'][0:7] ## degs
#bmElOff = wtHdu[1].data['BeamOff_EL'][0:7] ## degs
wtHdu = fits.open('/lustre/flag/' + projectID + '/BF/weight_files/' + wtFile)
bmXelOff = wtHdu[1].data['BeamOff_AZ'][0:7] ## degs
bmElOff = wtHdu[1].data['BeamOff_EL'][0:7] ## degs

## read in all beam offsets (only need one polarization)
matFile = scipy.io.loadmat('/lustre/flag/' + projectID + '/BF/mat/' + matFileName)


xelArr = matFile['AZ']
elArr = matFile['EL']


for i in range(0, len(bmXelOff)): ## loop through beam offset pairs from weight files
	bmXelVal = bmXelOff[i] 
	bmElVal = bmElOff[i]
	print('BYU Weight Beam: ' + np.str(i))
	print('Beam XEL Value [degs]: ' + np.str(bmXelVal))
	print('Beam EL Value [degs]: ' + np.str(bmElVal))
	distArr = np.zeros([len(elArr)]) ## create array to store distance values

	for j in range(0, len(xelArr)): ## loop through aggregated offset pairs
		aggXel = xelArr[j]
		aggEl = elArr[j]
		d = np.sqrt((aggXel - bmXelVal)**2 + (aggEl - bmElVal)**2)
		distArr[j] = d
	
	## compute minimum distance and beam index in agg array
	minDist = np.min(distArr)
	minInd = np.where(distArr == minDist)

	## inform user
	print('Minimum distace: ' + np.str(np.min(distArr)))
	print('Aggregated Beam Index: ' + np.str(np.int(minInd[0])))




