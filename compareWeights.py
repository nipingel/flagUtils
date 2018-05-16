"""
4/9/18
Plots the overall phase and amplitude difference of all of the complex weights collected with FLAG
on the GBT utilizing the BYU-WVU backend. No user inputs required 
Usage: 
ipython compareWeights.py
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""

## imports
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.rcsetup
import aplpy as apl
from astropy.io import fits
from numpy import linalg as la
import glob as glob
import sys
import os
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

font = {'family' : 'normal',
        'weight' : 'bold'}

matplotlib.rc('font', **font)
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('text', usetex='false')

# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

## dictionaries that store the scan numbers of a given calibration scan type per session
AGBT16B_400_03_ScanNums = {'grid':np.arange(12,56)}
AGBT16B_400_09_ScanNums = {'grid':np.arange(5,31)} 
AGBT16B_400_12_ScanNums = {'grid':np.arange(32,86), '7Pt-Cal':np.arange(130, 139)}
AGBT16B_400_13_ScanNums = {'grid':np.arange(9, 60), '7Pt-Cal':np.arange(101, 110)}
AGBT16B_400_14_ScanNums = {'grid':np.arange(25, 72), '7Pt-Cal':np.arange(16, 25)}
AGBT17B_360_01_ScanNums = {'7Pt-Cal1':np.arange(11,20), '7Pt-Cal2':np.arange(59, 68), 'grid':np.arange(19, 59)}
AGBT17B_360_02_ScanNums = {'7Pt-Cal1':np.arange(12,21), '7Pt-Cal2':np.arange(79, 88), 'grid':np.arange(24, 64)}
AGBT17B_360_03_ScanNums = {'7Pt-Cal':np.arange(5,14), 'grid':np.arange(14,54)}
AGBT17B_360_04_ScanNums = {'7Pt-Cal':np.arange(71, 80), 'grid':np.arange(31, 71)}
AGBT17B_360_05_ScanNums = {'7Pt-Cal1':np.arange(2,11), '7Pt-Cal2':np.arange(44,53)}
AGBT17B_360_06_ScanNums = {'grid':np.arange(72,112)}
AGBT17B_360_07_ScanNums = {'7Pt-Cal':np.arange(187, 196), 'grid':np.arange(5, 45)}
AGBT17B_455_01_ScanNums = {'7Pt-Cal1':np.arange(6, 15), '7Pt-Cal2':np.arange(111, 120)}
AGBT18A_443_01_ScanNums = {'7Pt-Cal1':np.arange(1, 10), '7Pt-Cal2':np.arange(181, 190)}

"""
Function to return the time stamp from provided project path and colleection of scan numbers:
"""
def getTimeStamp(projectName, scanArr, archFlag):
	if archFlag == 0:
		scanFitsPath = '/home/gbtdata/' + projectName + '/ScanLog.fits'
	else:
		scanFitsPath = '/home/archive/science-data/16B/' + projectName + '/ScanLog.fits'
	scanHdu = fits.open(scanFitsPath)
	scanData = scanHdu[1].data
	"""
	The dimensions of the binary data table of the ScanLog.fits file are (4*total scans)x3, 
	where the three columns are (0) time stamp (1) Scan number (3) Antenna/GO/DMJD Scan Start/ DMJD Scan Finish.
	Grab the correct element to grab is the 0th column of the 4*Scan row
	""" 
	timeList = [] ## define list to populate
	for i in range(0, len(scanArr)):
		scanNum = scanArr[i]
		timeStampElem = ((scanNum -1) * 4) 
		timeVal = scanData[timeStampElem][0]
		timeList.append(timeVal)
	timeList = [w.replace('T', '_') for w in timeList]
	timeList = [w.replace('-', '_') for w in timeList]
	timeList = [w + '.fits' for w in timeList]
	return timeList


"""
loop over time stamp list, read in DMJD values, and average
"""
def computeMean(projStr, seshNum, timeStampList, archFlag):
	dmjdList = []
	for timeStamp in timeStampList:
		if archFlag == 0:
			hdu = fits.open('/home/gbtdata/' + projStr + '_' + seshNum +'/Antenna/' + timeStamp)
		else:
			hdu = fits.open('/home/archive/science-data/16B/' + projStr + '_' + seshNum +'/Antenna/' + timeStamp)
		dmjdList.extend(hdu[2].data['DMJD'])
	return np.mean(dmjdList)


projectArr = ['AGBT16B_400', 'AGBT17B_360', 'AGBT17B_455', 'AGBT18A_443']
dataPath = '/home/scratch/npingel/FLAG/data/'

## initialize lists
PhaseDiffList_16B_XX = []
AmpDiffList_16B_XX = []
PhaseDiffList_16B_YY = []
AmpDiffList_16B_YY = []
meanDMJDList_16B = []

PhaseDiffList_17B_XX = []
AmpDiffList_17B_XX = []
PhaseDiffList_17B_YY = []
AmpDiffList_17B_YY = []
meanDMJDList_17B = []

"""
loop through the project and associated sessions to collate the weight vectors
"""
#for project in projectArr:
for projectIdx in range(0, 4):
	project = projectArr[projectIdx]
	sessPhaseList = []
	if project == 'AGBT16B_400':
		sessionArr = ['09', '12', '13', '14']
	elif project == 'AGBT17B_360':
		sessionArr = ['01', '02', '03', '04', '05', '06', '07']
	elif project == 'AGBT17B_455':
		sessionArr = ['01']
	elif project == 'AGBT18A_443':
		sessionArr = ['01']
	
	"""
	start loop over sessions
	"""
	for i in range(0, len(sessionArr)):
		sessNum = sessionArr[i]
		if project + '_' + sessNum == 'AGBT16B_400_03':
			weightLabel1 = 'weights'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT16B_400_03_ScanNums['grid']]
			archiveFlag = 1
		elif project + '_' + sessNum == 'AGBT16B_400_09':
			weightLabel1 = 'weights'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT16B_400_09_ScanNums['grid']]
			archiveFlag = 1
		elif project + '_' + sessNum == 'AGBT16B_400_12':
			weightLabel1 = 'FullGrid'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT16B_400_12_ScanNums['grid']]
			archiveFlag = 1
		elif project + '_' + sessNum == 'AGBT16B_400_13':
			weightLabel1 = 'FullGrid'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT16B_400_13_ScanNums['grid']]
			archiveFlag = 1
		elif project + '_' + sessNum == 'AGBT16B_400_14':
			weightLabel1 = 'w_14*'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT16B_400_14_ScanNums['7Pt-Cal']]
			archiveFlag = 1
		elif project + '_' + sessNum == 'AGBT17B_360_01':
			weightLabel1 = 'grid'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT17B_360_01_ScanNums['grid']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_360_02':
			weightLabel1 = '1st_seven'
			weightLabel2 = 'grid'
			weightLabel3 = '2nd_seven'
			labelArr = [weightLabel1, weightLabel2, weightLabel3]
			weightScanNumList = [AGBT17B_360_02_ScanNums['7Pt-Cal1'], AGBT17B_360_02_ScanNums['grid'], AGBT17B_360_02_ScanNums['7Pt-Cal2']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_360_03':
			weightLabel1 = 'seven'
			weightLabel2 = 'grid'
			labelArr = [weightLabel1, weightLabel2]
			weightScanNumList = [AGBT17B_360_03_ScanNums['7Pt-Cal'], AGBT17B_360_03_ScanNums['grid']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_360_04':
			weightLabel1 = 'grid'
			weightLabel2 = 'seven'
			labelArr = [weightLabel1, weightLabel2]
			weightScanNumList = [AGBT17B_360_04_ScanNums['grid'], AGBT17B_360_04_ScanNums['7Pt-Cal']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_360_05':
			weightLabel1 = '1st_seven'
			weightLabel2 = '2nd_seven'
			labelArr = [weightLabel1, weightLabel2]
			weightScanNumList = [AGBT17B_360_05_ScanNums['7Pt-Cal1'], AGBT17B_360_05_ScanNums['7Pt-Cal2']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_360_06':
			weightLabel1 = 'grid'
			labelArr = [weightLabel1]
			weightScanNumList = [AGBT17B_360_06_ScanNums['grid']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_360_07':
			weightLabel1 = 'seven'
			weightLabel2 = 'grid'
			labelArr = [weightLabel1, weightLabel2]
			weightScanNumList = [AGBT17B_360_07_ScanNums['7Pt-Cal'], AGBT17B_360_07_ScanNums['grid']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT17B_455_01':
			weightLabel1 = '1st_seven'
			weightLabel2 = '2nd_seven'
			labelArr = [weightLabel1, weightLabel2]
			weightScanNumList = [AGBT17B_455_01_ScanNums['7Pt-Cal1'], AGBT17B_455_01_ScanNums['7Pt-Cal2']]
			archiveFlag = 0
		elif project + '_' + sessNum == 'AGBT18A_455_01':
			weightLabel1 = '1st_seven'
			weightLabel2 = '2nd_seven'
			labelArr = [weightLabel1, weightLabel2]
			weightScanNumList = [AGBT18A_443_01_ScanNums['7Pt-Cal1'], AGBT18A_443_01_ScanNums['7Pt-Cal1']]
			archiveFlag = 0

			"""
			loop through calibration sessions for this project/session pair
			"""
		print('Collating data from ' + project + '_' + sessNum)
		for j in range(0, len(labelArr)):
			weightLabel = labelArr[j]
			scanNumArr = weightScanNumList[j]
			"""
			Now, read in the weight FITS files and select the data from XID = 10 (where the channel corresonding to the TOPO HI line is)
			"""
			weightList = glob.glob(dataPath + project + '/' + project + '_' + sessNum +'/weight_files/fits_files/' + '*' + weightLabel + '*')
			for wtFile in weightList:
				hdu = fits.open(wtFile)
				xidNum = hdu[0].header['XENGINE']
				if xidNum == 10:
					break
			"""
			select the 128 chunk of the 3200 length vector for the weights for the correct freq channel and boresight beam
			"""
			dataXX = hdu[1].data['Beam3X']
			dataYY = hdu[1].data['Beam3Y']
			vectorXX = dataXX[7*128:7*128 + 128]
			vectorYY = dataYY[7*128:7*128 + 128]
			if i == 0:
				initVecXX = np.zeros([19], dtype='complex64')
				initVecYY = np.zeros([19], dtype='complex64')
				initVecXX.real = vectorXX[0:38:2]
				initVecXX.imag = vectorXX[1:38:2]
				
				initVecYY.real = vectorYY[0:38:2]
				initVecYY.imag = vectorYY[1:38:2]

				## unit normalize			
				initVecMagXX = la.norm(initVecXX)
				initVec_NormXX = initVecXX / initVecMagXX

				initVecMagYY = la.norm(initVecYY)
				initVec_NormYY = initVecYY / initVecMagYY

				## make first element real
				initPhiXX = np.angle(initVec_NormXX)
				initVec_NormRealXX = np.exp(-1j*initPhiXX) * initVec_NormXX
				initVec_NormRealXX = np.matrix(initVec_NormRealXX)

				initPhiYY = np.angle(initVec_NormYY)
				initVec_NormRealYY = np.exp(-1j*initPhiYY) * initVec_NormYY
				initVec_NormRealYY = np.matrix(initVec_NormRealYY)
				"""
				acquire mean DMJD of calibration scan
				"""
				## get timestamps
				timeStampList = getTimeStamp(project + '_' + sessNum, scanNumArr, archiveFlag)
				## compute mean value
				initMeanDMJDVal = computeMean(project, sessNum, timeStampList, archiveFlag)
			if i > 0:
				vecXX = np.zeros([19], dtype = 'complex64')
				vecYY = np.zeros([19], dtype = 'complex64')

				vecXX.real = vectorXX[0:38:2]
				vecXX.imag = vectorXX[1:38:2]
				vecMagXX = la.norm(vecXX)

				vecYY.real = vectorYY[0:38:2]
				vecYY.imag = vectorYY[1:38:2]
				vecMagYY = la.norm(vecYY)

				## unit normalize 
				vec_NormXX = vecXX / vecMagXX
				vec_NormYY = vecYY / vecMagYY

				## compute correction angle
				phiHatXX = np.angle(np.matrix.transpose(initVec_NormRealXX) * np.matrix(vec_NormXX))
				phiHatYY = np.angle(np.matrix.transpose(initVec_NormRealYY) * np.matrix(vec_NormYY))

				## apply correction angle
				vec_NormCorrXX = np.exp(1j * phiHatXX)
				vec_NormCorrYY = np.exp(1j * phiHatYY)

				## compute vector difference
				diffVecXX = initVec_NormXX - vec_NormCorrXX
				diffVecMagXX = la.norm(diffVecXX)
				diffVecMagXX = la.norm(diffVecXX)

				diffVecYY = initVec_NormYY - vec_NormCorrYY
				diffVecMagYY = la.norm(diffVecYY)
				
				## compute explicit difference in phase
				#diffVecXX_Angle = np.angle(initVec_NormXX) - np.angle(vec_NormCorrXX)
				#print(diffVecXX_Angle)
				#diffVecMagXX = diffVecXX_Angle

				#diffVecYY_Angle = np.angle(initVec_NormYY) - np.angle(vec_NormCorrYY)
                                #diffVecMagYY = diffVecYY_Angle

				## compute amplitude difference
				scaleXX = vecMagXX / initVecMagXX
				scaleYY = vecMagYY / initVecMagYY

				## add to lists based on project
				if project in {'AGBT17B_360', 'AGBT17B_455', 'AGBT18A_443'}:
					PhaseDiffList_17B_XX.append(diffVecMagXX)
					AmpDiffList_17B_XX.append(scaleXX)	
					PhaseDiffList_17B_YY.append(diffVecMagYY)
					AmpDiffList_17B_YY.append(scaleYY)	

					## get timestamps
					timeStampList = getTimeStamp(project + '_' + sessNum, scanNumArr, archiveFlag)
					## compute mean value
					meanDMJDVal = computeMean(project, sessNum, timeStampList, archiveFlag)
					meanDMJDList_17B.append(meanDMJDVal - initMeanDMJDVal)
				elif project == 'AGBT16B_400':
					PhaseDiffList_16B_XX.append(diffVecMagXX)
					AmpDiffList_16B_XX.append(scaleXX)

					PhaseDiffList_16B_YY.append(diffVecMagYY)
					AmpDiffList_16B_YY.append(scaleYY)

					## get timestamps
					timeStampList = getTimeStamp(project + '_' + sessNum, scanNumArr, archiveFlag)
					## compute mean value
					meanDMJDVal = computeMean(project, sessNum, timeStampList, archiveFlag)
					meanDMJDList_16B.append(meanDMJDVal - initMeanDMJDVal)

majorYLocator = MultipleLocator(.1)
majorYFormatter = FormatStrFormatter('%.1f')
minorYLocator = MultipleLocator(.01)

majorXLocator = MultipleLocator(1)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(0.25)

print(np.std(np.array(PhaseDiffList_17B_XX)))

## plot XX Pol 
fig, ax = pyplot.subplots(figsize = (8,8))
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)
ax.yaxis.set_minor_locator(minorYLocator)
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
ax.xaxis.set_minor_locator(minorXLocator)
ax.tick_params(axis = 'both', which='both', width=2)
pyplot.xlabel('Days', fontsize = 18)
pyplot.ylim(0.85, 1.15)
pyplot.axhline(1.0, linestyle = '--', linewidth = 2, color = 'black')
pyplot.ylabel(r'Normalized $d_1$', fontsize = 18)
pyplot.plot(meanDMJDList_16B, PhaseDiffList_16B_XX[:] / PhaseDiffList_16B_XX[0], marker = 'o', linestyle = 'None', markersize=8, color =tableau20[0], label = '16B_400 (XX)')
pyplot.plot(meanDMJDList_17B, PhaseDiffList_17B_XX[:] / PhaseDiffList_17B_XX[0], marker = 'o', linestyle = 'None', markersize=8, color =tableau20[2], label = '17B_360 (XX)')
pyplot.legend(loc=0,prop={'size':14},fontsize='large')
#pyplot.savefig('Weight_Phase_Comparison_Boresight_XX.pdf')
pyplot.show()
pyplot.clf()
pyplot.close()

## plot YY Pol 
fig, ax = pyplot.subplots(figsize = (8,8))
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)
ax.yaxis.set_minor_locator(minorYLocator)
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
ax.xaxis.set_minor_locator(minorXLocator)
ax.tick_params(axis = 'both', which='both', width=2)
pyplot.xlabel('Days', fontsize = 18)
pyplot.ylabel(r'Difference', fontsize = 18)
pyplot.plot(meanDMJDList_16B, PhaseDiffList_16B_YY, marker = 'o', linestyle = 'None', markersize=8, color =tableau20[0], label = '16B_400 (YY)')
pyplot.plot(meanDMJDList_17B, PhaseDiffList_17B_YY, marker = 'o', linestyle = 'None', markersize=8, color =tableau20[2], label = '17B_360 (YY)')
pyplot.legend(loc=0,prop={'size':14},fontsize='large')
pyplot.axhline(1.0, linestyle = '--', linewidth = 2, color = 'black')
#pyplot.savefig('Weight_Phase_Comparison_Boresight_YY.pdf')
#pyplot.show()
pyplot.clf()
pyplot.close()

majorYLocator = MultipleLocator(5)
majorYFormatter = FormatStrFormatter('%.1f')
minorYLocator = MultipleLocator(1)

majorXLocator = MultipleLocator(1)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(0.25)

## plot XX pol
fig, ax = pyplot.subplots(figsize = (8,8))
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)
ax.yaxis.set_minor_locator(minorYLocator)
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
ax.xaxis.set_minor_locator(minorXLocator)
pyplot.axhline(1.0, linestyle = '--', linewidth = 2, color = 'black')
ax.tick_params(axis = 'both', which='both', width=2)
pyplot.xlabel('Days', fontsize = 18)
pyplot.ylabel(r'Amplitude', fontsize = 18)
pyplot.plot(meanDMJDList_16B, AmpDiffList_16B_XX, marker = 'o', linestyle = 'None', linewidth = 2, markersize=8, color =tableau20[0], label = '16B_400 (XX)')
pyplot.plot(meanDMJDList_17B, AmpDiffList_17B_XX, marker = 'o', linestyle = 'None', linewidth = 2, markersize=8, color =tableau20[2], label = '17B_360 (XX)')
pyplot.legend(loc=0,prop={'size':14},fontsize='large')
#pyplot.savefig('Weight_Scale_Comparison_Boresight_XX.pdf')
pyplot.show()
pyplot.clf()
pyplot.close()

## plot YY pol
fig, ax = pyplot.subplots(figsize = (8,8))
ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)
ax.yaxis.set_minor_locator(minorYLocator)
ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
ax.xaxis.set_minor_locator(minorXLocator)
ax.tick_params(axis = 'both', which='both', width=2)
pyplot.xlabel('Days', fontsize = 18)
pyplot.ylabel(r'Amplitude', fontsize = 18)
pyplot.plot(meanDMJDList_16B, AmpDiffList_16B_YY, marker = 'o', linestyle = 'None', linewidth = 2, markersize=8, color =tableau20[0], label = '16B_400 (YY)')
pyplot.plot(meanDMJDList_17B, AmpDiffList_17B_YY, marker = 'o', linestyle = 'None', linewidth = 2, markersize=8, color =tableau20[2], label = '17B_360 (YY)')
pyplot.legend(loc=0,prop={'size':14},fontsize='large')
#pyplot.savefig('Weight_Scale_Comparison_Boresight_YY.pdf')
#pyplot.show()
pyplot.clf()
pyplot.close()
