"""
Script to generate Tsys/eta as a function of DMJD from AGBT17B_360 (plus the single sessions from AGBT17B_455 and AGBT18A_443)
"""

from astropy.io import fits 
import numpy as np
import glob
import sys
import pickle
import matplotlib.pyplot as pyplot
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Helvetica Neue') 
matplotlib.rc('text', usetex='false') 
matplotlib.rc('xtick.major.width')
matplotlib.rcParams.update({'font.size': 14})

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

 ## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {0:1, 1:2, 2:6, 3:0, 4:3, 5:5,6:4}
byuBeamDict = {1:0, 2:1, 6:2, 0:3, 3:4, 5:5, 4:6}

## define dataDir
dataDir = '/home/scratch/npingel/FLAG/data'
projStr = 'AGBT16B_400'
sessionArr = ['02', '03', '09', '12', '13', '14', '01', '02', '03', '04', '05', '06', '07', '01']

## define lists to store data (7Pt-Cal1)
Tsys_Y_7Pt1 = []
Tsys_X_7Pt1 = []
TsysErr_Y_7Pt1 = []
TsysErr_X_7Pt1 = []

## define lists to store data (7Pt-Cal1)
Tsys_Y_7Pt2 = []
Tsys_X_7Pt2 = []
TsysErr_Y_7Pt2 = []
TsysErr_X_7Pt2 = []

## define lists to store data (7Pt-Cal)
Tsys_Y_7Pt = []
Tsys_X_7Pt = []
TsysErr_Y_7Pt = []
TsysErr_X_7Pt = []

## define lists to store data (grids)
Tsys_Y_Grid = []
Tsys_X_Grid = []
TsysErr_Y_Grid = []
TsysErr_X_Grid = []

## define lists to hold meanDMJDValues
dmjdMean_7PtCal1 = []
dmjdMean_7PtCal2 = []
dmjdMean_7PtCal = []
dmjdMean_grid = []

## initialize counters
PtCal1Cnt = 0
PtCal2Cnt = 0
PtCalCnt = 0
gridCnt = 0

## dictionaries that store the scan numbers of a given calibration scan type per session
AGBT16B_400_02_ScanNums = {'grid':np.arange(13,57)}
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
def getTimeStamp(projectName, scanArr, cnt):
	if cnt > 5:
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
def computeMean(projStr, seshNum, timeStampList, cnt):
	dmjdList = []
	for timeStamp in timeStampList:
		if cnt > 5:
			hdu = fits.open('/home/gbtdata/' + projStr + '_' + seshNum +'/Antenna/' + timeStamp)
		else:
			hdu = fits.open('/home/archive/science-data/16B/' + projStr + '_' + seshNum +'/Antenna/' + timeStamp)
		dmjdList.extend(hdu[2].data['DMJD'])
	return np.mean(dmjdList)

"""
loop to collate Tsys data
"""
cnt=0
mapVal = False
tMap = 31*0.5*71
tOn = tMap/115.38 ## GBT BEAMS PER 3 sq. deg map
tOff = 0.5 * 8
teff = tOn*tOff/(tOn + tOff)
for i in range(0, len(sessionArr)):
	seshNum = sessionArr[i]
	if i == len(sessionArr) - 1:
		projStr = 'AGBT17B_455'
	#elif i == 8:
	#	projStr = 'AGBT18A_443'
	elif i > 5:
		projStr = 'AGBT17B_360'

	fullProjName = dataDir + '/' + projStr + '/' + projStr + '_' + seshNum 
	print('Collating data from ' + fullProjName)
	txtList = glob.glob(fullProjName + '/Tsys/*.txt')
	"""
	loop to extract data from txt files and put in either grid or 7pt-Cal lists
	"""
	for j in range(0, len(txtList)):
		txtFile = txtList[j].split('/')[-1]
		splitVal = txtFile.split('_')
		pol = splitVal[1]
		data = np.loadtxt(txtList[j])
		calType = splitVal[-1]
		calType = calType[:-4]
		## based on calType, place in correct holder list
		if calType == '7Pt-Cal1':
			PtCal1Cnt += 1
			if pol == 'X':
				Tsys_X_7Pt1.extend(data[:, 0])
				TsysErr_X_7Pt1.extend(data[:, 1])
				Tsys_XArr = data[:, 0]
				TsysErrArr = data[:, 1]
			else:
				Tsys_Y_7Pt1.extend(data[:, 0])
				TsysErr_Y_7Pt1.extend(data[:, 1])
			if PtCal1Cnt % 2 == 0:
				## now, determine which dictionary to use and grab scan numbers
				if projStr + '_' + seshNum == 'AGBT17B_360_01':
					scanNumArr = AGBT17B_360_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr, i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal1.append(meanVal)

				elif projStr + '_' + seshNum == 'AGBT17B_360_02':
					scanNumArr = AGBT17B_360_02_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr, i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal1.append(meanVal)

				elif projStr + '_' + seshNum == 'AGBT17B_360_05':
					scanNumArr = AGBT17B_360_05_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal1.append(meanVal)
					mapVal = True
					numMaps = 1
				elif projStr + '_' + seshNum == 'AGBT17B_455_01':
					scanNumArr = AGBT17B_455_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal1.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT18A_443_01':
					scanNumArr = AGBT18A_443_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal1.append(meanVal)
			"""
			Compute the theoretical Noise for each beam if we made a map during that session
			"""
			if mapVal == True:
				sigArr = np.zeros([7]) 
				errArr = np.zeros([7])
				for z in range(0, 7):
					#print('The theoretical noise for Beam ' + np.str(z) + ': ' + np.str(Tsys_XArr[z]/np.sqrt(24.414*1000*teff*numMaps)) + '+/-' + np.str(TsysErrArr[z]/np.sqrt(24.414*1000*teff*numMaps)))
					sigArr[z] = Tsys_XArr[z]/np.sqrt(24.414*1000*teff*numMaps)
					errArr[z] = TsysErrArr[z]/np.sqrt(24.414*1000*teff*numMaps)
				tSysComb = np.mean(Tsys_XArr)/np.sqrt(7*24.414*1000*teff*numMaps)
				tSysCombErr = np.mean(errArr)/np.sqrt(7)
				#print('The theoretical noise for combined map: ' + np.str(tSysComb) + '+/-' + np.str(tSysCombErr))
				mapVal = False
		elif calType == '7Pt-Cal2':
			PtCal2Cnt += 1
			if pol == 'X':
				Tsys_X_7Pt2.extend(data[:, 0])
				TsysErr_X_7Pt2.extend(data[:, 1])
			else:
				Tsys_Y_7Pt2.extend(data[:, 0])
				TsysErr_Y_7Pt2.extend(data[:, 1])
			if PtCal2Cnt % 2 == 0:
				## now, determine which dictionary to use and grab scan numbers
				if projStr + '_' + seshNum == 'AGBT17B_360_01':
					scanNumArr = AGBT17B_360_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal2.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_02':
					scanNumArr = AGBT17B_360_02_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal2.append(meanVal)

				elif projStr + '_' + seshNum == 'AGBT17B_360_05':
					scanNumArr = AGBT17B_360_05_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal2.append(meanVal)

				elif projStr + '_' + seshNum == 'AGBT17B_455_01':
					scanNumArr = AGBT17B_455_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal2.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT18A_443_01':
					scanNumArr = AGBT18A_443_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal2.append(meanVal)
		elif calType == '7Pt-Cal':
			PtCalCnt += 1
			if pol == 'X':
				Tsys_X_7Pt.extend(data[:, 0])
				TsysErr_X_7Pt.extend(data[:, 1])
			else:
				Tsys_Y_7Pt.extend(data[:, 0])
				TsysErr_Y_7Pt.extend(data[:, 1])
			if PtCalCnt % 2 == 0:
				## now, determine which dictionary to use and grab scan numbers
				if 	projStr + '_' + seshNum == 'AGBT16B_400_12':
					scanNumArr = AGBT16B_400_12_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT16B_400_13':
					scanNumArr = AGBT16B_400_13_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT16B_400_14':
					scanNumArr = AGBT16B_400_14_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_03':
					scanNumArr = AGBT17B_360_03_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_04':
					scanNumArr = AGBT17B_360_04_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_07':
					scanNumArr = AGBT17B_360_04_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_7PtCal.append(meanVal)

		elif calType == 'grid':
			gridCnt += 1
			if pol == 'X':
				Tsys_X_Grid.extend(data[:, 0])
				TsysErr_X_Grid.extend(data[:, 1])
				Tsys_XArr = data[:, 0]
				TsysErrArr = data[:, 1]
			else:
				Tsys_Y_Grid.extend(data[:, 0])
				TsysErr_Y_Grid.extend(data[:, 1])	

			if gridCnt % 2 == 0:
				## now, determine which dictionary to use and grab scan numbers
				if projStr + '_' + seshNum == 'AGBT16B_400_02':
					scanNumArr = AGBT16B_400_02_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 1
				elif projStr + '_' + seshNum == 'AGBT16B_400_03':
					scanNumArr = AGBT16B_400_03_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 1
				elif projStr + '_' + seshNum == 'AGBT16B_400_09':
					scanNumArr = AGBT16B_400_09_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 1
				elif projStr + '_' + seshNum == 'AGBT16B_400_12':
					scanNumArr = AGBT16B_400_12_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 1
				elif projStr + '_' + seshNum == 'AGBT16B_400_13':
					scanNumArr = AGBT16B_400_13_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 1
				elif projStr + '_' + seshNum == 'AGBT16B_400_14':
					scanNumArr = AGBT16B_400_14_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_01':
					scanNumArr = AGBT17B_360_01_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_02':
					scanNumArr = AGBT17B_360_02_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_03':
					scanNumArr = AGBT17B_360_03_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 3
				elif projStr + '_' + seshNum == 'AGBT17B_360_04':
					scanNumArr = AGBT17B_360_04_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
					mapVal = True
					numMaps = 4
				elif projStr + '_' + seshNum == 'AGBT17B_360_06':
					scanNumArr = AGBT17B_360_06_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)
				elif projStr + '_' + seshNum == 'AGBT17B_360_07':
					scanNumArr = AGBT17B_360_07_ScanNums[calType]
					timeStampList = getTimeStamp(projStr + '_' + seshNum, scanNumArr,i)
					meanVal = computeMean(projStr, seshNum, timeStampList, i)
					dmjdMean_grid.append(meanVal)	
					mapVal = True
					numMaps = 4.5
			"""
			Compute the theoretical Noise for each beam if we made a map during that session
			"""
			if mapVal == True:
				sigArr = np.zeros([7]) 
				errArr = np.zeros([7])
				for z in range(0, 7):
					print('The theoretical noise for Beam ' + np.str(z) + ': ' + np.str(Tsys_XArr[z]/np.sqrt(24.414*1000*teff*numMaps)) + '+/-' + np.str(TsysErrArr[z]/np.sqrt(24.414*1000*teff*numMaps)))
					sigArr[z] = Tsys_XArr[z]/np.sqrt(24.414*1000*teff*numMaps)
					errArr[z] = TsysErrArr[z]/np.sqrt(24.414*1000*teff*numMaps)
				tSysComb = np.mean(Tsys_XArr)/np.sqrt(7*24.414*1000*teff*numMaps)
				tSysCombErr = np.mean(errArr)/np.sqrt(7)
				print('The theoretical noise for combined map: ' + np.str(tSysComb) + '+/-' + np.str(tSysCombErr))
				mapVal = False

## replace zero values with nans
Tsys_X_7Pt1 = [np.nan if x == 0 else x for x in Tsys_X_7Pt1]
Tsys_Y_7Pt1 = [np.nan if x == 0 else x for x in Tsys_Y_7Pt1]
Tsys_X_7Pt2 = [np.nan if x == 0 else x for x in Tsys_X_7Pt2]
Tsys_Y_7Pt2 = [np.nan if x == 0 else x for x in Tsys_Y_7Pt2]
Tsys_X_7Pt = [np.nan if x == 0 else x for x in Tsys_X_7Pt]
Tsys_Y_7Pt = [np.nan if x == 0 else x for x in Tsys_Y_7Pt]
Tsys_X_Grid = [np.nan if x == 0 else x for x in Tsys_X_Grid]
Tsys_Y_Grid = [np.nan if x == 0 else x for x in Tsys_Y_Grid]


## now, combine 7Pt scans
Tsys_X_Full7Pt = Tsys_X_7Pt + Tsys_X_7Pt1 + Tsys_X_7Pt2
Tsys_Y_Full7Pt = Tsys_Y_7Pt + Tsys_Y_7Pt1 + Tsys_Y_7Pt2
TsysErr_X_Full7Pt = TsysErr_X_7Pt + TsysErr_X_7Pt1 + TsysErr_X_7Pt2
TsysErr_Y_Full7Pt = TsysErr_Y_7Pt + TsysErr_Y_7Pt1 + TsysErr_Y_7Pt2

dmjdMean_7PtCal[0:3] = dmjdMean_7PtCal[0:3] - dmjdMean_7PtCal[0]
dmjdMean_7PtCal[3:] = dmjdMean_7PtCal[3:] - dmjdMean_7PtCal[3]
dmjdMean_7PtCal1[:] = dmjdMean_7PtCal1[:] - dmjdMean_7PtCal1[0]
dmjdMean_7PtCal2[:] = dmjdMean_7PtCal2[:] - dmjdMean_7PtCal2[0]


dmjdMean_Full7Pt = dmjdMean_7PtCal + dmjdMean_7PtCal1 + dmjdMean_7PtCal2


#dmjdMean_Full7Pt[0:3] = dmjdMean_Full7Pt[0:3] - dmjdMean_Full7Pt[0]
#dmjdMean_Full7Pt[3:] = dmjdMean_Full7Pt[3:] - dmjdMean_Full7Pt[3]
dmjdMean_grid[:2] = dmjdMean_grid[:2] - dmjdMean_grid[0]
dmjdMean_grid[2:6] = dmjdMean_grid[2:6] - dmjdMean_grid[2]
dmjdMean_grid[6:] = dmjdMean_grid[6:] - dmjdMean_grid[6] 

Tsys_X_ALL = Tsys_X_7Pt1 + Tsys_X_7Pt2 + Tsys_X_7Pt + Tsys_X_Grid
Tsys_Y_ALL = Tsys_Y_7Pt1 + Tsys_Y_7Pt2 + Tsys_Y_7Pt + Tsys_Y_Grid
TsysErr_X_ALL = TsysErr_X_7Pt1 + TsysErr_X_7Pt2 + TsysErr_X_7Pt + TsysErr_X_Grid
TsysErr_Y_ALL = TsysErr_Y_7Pt1 + TsysErr_Y_7Pt2 + TsysErr_Y_7Pt + TsysErr_Y_Grid
dmjd_ALL = dmjdMean_7PtCal1 + dmjdMean_7PtCal2 + dmjdMean_7PtCal + dmjdMean_grid




majorYLocator = MultipleLocator(5)
majorYFormatter = FormatStrFormatter('%.1f')
minorYLocator = MultipleLocator(1)

majorXLocator = MultipleLocator(1)
majorXFormatter = FormatStrFormatter('%d')
minorXLocator = MultipleLocator(0.25)



"""
Now plot all the computed values. A total of nine plots will be produced:
- Tsys/eta vs. DMJD (of both polarizations) for each inididual beam 
- Tsys/eta vs. DMJD for all beams in XX pol
- Tsys/eta vs. DMJD for all beams in YY pol
"""

"""
Plot in a series of 4x2 (columns x rows) subplot. Leave last empty
"""
for bm in range(0, 9):
	if bm < 7:
		fig, ax = pyplot.subplots(figsize=(8,8))
		ax.yaxis.set_major_locator(majorYLocator)
		ax.yaxis.set_major_formatter(majorYFormatter)
		ax.yaxis.set_minor_locator(minorYLocator)
		ax.xaxis.set_major_locator(majorXLocator)
		ax.xaxis.set_major_formatter(majorXFormatter)
		ax.xaxis.set_minor_locator(minorXLocator)
		ax.tick_params(axis = 'both', which='both', width=2)
		pyplot.xlim(-0.25, np.max(dmjdMean_Full7Pt - dmjdMean_Full7Pt[0])+0.5)
		pyplot.ylim(20, 70)
		bmName = np.str(bm)
		pyplot.xlabel('Days')
		pyplot.ylabel(r'Tsys [K]/$\eta$')
		pyplot.title(r'Beam ' + np.str(bmName))
		pyplot.errorbar(dmjdMean_Full7Pt, np.array(Tsys_X_Full7Pt[bm::7]) / 0.6, np.array(TsysErr_X_Full7Pt[bm::7]) / 0.6, color=tableau20[0], fmt='o', label = 'XX (7-Pt)', markersize=8)
		pyplot.errorbar(dmjdMean_Full7Pt, np.array(Tsys_Y_Full7Pt[bm::7]) / 0.6, np.array(TsysErr_Y_Full7Pt[bm::7]) / 0.6, color=tableau20[2], fmt='o', label = 'YY (7-Pt)', markersize=8)
		pyplot.errorbar(dmjdMean_grid, np.array(Tsys_X_Grid[bm::7]) / 0.6, np.array(TsysErr_X_Grid[bm::7]) / 0.6, color=tableau20[0], fmt='s', label = 'XX (grid)', markersize=8)
		pyplot.errorbar(dmjdMean_grid, np.array(Tsys_Y_Grid[bm::7]) / 0.6, np.array(TsysErr_Y_Grid[bm::7]) / 0.6, color=tableau20[2], fmt='s', label = 'YY (grid)', markersize=8)
		pyplot.axhline(18/0.65, linestyle = '--', color = 'black', linewidth = 2, label = 'Single Pixel')
                pyplot.legend(loc=0,prop={'size':14},fontsize='large')
		pyplot.tight_layout()
		pyplot.savefig('Tsys_eta_vs_Time_BEAM_' + bmName + '.pdf')
		pyplot.show()
		pyplot.clf()
		pyplot.close()

		"""
		compute mean Tsys value while being careful to 
		correctly handle nan's and associated uncertainty 
		"""
		xNanInds = np.where(np.isnan(np.array(Tsys_X_ALL[bm::7])))
		yNanInds = np.where(np.isnan(np.array(Tsys_Y_ALL[bm::7])))
		TsysErr_X_ALL_Bm = np.array(TsysErr_X_ALL[bm::7])/0.6
		TsysErr_Y_ALL_Bm = np.array(TsysErr_Y_ALL[bm::7])/0.6
		TsysErr_X_ALL_Bm[xNanInds] = np.nan
		TsysErr_Y_ALL_Bm[yNanInds] = np.nan
		numXNans = len(xNanInds[0])
		numYNans = len(yNanInds[0])
		print('Mean Tsys/eta for XX, Beam ' + bmName + ': ' + np.str(np.nanmean(np.array(Tsys_X_ALL[bm::7])/0.6)) + '+/-' + np.str(np.sqrt(np.nansum((np.array(TsysErr_X_ALL_Bm)/.6)**2)/(len(TsysErr_X_ALL_Bm) - numXNans))))
		print('Mean Tsys/eta for YY, Beam ' + bmName + ': ' + np.str(np.nanmean(np.array(Tsys_Y_ALL[bm::7])/0.6)) + '+/-' + np.str(np.sqrt(np.nansum((np.array(TsysErr_Y_ALL_Bm)/.6)**2)/(len(TsysErr_Y_ALL_Bm) - numYNans))))
	else:
		for i in range(0, 7):
			bmN = i
			if i == 0:
				fig, ax = pyplot.subplots(figsize = (8,8))
				ax.yaxis.set_major_locator(majorYLocator)
				ax.yaxis.set_major_formatter(majorYFormatter)
				ax.yaxis.set_minor_locator(minorYLocator)
				ax.xaxis.set_major_locator(majorXLocator)
				ax.xaxis.set_major_formatter(majorXFormatter)
				ax.xaxis.set_minor_locator(minorXLocator)
				ax.tick_params(axis = 'both', which='both', width=2)
				pyplot.xlim(-0.25, np.max(dmjd_ALL - dmjd_ALL[0])+0.5)
				pyplot.ylim(20, 80)
				pyplot.axhline(18/0.65, linestyle = '--', color = 'black', linewidth = 2, label = 'Single Pixel')
				pyplot.xlabel('Days')
				pyplot.ylabel(r'Tsys [K]/$\eta$')
			if bm == 7:
				polStr = 'XX'
				pyplot.errorbar(dmjd_ALL, np.array(Tsys_X_ALL[bmN::7])/0.6, np.array(TsysErr_X_ALL[bmN::7])/0.6, color=tableau20[i*2], fmt='o', label = 'Beam ' + np.str(i) + ' (XX)', markersize=8)
			else:
				polStr = 'YY'
				pyplot.errorbar(dmjd_ALL, np.array(Tsys_Y_ALL[bmN::7])/0.6, np.array(TsysErr_Y_ALL[bmN::7])/0.6, color=tableau20[i*2], fmt='o', label = 'Beam ' + np.str(i) + ' (YY)', markersize=8)
		pyplot.legend(loc=0,prop={'size':10},fontsize='xx-small')
		pyplot.savefig('Tsys_eta_vs_Time_ALL_BEAMS_' + polStr + '.pdf')
		pyplot.show()
		pyplot.clf()
		pyplot.close()




