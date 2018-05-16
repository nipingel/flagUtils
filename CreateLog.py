"""
Script to generate a .txt file that summarizes an observations session. The timestamp, procedure name, 
observed object, scan number, schedule block name, current sequence number, total sequence number, integration length,
and mode name are written out as respective rows to PROJECT_SESSION_LOG.txt

usage: ipython CreateLog.py PROJECT_SESSION
"""

from astropy.io import fits
from astropy.io import ascii
import numpy as np
import os
import glob
import sys

## read in project_session
projectSession = sys.argv[1]

if projectSession[:-3] == 'AGBT16B_400':
	dataPath = '/home/archive/science-data/16B/'
else:
	dataPath = '/home/gbtdata/'

## get scanLog.fits
scanHdu = fits.open(dataPath + projectSession + '/ScanLog.fits')
goHduList = glob.glob(dataPath + projectSession + '/GO/*.fits')
goHduList.sort()


cnt = 0
scanNums = scanHdu[1].data['SCAN']
scanNums = scanNums[::4]

## initialize lists to store table entries
procNameList = []
objList = []
obsTimeList = []
scanList = []
blockNameList = []
seqSizeList = []
seqNumList = []
intList = []
modeList = []


for i in range(0, len(goHduList)):
	goFitsFile = goHduList[i] 
	fileName = goFitsFile.split('/')[-1]
	## if file not found, skip. A GO file was written, but no data was taken
	try:
		dataHdu = fits.open('/lustre/projects/flag.old/' + projectSession + '/BF/' + fileName[:-5] + 'A.fits')
	except IOError:
		continue
		pass
	hdu = fits.open(goFitsFile)
	procNameList.append(hdu[0].header['PROCNAME'])
	objList.append(hdu[0].header['OBJECT'])
	obsTimeList.append(fileName[:-5])
	intList.append(dataHdu[0].header['REQSTI'])
	modeList.append(dataHdu[0].header['MODENAME'])
	scanList.append(np.str(scanNums[i]))
	blockNameList.append(hdu[0].header['BLOCK'])
	seqSizeList.append(np.str(hdu[0].header['PROCSIZE']))
	seqNumList.append(np.str(hdu[0].header['PROCSEQN']))

## create dictionary that will write out as table
data = {'Time':obsTimeList, 
	'ProcedureName': procNameList, 
	'Object': objList,
	'ScanNumber':scanList,
	'ScheduleBlock': blockNameList,
	'SequenceNumber': seqNumList,
	'TotalSequenceNumber': seqSizeList, 
	'IntLength': intList,
	'Mode': modeList
}

ascii.write(data, output = projectSession + '_LOG.txt', names=['Time', 'ProcedureName', 'Object', 'ScanNumber', 'ScheduleBlock', 'SequenceNumber', 'TotalSequenceNumber', 'IntLength', 'Mode'])
#thefile = open(projectSession + '.txt', 'w')
#with open(projectSession + '.txt', 'wb') as out:
#    for row in rowList:
#        out.write(' '.join(str(num) for num in row))
#        out.write('\n')



