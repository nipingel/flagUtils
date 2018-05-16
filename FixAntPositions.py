"""
3/2/18
Script to correct Antenna positions of PAF observations. Based on the input file and object, the associated Antenna 
and data FITS files will be read in, the Antenna position interpolated to the data (using the respective DMJD values), 
and beam offsets applied. Once the beam offsets are applied, the horiztonal coordinates are transformed to J2000 or Galactic
coordinates. A doppler correction is additionally performed. 
User Inputs:
fileName - input FITS file
observedObj - the observed object
Usage:
ipython FixAntPositions.py fileName observedObj
Eample:
ipython FixAntPositions.py AGBT17B_455_01_G353-4.0_Beam0_1st_seven.fits G353-4.0
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as pyplot
import pyslalib as pysla
import os 
import numpy as np
import datetime
import collections
import glob
import sys
sys.path.insert(0, '/home/scratch/npingel/FLAG/pros/SpectralFiller/modules/')
import RadVelCorr

"""
Function that generates and returns list of observed objects and the associated timestamp FITS files. 
The GO FITS files are searched and sorted to generate the required lists, which is why the path to these 
files is an input to this function. 
"""
def generateObjAndFitsLists(goFitsPath):        

    goFits = glob.glob(goFitsPath+'/*.fits') ## read in GO FITS files to
    goFits.sort() ## sort to get correct time stamp order
    genObjList = [] ## define empty list to be returned (after casting to list for indexing)
    """
    define empty list. Once we add new object, an empty list will be appended and
    associated FITS files added. For example, if the second object added is '3C147', 
    genFitsList[1] will be a list of all associated timestamp FITS files. 
    """
    genFitsList = [] 
    itr = 0 ## counter variable to check if first iteration
    ## iterate over all GO FITS files to check the observed object
    for goName in goFits:
       goHDU = fits.open(goName)
       if itr == 0: ## upon first iteration, append an empty list to the generated list of FITS files
           genObjList.append(goHDU[0].header['OBJECT'])
           genFitsList.append([])
           genFitsList[0].append(goName[-24:])
           itr+=1
       else: ## now processing for remaining 
           obj = goHDU[0].header['OBJECT'] 
           if not any(x == obj for x in genObjList): ## test if object already in list
               genObjList.append(obj) 
               genFitsList.append([]) ## append empty list
               genFitsList[-1].append(goName[-24:])
           else:
            ind = genObjList.index(obj)
            genFitsList[ind].append(goName[-24:])
    return genObjList, genFitsList

fileName = sys.argv[1]
weightPath = sys.argv[2]
observedObj = sys.argv[3]

hdu = fits.open(fileName, mode = 'update') ## read in FITS file
projectID = fileName[:14] ## get project name
weightFileList = glob.glob(weightPath)
weightHdu = fits.open(weightFileList[0]) ## get weight file


## read in ancilllary FITS files
scanFitsHdu = fits.open('/home/gbtdata/' + projectID + '/ScanLog.fits')
goFitsPath = '/home/gbtdata/' + projectID + '/GO'

## generate object and list of FITS files in case user wishes to process all object/timestamps
allObjList, allFitsList = generateObjAndFitsLists(goFitsPath)

print(observedObj)
objInd = allObjList.index(observedObj)
source = allObjList[objInd] ## get source
fileList = allFitsList[objInd] ## get FITS list for object (just fileneams; no path)

## define lists that will hold data coordinates
majList = []
minList = []
raList = []
decList = []
azList = []
elList = []
lstStartList = []
dataDMJDList = []
refractList = []

GBTLAT  = np.deg2rad(38.4331294)
GBTLONG = np.deg2rad(79.8398397)
GBTHGT  = 824.36                     # meters above the ellipsoid

def radVelCorrection(raArr, decArr):

	print('\n')
	print('Applying Doppler Correction to HELIOCENTRIC...')

	radvelcorrObj = RadVelCorr.RadVelCorr() ## initialize object

	## collect all relevant data for doppler calculation
	utDate_Time = hdu[1].data['DATE-OBS']
	cenFreqsArr = hdu[1].data['CRVAL1']
	restFreqArr = hdu[1].data['RESTFREQ']
	c = 299792458.0 ## m/s
	for velIter in range(0,len(raArr)): ## loop through to calculate each integration/polarization's correciton 
	  raVal = raArr[velIter] ## RA
	  decVal = decArr[velIter] ## Dec 
	  utDateTimeVal = utDate_Time[velIter] ## UT Time
	  cenFreqVal = cenFreqsArr[velIter] ## center frequency value
	  restFreqVal = restFreqArr[velIter] ## restfreq
	  ## put time in necessary format
	  t = Time(utDateTimeVal, format='isot')
	  utdate = t.iso[0:10]
	  uttime = t.iso[11:]
	  radVelCorrection = radvelcorrObj.correctVel(utdate,uttime,raVal,decVal) ## calculate correction
	  ## compute optical velocity of ref freq
	  vOpt = c*(1-cenFreqVal/restFreqVal)
	  ## add radial correction
	  newVOpt = vOpt + radVelCorrection
	  ## now convert back to frequency
	  newCenFreq = (1-newVOpt/c)*restFreqVal 
	  ## update reference frequency value
	  cenFreqsArr[velIter] = newCenFreq

	return cenFreqsArr

def az2ra(LST, Az, El, dec, coordSys):
	## convert all angles to radians
	azRad = np.deg2rad(Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
	posInds = np.where(Az > 360)
	negInds = np.where(np.logical_and(Az < 180, Az > 0)) ## need to rotate multiply by negative 1 while in 
	elRad = np.deg2rad(El)
	decRad = np.deg2rad(dec)
	ha = (-1) * np.arccos(1/np.cos(decRad)*(np.sin(elRad)*np.cos(GBTLAT) - np.cos(elRad)*np.cos(azRad)*np.sin(GBTLAT)))
	ha[posInds] = ha[posInds] * (-1)
	ha[negInds] = ha[negInds] * (-1)
	return LST + np.rad2deg(ha)
 
"""  
method to convert Elevation to Declination (geoapparent)
""" 
def el2dec(LST, Az, El):
	## convert all angles to radians
	azRad = np.deg2rad(Az) ## Azimuth of NCP is 0 for GBT. This code assumes 180. 
	elRad = np.deg2rad(El)
	el = np.rad2deg(np.arcsin(np.sin(elRad)*np.sin(GBTLAT) + np.cos(elRad)*np.cos(azRad)*np.cos(GBTLAT)))
	return el

"""
method to transform the beam offsets to Ra/Dec coordinates. The current state of the binary table HDU is passed in. 
That same binary table HDU with the updated coordinates is returned. 
"""
def offsetCorrection(majList, minList, azList, elList, raList, decList, dataDMJDList, refractList, lstStartSec, coordSys):
	print('\n')
	print('Correcting coordinates for beam offset...')
	## we have refract value for a single pol
	## extend by 2x to describe both XX and YY pol
	extRefractArr = np.zeros(len(refractList)*2)
	extMajArr = np.zeros(len(refractList)*2)
	extMinArr = np.zeros(len(refractList)*2)
	extAzArr = np.zeros(len(refractList)*2)
	extElArr = np.zeros(len(refractList)*2)
	extRaArr = np.zeros(len(refractList)*2)
	extDecArr = np.zeros(len(refractList)*2)
	extDataDMJDArr = np.zeros(len(refractList)*2)
	extLstArr = np.zeros(len(refractList)*2)
	dataDMJDArr = np.zeros(len(refractList)*2)
	beamXElArr = hdu[1].data['FEEDXOFF'] 
	beamElArr = hdu[1].data['FEEDEOFF']


	extRefractArr[0::2] = refractList
	extRefractArr[1::2] = refractList

	extMajArr[0::2] = majList
	extMajArr[1::2] = majList

	extMinArr[0::2] = minList
	extMinArr[1::2] = minList

	extAzArr[0::2] = azList
	extAzArr[1::2] = azList

	extElArr[0::2] = elList
	extElArr[1::2] = elList

	extRaArr[0::2] = raList
	extRaArr[1::2] = raList
	extDecArr[0::2] = decList
	extDecArr[1::2] = decList

	extLstArr[0::2] = lstStartSec
	extLstArr[1::2] = lstStartSec

	lstStartDegs = extLstArr/3600*15 ## convert to degs

	dataDMJDArr[0::2] = dataDMJDList
	dataDMJDArr[1::2] = dataDMJDList

	azOffVal = beamXElArr/np.cos(np.deg2rad(beamElArr)) ## put Cross-El back into Azimuth
	elOffVal = beamElArr
	newAzArr = extAzArr - azOffVal ## apply azimuth offset
	newElArr = extElArr - extRefractArr - elOffVal ## apply elevation offset (and atmospheric refraction)

	"""
	we have data DMJD for a single pol 
	extend by 2x to describe both XX and YY pol
	""" 
	dataDMJDArr = np.zeros(len(dataDMJDList)*2)
	dataDMJDArr[0::2] = dataDMJDList
	dataDMJDArr[1::2] = dataDMJDList 

	newDecArr = el2dec(lstStartDegs, newAzArr, newElArr) ## convert to Declination
	newRaArr = az2ra(lstStartDegs, newAzArr, newElArr, newDecArr, coordSys) ## convert to Right Ascension
	if any(t > 360 for t in newRaArr): ## wrap if any values > 2pi (360 deg)
		newRaArr = newRaArr - 360 

	## loop through each coordinate to precess from geoapparent RA/Dec to J2000. 
	for coordIdx in range(0, len(newRaArr)):
		raVal = newRaArr[coordIdx]
		decVal = newDecArr[coordIdx]
		dmjdVal = dataDMJDArr[coordIdx]
		## convert from geocentric apparant to mean place (J2000.0) 
		newRa, newDec = pysla.slalib.sla_amp(np.deg2rad(raVal), np.deg2rad(decVal), dmjdVal, 2000.0)
		## convert radians to degrees
		newRaArr[coordIdx] = np.rad2deg(newRa)
		newDecArr[coordIdx] = np.rad2deg(newDec)
	newMajArr = np.copy(newRaArr)
	newMinArr = np.copy(newDecArr)
	if coordSys == 'GALACTIC': ## if we need GALACTIC coordinates, convert using astropy
		c = SkyCoord(ra=newMajArr * u.degree, dec = newMinArr * u.degree, frame = 'fk5')
		lVal = c.galactic.l.deg
		bVal = c.galactic.b.deg

		newMajArr = lVal
		newMinArr = bVal

	return newMajArr, newMinArr, extMajArr, extMinArr, newAzArr, newElArr, extAzArr, extElArr, newRaArr, newDecArr, extRaArr, extDecArr, lstStartDegs
"""
begin looping over files to get Antenna and data positions, DMJD values, and perfrom correct coord
"""
fileList = fileList[1:]
for flName in fileList:
	print(flName)
	antHdu = fits.open('/home/gbtdata/' + projectID + '/Antenna/' + flName) 
	dataHdu = fits.open('/lustre/projects/flag.old/' + projectID + '/BF/' + flName[:-5] + 'A.fits') ## only need one correlator bank FITS file
	goHDU = fits.open(goFitsPath + '/' + flName)

	## extract pointing model information 
	smntAz = antHdu[0].header['SMNTC_AZ']
	sobscAz = antHdu[0].header['SOBSC_AZ']
	sobscEl = antHdu[0].header['SOBSC_EL']
	smntEl = antHdu[0].header['SMNTC_EL']
	  
	## calculate pointing model
	azPt = smntAz - sobscAz
	elPt = smntEl - sobscEl


	## get DMJD values
	antDMJD = antHdu[2].data['DMJD']
	corrDMJD = dataHdu[1].data['DMJD']

	##get Antenna MAJOR and MINOR axis and refraction values for interpolation
	antMaj = antHdu[2].data['MAJOR']
	antMin = antHdu[2].data['MINOR']

	antRa = antHdu[2].data['RAJ2000']
	antDec = antHdu[2].data['DECJ2000']

	## get refraction correction values
	antRefract = antHdu[2].data['REFRACT']
	refractInterp = np.interp(corrDMJD, antDMJD, antRefract)

	## get Azimuth/Elevation and subtract off pointing model
	az = antHdu[2].data['MNT_AZ']- azPt
	el = antHdu[2].data['MNT_EL'] - elPt
	#az = antHdu[2].data['OBSC_AZ']# - azPt
	#el = antHdu[2].data['OBSC_EL']# - elPt
	#az = antHdu[2].data['MNT_AZ']- antHdu[2].data['OBSC_AZ'] - azPt
	#el = antHdu[2].data['MNT_EL']- antHdu[2].data['OBSC_EL'] - elPt
	azInterp = np.interp(corrDMJD, antDMJD, az)
	elInterp = np.interp(corrDMJD, antDMJD, el)

	lstStart = antHdu[0].header['LSTSTART']
	intLen = np.float(dataHdu[0].header['ACTSTI'])
	lstStartArr = np.zeros([len(azInterp)])
	for idx in range(0, len(azInterp)): ## calculate LST 
  		val = lstStart + (idx*intLen)
  		lstStartArr[idx] = val
	
	## interpolate antenna values to data timestamps
	corrMaj = np.interp(corrDMJD, antDMJD, antMaj)
	corrMin = np.interp(corrDMJD, antDMJD, antMin)

	## interpolate atnenna values to data timestamps (RA/Dec)
	raInterp = np.interp(corrDMJD, antDMJD, antRa)
	decInterp = np.interp(corrDMJD, antDMJD, antDec)


	## place in lists
	majList.extend(corrMaj)
	minList.extend(corrMin)
	raList.extend(raInterp)
	decList.extend(decInterp)
	azList.extend(azInterp)
	elList.extend(elInterp)
	lstStartList.extend(lstStartArr)
	dataDMJDList.extend(corrDMJD)
	refractList.extend(refractInterp)

	coordSys = goHDU[0].header['COORDSYS']


newMajArr, newMinArr, dataMajArr, dataMinArr, newAzArr, newElArr, extAzArr, extElArr, newRaArr, newDecArr, dataRaArr, dataDecArr, lstStartArr = offsetCorrection(majList, minList, azList, elList, raList, decList, dataDMJDList, refractList, lstStartList, coordSys)
## inform about percent difference
majPercDiffMaj = np.abs((newMajArr - dataMajArr)/dataMajArr*100)
minPercDiffMin = np.abs((newMinArr - dataMinArr)/dataMinArr*100)

##
print('Maximum Percent Difference in MAJOR axis:' + np.str(np.max(majPercDiffMaj)))
print('Maximum Percent Difference in MINOR axis: ' + np.str(np.max(minPercDiffMin)))

pyplot.figure()
pyplot.plot(newMinArr, label = 'Computed', linewidth = 2)
pyplot.plot(dataMinArr, label = 'Data', linewidth = 2)
pyplot.xlabel('Element')
pyplot.ylabel('Coordinate [deg]')
pyplot.title(fileName[:-5] + ' (' + observedObj + ')' + ' MINOR Coordinate', fontsize = 10)
pyplot.legend(loc=0, fontsize=14)
pyplot.savefig(fileName[:-5] + '_' + observedObj + '_MINOR_Coord.pdf')
pyplot.show()

pyplot.figure()
pyplot.plot(newMajArr, label = 'Computed', linewidth = 2)
pyplot.plot(dataMajArr, label = 'Data', linewidth = 2)
pyplot.title(fileName[:-5] + ' (' + observedObj + ')' + ' MAJOR Coordinate', fontsize = 10)
pyplot.xlabel('Element')
pyplot.ylabel('Coordinate [deg]')
pyplot.legend(loc=0, fontsize=14)
pyplot.savefig(fileName[:-5] + '_' + observedObj + '_MAJOR_Coord.pdf')
pyplot.show()

## update values for spatial coordinates
hdu[1].data['CRVAL2'] = newMajArr
hdu[1].data['CRVAL3'] = newMinArr
hdu[1].data['TRGTLONG'] = newMinArr
hdu[1].data['TRGTLAT'] = newMajArr
hdu[1].data['AZIMUTH'] = newAzArr
hdu[1].data['ELEVATIO'] = newElArr

cenFreqsArr = radVelCorrection(newRaArr, newDecArr)

## update values for spectral coordinates
hdu[1].data['CRVAL1'] = cenFreqsArr
hdu[1].data['OBSFREQ'] = cenFreqsArr
hdu[1].data['VELDEF'] = 'OPTI-HEL'

hdu.flush()


