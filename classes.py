import numpy as np
import os
import glob
import cPickle
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn, vstack, hstack

class LightCurve():
	"""
	Once fully initiated, instances of this class have the following important 
	properties
	band			(string) 'g','r','i' or 'z'
	mjd				(array) modified julian dates of observations
	flux			(array) the observed flux
	fluxErr			(array) the error in the observed flux
	shifted_mjd		(array) mjd shifted such that the peak has mjd = 0. To be 
							modified
	"""
	badCurve = False
	shifted_mjd = np.zeros(0)
	def __init__(self, band):
		self.band = band
		self.mjd = np.zeros(0)
		self.flux = np.zeros(0)
		self.fluxErr = np.zeros(0)

	@classmethod
	def data(band, mjd, flux, fluxErr):
		self.band = band
		self.mjd = mjd
		self.flux = flux
		self.fluxErr = fluxErr
		self.set_badCurve()

	def set_badCurve(self):	
		if len(self.flux) == 0:
			self.badCurve = True
		# else:
		# 	self.badCurve = False

	def addDataPoint(self, mjd, flux, fluxErr):
		"""
		Adds a data point to the light curve.
		"""
		self.mjd = np.append(self.mjd, mjd)
		self.flux = np.append(self.flux, flux)
		self.fluxErr = np.append(self.fluxErr, fluxErr)

	@property
	def make_shifted_mjd(self, distance):
		"""
		Construct shifted_mjd, by subtracting 'distance' from 'self.flux'
		"""
		self.shifted_mjd = self.mjd - distance
	
	@property
	def get_maxfluxIndex(self):
		"""
		Return the index of the maximum flux
		"""
		return np.argmax(self.flux)
	
	@property
	def get_max_fmfe_Index(self):
		"""
		Return the index of max (flux - fluxErr)
		"""
		difference = np.subtract(self.flux, self.fluxErr)
		
		return np.argmax(difference)
	
	@property
	def get_max_flux_p(self, p):
		"""
		Returns max (flux - p*fluxErr)
		"""
		return np.max(np.subtract(self.flux, p*self.fuxErr))
		



class Supernova():
	"""
	Has the following properties

	g				(LightCurve)
	r				(LightCurve)
	i				(LightCurve)
	z				(LightCurve)
	lightCurvesDict	(dictionary) 4 entries, g,r,i,z returning the corresponding 
								 LightCurves
	SNID			(int)	supernova ID
	SNTypeInt		(int)	supernova type integer (see relationship between 
													number and type)
	zSpec			(float)	If known via spectroscope, otherwise None
	hostGalaxyID	(int)	THe host galaxy ID (all supernovae (in +zPhotHost) 
												catalog have this)
	zPhotHost		(float)	The redshift of the host galaxy (all supernovae in 
														the catalog have this)
	zPhotHostErr	(float)	Error in zPhotHost
	"""

	def __init__(self, inFileName):
		"""
		Parses all the light curve data in inFileName into a Supernova object.
		"""

		inFile = file(inFileName, "r")
		lines = inFile.readlines()
		inFile.close()

		self.g = LightCurve("g")
		self.r = LightCurve("r")
		self.i = LightCurve("i")
		self.z = LightCurve("z")


		self.lightCurvesDict = {'g':self.g, 
								'r':self.r, 
								'i':self.i, 
								'z':self.z}
		
		for line in lines:
			if len(line) > 3 and line[0] != "#":
				
				tag = line.split(":")[0]
				data = line.split(":")[-1].split()
				
				if tag == "OBS":
					mjd = float(data[0])
					passband = data[1]
					flux = float(data[3])
					fluxErr = float(data[4])
					if fluxErr > 0:
						if passband == "g":
							self.g.addDataPoint(mjd, flux, fluxErr)
						elif passband == "r":
							self.r.addDataPoint(mjd, flux, fluxErr)
						elif passband == "i":
							self.i.addDataPoint(mjd, flux, fluxErr)
						elif passband == "z":
							self.z.addDataPoint(mjd, flux, fluxErr)
						else:
							print "Filter not recognized: {:<5}".format(passband)
				elif tag == "SNID":
					self.SNID = int(data[0])
				elif tag == "SNTYPE":
					self.SNTypeInt = int(data[0])                    
				elif tag == "RA":
					self.RADeg = float(data[0])
				elif tag == "DECL":
					self.decDeg = float(data[0])
				elif tag == "MWEBV":
					self.MWEBV = float(data[0])
				elif tag == "REDSHIFT_SPEC":
					if float(data[0]) == -9:
						self.zSpec = None
					else:
						self.zSpec = float(data[0])
				elif tag == "HOST_GALAXY_GALID":
					self.hostGalaxyID = int(data[0])
				elif tag == "HOST_GALAXY_PHOTO-Z":
					self.zPhotHost = float(data[0])
					self.zPhotHostErr = float(data[2])

		for b in self.lightCurvesDict.keys():
			self.lightCurvesDict[b].set_badCurve()

	def __cmp__(self, other):
		return 2*(self.zPhotHost - other.zPhotHost > 0) - 1 


class SupernovaFit():
	def __init__(self, SNID, 
		SNType=None, RADeg=None, decDeg=None, MWEBV=None, zSpec=None,
		hostGalaxyID=None, zPhotHost=None, zPhotHostErr=None):
		self.g = LightCurve("g")
		self.r = LightCurve("r")
		self.i = LightCurve("i")
		self.z = LightCurve("z")
		self.lightCurvesDict = {"g":self.g, 
								"r":self.r, 
								"i":self.i, 
								"z":self.z}
		self.SNID = SNID
		self.SNType = SNType
		self.SNType = SNType 
		self.RADeg = RADeg 
		self.decDeg = decDeg 
		self.MWEBV = MWEBV 
		self.zSpec = zSpec 
		self.hostGalaxyID = hostGalaxyID 
		self.zPhotHost = zPhotHost 
		self.zPhotHostErr = zPhotHostErr 

	def setLightCurve(self, band, mjd, flux, fluxErr):
		self.lightCurvesDict[band].mjd = mjd
		self.lightCurvesDict[band].flux = flux
		self.lightCurvesDict[band].fluxErr = fluxErr
		self.lightCurvesDict[band].set_badCurve()

	def setLCZeroPoints(self):
		try:
			mjd_rMax = self.r.mjd[self.r.flux == self.r.flux.max()][0]
			idx_rMax = np.where(self.r.flux == self.r.flux.max())[0][0]
			for b in self.lightCurvesDict.keys():
				if not self.lightCurvesDict[b].badCurve:
					self.lightCurvesDict[b].mjd -= mjd_rMax
		except:
			print "Problem detected!"

	def normalizeLC(self):
		for b in self.lightCurvesDict.keys():
			if not self.lightCurvesDict[b].badCurve:
				self.lightCurvesDict[b].flux /= self.lightCurvesDict[b].flux.max()

	def save_on_txt(self, fileName):
		t = Table(masked=True)
		colNames = [["MJD_r_band", "{0:5d}"],
					["FLUX_g", "{0:10.5f}"], ["FLUX_ERR_g", "{0:10.5f}"],
					["FLUX_r", "{0:10.5f}"], ["FLUX_ERR_r", "{0:10.5f}"],
					["FLUX_i", "{0:10.5f}"], ["FLUX_ERR_i", "{0:10.5f}"],
					["FLUX_z", "{0:10.5f}"], ["FLUX_ERR_z", "{0:10.5f}"]]

		t["MJD_r_band"] = self.r.mjd
		for b in self.lightCurvesDict.keys():
			colNameFlux = "FLUX_{:<1}".format(b)
			colNameErr = "FLUX_ERR_{:<1}".format(b)
			if self.lightCurvesDict[b].badCurve:
				colFlux = MaskedColumn(np.zeros(self.r.mjd.size),
					name=colNameFlux, 
					dtype=np.float, fill_value=-9,
					mask=np.ones(self.r.mjd.size))
				colErr = MaskedColumn(np.zeros(self.r.mjd.size), 
					name=colNameErr,
					dtype=np.float, fill_value=-9,
					mask=np.ones(self.r.mjd.size))
			else:
				colFlux = MaskedColumn(self.lightCurvesDict[b].flux,
					name=colNameFlux,
					dtype=np.float32, fill_value=-9)
				colErr = MaskedColumn(self.lightCurvesDict[b].fluxErr,
					name=colNameErr,
					dtype=np.float32, fill_value=-9)
			t.add_columns([colFlux, colErr])

		t.filled()
		fOut = open(fileName, 'w')
		fOut.write("# SNID: {:>10d}\n".format(self.SNID))
		ascii.write(t, output=fOut, delimiter='  ', 
			format='fixed_width_two_line')
		# t.write(output=fOut, delimiter='  ', format='ascii.fixed_width_two_line')
		fOut.close()
		# np.savetxt(fileName, (self.r.mjd,
		# 	self.g.flux, self.g.fluxErr,
		# 	self.r.flux, self.r.fluxErr,
		# 	self.i.flux, self.i.fluxErr,
		# 	self.z.flux, self.z.fluxErr), 
		# 	fmt="%5d, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f, \
		# 	%10.5f, %10.5f",
		# 	newline="\n", header="SNID: {:>10d}".format(self.SNID))


class SupernovaeCatalog():
	"""
	Class variables are
	sne		(object array)	array of Supernova objects
	zSpec	(float array) the spectroscopically observed redshift 
			(None if no spctrscp) 
	zPhotHost	(float array)	the redshift of the host galaxy
	SNID	(int array)	the IDs of the Supernovae
	SNType	(int array)	The types of the supernovae
	"""
	def __init__(self, dataDir, load_all):
		"""
		This class loads in all the light curve data under dataDir into a big 
		list, and creates a series of top level arrays that we can use to cut 
		the catalog by z, type etc. load_all: Do you want to
		load all of the data, or will your computer then crash?
		"""
		
		print ">>> Loading data from text files ..."
		inFileNames = glob.glob(dataDir+os.path.sep+"DES_SN*.DAT")
		self.sne = np.zeros(0, dtype=np.object)
		self.zPhotHost = np.zeros(0, dtype=np.float32)
		self.SNID = np.zeros(0, dtype=np.int)
		self.SNType = np.zeros(0, dtype=np.int)
		count = 0
		for inFileName in inFileNames:
			count += 1
			# Progress update
			tenPercent = len(inFileNames)/10
			for j in range(0,11):
				if count == j*tenPercent:
					print "... "+str(j*10)+"% complete ..."
			
			
			inFile = file(inFileName, "r")	
			for i in range(4):
				line = inFile.readline()
			inFile.close()
			
			SNTYPE = line.split(":")[-1].split()[0]
			if SNTYPE == '-9' and load_all == False:
				pass		

			else:	
				sn = Supernova(inFileName)            
				self.SNType = np.append(self.SNType, sn.SNTypeInt)
				self.SNID   = np.append(self.SNID, sn.SNID)
				self.sne    = np.append(self.sne, sn)

		self.zPhotHost = np.nan_to_num(self.zPhotHost)
			
	def findSupernovae(self, SNType, zSpecLow, zSpecHigh):
		"""
		Given a SNType code and a redshift range, return a list of matching 
		supernovae in the catalog.		
		"""
		typeMask = np.equal(self.SNType, SNType)
		zSpecMask = np.logical_and(np.greater(self.zSpec, zSpecLow), 
					   			   np.less(self.zSpec, zSpecHigh))
		mask = np.logical_and(typeMask, zSpecMask)
		foundSupernovae = self.sne[mask]
		
		return foundSupernovae

	def getSNTypeStr(self, SNTypeInt):
		"""Given a SNTypeInt, returns a string, e.g. 'Ia'
		
		Mapping to type names (from DES_BLIND+HOSTZ.README):
		
			1  (Ia)
			2  (II)
			3  (Ib/c)
			11  (pec. Ia)
			66  (other)
			-1  (rejected)
			-9 (unknown)
			
		"""
		
		if SNTypeInt == 1:
			SNTypeStr="Ia"
		elif SNTypeInt == 2:
			SNTypeStr="II"
		elif SNTypeInt == 3:
			SNTypeStr="Ib/c"
		elif SNTypeInt == 11:
			SNTypeStr="pec. Ia"
		elif SNTypeInt == 66:
			SNTypeStr="other"
		elif SNTypeInt == -1:
			SNTypeStr="rejected"
		elif SNTypeInt == -9:
			SNTypeStr="unclassified"
		
		return SNTypeStr

	def getSNTypeInt(self, SNTypeStr):
		"""Given a SNTypeStr, returns the corresponding int, e.g. 1
		
		Mapping to type names (from DES_BLIND+HOSTZ.README):
		
			1  (Ia)
			2  (II)
			3  (Ib/c)
			11  (pec. Ia)
			66  (other)
			-1  (rejected)
			-9 (unknown)
			
		"""
			
		if SNTypeStr == "Ia":
			SNTypeInt = 1
		elif SNTypeStr == "II":
			SNTypeInt = 2
		elif SNTypeStr == "Ib/c":
			SNTypeInt = 3
		elif SNTypeStr == "pec. Ia":
			SNTypeInt = 11
		elif SNTypeStr == "other":
			SNTypeInt = 66
		elif SNTypeStr == "rejected":
			SNTypeInt = -1
		elif SNTypeStr == "unclassified":
			SNTypeInt = -9
		
		return SNTypeInt