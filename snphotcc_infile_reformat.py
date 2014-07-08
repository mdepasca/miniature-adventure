'''
MDEPACSA 31/01/2014 created
-----------------------------------------------------------------------------
Reformats ASCII files from SNPhotCC (Kessler+ 2010) into ASCII files
readable by Faraway+ 2013 R code

SNPhotCC ASCII files format
----------------------------------------------- 
SNID:   41                                     <<< to R readable file       1  
IAUC:    UNKNOWN                                                            2
SNTYPE:  -9                                    <<< to R readable file       3
FILTERS: griz                                                               4
RA:      34.500000  deg                        <<< to R readable file       5
DECL:    -5.500000  deg                        <<< to R readable file       6
FAKE:    3   (=> BLIND-TEST simulation)                                     7
MWEBV:   0.0227    MW E(B-V)                                                8
REDSHIFT_SPEC:    -9.0000 +-  9.0000                                        9
                                                                           10
HOST_GALAXY_GALID:   19513                                                 11
HOST_GALAXY_PHOTO-Z:   0.8642  +- 0.0171       <<< needed for K-correction 12
                                                                           13 
                                                                           14
                                                                           15
                                                                           16
# ============================================                             17
# TERSE LIGHT CURVE OUTPUT:                                                18
#                                                                          19
NOBS: 41                                                                   20
NVAR: 5                                                                    21
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR                            22
OBS:  56267.156  i NULL   2.217e+00   2.211e+00 <<< to R readable file    ...

-----------------------------------------------

R readable file (1 per filter) has to contain
SNID and SNTYPE
RA
DECL
MJD FLUXCAL FLUXCALERR

-----------------------------------------------

The script uses AstrPy table to read in a file containing
light curves file names.

??? The same package is used to write the R readable files

Constants have a _ before name

Rregarding flux from SNPhotCC (Kessler+ 2010)
"...The sky-noise, point-spread
function and atmospheric transparency are evaluated in
each filter based on years of observational data from the
ESSENCE project at the Cerro Tololo Inter-American
Observatory (CTIO)..."
'''

'''
PACKAGES VER. 
NumPy    1.7
astropy  0.3 
'''

from astropy.io import ascii
import numpy as np
import numpy.ma as ma                     # Masked arrays
import matplotlib.pyplot as plt
'''
Converting flux arrays into mgnitude arrays.
It considers negative fluxes as magnitude limited measures.
'''
def flux_to_mag(flux, limFlux, isFlux=True):
    # flux = 10^(-0.4 * m + 11) => m = -(log10(flux) - 11) / 0.4
    
    # if statement to convert limiting flux into limiting magnitude
    if isFlux:
        if limFlux <= 0:
            print "Limiting flux was <= 0"
        limMag = -2.5 * (-11 + np.log10(limFlux))
        intLimFlux = limFlux
    else:
        limMag = limFlux
        intLimFlux = np.power(10, -0.4 * limMag + 11)

    # applying the mask to detection below the limiting flux
    maFlux = ma.masked_where(flux < intLimFlux, flux)

    fluxMask = maFlux.mask                # to avoid warnings due to values passed to np.log10
    maMag = -2.5 * (-11.0 + np.log10(ma.filled(maFlux,1)))
    maMag = ma.array(maMag, mask = fluxMask)
    mag = ma.filled(maMag, limMag)
    return mag

def error_to_mag(flux, fluxErr):
    # error prop. zErr = abs(dF(x)/dx) xErr
    magErr = np.multiply(np.abs(np.divide(-2.5, np.log(10) * flux)), np.abs(fluxErr)) 
    return magErr

def main():
    inDir = "../DES_BLIND+HOSTZ"
    outDir = inDir + "_Rfiles"
    fNameFile = "DES_BLIND+HOSTZ.LIST"
    
    # Reading list on light curves files
    fName = ascii.read (inDir + '/' + fNameFile, format = "no_header")
    _inFileIdx = 0
    
    # DES limiting magnitudes (Bernstein+ 2012 ApJ.753)
    gLimMag = 25.2
    rLimMag = 25.4
    iLimMag = 25.1
    zLimMag = 24.9

    # DES limiting fluxes
    gLimFlux = np.power(10, -0.4 * gLimMag + 11)
    rLimFlux = np.power(10, -0.4 * rLimMag + 11)
    iLimFlux = np.power(10, -0.4 * iLimMag + 11)
    zLimFlux = np.power(10, -0.4 * zLimMag + 11)

    count = 0
    tenPercent = len(fName) * 4 / 10

    print "File creation started ..."
    for i in range(len(fName) - 1):
        # Reading light curve file
        inFile = file(inDir + '/' + fName[i][_inFileIdx], 'r')
        lines = inFile.readlines()
        inFile.close()
    
        '''
        Observations per band are made in differents MJDs
        One vector of observations per band
        One file of observation per band, since the R structure
        is per band
        '''
        
        gLightCurve = np.zeros(0)
        rLightCurve = np.zeros(0)
        iLightCurve = np.zeros(0)
        zLightCurve = np.zeros(0)
        
        gLightCurveErr = np.zeros(0)
        rLightCurveErr = np.zeros(0)
        iLightCurveErr = np.zeros(0)
        zLightCurveErr = np.zeros(0)    

        gMjd = np.zeros(0)
        rMjd = np.zeros(0)
        iMjd = np.zeros(0)
        zMjd = np.zeros(0)

        snId = ''
        snType = 0
        snRa = 0.0
        snDecl= 0.0
        snZ = 0.0
        snZErr = 0.0

        hostZ = 0.0
        hostZErr = 0.0

        ebvMW = 0.0
        for line in lines:
            if len(line) > 3 and line[0] != '#':
                lineTag = line.split(':')[0]
                lineData = line.split(':')[-1].split()

                if lineTag == "SNID":
                    snId = lineData[0]

                if lineTag == "SNTYPE":
                    snType = float(lineData[0])

                if lineTag == "RA":
                    snRa = float(lineData[0])

                if lineTag == "DECL":
                    snDecl = float(lineData[0])

                if lineTag == "REDSHIFT_SPEC":
                    if float(lineData[0]) == -9:
                        snZ    = -1
                        snZErr = -1
                    else:
                        snZ    = float(lineData[0])
                        snZErr = float(lineData[2])

                if lineTag == "HOST_GALAXY_PHOTO-Z":
                    hostZ    = float(lineData[0])
                    hostZErr = float(lineData[2])

                if lineTag == "MWEBV":
                    ebvMW = float(lineData[0])
                
                # Filling observation arrays
                if lineTag == "OBS":
                    if lineData[1] == 'g':
                        gMjd           = np.append(gMjd, float(lineData[0]))
                        gLightCurve    = np.append(gLightCurve, float(lineData[3]))
                        gLightCurveErr = np.append(gLightCurveErr, float(lineData[4]))
                    elif lineData[1] == 'r':
                        rMjd           = np.append(rMjd, float(lineData[0]))
                        rLightCurve    = np.append(rLightCurve, float(lineData[3]))
                        rLightCurveErr = np.append(rLightCurveErr, float(lineData[4]))
                    elif lineData[1] == 'i':
                        iMjd           = np.append(iMjd, float(lineData[0]))
                        iLightCurve    = np.append(iLightCurve, float(lineData[3]))
                        iLightCurveErr = np.append(iLightCurveErr, float(lineData[4]))
                    elif lineData[1] == 'z':
                        zMjd           = np.append(zMjd, float(lineData[0]))
                        zLightCurve    = np.append(zLightCurve, float(lineData[3]))
                        zLightCurveErr = np.append(zLightCurveErr, float(lineData[4]))
                    else:
                        raise IOError('Wrong filter!')
                

            # end of loop on file's lines
        gMag = flux_to_mag(gLightCurve, gLimMag, False)
        gMagErr = error_to_mag(gLightCurve, gLightCurveErr)
        rMag = flux_to_mag(rLightCurve, rLimMag, False)
        rMagErr = error_to_mag(rLightCurve, rLightCurveErr)
        iMag = flux_to_mag(iLightCurve, iLimMag, False)
        iMagErr = error_to_mag(iLightCurve, iLightCurveErr)
        zMag = flux_to_mag(zLightCurve, zLimMag, False)
        zMagErr = error_to_mag(zLightCurve, zLightCurveErr)
        # wrinting on new file
        if False:
            outName = 'DES_BLIND+HOSTZ.SUMMARY.Rdat'
            try:
                with open(outDir + '/' + outName):
                    fileFlag = 'a'
            except IOError:
                fileFlag = 'w'
            outFile = file(outDir + '/' + outName, fileFlag) 
            sep = '  '

            # writing header
            if fileFlag == 'w':
                header = "%6s" % "SNID" + sep + \
                    "%4s" % "TYPE" + sep + \
                    "%7s" % "SUBTYPE" + sep + \
                    "%8s" % "SN_Z" + sep + \
                    "%8s" % "SN_Z_ERR" + sep + \
                    "%8s" % "HOST_Z" + sep + \
                    "%10s" % "HOST_Z_ERR" + sep + \
                    "%6s" % "MW_EBV" + sep + \
                    "%5s" % "gNOBS" + sep + \
                    "%5s" % "rNOBS" + sep + \
                    "%5s" % "iNOBS" + sep + \
                    "%5s" % "zNOBS" + '\n'
            
                outFile.write(header)
                # header = "# %6s" % "------" + sep + \
                    #   "%4s" % "----" + sep + \
                    #   "%7s" % "-------" + sep + \
                    #   "%8s" % "--------" + sep + \
                    #   "%8s" % "--------" + sep + \
                    #   "%8s" % "--------" + sep + \
                    #   "%10s" % "----------" + sep + \
                    #   "%6s" % "------" + sep + \
                    #   "%5s" % "----" + sep + \
                    #   "%5s" % "----" + sep + \
                    #   "%5s" % "----" + sep + \
                    #   "%5s" % "----" + '\n'
                    # outFile.write(header)

                outLine = "  %6s" % str(snId) + sep + \
                    "%4s" % "SN" + sep + \
                    "%7i" % snType + sep + \
                    "%8.4f" % snZ + sep + \
                    "%8.4f" % snZErr + sep + \
                    "%8.4f" % hostZ + sep + \
                    "%10.4f" % hostZErr + sep + \
                    "%6.4f" % ebvMW + sep + \
                    "%5i" % (len(gMjd)) + sep + \
                    "%5i" % (len(rMjd)) + sep + \
                    "%5i" % (len(iMjd)) + sep + \
                    "%5i" % (len(zMjd)) + '\n'

            outFile.write(outLine)
            outLine = ''
            header = ''
            outFile.close()                       # line in summary written
    
        # Writing on SN specific file (one file per band)
        for b in 'griz':
            outName = fName[i][_inFileIdx].split('.')[0] + '_' + b + '.Rdat'
            fileFlag = 'w'
            outFile = file(outDir + '/' + outName, fileFlag)
            
            sep = '  '
        
            header = "%14s" % "RA" + sep + \
                "%14s" % "DECL" + sep + \
                "%9s" % "MJD" + sep + \
                "%7s" % "MAG" + sep + \
                "%7s" % "MAG_ERR" + sep + \
                "%10s" % "FLUX" + sep + \
                "%10s" % "FLUX_ERR" + '\n'
          
            outFile.write(header)
            # header = "# %14s" % "--------------" + sep + \
                #   "%14s" % "--------------" + sep + \
                #   "%9s" % "---------" + sep + \
                #   "%7s" % "-------" + sep + \
                #   "%7s" % "-------" + sep + \
                #   "%10s" % "----------" + sep + \
                #   "%10s" % "----------" + '\n'
            
            # outFile.write(header)

            if b == 'g':
                jMax = len(gMjd) - 1
                mjd = np.copy(gMjd)
                flux = np.copy(gLightCurve)
                fluxErr = np.copy(gLightCurveErr)
                mag = np.copy(gMag)
                magErr = np.copy(gMagErr)
            elif b == 'r':
                jMax = len(rMjd) - 1
                mjd = np.copy(rMjd)
                flux = np.copy(rLightCurve)
                fluxErr = np.copy(rLightCurveErr)
                mag = np.copy(rMag)
                magErr = np.copy(rMagErr)
            elif b == 'i':
                jMax = len(iMjd) - 1
                mjd = np.copy(iMjd)
                flux = np.copy(iLightCurve)
                fluxErr = np.copy(iLightCurveErr)
                mag = np.copy(iMag)
                magErr = np.copy(iMagErr)
            elif b == 'z':
                jMax = len(zMjd) - 1
                mjd = np.copy(zMjd)
                flux = np.copy(zLightCurve)
                fluxErr = np.copy(zLightCurveErr)
                mag = np.copy(zMag)
                magErr = np.copy(zMagErr)
            
            for j in range(jMax):
                outLine = "  %+14.9f" % snRa + sep + \
                    "%+14.9f" % snDecl + sep + \
                    "%9.3f" % mjd[j] + sep + \
                    "%7.4f" % mag[j] + sep + \
                    "%7.4f" % magErr[j] + sep + \
                    "%+5.3e" % flux[j] + sep + \
                    "%+5.3e" % fluxErr[j] + '\n'

                outFile.write(outLine)
                outLine = ''
            outFile.close()                   # specific SN band file written
        
            count = count + 1
            # Progress update
        
            for k in range(0, 11):
                if count == k * tenPercent:
                    print "... " + str(k * 10) + "% complete ..."	
 
        
    print "Done!"
