import argparse
import os
from  os import path
import subprocess
import sys
import socket
import time
import warnings
from math import floor
import gc # garbage collector
import smtplib
import numpy as np
from scipy import signal, linalg
from matplotlib import pyplot as plt
import GPy

import classes as cls
import utilities as util
from utilities import bcolors

# import rpy2.robjects as ro
# from rpy2.robjects.packages import importr
# from rpy2.robjects.numpy2ri import numpy2ri
# # Activate automatic conversion of ndarray to R objects
# ro.conversion.py2ri = numpy2ri

from progressbar import ProgressBar, SimpleProgress, ETA, Percentage, Bar, \
                        AnimatedMarker, Timer, Counter


if __name__ == "__main__":
    # gc.set_debug(gc.DEBUG_LEAK)
    # Parsing input from command line
    parser = argparse.ArgumentParser(
        description = "SN lightcurve fitter and classifier.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    actionGroup = parser.add_argument_group('ACTION')

    inputGroup = parser.add_argument_group('INPUT')

    """

    ACTION OPTIONS
    ----------------------------------------------------------------------------
    """

    actionGroup.add_argument(
        "--fit", dest="fit",
        action="store_true",
        help="Fit lightcurves with Gaussian processes method."
        )

    actionGroup.add_argument(
        '--prior', dest='prior',
        action='store_true', help='Use priors in GP regression.'
        )

    actionGroup.add_argument(
        '--length', dest='testLength',
        action='store_true',
        help='Set length scale hyper parameter to random value to ease \
            optimization.'
        )

    actionGroup.add_argument(
        "--cross-correlation", dest="crossCor",
        action="store_true",
        help="Performs cross correlation between non peaked lcs (with maximum in \
            r-band at one of the MJD extremes) and all the peaked lcs. Produces \
            an estimate for maximum in r-band. VERY TIME CONSUMING."
        )

    actionGroup.add_argument(
        "--distance-matrix", dest="distMatrix",
        action="store_true",
        help="Calculate distance between fitted lightcurves in same band. \
        It is use to build a diffusion map (see Coifman & Lafon (2006) \
        and Lafon & Lee (2006)).")

    actionGroup.add_argument(
        "--diffuse", dest="diffuse",
        action="store_true",
        help="Computes the diffusion map coefficients. Run together or after \
        --distance-matrix option. Uses `diffusionMap` R package developed \
        by Joseph Richards.")

    actionGroup.add_argument(
        "--train", dest="train",
        action="store_true",
        help="Train the classifier - Random Forest. Uses `randomForest` R \
        package.")

    actionGroup.add_argument(
        "--classify", dest="classify",
        action="store_true")

    actionGroup.add_argument(
        "--plot", dest="plot",
        action="store_true",
        help="Save on `pdf` file the plot of fitting curve over data.")

    actionGroup.add_argument(
        '--nice-plots', dest='nicePlots',
        action='store_true',
        help='Produces plot suitable for publication (pdf, 300dpi).'
        )

    """-------------------------------------------------------------------------

    INPUT OPTIONS
    ----------------------------------------------------------------------------
    """

    inputGroup.add_argument(
        "--data-directory", dest="dirData",
        default="train_data" + os.sep + "SIMGEN_PUBLIC_DES",
        help="Path to directory containing training data.")

    inputGroup.add_argument(
        "--fit-directory", dest="dirFit",
        default="results" + os.sep + "FIT",
        help="Path to directory containing fitted data.")

    # the use of this keyword is developed in dev_magnitudes branch
    inputGroup.add_argument(
        "--mag", dest="mag",
        action="store_true",
        help="Reads in magnitudes from file."
        )

    inputGroup.add_argument(
        "--fit-file", dest="fitFile",
        help="Path to file in which to dump fitting results.")

    inputGroup.add_argument(
        "-f", "--file",
        help="")

    inputGroup.add_argument(
        "-c", "--candidate", dest="cand",
        default=-1, type=int,
        help="ID of a candidate."
        )

    inputGroup.add_argument(
        "--all-bands", dest="allBands",
        action="store_true",
        help="Plot all bands --nice-plots option."
        )

    inputGroup.add_argument(
        "-b", "--band", dest="band", default='r',
        help="Which band to plot with --nice-plots.")

    inputGroup.add_argument(
        "--nBands", dest="nBands",
        default=-1, type=int,
        help="Number of bands to plot with --nice-plots.")

    inputGroup.add_argument(
        '--limits', nargs=2, dest='limits',
        default=[0, 5], type=int,
        help='Starting ending indeces for fitting and cross-correlation.'
        )

    inputGroup.add_argument(
        '--offset', '-o', dest='offset',
        default=0, type=int,
        help='Offset for columns WRT limits (which are referred to rows).'
    )

    inputGroup.add_argument(
        '--plot-offset', dest='plotOffset',
        default=-1, type=int,
        help='Offset in index to begin light curves plotting from.'
    )
    """-------------------------------------------------------------------------
    """

    args = parser.parse_args()

    bands = ['g', 'r', 'i', 'z']
else:
    pass

if __name__ == "__main__":
    # os.system("clear")
    fromAddress = 'mothra@oapd.inaf.it'
    toAddress = 'marco.depa@gmail.com'
    sent = False

    indent = "          "
    resDir = "results"+os.sep

    peakIdx = np.empty(0)
    nopeakIdx = np.empty(0)

    print bcolors.bldpur
    print indent + "* * * * * * * * * * * * * * *"
    print indent + "*    Miniature Adventure    *"
    print indent + "*    -------------------    *"
    print indent + "*    lightcurves fitting    *"
    print indent + "*             and           *"
    print indent + "*      SN classification    *"
    print indent + "* * * * * * * * * * * * * * *"
    print bcolors.txtrst

    if args.dirFit == 'results/FIT':
        yesno = str(raw_input(indent + 'Set fit directory other then default (' + \
            parser.get_default('dirFit') + ')? (y/n)'))
        if yesno == 'y':
            args.dirFit = str(raw_input(indent + 'Specify new directory '\
                +'for fit: '))

    if args.dirData[-1] != os.sep:
        args.dirData += os.sep

    if args.dirFit[-1] != os.sep:
        args.dirFit += os.sep

    print indent + 'Fit directory will be: ' + path.abspath(args.dirFit)

    if not os.path.exists(path.abspath(args.dirFit)):
        os.makedirs(path.abspath(args.dirFit))

    start_time = time.time()

    """
    Get list of files in data directory and fit directory
    ----------------------------------------------------------------------------
    """
    p = subprocess.Popen("ls *SN*.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirData)
    lsDirData = p.stdout.read()
    lsDirData = lsDirData.split('\n')
    lsDirData.sort()
    lsDirData.remove('')

    p = subprocess.Popen("ls *SN*.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirFit)
    lsDirFit = p.stdout.read()
    lsDirFit = lsDirFit.split('\n')
    lsDirFit.sort()
    lsDirFit.remove('')

    """-------------------------------------------------------------------------
    """

    """

    PERFORMS LCs FITTING

    """
    if args.fit:
        if args.limits[1] > len(lsDirData):
            print indent + \
                "WARNING: upper limit > than the number of files. Corrected.\n"
            args.limits[1] = len(lsDirData)
        filePath = args.dirFit + 'PEAKED_{:<}_{:<5.3f}.LIST'.format(
            socket.gethostname(), time.time()
            )

        fPeaked = open(filePath, 'w')


        filePath = args.dirFit + 'NOPEAKED_{:<}_{:<5.3f}.LIST'.format(
            socket.gethostname(), time.time()
            )

        fNopeaked = open(filePath, 'w')

        # Relevant input data
        print "\n" + indent + "[1] * Fit lightcurves ..."

        print "\n" + indent + "Index interval [{:<},{:<})".format(
                args.limits[0], args.limits[1]
                )

        print "\n" + indent + \
            "Data directory: " + os.curdir + args.dirData

        print "\n" + indent \
            + "Number of candidates = {:<d}".format(len(lsDirData))

        """


        GP kernel specification
        ------------------------------------------------------------------------

        """
        # kern = GPy.kern.RatQuad(1)
        kern = GPy.kern.RBF(1)
        # kern = GPy.kern.Matern32(1)
        # kern = GPy.kern.Matern52(1)

        """---------------------------------------------------------------------
        """

        print "\n" + indent \
            + "Data will be smoothed using GP kernel " + kern.name.upper()

        print '\n' + indent + \
                    "INDEX | SN ID | BAND"
        for i in range(args.limits[0], args.limits[1]):
            filePath = path.splitext(lsDirData[i])[0] + "_FIT.DAT"

            """
            Check if file with fit results already exits. If positive skip
            to next loop iteration.
            """
            if filePath in lsDirFit:
                continue

            candidate = util.get_sn_from_file(
                args.dirData + lsDirData[i],
                args.mag
                )

            # Creating SupernovaFit object
            candidateFit = cls.SupernovaFit(candidate, kern.name)

            for b in candidate.lcsDict.keys():
                # Correcting for time dilution
                epoch = util.time_correct(
                    candidate.lcsDict[b].mjd,
                    candidate.zSpec if candidate.zSpec else candidate.zPhotHost
                    )

                # Correcting for absorption
                flux = util.correct_for_absorption(
                    candidate.lcsDict[b].flux,
                    candidate.MWEBV, b
                    )

                errFlux = candidate.lcsDict[b].fluxErr

                if (candidate.lcsDict[b].badCurve) or (len(flux) <= 3):
                    candidateFit.lcsDict[b].badCurve = True
                    print indent + bcolors.FAIL + \
                        "{:<}   {:<}   {:<} Bad Curve".format(i, candidate.SNID, b) + \
                        bcolors.txtrst
                    """
                    >>> if 'break' instead of 'continue' the candidate would not be
                    >>> processed and the further code would be easier (no double
                    >>> checks both on data and fit).
                    """
                    continue

                """

                Fitting Lightcurve
                ----------------------------------------------------------------

                """
                try:
                    predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                    epoch, flux, errFlux,
                                                    kern, n_restarts=10,
                                                    parallel=False,
                                                    test_length=args.testLength,
                                                    test_prior=args.prior)
                except linalg.LinAlgError as e:
                    if sent == False:
                        server = smtplib.SMTP('mailauth.oapd.inaf.it',587)
                        server.starttls()
                        server.login('marco.depascale', 'M@p3d_8$')
                        msg = 'Subject: LinAlgError\n\n' + \
                            'index = {:<d}, SNID = {:<d}'.format(i, candidate.SNID)
                        server.sendmail(fromAddress, toAddress, msg)
                        server.close()
                        sent = True
                    """
                    if LinAlgError light curve won't be saved.
                    """
                    print indent + \
                        "{:>5d}   {:>5d}   {:>4s}  >  FAIL".format(
                                                        i, candidate.SNID, b
                            ) + bcolors.FAIL + ' LinAlgError' + bcolors.txtrst
                    candidateFit.r.badCurve = True
                    raise ValueError(
                        'LinAlgError from GPy. Mail sent to {:s}'.format(
                                                                    toAddress
                            )
                        )
                else:
                    candidateFit.set_lightcurve(b, predMjd, predFlux, predErr)

                    print indent + bcolors.OKGREEN + \
                        "{:>5d}   {:>5d}   {:>4s}  >  DONE".format(
                                                        i, candidate.SNID, b
                            ) + bcolors.txtrst
                """-------------------------------------------------------------
                """
            else:
                """
                Saving fit results on file
                ----------------------------------------------------------------
                """
                if (candidateFit.r.badCurve == False):
                    filePath = args.dirFit + \
                        path.splitext(lsDirData[i])[0] + "_FIT.DAT"

                    candidateFit.save_on_txt(filePath)
                    print indent + 'file saved!'

                    if candidateFit.peaked:
                        peakIdx = np.append(peakIdx, i)
                        fPeaked.write('{:<}\n'.format(filePath))
                    else:
                        nopeakIdx = np.append(nopeakIdx, i)
                        fNopeaked.write('{:<}\n'.format(filePath))
                """-------------------------------------------------------------
                """
                gc.collect()
            # free memory
            gc.collect()

        fPeaked.close()
        fNopeaked.close()

        filePath = 'peaked_{:<}_{:<5.3f}.dat'.format(
            socket.gethostname(), time.time()
            )

        np.savetxt(args.dirFit + filePath, peakIdx,
            header='Indexes of fitted LCs with r maximum.', fmt='%d')


        filePath = args.dirFit + 'nopeaked_{:<}_{:<5.3f}.dat'.format(
            socket.gethostname(), time.time()
            )

        np.savetxt(filePath, nopeakIdx,
            header='Indexes of fitted LCs without an r maximum.', fmt='%d')

        gc.collect()

    """#########################################################################
    ############################################################################
    PERFORMING CROSS-CORRELATION
    ############################################################################
    ############################################################################
    """


    if args.crossCor:
        """
        File are sorted by SNID.
        In the following peakIdx and nopeakIdx contain index referring to the
        full list of files. For this reason the list of files it is queried on
        dirData. It is then filtered using the above variables.
        """

        print "\n" + indent + bcolors.undwht + \
            "(*) Calculate cross-correlation of not peaked- with " + \
            "peaked-lcs ..." +  bcolors.txtrst

        print "\n" + indent + "Interval [{:<},{:<})".format(args.limits[0], args.limits[1])

        filePath = args.dirFit +  'PEAKED.LIST'
        if path.exists(filePath) == False:
            # create the file concatenating existing partial files
            print '{:<s} created!'.format(filePath)
            peakedFileList = util.list_files(args.dirFit+'PEAKED*.LIST')
            util.concat_files(peakedFileList, filePath)
        peakList = np.loadtxt(filePath, dtype=np.str)

        filePath = args.dirFit + 'NOPEAKED.LIST'
        if path.exists(filePath) == False:
            # create the file from existing partial files
            print '{:<s} created!'.format(filePath)
            noPeakedFileList = util.list_files(args.dirFit+'NOPEAKED*.LIST')
            util.concat_files(noPeakedFileList, filePath)
        tmp = np.loadtxt(filePath, dtype=np.str)
        if tmp.size == 1:
            nopeakList = np.asarray([tmp])
        else:
            nopeakList = np.asarray(tmp)

        if args.limits[1] > len(nopeakList):
            args.limits[1] = len(nopeakList)
        #
        # filePath = 'repeats.txt'
        # repeats = np.loadtxt(args.dirFit + filePath, dtype=np.str)

        filePath = 'cross_correlated_files_{:<5.3f}.dat'.format(time.time())
        reWrite = open(args.dirFit + filePath, 'w')
        prog = 0

        for i in nopeakList[args.limits[0]:args.limits[1]]:

            z = 0 # goes on peakIdx to index the progress bar

            """
            READ DATA FROM NOT-PEAKED FILE
            creates a Supernova object
            """
            filePath = i
            try:
                tmpSN = util.get_sn_from_file(filePath)

                print "Progress: {:<d} -- {:<}".format(prog, filePath)
                prog += 1

                ccIndent = "ID:{: ^7d}".format(tmpSN.SNID)

                widgets = [ccIndent, Percentage(), ' ',
                   Bar(marker='#',left='[',right=']'),
                   ' ', ETA()]
                pbar = ProgressBar(widgets=widgets, maxval=len(peakList)).start()
            except IOError:
                print "IOError: {:<}".format(filePath)
                continue

            if tmpSN.r.badCurve:
                print "IOError (BAD r curve): {:<}".format(filePath)
                continue
            """
            create SupernovaFit object
            """
            notPeaked = cls.SupernovaFit(tmpSN)
            for l in tmpSN.lcsDict.keys():
                notPeaked.set_lightcurve(l,
                    tmpSN.lcsDict[l].mjd,
                    tmpSN.lcsDict[l].flux,
                    tmpSN.lcsDict[l].fluxErr
                    )
            """
            Shifting mjds in not-peaked
            """
            notPeaked.shift_mjds()

            ccMax = list()#np.zeros(peakIdx.size)
            k = 0 # goes on ccMax
            # for j in peakIdx:
            for j in peakList:
                """
                READ DATA FROM PEAKED FILE
                """
                # if j in repeats:
                #     print indent + bcolors.WARNING + \
                #         'File appears also in unpeaked list: ignoring it.' + \
                #             bcolors.txtrst
                #     continue
                filePath = j#args.dirFit + lsDirData[j][0:12] + '_FIT.DAT'
                try:
                    tmpSN = util.get_sn_from_file(filePath)
                except IOError:
                    print indent + bcolors.WARNING + \
                    'File appears also in peaked list but it does not exists: ignoring it.' + \
                    bcolors.txtrst
                    continue

                if tmpSN.r.badCurve:
                    print indent + bcolors.WARNING + \
                    'Peaked file has bad r curve: ignoring it.' + \
                    bcolors.txtrst
                    continue
                peaked = cls.SupernovaFit(tmpSN)
                for l in tmpSN.lcsDict.keys():
                    peaked.set_lightcurve(l,
                        tmpSN.lcsDict[l].mjd,
                        tmpSN.lcsDict[l].flux,
                        tmpSN.lcsDict[l].fluxErr
                        )
                """
                Shifting mjds in peaked
                """
                peaked.shift_mjds()

                """
                Performing cross-correlation
                """
                ycorr = signal.correlate(
                    notPeaked.normalized_flux('r'),
                    peaked.normalized_flux('r')
                    )
                xcorr = np.arange(ycorr.size)
                lags = xcorr - (
                    len(notPeaked.normalized_flux('r'))-1
                    )
                distancePerLag = (
                    notPeaked.r.shiftedMjd[-1] - \
                    notPeaked.r.shiftedMjd[0])/float(
                                        len(notPeaked.r.shiftedMjd)
                                        )
                offsets = -lags*distancePerLag

                # ccMax[k] = offsets[np.argmax(ycorr)]
                ccMax.append(offsets[np.argmax(ycorr)])
                # k += 1

                pbar.update(z+1)
                z += 1

                # gc.collect()
            notPeaked.ccMjdMaxFlux = np.mean(ccMax)#ccMax.mean()
            """
            re-writing file of not peaked lc to include information on maximum
            position from CC.
            """
            filePath = i#args.dirFit + lsDirData[i][0:12] + '_FIT.DAT'
            notPeaked.save_on_txt(filePath)

            reWrite.write(filePath+'\n')
            pbar.finish()

            # gc.collect()
        reWrite.close()
        print 'CC ended!'
        gc.collect()

    """

    CALCULATING DISTANCE MATRIX



    needs:
    - args.distMatrix
    - args.limits
    - args.offset
    - args.dirFit
    """
    if args.distMatrix:
        if not os.path.exists(path.abspath(args.dirFit + 'distance_matrix' + os.sep)):
            os.makedirs(path.abspath(args.dirFit + 'distance_matrix' + os.sep))

        """
        Calculate distance between fitted lightcurves.
        Distance values are saved in a R matrix. This will be used by the R
        package `diffusionMap` through rpy2 Python package.
        """
        j_offset = args.offset
        i_start = args.limits[0]
        i_end = args.limits[1]
        j_start = i_start + j_offset
        j_end = (i_end + j_offset) if (i_end+j_offset<=len(lsDirFit)) else len(lsDirFit)
        print "\n" + indent + bcolors.undwht + \
            "(*) Calculate distances between lightcurves ..." + \
            bcolors.txtrst
        print indent + "Rows in [{:<d}, {:<d})".format(i_start, i_end)
        print indent + "Cols in [{:<d}, {:<d})".format(j_start, j_end)

        """
        setting value for big distance
        """
        distFlag = 5
        missColCount = 0
        missRowlist = list()
        bandDict = {
            'g':0,
            'r':1,
            'i':2,
            'z':3
            }
        widgets = [indent, 'Processing:', ' ', Counter(), ' ',
            AnimatedMarker(), indent, Timer()]

        # creating list of 4 lists
        distList = list([[], [], [], []])
        nCols = 0
        # distList = np.zeros((4,
        #     len(lsDirFit[i_start:i_end]), len(lsDirFit[i_start:i_end])),
        #     dtype=float
        #     )


        pbar = ProgressBar(widgets=widgets, maxval=(i_end-i_start)).start()

        for i in range(i_start, i_end):
            missColCount = 0
            """
            Reading in i-candidate
            """
            tmpSN = util.get_sn_from_file(
                args.dirFit+lsDirFit[i]
                )
            if tmpSN.r.badCurve:
                # nothing has to be added to the distance matrix. Print and
                #
                # continue to nex object

                # print "{:<} Has bad curve in r band - ".format(lsDirFit[i]) + \
                #     "THE FILE HAS TO BE DELETED" +\
                #     " indices {:<d}".format(i)
                missRowlist.append(i)
                continue

            iCandidate = cls.SupernovaFit(tmpSN)

            for b in tmpSN.lcsDict.keys():
                # set_lightcurve set also if the lc is peaked or not
                iCandidate.set_lightcurve(b,
                    tmpSN.lcsDict[b].mjd,
                    tmpSN.lcsDict[b].flux,
                    tmpSN.lcsDict[b].fluxErr
                    )

            """
            Shifting mjds in i-candidate
            """
            iCandidate.shift_mjds()
            if iCandidate.peaked == False:
                # print i, iCandidate.SNID
                """
                keeping to perform check with other non peaked LC
                """
                iElMax = iCandidate.r.shiftedMjd.index(0.)
                """
                correcting using CC results
                """
                for b in bands:
                    iCandidate.lcsDict[b].shiftedMjd = [
                        iCandidate.lcsDict[b].shiftedMjd[l] +
                        iCandidate.ccMjdMaxFlux for l in range(len(
                            iCandidate.lcsDict[b].shiftedMjd
                            ))
                    ]
            iElSize = iCandidate.r.size
            iPeaked = iCandidate.peaked

            for j in range(j_start, j_end):
                """
                if this SN has badCurve in this band it will be far from all
                the others by default.
                here will save time from not opening all the other files
                to create new SupernovaFit objcets.
                """

                if j == i:
                    # filling elements on the distance matrix diagonal
                    for b in bands:
                        # adding one element to each sub list in distList
                        distList[bandDict[b]].append(0.)
                        # distList[bandDict[b], i-i_start, j-j_start] = 0.
                    continue

                if j < i:
                    # filling matrix elements below the diagonal
                    if j in missRowlist:
                        missColCount += 1
                        continue
                    for b in bands:
                        # appending the symmetric element in the list: i-i_start
                        distList[bandDict[b]].append(
                            distList[bandDict[b]][
                                (j-j_start-missColCount)*nCols+\
                                    i-i_start-len(missRowlist)
                                ])
                        # distList[bandDict[b], i-i_start, j-j_start] = \
                        #             distList[bandDict[b], j-j_start, i-i_start]
                    continue # jump to the next iteration of the loop

                """
                Reading in j-candidate
                """
                try:
                    tmpSN = util.get_sn_from_file(
                        args.dirFit+lsDirFit[j]
                    )
                except IndexError:
                    print j, len(lsDirFit)
                    raise IndexError("list index out of range")
                if tmpSN.r.badCurve:
                    # nothing has to be added to the distance matrix. Print and
                    #
                    # continue to nex object

                    # print "{:<} Has bad curve in r band -".format(lsDirFit[j])+\
                    #     " THE FILE HAS TO BE DELETED:" +\
                    #     " indices {:<d}, {:<d}".format(i, j)
                    continue

                jCandidate = cls.SupernovaFit(tmpSN)
                for b in tmpSN.lcsDict.keys():
                    jCandidate.set_lightcurve(b,
                        tmpSN.lcsDict[b].mjd,
                        tmpSN.lcsDict[b].flux,
                        tmpSN.lcsDict[b].fluxErr
                        )

                """
                Shifting mjds in j-candidate
                """
                jCandidate.shift_mjds()
                if jCandidate.peaked == False:
                    """
                    keeping to perform check with other non peaked LC
                    """
                    jElMax = jCandidate.r.shiftedMjd.index(0.)
                    """
                    correcting using CC results
                    """
                    for b in bands:
                        jCandidate.lcsDict[b].shiftedMjd = [
                            jCandidate.lcsDict[b].shiftedMjd[l] +
                            jCandidate.ccMjdMaxFlux for l in range(len(
                                jCandidate.lcsDict[b].shiftedMjd
                                ))
                        ]

                jElSize = jCandidate.r.size

                for b in bands:
                    if not jCandidate.lcsDict[b].badCurve \
                    and not iCandidate.lcsDict[b].badCurve:
                        distList[bandDict[b]].append(
                            iCandidate.get_distance(jCandidate, b)
                            )
                        # distList[bandDict[b], i-i_start, j-j_start] = \
                        #     iCandidate.get_distance(jCandidate, b)
                    else:
                        # in case of bad curve
                        """
                        This works like a flag. These elements will be set
                        equal to a neutral value (the mean of the other)
                        """
                        distList[bandDict[b]].append(distFlag)
                        # distList[bandDict[b], i-i_start, j-j_start] = distFlag

            """
            # >>> !! Checking for i being equal to its beginning value in the loop
            does not take into account the
            possibility of the first SN having a bad r curve, in which case
            the loop will never arrive here, since it is reset by a continue.
            Checking on nCols being still equal to zero is much better, since is
            the only way to verify if the first loop has been completed.
            """
            # if (i == i_start):
            if (nCols == 0):
                nCols = len(distList[0])
                print 'nCols updated! {:<d}'.format(nCols)
            pbar.update(i-i_start+1)



        pbar.finish()

        # del iCandidate
        # del jCandidate
        # del tmpSN

        gc.collect()

        distMatrix = np.zeros((4,
            len(distList[0])/nCols, nCols),
            dtype=float
            )

        for b in bands:
            distMatrix[bandDict[b]] = np.reshape(
                distList[bandDict[b]], (len(distList[bandDict[b]])/nCols, nCols)
                )

        """
        distList is no more used from now on. I delete it to save memory
        """
        del distList

        gc.collect()

        # fixing flagged elements
        # raise SystemExit
        if distMatrix[0, distMatrix[0] == distFlag].size > 0:
            ind = np.where(distMatrix[0] == distFlag)
            distMatrix[0, ind[0], ind[1]] = np.add(
                np.add(
                    distMatrix[1, ind[0], ind[1]],
                    distMatrix[2, ind[0], ind[1]]
                    ),
                distMatrix[3, ind[0], ind[1]]
                )/3.


        if distMatrix[1, distMatrix[1] == distFlag].size > 0:
            ind = np.where(distMatrix[1] == distFlag)
            # distMatrix[1, ind[0], ind[1]] = distMatrix[1,:,:].max()
            distMatrix[1, ind[0], ind[1]] = np.add(
                np.add(
                    distMatrix[0, ind[0], ind[1]],
                    distMatrix[2, ind[0], ind[1]]
                    ),
                distMatrix[3, ind[0], ind[1]]
                )/3.

        if distMatrix[2, distMatrix[2] == distFlag].size > 0:
            ind = np.where(distMatrix[2] == distFlag)
            # distMatrix[2, ind[0], ind[1]] = distMatrix[2].max()
            distMatrix[2, ind[0], ind[1]] = np.add(
                np.add(
                    distMatrix[0, ind[0], ind[1]],
                    distMatrix[1, ind[0], ind[1]]
                    ),
                distMatrix[3, ind[0], ind[1]]
                )/3.

        if distMatrix[3, distMatrix[3] == distFlag].size > 0:
            ind = np.where(distMatrix[3] == distFlag)
            # distMatrix[3, ind[0], ind[1]] = distMatrix[3].max()
            distMatrix[3, ind[0], ind[1]] = np.add(
                np.add(
                    distMatrix[0, ind[0], ind[1]],
                    distMatrix[1, ind[0], ind[1]]
                    ),
                distMatrix[2, ind[0], ind[1]]
                )/3.

        distMatrixSum = np.sum(distMatrix, 0)
        """
        Saving on text files
        """

        fileHeader = "distMatrix[{:<d}:{:<d},{:<d}:{:<d}] --- ".format(
            i_start, i_end, j_start, j_end
            ) + \
            "Created by {:<}".format(socket.gethostname())

        filePath = args.dirFit + 'distance_matrix' + os.sep + \
            'dist_matrix_Sum_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, distMatrixSum, fmt='%6.4f', header=fileHeader)

        del distMatrixSum
        gc.collect()

        filePath = args.dirFit + 'distance_matrix' + os.sep + \
            'dist_matrix_g_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, distMatrix[0], fmt='%6.4f', header=fileHeader)

        filePath = args.dirFit + 'distance_matrix' + os.sep + \
            'dist_matrix_r_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, distMatrix[1], fmt='%6.4f', header=fileHeader)

        filePath = args.dirFit + 'distance_matrix' + os.sep + \
            'dist_matrix_i_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, distMatrix[2], fmt='%6.4f', header=fileHeader)

        filePath = args.dirFit + 'distance_matrix' + os.sep + \
            'dist_matrix_z_{:<}_{:<5.3f}.txt'.format(
                socket.gethostname(), time.time()
            )
        np.savetxt(filePath, distMatrix[3], fmt='%6.4f', header=fileHeader)

        del distMatrix
        gc.collect()


    """

    CALCULATING DIFFUSION MAP

    """


    if args.diffuse:
        if 'diffusionMap' not in globals():
            diffusionMap = importr('diffusionMap')

        ndim = ro.r.attributes(Rmatrix)[0][0]
        dmap = diffusionMap.diffuse(Rmatrix, neigen=5)
        util.dump_pkl('diffusion_map.pkl', dmap)


    """

    TRAINING RANDOM FOREST CLASSIFIER

    """

    if args.train:
        randomForest = importr('randomForest')
        if 'dmap' not in globals():
            print indent + 'Loading catalog from dump file ...'
            dmap = util.open_pkl('tmp_diffusion_map.pkl')

        dmap_rf = randomForest.randomForest(dmap)


    """

    PLOT OBSERVATION AND FIT
    --plot
    """

    if args.plot:
        timeMark = time.time()
        """
        getting file list from directory
        File will be sorted by SNID
        """

        print indent + 'Plotting ...'
        '''
        Column index is always increasing, no check on its value.
        '''
        nrows = 5
        ncols = 5
        """
        If plotOffset is to specified, get a proper random value
        """
        if (args.plotOffset == -1):
            np.random.RandomState
            offset = int(np.random.uniform(low=0, high=len(lsDirFit)-nrows*ncols))
        else:
            offset = args.plotOffset

        fig_g, ax_g = plt.subplots(nrows=nrows, ncols=ncols,
                    figsize=(16.5, 11.7)#,
                    #tight_layout=True
                    )
        fig_r, ax_r = plt.subplots(nrows=nrows, ncols=ncols,
                    figsize=(16.5, 11.7)#,
                    #tight_layout=True
                    )
        fig_i, ax_i = plt.subplots(nrows=nrows, ncols=ncols,
                    figsize=(16.5, 11.7)#,
                    #tight_layout=True
                    )
        fig_z, ax_z = plt.subplots(nrows=nrows, ncols=ncols,
                    figsize=(16.5, 11.7)#,
                    # tight_layout=True
                    )

        dictFig = {'g':fig_g,
                   'r':fig_r,
                   'i':fig_i,
                   'z':fig_z}

        dictAx = {'g':ax_g,
                  'r':ax_r,
                  'i':ax_i,
                  'z':ax_z}

        r = {'g':0,
             'r':0,
             'i':0,
             'z':0}

        c = {'g':0,
             'r':0,
             'i':0,
             'z':0}

        """
        Adjust subplot margins and title
        """
        for b in dictFig.keys():
            dictFig[b].subplots_adjust(
                top=0.96, right=0.99, bottom=0.03, left=0.02,
                wspace=0.08, hspace=0.13
                )
            dictFig[b].suptitle('band {:<1} - offset {:<d}'.format(b, offset))

        GPkern = ''
        for i in range(nrows*ncols):
            """
            Getting the observational data from file
            """
            candidate = util.get_sn_from_file(
                args.dirData + lsDirData[i+offset]#candidateIdx]
                )

            """
            Reading fit data from file
            """
            try:
                tmpSN = util.get_sn_from_file(
                            args.dirFit+lsDirFit[i+offset],
                            magFlag=args.mag,
                        )
            except IndexError:
                warnStr = 'IndexError: list index out of range. '+\
                    'i={:<d}.'.format(i+offset)
                print warnings.warn(warnStr)

                print '\n'+indent+'Saving files as they are and stopping.'
            else:
                """
                Initializing SupernovaFit object
                """
                fit = cls.SupernovaFit(tmpSN,
                        tmpSN.kern if hasattr(tmpSN, 'kern') else None)
                if (i == 0) and hasattr(tmpSN, 'kern'):
                    GPkern = tmpSN.kern
                for b in tmpSN.lcsDict.keys():
                    fit.set_lightcurve(b,
                        tmpSN.lcsDict[b].mjd,
                        tmpSN.lcsDict[b].flux,
                        tmpSN.lcsDict[b].fluxErr,
                        magFlag=args.mag
                        )
                if fit.r.badCurve:
                    print 'SN ID{:>06d} has bad r band light curve!'.format(
                    fit.SNID)
                    # continue
                else:
                    """
                    Shift fit mjd to have 0 at r band maximum
                    """
                    fit.shift_mjds()
                """
                Fixing shiftedMjd for not-peaked LCs
                """
                if (fit.peaked == False) and (fit.r.badCurve == False) :
                    """
                    correcting using CC results
                    """
                    for b in bands:
                        fit.lcsDict[b].shiftedMjd = [
                        el + fit.ccMjdMaxFlux for el in fit.lcsDict[b].shiftedMjd
                        ]


                for b in dictAx.keys():
                    """
                    variable `data` initialized as light curve in band b for
                    cleaner code.
                    """
                    data = candidate.lcsDict[b]

                    fit_b = fit.lcsDict[b]

                    fit_r = fit.lcsDict['r']

                    if c[b] > nrows-1:
                        c[b] = 0
                        r[b] += 1

                    xlim = dictAx[b][r[b], c[b]].get_xlim()
                    ylim = dictAx[b][r[b], c[b]].get_ylim()

                    dictAx[b][r[b], c[b]].set_xticks([0])
                    dictAx[b][r[b], c[b]].set_yticks([0])

                    dictAx[b][r[b], c[b]].set_xticklabels(['0'])
                    dictAx[b][r[b], c[b]].set_yticklabels(['0'])

                    if (data.badCurve == False) and (fit_b.badCurve == False) and (fit.r.badCurve == False):
                        epoch = util.time_correct(data.mjd,
                            candidate.zSpec if candidate.zSpec else candidate.zPhotHost)

                        epoch = [val-fit_r.mjd[fit_r.max_flux_index] for val in epoch]
                        if fit.peaked == False:
                            epoch = [val+fit.ccMjdMaxFlux for val in epoch]

                        flux = util.correct_for_absorption(data.flux,
                            candidate.MWEBV, b)

                        """
                        Setting limits for plot axes
                        """
                        if min(fit_b.flux) < min(flux):
                            y_min = min(fit_b.flux) - 3*max(fit_b.fluxErr)
                        else:
                            y_min = min(flux) - np.median(data.fluxErr)

                        if max(fit_b.flux) > max(flux):
                            y_max = max(fit_b.flux) + 3*max(fit_b.fluxErr)
                        else:
                            y_max = max(flux) + np.median(data.fluxErr)

                        dictAx[b][r[b], c[b]].set_ylim(y_min, y_max)

                        """
                        Setting limits for fill_between
                        """
                        fluxUpLim = [val for val in [
                            fit_b.flux[el] + fit_b.fluxErr[el]
                                for el in range(len(fit_b.flux))
                            ]]
                        fluxLowLim = [val for val in [
                            fit_b.flux[el] - fit_b.fluxErr[el]
                                for el in range(len(fit_b.flux))
                            ]]

                        dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd,
                            fluxUpLim, fluxLowLim,
                            facecolor='red', alpha=0.4, linewidth=0.5)

                        """
                        Setting limits for fill_between
                        """
                        fluxUpLim = [val for val in [
                            fit_b.flux[el] + 2*fit_b.fluxErr[el]
                                for el in range(len(fit_b.flux))
                            ]]
                        fluxLowLim = [val for val in [
                            fit_b.flux[el] - 2*fit_b.fluxErr[el]
                                for el in range(len(fit_b.flux))
                            ]]

                        dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd,
                            fluxUpLim, fluxLowLim,
                            facecolor='red', alpha=0.2, linewidth=0.5)

                        """
                        Setting limits for fill_between
                        """
                        fluxUpLim = [val for val in [
                            fit_b.flux[el] + 3*fit_b.fluxErr[el]
                                for el in range(len(fit_b.flux))
                            ]]
                        fluxLowLim = [val for val in [
                            fit_b.flux[el] - 3*fit_b.fluxErr[el]
                                for el in range(len(fit_b.flux))
                            ]]

                        dictAx[b][r[b], c[b]].fill_between(fit_b.shiftedMjd,
                            fluxUpLim, fluxLowLim,
                            facecolor='red', alpha=0.1, linewidth=0.5)


                        dictAx[b][r[b], c[b]].plot(fit_b.shiftedMjd, fit_b.flux,
                            color='#7f0000',
                            linewidth=2)

                        scatterLab = 'SN ID {:<d}'.format(candidate.SNID)
                        dictAx[b][r[b], c[b]].scatter(epoch, flux,
                            s=10, label=scatterLab, c='black', marker='x')

                        dictAx[b][r[b], c[b]].errorbar(epoch, flux,
                            data.fluxErr, fmt=None, color='black', ecolor='black')


                        if not fit.peaked:
                            pass

                        dictAx[b][r[b], c[b]].legend(
                            loc='best', framealpha=0.3, fontsize='10')
                    else:
                        label = str(candidate.SNID)+" BAD CURVE"
                        dictAx[b][r[b], c[b]].plot([0, 1], [0, 1], color='red',
                            label=label)
                        dictAx[b][r[b], c[b]].plot([0, 1], [1, 0], color='red')
                        dictAx[b][r[b], c[b]].legend(
                            loc='best', framealpha=0.3, fontsize='10')
                    c[b] += 1


        print indent + "Plots saved in files:"
        if not os.path.exists(path.abspath(args.dirFit + "plots" + os.sep)):
            os.makedirs(args.dirFit + "plots")
        for b in dictFig.keys():
            dictFig[b].savefig(
                args.dirFit + "plots"+ os.sep + GPkern + \
                "_band_{:<1}_{:<f}.png".format(b,timeMark),
                dpi=300
                )
            print indent + " - " + args.dirFit + "plots" + os.sep + \
                GPkern + "_band_{:<1}_{:<f}.png".format(b,timeMark)

        plt.close('all')


    """

    PLOT OBSERVATION AND FIT (publication style)
    --nice-plots
    """

    if args.nicePlots:
        """
        1 candidate
        choose how many bands
        make the plot with confidence regions
        """
        # if args.nBands != 1 or args.nBands != 4:
            # args.nBands = 1

        if args.cand == -1:
            args.cand = np.random.random_integers(
                low=0, high=len(lsDirData))

        fname = 'DES_SN{:0>6d}.DAT'.format(args.cand)
        candidate = util.get_sn_from_file(
                args.dirData+fname
                )

        fname = 'DES_SN{:0>6d}_FIT.DAT'.format(args.cand)
        tmpSN = util.get_sn_from_file(
            args.dirFit+fname,
            magFlag=args.mag,
                )
        """
        Initializing SupernovaFit object
        """
        fit = cls.SupernovaFit(tmpSN, tmpSN.kern if hasattr(tmpSN, 'kern') else None)

        for b in tmpSN.lcsDict.keys():
            fit.set_lightcurve(b,
                tmpSN.lcsDict[b].mjd,
                tmpSN.lcsDict[b].flux,
                tmpSN.lcsDict[b].fluxErr,
                magFlag=args.mag
                )

        if fit.r.badCurve:
            raise SystemExit('Bad r curve!')

        fit.shift_mjds()
        """
        Fixing shiftedMjd for not-peaked LCs
        """
        if fit.peaked == False:
            """
            correcting using CC results
            """
            for b in candidate.lcsDict.keys():
                fit.lcsDict[b].shiftedMjd = [el + fit.ccMjdMaxFlux
                                        for el in fit.lcsDict[b].shiftedMjd]


        bands = candidate.lcsDict.keys() if args.allBands else args.band
        """
        Pre-process data so to be compared with fit (made from
        pre-precessed data)
        """
        for b in bands:
            if (not candidate.lcsDict[b].badCurve) and (not fit.lcsDict[b].badCurve):

                candidate = util.pre_process(candidate, b)

                candidate.lcsDict[b].mjd = [el - fit.r.mjd[fit.r.max_flux_index]
                        for el in candidate.lcsDict[b].mjd]
                if fit.peaked == False:
                    candidate.lcsDict[b].mjd = [el + fit.ccMjdMaxFlux
                            for el in candidate.lcsDict[b].mjd]

            else:
                raise SystemExit('Bad {:1s} curve!'.format(b))

        if args.allBands:

            fig, ax = plt.subplots(nrows=2, ncols=2,
                    # figsize=(16.5, 11.7),
                    tight_layout=False
                    )

            axDict = {
            'g':ax[0,0],
            'r':ax[0,1],
            'i':ax[1,0],
            'z':ax[1,1]
            }
            # fig.subplots_adjust(left=0.05, right=0.97, top=0.94, wspace=0.29)
        else:
            fig = plt.figure()
            xlim = [-35,12]
            ylim = [-10,10]
            # fig, ax = plt.subplots(nrows=2, ncols=1,
            #         # figsize=(16.5, 11.7),
            #         tight_layout=False
            #         )

            # axDict = {
            # 'g':ax[0,0],
            # 'r':ax[0,1],
            # 'i':ax[1,0],
            # 'z':ax[1,1]
            # }


        if not args.allBands:
            fit_b = fit.lcsDict[args.band]
            data = candidate.lcsDict[args.band]
            if not data.badCurve and not fit_b.badCurve:
                epoch = data.mjd

                flux = data.flux
            """
            Setting limits for fill_between
            """
            fluxUpLim = [el for el in [
                fit_b.flux[i] + fit_b.fluxErr[i]
                    for i in range(len(fit_b.flux))
                ]]
            fluxLowLim = [el for el in [
                fit_b.flux[i] - fit_b.fluxErr[i]
                    for i in range(len(fit_b.flux))
                ]]

            plt.fill_between(fit_b.shiftedMjd,
                fluxUpLim, fluxLowLim,
                facecolor='red', alpha=0.4, linewidth=0.5)
            # axDict[b].fill_between(fit_b.shiftedMjd,
            #     fluxUpLim, fluxLowLim,
            #     facecolor='red', alpha=0.4, linewidth=0.5)

            """
            Setting limits for fill_between
            """
            fluxUpLim = [el for el in [
                fit_b.flux[i] + 2*fit_b.fluxErr[i]
                    for i in range(len(fit_b.flux))
                ]]
            fluxLowLim = [el for el in [
                fit_b.flux[i] - 2*fit_b.fluxErr[i]
                    for i in range(len(fit_b.flux))
                ]]

            plt.fill_between(fit_b.shiftedMjd,
                fluxUpLim, fluxLowLim,
                facecolor='red', alpha=0.2, linewidth=0.5)
            # axDict[b].fill_between(fit_b.shiftedMjd,
            #     fluxUpLim, fluxLowLim,
            #     facecolor='red', alpha=0.2, linewidth=0.5)

            """
            Setting limits for fill_between
            """
            fluxUpLim = [el for el in [
                fit_b.flux[i] + 3*fit_b.fluxErr[i]
                    for i in range(len(fit_b.flux))
                ]]
            fluxLowLim = [el for el in [
                fit_b.flux[i] - 3*fit_b.fluxErr[i]
                    for i in range(len(fit_b.flux))
                ]]

            plt.fill_between(fit_b.shiftedMjd,
                fluxUpLim, fluxLowLim,
                facecolor='red', alpha=0.1, linewidth=0.5)
            # axDict[b].fill_between(fit_b.shiftedMjd,
            #     fluxUpLim, fluxLowLim,
            #     facecolor='red', alpha=0.1, linewidth=0.5)

            plt.plot(fit_b.shiftedMjd, fit_b.flux,
                color='#7f0000',
                linewidth=2,
                label='GP fit')
            # axDict[b].plot(fit_b.shiftedMjd, fit_b.flux,
            #     color='#7f0000',
            #     linewidth=2)

            plt.scatter(epoch, flux,
                s=30, label='data', c='black', marker='x')

            # axDict[b].scatter(epoch, flux,
            #     s=10, label=str(candidate.SNID), c='black', marker='x')

            plt.errorbar(epoch, flux,
                data.fluxErr, fmt=None, color='black', ecolor='black')

            # plt.xlim(xlim)
            plt.ylim(ylim)

            title = 'SN ID {:d} - Band {:s}'.format(candidate.SNID, args.band)

            plt.title(title)
            plt.xlabel('Epoch [mjd]')
            plt.ylabel('Flux [adu]')
            plt.legend(loc='upper right', scatterpoints=1)
            # axDict[b].errorbar(epoch, flux,
            #     data.fluxErr, fmt=None, color='black', ecolor='black')




    print "\n" + indent \
        + "The process took {:5.3f} secs.".format(time.time()-start_time)
