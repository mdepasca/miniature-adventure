import utilities as util
import argparse
import os
import sys
import numpy as np
import GPy
import time
import classes as cls
from matplotlib import pyplot as plt
# from progressbar import ProgressBar, SimpleProgress

if __name__ == "__main__":
    # Parsing input from command line
    parser = argparse.ArgumentParser(
        description = "SN lightcurve fitter and classifier.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-f", "--fitting", dest="fitting",
        action="store_true",
        help="Fit lightcurves")

    parser.add_argument(
        "-z", "--zero-point", dest="zeroPoint",
        action="store_true",
        help="Set zero point in time. For each object it is the time, in MJD, \
        of maximum observed in r-band (tailored on simulated data from SNANA).")

    parser.add_argument(
        "-d", "--distance-metric", dest="distMetric",
        action="store_true",
        help="Calculate distance between fitted lightcurves in same band. \
        It is use to build a diffusion map (see Coifman & Lafon (2006) \
        and Lafon & Lee (2006)).")

    parser.add_argument(
        "--training-directory", dest="args.dirData",
        default="train_data" + os.sep + "DES_BLIND+HOSTZ",
        help="Path to directory containing training data.")

    parser.add_argument(
        "--fitting-directory", dest="args.dirFit",
        default="fit_data" + os.sep,
        help="Path to directory in which to save fitting results.")
    args = parser.parse_args()
else:
    pass
    
if __name__ == "__main__":
    os.system("clear")
    indent = "          "
    KERN_RATQUAD = "RatQuad"
    # should be possible to change the next two variables
    # args.dirData = "train_data" + os.sep + "DES_BLIND+HOSTZ"
    # args.dirFit = "fit_data" + os.sep
    fNameCandidatesList = "DES_BLIND+HOSTZ.LIST"

    print indent + "* * * * * * * * * * * * * * *"
    print indent + "*    Miniature Adventure    *"
    print indent + "*    -------------------    *"
    print indent + "*    lightcurves fitting    *"
    print indent + "*             and           *"
    print indent + "*      SN classification    *"
    print indent + "* * * * * * * * * * * * * * *" 
    
    start_time = time.time()
    if args.fitting:
        # Perform fitting

        # Relevant input data
        print "\n" + indent + "[1] * Fit lightcurves"
        print "\n" + indent + "Data directory: " + os.curdir + args.dirData + os.sep
        print "\n" + indent + "List of candidates contained in: " \
            + os.curdir + args.dirData + os.sep + fNameCandidatesList

        vecCandidates = np.genfromtxt(
                args.dirData+os.sep+fNameCandidatesList, dtype=None)
        # tenPercent = vecCandidates.size / 10
        # const = 0
        print "\n" + indent \
            + "Number of candidates = {:<d}".format(vecCandidates.size)

        print "\n" + indent \
            + "Data are fitted using GP with Radial Basis Function kernel."

        # kern = GPy.kern.RatQuad(1)
        kern = GPy.kern.RBF(1)
    
        # Redirecting stderr output to file
        saveErr = sys.stderr
        ferr = open('error.log', 'w')
        sys.stderr = ferr

        # Setting up Progress bar using progressbar module
        # pbar = ProgressBar(
        #     widgets=[SimpleProgress()], 
        #     maxval=vecCandidates.size).start()

        # Fitting single lightcurves 
        #
        # THIS PIECE NEEDS TO BE PARALLELIZED
        # 
        # optimize_restarts parallel using multiprocessing
        candidatesCat = cls.CandidatesCatalog()
        for i in range(vecCandidates.size):
            candidate = util.get_sn_from_file(args.dirData+os.sep+vecCandidates[i])
            candidateFit = cls.SupernovaFit(candidate.SNID)

            for b in candidate.lightCurvesDict.keys():

                phase = candidate.lightCurvesDict[b].mjd
                flux = candidate.lightCurvesDict[b].flux
                errFlux = candidate.lightCurvesDict[b].fluxErr

                # test_prior should be deleted as option. Prior too weak.
                # 
                # Fitting Lightcurve
                if (not candidate.lightCurvesDict[b].badCurve) and \
                    (flux.size >= 3):
                    saveOut = sys.stdout
                    fout = open('out.log', 'w')
                    # fout = open('/dev/null', 'w')
                    sys.stdout = fout

                    
                    predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                    phase, flux, errFlux, 
                                                    kern, n_restarts=10, 
                                                    test_length=True)
                    sys.stdout = saveOut
                    fout.close()

                    candidateFit.set_lightcurve(b, 
                        predMjd.reshape(predMjd.size),
                        predFlux.reshape(predFlux.size), 
                        predErr.reshape(predErr.size))
                    
                    print indent + \
                        "{:<} {:<} {:<}".format(i, candidate.SNID, b)
                else:
                    candidateFit.lightCurvesDict[b].badCurve = True
                    print indent + util.bcolors.FAIL + \
                        "{:<} {:<} {:<}".format(i, candidate.SNID, b) + \
                        util.bcolors.ENDC

                    # print indent + \
                    #     "Candidate {:<d} has ".format(candidate.SNID) + \
                    #     util.bcolors.FAIL + "BAD " + util.bcolors.ENDC + \
                    #     "{:<1} lightcurve \n".format(b)
                                
            # Setting phase 0 point to phase or r maximum
            candidateFit.set_lc_zero_points()    
            candidatesCat.add_candidate(candidateFit)
            
            # candidateFit.save_on_txt(
            #         args.dirFit+"DES_FIT_{:0>6d}.dat".format(candidate.SNID))
            # pbar.update(i + 1)
            if i == 10:
                break
            # if i > const * tenPercent:
            #     pg + 1
            #     print pg
            #     # print "\n" + indent + util.bcolors.OKGREEN + \
            #     #     "{:<d}% Completed".format(const*10) + util.bcolors.ENDC
            #     const += 1
        # pbar.finish()
        sys.stderr = saveErr
        ferr.close()
    if args.zeroPoint:
        
        # 
        
        
        # Set zero point in time 
        pass
    if args.distMetric:
        print"qq"
        # Calculate distance between fitted lightcurves
        pass

    print "\n" + indent \
        + "The process took {:5.3f} secs.".format(time.time()-start_time)