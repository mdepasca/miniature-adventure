import lightcurve
import supernova
import supernova_fit
import numpy as np
import os
import sys
import time
import argparse

import subprocess

warnings.filterwarnings(
    'error',
    message=".*divide by zero encountered in double_scalars.*",
    category=RuntimeWarning
    )
from math import sqrt

if __name__ == '__main__':
    import utilities as util
    import GPy
    from matplotlib import pyplot as plt



    parser = argparse.ArgumentParser(
        description="Test of general functions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    actionGroup = parser.add_argument_group('ACTION')
    inputGroup = parser.add_argument_group('INPUT')

    """

    ACTION OPTIONS

    """
    actionGroup.add_argument(
        "--k-correction", dest='kcor',
        action='store_true', help='Switch on k correction.'
        )

    actionGroup.add_argument(
        '--distance', dest='distance',
        action='store_true', help='Calculate distance between fitted lightcurves \
        in same band.'
        )

    actionGroup.add_argument(
        '--test-prior', dest='testPrior',
        action='store_true', help='Test prior in GP regression.'
        )

    """

    INPUT OPTIONS

    """
    inputGroup.add_argument(
        "--data-directory", dest="dirData",
        default="train_data" + os.sep + "SIMGEN_PUBLIC_DES",
        help="Path to directory containing training data.")

    inputGroup.add_argument(
        "-b", "--band", dest="band", default='r',
        help="Photometric band.")

    inputGroup.add_argument(
        "-c1", "--candidate1", dest="candidate1",
        type=np.int32, default=None,
        help="First candidate idx")

    inputGroup.add_argument(
        "-c2", "--candidate2", dest="candidate2",
        type=np.int32, default=None,
        help="Second candidate idx")

    inputGroup.add_argument(
        "--mag", dest="mag",
        action="store_true",
        help="Reads in magnitudes from file."
        )

    args = parser.parse_args()

if __name__ == '__main__':
    indent = "          "
    lambda_obs = [479.66, 638.26, 776.90, 910.82]
    limMagDict = {
        'g': 25.2,
        'r': 25.4,
        'i': 25.1,
        'z': 24.9
    }
    p = subprocess.Popen("ls *SN*.DAT", shell=True, stdout=subprocess.PIPE,
            cwd=args.dirData+os.sep)
    lsDirData = p.stdout.read()
    lsDirData = lsDirData.split('\n')
    lsDirData.sort()
    lsDirData.remove('')
    # fCandidatesList = "DES_BLIND+HOSTZ.LIST"
    # candidatesFileList = np.genfromtxt(dirData+os.sep+fCandidatesList, dtype=None)\

    """

    KERNEL SPECIFICATION

    """
    # kern = GPy.kern.RBF(1)
    kern = GPy.kern.RatQuad(1)

    """
    ----------------------
    """

    if args.band not in ['g', 'r', 'i', 'z']:
        print 'Band {:<} not recognised! Changing to r'.format(args.band)
        args.band = 'r'

    if args.candidate1 is None:
        args.candidate1 = np.random.random_integers(
                low=0, high=len(lsDirData))

    if args.candidate2 is None:
        args.candidate2 = np.random.random_integers(
                low=0, high=len(lsDirData))

    while args.candidate2 == args.candidate1:
        args.candidate2 = np.random.random_integers(
            low=0, high=len(lsDirData))

    print args.candidate1
    print args.candidate2

    candidates = list()
    fit = list()

    """
    Getting observation's data
    """

    candidates.append(Supernova(
        args.dirData+os.sep+lsDirData[args.candidate1], args.mag))

    candidates.append(Supernova(
        args.dirData+os.sep+lsDirData[args.candidate2], args.mag))


    for candidate in candidates:
        """
        Setting limits in magnitudes
        """
        for b in candidate.lcsDict.keys():
            if args.mag:
                candidate.lcsDict[b].lim = limMagDict[b]


        print 'candidate z = {:>6.4f}'.format(
            candidate.zSpec if candidate.zSpec else candidate.zPhotHost)

        if args.kcor:
            lambda_rf = [el/(
                1+(candidate.zSpec if candidate.zSpec else candidate.zPhotHost)
                ) for el in lambda_obs]

            plt.figure()
            int_mjd_g = [int(el) for el in candidate.g.mjd]
            int_mjd_r = [int(el) for el in candidate.r.mjd]
            int_mjd_i = [int(el) for el in candidate.i.mjd]
            int_mjd_z = [int(el) for el in candidate.z.mjd]
            mjd = [el for el in int_mjd_g if el in int_mjd_r
                        and el in int_mjd_i and el in int_mjd_z]
            jd = mjd[0]

            flux = [
                candidate.g.flux[int_mjd_g.index(jd)],
                candidate.r.flux[int_mjd_r.index(jd)],
                candidate.i.flux[int_mjd_i.index(jd)],
                candidate.z.flux[int_mjd_z.index(jd)]
                ]
            fluxErr = [
                candidate.g.fluxErr[int_mjd_g.index(jd)],
                candidate.r.fluxErr[int_mjd_r.index(jd)],
                candidate.i.fluxErr[int_mjd_i.index(jd)],
                candidate.z.fluxErr[int_mjd_z.index(jd)]
                ]
            plt.errorbar(lambda_obs, flux, yerr=fluxErr, fmt='k--', ecolor='black')
            # if args.mag:
            #     plt.xlim([plt.ylim()[1]], plt.ylim()[0]])
            plt.scatter(lambda_obs, flux, color='black')
            plt.errorbar(lambda_rf, flux, yerr=fluxErr, fmt='b--', ecolor='blue')
            plt.scatter(lambda_rf, flux, color='blue')
            plt.show()
            # mjd_k_corr, k_correct_flux = candidate.k_correct_flux()


            # (a, b) = np.polyfit(
            #     # [lambda_rf[0], lambda_rf[1]],
            #     [
            #     lambda_rf[0],
            #     lambda_rf[1]#,
            #     # lambda_rf[2],
            #     # lambda_rf[3],
            #     ],
            #     [
            #     candidate.g.flux[int_mjd_g.index(jd)],
            #     candidate.r.flux[int_mjd_r.index(jd)]#,
            #     # candidate.i.flux[int_mjd_i.index(jd)],
            #     # candidate.z.flux[int_mjd_z.index(jd)]
            #     ], deg = 1,
            #     w = [
            #         1/candidate.g.fluxErr[int_mjd_g.index(jd)],
            #         1/candidate.r.fluxErr[int_mjd_r.index(jd)]]#,
            #     # cov=True
            #     )
            # raise ValueError

        """
        Create SupernovaFit objects
        """
        candidateFit = SupernovaFit(candidate, kern.name)

        """
        Looping in bands and fit of flux
        """
        for b in candidate.lcsDict.keys():

            phase = candidate.lcsDict[b].mjd
            flux = candidate.lcsDict[b].flux

            """
            Clipping to limiting magnitudes when flux is above 90th mag
            """
            if args.mag :
                flux = [candidate.lcsDict[b].lim if \
                        (el>90) else el for el in flux]

            errFlux = candidate.lcsDict[b].fluxErr

            # Fitting Lightcurve
            if (not candidate.lcsDict[b].badCurve) and (len(flux) >= 3):

                start_time = time.time()
                predMjd, predFlux, predErr, GPModel = util.gp_fit(
                                                phase, flux, errFlux,
                                                kern, n_restarts=10,
                                                parallel=False,
                                                test_length=True,
                                                test_prior=args.testPrior)
                print "\n" + indent \
                    + "The process took {:5.3f} secs.".format(time.time()-start_time)
                print GPModel

                candidateFit.set_lightcurve(
                    b, predMjd, predFlux, predErr, args.mag
                    )

                print indent + \
                    "{:<} {:<}".format(candidate.SNID, b)
            else:
                candidateFit.lcsDict[b].badCurve = True
                print indent + util.bcolors.FAIL + \
                    "{:<} {:<}".format(candidate.SNID, b) + \
                    util.bcolors.ENDC

        candidateFit.shift_mjds()
        fit.append(candidateFit)


    if args.distance:
        if fit[0].peaked and fit[1].peaked:
            print 'Distance between the 2 normalized lcs in ' + \
            '{:<} band = {:<2.4f}'.format(args.band,
                fit[0].get_distance(fit[1],
                args.band))

            # if plt.get_fignums():
            #     figNum = plt.get_fignums()[-1]+1
            # else:
            #     figNum = 1

            # plt.figure(figNum)
        else:
            print 'One of the 2 candidate has not r-band peak: '
            print 'Candidate 1 - {:<}'.format(fit[0].peaked)
            print 'Candidate 2 - {:<}'.format(fit[1].peaked)

    nrows = 2
    ncols = 4

    colorList = ['blue', 'orange']

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols,
                    figsize=(16.5, 11.7),
                    tight_layout=False
                    )

    axDict = {
    'g':[ax[0,0], ax[0,2]],
    'r':[ax[0,1], ax[0,3]],
    'i':[ax[1,0], ax[1,2]],
    'z':[ax[1,1], ax[1,3]]
    }

    fig.subplots_adjust(left=0.05, right=0.97, top=0.94, wspace=0.29)
    fig.suptitle(eval('\'Band \' + args.band + (\' with Prior Test\' ' +
        'if args.testPrior else \'Band \' + args.band + \' No prior\') + \' -- \'' +
        '+ \'Kernel: \' + kern.name'))

    for j in [0,1]:
        for b in axDict.keys():
            if args.mag:
                upperlimits = [
                    0 if el < 90 else 1 for el in candidates[j].lcsDict[b].fluxErr
                ]
                axDict[b][j].set_ylim(candidates[j].lcsDict[b].lim+2, 22)
                lowerlimits = False
            else:
                upperlimits = False
                lowerlimits = [0 if el > 0 else 1 for el in candidates[j].lcsDict[b].flux]
            axDict[b][j].plot([min(candidates[j].lcsDict[b].mjd),
                max(candidates[j].lcsDict[b].mjd)],
                [candidates[j].lcsDict[b].lim]*2,
                c='k')

            fluxUpLim = [val for val in [
                        fit[j].lcsDict[b].flux[i] +
                        2*fit[j].lcsDict[b].fluxErr[i]
                            for i in range(len(fit[j].lcsDict[b].flux))
                        ]]
            fluxLowLim = [val for val in [
                fit[j].lcsDict[b].flux[i] -
                2*fit[j].lcsDict[b].fluxErr[i]
                    for i in range(len(fit[j].lcsDict[b].flux))
                        ]]
            axDict[b][j].fill_between(fit[j].lcsDict[b].mjd,
                fluxUpLim, fluxLowLim,
                alpha=0.2, linewidth=0.5)

            axDict[b][j].plot(
                fit[j].lcsDict[b].mjd,
                fit[j].lcsDict[b].flux, c=colorList[j],
                )
            axDict[b][j].errorbar(
                candidates[j].lcsDict[b].mjd,
                candidates[j].lcsDict[b].flux,
                candidates[j].lcsDict[b].fluxErr,
                uplims=upperlimits, lolims=lowerlimits, ecolor=colorList[j],
                fmt=None
                )
            axDict[b][j].scatter(
                candidates[j].lcsDict[b].mjd,
                candidates[j].lcsDict[b].flux,
                c=colorList[j],
                label='Band {:>s} | IDX {:>5d} | SNID {:>5d}'.format(b,
                    eval('args.candidate1 if j == 0 else args.candidate2'),
                    candidates[j].SNID)
                )
            axDict[b][j].legend(
                loc='best', framealpha=0.3, fontsize='10'
                )

            axDict[b][j].set_xlabel('epoch [MJD]')
            if args.mag:
                axDict[b][j].set_ylabel('flux [mag]')
            else:
                axDict[b][j].set_ylabel('flux [adu]')


    plt.show()
