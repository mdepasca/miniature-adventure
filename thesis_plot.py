"""
23/05/2015

Plots to compare results from different Gaussian processes. Using data in
band 'r'.
"""

import classes as cl
import util
import matplotlib.pyplot as plt
import utilities as ut

dictR = dict([('g', 3.237), ('r', 2.176), ('i', 1.595), ('z', 1.217)])
band = 'i'# 'r'
snid = 2896 # 3644 # 3776
xlim = [-36, 46]# [-30, 70] # [20.1, 97]
ylim = [-10, 30]# [-12, 180] # [-38.6,50.5]
fitRBFl = util.IO.get_fit_from_file(
        'results/RBF_test-length/DES_SN{:>06d}_FIT.DAT'.format(snid)
        )
fitRBFp = util.IO.get_fit_from_file(
        'results/RBF-with_prior/DES_SN{:>06d}_FIT.DAT'.format(snid)
        )
fitRATQUADp = util.IO.get_fit_from_file(
        'results/RATQUAD-with_prior/DES_SN{:>06d}_FIT.DAT'.format(snid)
        )

sn = util.IO.get_sn_from_file(
        'train_data/SIMGEN_PUBLIC_DES/DES_SN{:>06d}.DAT'.format(snid)
        )

print 'FIT'
print 'flux_min = {:<f}'.format(min(fitRBFl.lcsDict[band].flux))
print 'flux_max = {:<f}'.format(max(fitRBFl.lcsDict[band].flux))
print 'flux_min = {:<f}'.format(min(fitRBFp.lcsDict[band].flux))
print 'flux_max = {:<f}'.format(max(fitRBFp.lcsDict[band].flux))
print 'flux_min = {:<f}'.format(min(fitRATQUADp.lcsDict[band].flux))
print 'flux_max = {:<f}'.format(max(fitRATQUADp .lcsDict[band].flux))

print 'mjd_min = {:<5.2f}'.format(min(fitRBFl.lcsDict[band].shiftedMjd))
print 'mjd_max = {:<5.2f}'.format(max(fitRBFl.lcsDict[band].shiftedMjd))
print 'mjd_min = {:<5.2f}'.format(min(fitRBFp.lcsDict[band].shiftedMjd))
print 'mjd_max = {:<5.2f}'.format(max(fitRBFp.lcsDict[band].shiftedMjd))
print 'mjd_min = {:<5.2f}'.format(min(fitRATQUADp.lcsDict[band].shiftedMjd))
print 'mjd_max = {:<5.2f}'.format(max(fitRATQUADp.lcsDict[band].shiftedMjd))

print '\nDATA'
print 'flux_min = {:<f}'.format(min(sn.lcsDict[band].flux))
print 'flux_max = {:<f}'.format(max(sn.lcsDict[band].flux))


fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
fig.subplots_adjust(top=0.97, right=0.98, left=0.11, hspace=0)

snFlux = sn.lcsDict[band].calc_dereddened_flux(dictR['r'], sn.MWEBV)


sn.lcsDict[band].set_shifted_mjd_2(
            fitRBFl.lcsDict[band], fitRBFl.r.mjd[fitRBFl.r.max_flux_index],
            fitRBFl.zSpec if fitRBFl.zSpec else fitRBFl.zPhotHost,
            fitRBFl.ccMjdMaxFlux
            )
# lineRBFl = util.plot.plot_lc_fit(
#             fitRBFl.lcsDict[band].shiftedMjd, fitRBFl.lcsDict[band].flux, fitRBFl.lcsDict[band].fluxErr,
#             pltObj=ax[0]
#             )
scatterRBFl = util.plot.plot_lc_data(sn.lcsDict[band].shiftedMjd, snFlux, sn.lcsDict[band].fluxErr,
            pltObj=ax[0])
# lineRBFl.set_label('GP fit')
scatterRBFl.set_label('data')
ax[0].legend(scatterpoints=1, fontsize=10)
ax[0].set_xlabel('Epoch [mjd]')
ax[0].set_ylabel('Flux [adu]')


# sn.lcsDict[band].set_shifted_mjd_2(
#             fitRBFp.lcsDict[band], fitRBFp.r.mjd[fitRBFp.r.max_flux_index],
#             fitRBFp.zSpec if fitRBFp.zSpec else fitRBFp.zPhotHost,
#             fitRBFp.ccMjdMaxFlux
#             )
# lineRBFp = util.plot.plot_lc_fit(
#             fitRBFp.lcsDict[band].shiftedMjd, fitRBFp.lcsDict[band].flux, fitRBFp.lcsDict[band].fluxErr,
#             pltObj=ax[1]
#             )
# scatterRBFp = util.plot.plot_lc_data(sn.lcsDict[band].shiftedMjd, snFlux, sn.lcsDict[band].fluxErr,
#             pltObj=ax[1])
# lineRBFp.set_label('GP fit, RBF kernel')
# scatterRBFp.set_label('data')
# # ax[1].legend(scatterpoints=1)
#
#
# sn.lcsDict[band].set_shifted_mjd_2(
#             fitRATQUADp.lcsDict[band], fitRATQUADp.r.mjd[fitRATQUADp.r.max_flux_index],
#             fitRATQUADp.zSpec if fitRATQUADp.zSpec else fitRATQUADp.zPhotHost,
#             fitRATQUADp.ccMjdMaxFlux
#             )
# lineRATQUADp = util.plot.plot_lc_fit(
#             fitRATQUADp.lcsDict[band].shiftedMjd, fitRATQUADp.lcsDict[band].flux, fitRATQUADp.lcsDict[band].fluxErr,
#             pltObj=ax[2]
#             )
# scatterRATQUADp = util.plot.plot_lc_data(sn.lcsDict[band].shiftedMjd, snFlux, sn.lcsDict[band].fluxErr,
#             pltObj=ax[2])
# lineRATQUADp.set_label('GP fit, RATQUAD kernel')
# scatterRATQUADp.set_label('data')
# # ax[2].legend(scatterpoints=1)
#
# for i in range(len(ax)):
#     ax[i].set_xlim(xlim)
#     ax[i].set_ylim(ylim)
# ax[2].set_xlabel('Epoch [mjd]')
# ax[1].set_ylabel('Flux [adu]')

plt.show()
