from matplotlib import pyplot as plt
import numpy as np
# import supernova
# import supernova_fit
# import utilities as util
# import classes
def plot_k(snObj, save=False, filename='plot_k.png'):
    lambda_obs = [479.66, 638.26, 776.90, 910.82]

    z = 0.10#(snObj.zSpec if snObj.zSpec else snObj.zPhotHost)
    lambda_rest = [el/(1+z) for el in lambda_obs]

    snrSortIdx = np.argsort(snObj.r.snr)

    bandList = snObj.lcsDict.keys()
    print snrSortIdx
    print [snObj.lcsDict[el].mjd[snrSortIdx[-1]] for el in bandList]
    flux = [snObj.lcsDict[el].flux[snrSortIdx[-1]] for el in bandList]

    gKcor = np.interp(lambda_obs[0], lambda_rest[0:2], flux[0:2])

    fig = plt.figure()
    plt.subplots_adjust(top=0.93, right=0.96, left=0.08)
    plt.suptitle('Source at z={:<3.2f}'.format(z))
    plt.plot([0, lambda_obs[0], lambda_obs[0]], [gKcor, gKcor, 0], figure=fig,
                color='black', linestyle='-.')
    plt.plot(lambda_obs, flux, color='red', linestyle='-',
                marker='o', markeredgecolor='red', markerfacecolor='red',
                label='Obs. rest frame')
    plt.plot(lambda_rest, flux, color='blue', linestyle='--',
                marker='o', markeredgecolor='blue', markerfacecolor='blue',
                figure=fig, label='SN rest frame')
    plt.scatter(lambda_obs[0], gKcor, marker='o',
                edgecolor='black', facecolor='black')
    plt.xlim([np.floor(lambda_obs[0]/100.)*100, np.ceil(lambda_obs[-1]/100.)*100])
    plt.ylim([min(flux)-2, max(flux)+2])

    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Flux (at max S/N)')

    plt.legend(loc='best', numpoints=2, handlelength=3)
    if save:
        plt.savefig(filename, dpi=300, format='png')
    plt.show()

"""
either get the observation file or get a supernova object...

- get a supernova object (independt on storing format!)
- needs:
    .filter's central wavelength
    .to check for observations in all bands at the same *integer* mjd
- blueshifts filters to restfram wavelength
- plot same mjd flux measures (either max flux or max SNR)
- plot lines with dots
- plot observer rest frame central wavelength of each filter
- plot shift direction arrow
- title containing source redshift
- axes labels
"""

if __name__ == "__main__":
    import utilities as util

    snObj = util.get_sn_from_file('train_data/SIMGEN_PUBLIC_DES/DES_SN000017.DAT')
    plot_k(snObj, save=True)
