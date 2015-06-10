try:
    from matplotlib import pyplot as plt
except:
    pass

def plot_lc_fit(epoch, flux, fluxErr, pltObj=None):
    """
    Plots interpolated lightcurve including shaded regions up to 3 sigma

    INPUT
    epoch   -   time array
    flux    -   flux array
    fluxErr -   error on flux array
    pltObj  -   matplotlib object where the plot has to be done, can be
                an `axis` or `None`. If `None` pyplot methods are called
                directly from pyplot
    """

    shadeLevel = [1, 2, 3]
    shade = [0.4, 0.2, 0.1]

    for i in range(len(shadeLevel)):

        fluxUpLim = [val for val in [flux[el] + shadeLevel[i]*fluxErr[el]
                for el in range(len(flux))
            ]]
        fluxLowLim = [val for val in [flux[el] - shadeLevel[i]*fluxErr[el]
                for el in range(len(flux))
            ]]

        if (pltObj == None):
            plt.fill_between(epoch, fluxUpLim, fluxLowLim,
                facecolor='red', alpha=shade[i], linewidth=0.5)
        else:
            pltObj.fill_between(epoch, fluxUpLim, fluxLowLim,
                facecolor='red', alpha=shade[i], linewidth=0.5)

        if (pltObj == None):
            line, = plt.plot(epoch, flux, color='#7f0000', linewidth=2, zorder=1)
        else:
            line, = pltObj.plot(epoch, flux, color='#7f0000', linewidth=2, zorder=1)

    return line

def plot_lc_data(epoch, flux, fluxErr, pltObj=None):
    """
    Scatter plot of light curve observation, with error bar on top.

    INPUT
    epoch   -   time array
    flux    -   flux array
    fluxErr -   error on flux array
    pltObj  -   matplotlib object where the plot has to be done, can be
                an `axis` or `None`. If `None` pyplot methods are called
                directly from pyplot
    """

    if (pltObj == None):
        scatter = plt.scatter(epoch, flux, s=30, c='black', marker='x', zorder=2)
    else:
        scatter = pltObj.scatter(epoch, flux, s=30, c='black', marker='x', zorder=2)

    if (pltObj == None):
        plt.errorbar(epoch, flux,
            fluxErr, fmt=None, color='black', ecolor='black', zorder=2)
    else:
        pltObj.errorbar(epoch, flux,
            fluxErr, fmt=None, color='black', ecolor='black', zorder=2)

    return scatter
