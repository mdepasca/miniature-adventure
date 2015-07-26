import numpy as np

class LightCurve():
    """
    Once fully initiated, instances of this class have the following important
    properties
    band            (string) 'g','r','i' or 'z'
    lim             (float) brightness threshold for the band (in correct units
                            flux or mag)
    mjd             (array) modified julian dates of observations
    flux            (array) the observed flux
    fluxErr         (array) the error in the observed flux
    shifted_mjd     (array) mjd shifted such that the peak has mjd = 0. To be
                            modified
    magFlag         (boolean) if the flux is expressed in magnitudes
    """
    badCurve = False
    shifted_mjd = np.zeros(0)
    normFlux = list()
    normErr = list()
    """
    TRY TO USE __slots__
    """
    __slots__ = ['band', 'lim', 'mjd', 'shiftedMjd', 'flux', 'fluxErr', 'badCurve',
        'shifted_mjd', 'normFlux', 'normErr', 'magFlag', 'SNR']
    def __init__(self, band, magFlag=False, lim=0):
        self.band = band
        self.mjd = list()#np.zeros(0, dtype=float)
        self.shiftedMjd = list()#np.zeros(0, dtype=float)
        self.flux = list()#np.zeros(0, dtype=float)
        self.fluxErr = list()#np.zeros(0, dtype=float)
        self.SNR = list()
        self.magFlag = magFlag
        self.lim = lim

    @classmethod
    def data(band, mjd, flux, fluxErr):
        self.band = band
        self.mjd = mjd
        self.flux = flux
        self.fluxErr = fluxErr
        self.set_badCurve()

    def set_badCurve(self):

        # this could be generalised as "if min(self.flux) <= self.lim:"
        if len(self.flux) == 0:
            self.badCurve = True
        elif self.magFlag:
            if min(self.flux) <= self.lim:
                self.badCurve = True
        elif max(self.flux) == 0:
                self.badCurve = True

    def set_shifted_mjd(self, distance):
        """
        Construct shiftedMjd, by subtracting 'distance' from 'self.flux'
        """
        self.shiftedMjd = [self.mjd[i]-distance for i in range(len(self.mjd))]

    def set_shifted_mjd_2(self, supernovaFitLC, max_flux_epoch, redshift=-1, ccMjdMaxFlux=0):
        """Initilise `shiftedMjd` on the basis of the fitted light curve
        `SupernovaFitLC`. If redshift is provided shiftedMjd will be also
        corrected from time dilation.
        """

        if supernovaFitLC.band != self.band:
            raise TypeError('Light curves of different bands.')

        if self.badCurve:
            raise TypeError('Operating on a bad light curve.')

        if supernovaFitLC.badCurve:
            raise TypeError('Bad input light curve.')

        if redshift != -1:
            self.shiftedMjd = self.calc_destretched_time(redshift)

        print supernovaFitLC.mjd[supernovaFitLC.max_flux_index]
        self.shiftedMjd = [el - max_flux_epoch# supernovaFitLC.mjd[supernovaFitLC.max_flux_index]
                            for el in self.shiftedMjd]

        if ccMjdMaxFlux != 0:
            epochs = [el + ccMjdMaxFlux for el in self.shiftedMjd]
            self.shiftedMjd = epochs
            del epochs

    def calc_destretched_time(self, redshift):
        """Perform correction for time dialtion. Returns an array of time in
        mjd.

        Keywords arguments
        redshift    -   from a object of type Supernova.
        """
        epochs = [el/(1.+redshift) for el in self.mjd]

        return epochs

    def calc_dereddened_flux(self, R, ebv):
        """Perform correction for MW dust reddening

        Keyword arguments:
        R   -   depends on the band
        ebv -   colour excess, from Supernova object
        """

        a_mag = R * ebv
        a_flux = 10**(0.4*a_mag)

        return [el*a_flux for el in self.flux]

    @property
    def max_flux(self):
        if not self.badCurve:
            result = min(self.flux) if self.magFlag else max(self.flux)
        else:
            result = self.lim

        return result

    @property
    def max_error(self):
        if not self.badCurve:
            result = max(self.fluxErr)
        else:
            result = 0

        return result

    @property
    def max_flux_index(self):
        """
        Return the index of the maximum flux
        """
        # return np.argmax(self.flux)
        return self.flux.index(min(self.flux)) if self.magFlag else self.flux.index(max(self.flux))

    @property
    def size(self):
        # return self.mjd.size
        return len(self.mjd)
