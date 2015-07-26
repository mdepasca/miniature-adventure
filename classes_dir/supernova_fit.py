import lightcurve
import supernova
class SupernovaFit():
    ccMjdMaxFlux = 0

    def __init__(self, supernova, GPkern=''):
        # check on supernova type HAS TO BE ADDED
        self.g = LightCurve("g")
        self.r = LightCurve("r")
        self.i = LightCurve("i")
        self.z = LightCurve("z")
        self.lcsDict = {"g":self.g,
                        "r":self.r,
                        "i":self.i,
                        "z":self.z}
        self.peaked = False

        if hasattr(supernova, 'kern'):
            self.kern = supernova.kern
        elif GPkern != '':
            self.kern = GPkern

        if hasattr(supernova, 'ccMjdMaxFlux'):
            self.ccMjdMaxFlux = supernova.ccMjdMaxFlux
        self.SNID = supernova.SNID
        self.SNTypeInt = supernova.SNTypeInt
        self.RADeg = supernova.RADeg
        self.decDeg = supernova.decDeg
        self.MWEBV = supernova.MWEBV
        if hasattr(supernova, 'zSpec'):
            self.zSpec = supernova.zSpec
        if hasattr(supernova, 'zSpecErr'):
            self.zSpecErr = supernova.zSpecErr
        else:
            self.zSpecErr = None
        self.hostGalaxyID = supernova.hostGalaxyID
        self.zPhotHost = supernova.zPhotHost
        self.zPhotHostErr = supernova.zPhotHostErr

    def set_lightcurve(self, band, mjd, flux, fluxErr, magFlag=False):
        self.lcsDict[band].mjd = mjd
        self.lcsDict[band].flux = flux
        self.lcsDict[band].fluxErr = fluxErr
        self.lcsDict[band].magFlag = magFlag

        self.lcsDict[band].set_badCurve()

        if not self.lcsDict[band].badCurve:
            if (band == 'r') \
                and self.r.max_flux_index not in set([0, self.r.size-1]):
                    self.peaked = True


    def shift_mjds(self):
        """ Shifts mjd attribute of each lc according to flux maximum
        in r band for peaked lcsself.
        """

        mjdrMax = self.r.mjd[self.r.max_flux_index]

        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                continue
            self.lcsDict[b].set_shifted_mjd(mjdrMax)


    def normalized_flux(self, band):
        """Normalizes the light curve in `band` using the sum of the maximums
        in all that band.
        """

        if not self.lcsDict[band].normFlux:
            den = self.g.max_flux + \
                self.r.max_flux + \
                self.i.max_flux + \
                self.z.max_flux

            flux = self.lcsDict[band].flux
        # result = [flux[idx]/den for idx in range(len(flux))]
            self.lcsDict[band].normFlux = [
                flux[idx]/den for idx in range(len(flux))
                ]
        # return result
        return self.lcsDict[band].normFlux

    def normalized_error(self, band):
        """Normalizes the light curve in band b using the maximum in that band.
        s is a slice on the array.
        """

        if not self.lcsDict[band].normErr:
            den = self.g.max_flux + \
                self.r.max_flux + \
                self.i.max_flux + \
                self.z.max_flux

            fluxErr = self.lcsDict[band].fluxErr
        # result = [fluxErr[idx]/den for idx in range(len(fluxErr))]
            self.lcsDict[band].normErr = [
                fluxErr[idx]/den for idx in range(len(fluxErr))
                ]
        # return result
        return self.lcsDict[band].normErr

    def set_peaked(self):
        self.peaked = True

    @property
    def peaked(self):
        return self.peaked

    def get_distance(self, candidate, band):
        """Calculate difference (aka distance) between two
        interpolated light curves.
        """
        if type(band) is not str:
            raise TypeError("variable `band` is not of type string")
        distFlag = 5

        sizeSelf = self.lcsDict[band].size
        sizeCandidate = candidate.lcsDict[band].size


        if sizeSelf >= sizeCandidate:
            mjd1 = [round(el) for el in self.lcsDict[band].shiftedMjd]
            flux1 = self.normalized_flux(band)
            fluxErr1 = self.normalized_error(band)
            mjd2 = [round(el) for el in candidate.lcsDict[band].shiftedMjd]
            flux2 = candidate.normalized_flux(band)
            fluxErr2 = candidate.normalized_error(band)
        else:
            mjd1 = [round(el) for el in candidate.lcsDict[band].shiftedMjd]
            flux1 = candidate.normalized_flux(band)
            fluxErr1 = candidate.normalized_error(band)
            mjd2 = [round(el) for el in self.lcsDict[band].shiftedMjd]
            flux2 = self.normalized_flux(band)
            fluxErr2 = self.normalized_error(band)

        mjdIntersection = [el for el in mjd1 if el in mjd2]

        if len(mjdIntersection) < 2:
            distance = distFlag
        else:
            flux1Int = [
                flux1[i] for i in [mjd1.index(el) for el in mjdIntersection]
                ]


            flux2Int = [
                flux2[i] for i in [mjd2.index(el) for el in mjdIntersection]
                ]


            fluxErr1Int = [
                fluxErr1[i] for i in [mjd1.index(el) for el in mjdIntersection]
                ]


            fluxErr2Int = [
                fluxErr2[i] for i in [mjd2.index(el) for el in mjdIntersection]
                ]


            num = [
                (el)**2 for el in [
                    flux1Int[i] - flux2Int[i] for i in range(
                        len(mjdIntersection)
                        )
                    ]
                ]

            den = [(fluxErr1Int[i])**2 + (fluxErr2Int[i])**2 for i in range(
                    len(mjdIntersection)
                    )
                ]

            try:
                distance = sqrt(
                    sum([r for r in [num[i]/den[i] for i in range(
                            len(mjdIntersection)
                            )]]
                    )
                    )/(max(mjdIntersection) - min(mjdIntersection))

            except RuntimeWarning:
                print "selfID: {:<d} -- CandidateID {:<d}".format(self.SNID,
                    candidate.SNID)
                print "1: {:<d}".format(id1)
                print "len(num) {:<d}".format(len(num))
                print "len(den) {:<d}".format(len(den))
                print den.index(0)
                print fluxErr1Int[den.index(0)], fluxErr2Int[den.index(0)]
                print fluxErr1.index(0), fluxErr2.index(0)
                print "len(mjdIntersection) {:<d}".format(len(mjdIntersection))
                print distance
                print '--------------------------------------------------------'

        return distance


    def save_on_txt(self, fileName, survey="DES"):
        t = Table(masked=True)
        # OBS is used to reproduce original SNPhotCC files can be deleted
        #
        # provided to change init method of Supernova class
        colNames = [
                    ["OBS", "{:4s}"],
                    ["MJD", "{0:9.3f}"], ["BAND", "{:s}"], ["FIELD", "{:6s}"]
                    ]
        if self.r.magFlag:
            colNames.extend([["MAG", "{0:7.3f}"], ["MAG_ERR", "{0:7.3f}"]])
        else:
            colNames.extend([["FLUX", "{0:10.5f}"], ["FLUX_ERR", "{0:10.5f}"]])

        bandArr = list()
        mjd = list()
        flux = list()
        fluxErr = list()
        mjdArgsort = list()

        for b in self.lcsDict.keys():
            if self.lcsDict[b].badCurve:
                continue

            if len(bandArr) == 0:
                bandArr = [b]*len(self.lcsDict[b].mjd)
                # bandArr = np.empty(len(self.lcsDict[b].mjd), dtype=np.str)
                # bandArr[:] = b
                mjd = self.lcsDict[b].mjd
                flux = self.lcsDict[b].flux
                fluxErr = self.lcsDict[b].fluxErr
            else:
                tmp = [b]*len(self.lcsDict[b].mjd)
                # tmp = np.empty(len(self.lcsDict[b].mjd), dtype=np.str)
                # tmp[:] = b
                bandArr.extend(tmp)
                # bandArr = np.concatenate((bandArr, tmp))
                mjd.extend(self.lcsDict[b].mjd)
                flux.extend(self.lcsDict[b].flux)
                fluxErr.extend(self.lcsDict[b].fluxErr)
                # mjd = np.concatenate((mjd, self.lcsDict[b].mjd))
                # flux = np.concatenate((flux, self.lcsDict[b].flux))
                # fluxErr = np.concatenate((fluxErr, self.lcsDict[b].fluxErr))

        mjdArgsort = np.argsort(mjd)

        mjd = [mjd[i] for i in mjdArgsort]
        flux = [flux[i] for i in mjdArgsort]
        fluxErr = [fluxErr[i] for i in mjdArgsort]
        bandArr = [bandArr[i] for i in mjdArgsort]

        """
        Adding and setting column to Table
        """

        for c in range(len(colNames)):
            if colNames[c][0] == "BAND" \
            or colNames[c][0] == "OBS" \
            or colNames[c][0] == "FIELD":
                col = MaskedColumn(np.zeros(len(mjd)),
                    name=colNames[c][0],
                    format=colNames[c][1],
                    dtype=np.str, fill_value='-',
                    mask=np.zeros(len(mjd))
                    )
            else:
                col = MaskedColumn(np.zeros(len(mjd)),
                    name=colNames[c][0],
                    format=colNames[c][1],
                    dtype=np.float, fill_value=-9,
                    mask=np.zeros(len(mjd)))
            t.add_column(col)

        """
        Initializing columns
        """

        t["OBS"] = np.empty(len(mjd), dtype=np.str)
        t["OBS"][:] = "OBS:"
        t["FIELD"] = np.empty(len(mjd), dtype=np.str)
        t["FIELD"][:] = "NULL"
        t["MJD"] = mjd
        t["BAND"] = bandArr
        t["FLUX"] = flux
        t["FLUX_ERR"] = fluxErr

        t.filled()

        fOut = open(fileName, 'w')
        fOut.write("# File produced by Miniature Adventure on " + \
            "{:<02d}/{:<02d}/{:<4d} at {:<02d}:{:<02d}:{:<02d} GMT\n".format(
                time.gmtime().tm_mday, time.gmtime().tm_mon,
                time.gmtime().tm_year,
                time.gmtime().tm_hour, time.gmtime().tm_min,
                time.gmtime().tm_sec))
        fOut.write("SURVEY:  {:<}\n".format(survey))
        fOut.write("SNID:  {:<d}\n".format(self.SNID))
        fOut.write("GPKERNEL: {:<}\n".format(self.kern))
        # if self.SNTypeInt :
        fOut.write("SNTYPE: {:>d}\n".format(self.SNTypeInt))
        # if self.RADeg :
        fOut.write("RA:     {:>9.6f} deg\n".format(self.RADeg))
        # if self.decDeg :
        fOut.write("DECL:   {:>9.6f} deg\n".format(self.decDeg))
        # if self.MWEBV :
        fOut.write("MWEBV:  {:>6.4f}\n".format(self.MWEBV))
        if hasattr(self, "zSpec"):
            if self.zSpec:
                fOut.write("REDSHIFT_SPEC:  {:>6.4f} +- {:>6.4f}\n".format(
                    self.zSpec, self.zSpecErr
                    ))
            else:
                fOut.write("REDSHIFT_SPEC: -9.0000 +- 9.0000\n")
        if hasattr(self, "hostGalaxyID"):
            fOut.write("HOST_GALAXY_GALID: {:>d}\n".format(self.hostGalaxyID))
        if hasattr(self, "zPhotHost"):
            fOut.write("HOST_GALAXY_PHOTO-Z:  {:>6.4f} +- {:>6.4f}\n".format(
                self.zPhotHost, self.zPhotHostErr
                ))
        # if self.ccMjdMaxFlux != 0:
        fOut.write("MJD_MAX_FLUX-CCF:  {:>9.3f}\n".format(self.ccMjdMaxFlux))
        fOut.write("\n\n\n")
        fOut.write("# ======================================\n")
        fOut.write("# LIGHT CURVE FIT USING GAUSSIAN PROCESS\n")
        fOut.write("#\n")
        fOut.write("# NOBS: {:<}\n".format(len(mjd)))

        ascii.write(t, output=fOut, delimiter='  ',
            format='fixed_width_two_line')

        fOut.close()
