class Supernova():
    """
    Has the following properties

    g               (LightCurve)
    r               (LightCurve)
    i               (LightCurve)
    z               (LightCurve)
    lightCurvesDict (dictionary) 4 entries, g,r,i,z returning the corresponding
                                 LightCurves
    SNID            (int)   supernova ID
    SNTypeInt       (int)   supernova type integer (see relationship between
                                                    number and type)
    zSpec           (float) If known via spectroscope, otherwise None
    hostGalaxyID    (int)   THe host galaxy ID (all supernovae (in +zPhotHost)
                                                catalog have this)
    zPhotHost       (float) The redshift of the host galaxy (all supernovae in
                                                        the catalog have this)
    zPhotHostErr    (float) Error in zPhotHost
    """
    """
    TRY TO USE __slots__
    """
    def __init__(self, inFileName, magFlag=False):
        """
        Parses all the light curve data in inFileName into a Supernova object.
        """

        inFile = file(inFileName, "r")
        lines = inFile.readlines()
        inFile.close()

        self.g = LightCurve("g", magFlag=magFlag)
        self.r = LightCurve("r", magFlag=magFlag)
        self.i = LightCurve("i", magFlag=magFlag)
        self.z = LightCurve("z", magFlag=magFlag)

        self.lcsDict = {'g':self.g,
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
                    if magFlag:
                        flux = float(data[6])
                        fluxErr = float(data[7])
                    else:
                        flux = float(data[3])
                        fluxErr = float(data[4])
                    if fluxErr > 0:
                        if passband == "g":
                            self.g.mjd.append(mjd)
                            self.g.flux.append(flux)
                            self.g.fluxErr.append(fluxErr)
                            # self.g.add_data_point(mjd, flux, fluxErr)
                        elif passband == "r":
                            self.r.mjd.append(mjd)
                            self.r.flux.append(flux)
                            self.r.fluxErr.append(fluxErr)
                            # self.r.add_data_point(mjd, flux, fluxErr)
                        elif passband == "i":
                            self.i.mjd.append(mjd)
                            self.i.flux.append(flux)
                            self.i.fluxErr.append(fluxErr)
                            # self.i.add_data_point(mjd, flux, fluxErr)
                        elif passband == "z":
                            self.z.mjd.append(mjd)
                            self.z.flux.append(flux)
                            self.z.fluxErr.append(fluxErr)
                            # self.z.add_data_point(mjd, flux, fluxErr)
                        else:
                            print "Filter not recognized: {:<}".format(passband)
                elif tag == "SNID":
                    self.SNID = int(data[0])
                # SupernovaFit attribute
                elif tag == "GPKERNEL":
                    self.kern = data[0]
                elif tag == "SNTYPE":
                    self.SNTypeInt = int(data[0])
                elif tag == "RA":
                    self.RADeg = float(data[0])
                elif tag == "DECL":
                    self.decDeg = float(data[0])
                elif tag == "MWEBV":
                    self.MWEBV = float(data[0])
                elif (tag == "REDSHIFT_SPEC") or (tag == "REDSHIFT_FINAL"):
                    if float(data[0]) == -9:
                        self.zSpec = None
                    else:
                        self.zSpec = float(data[0])
                        self.zSpecErr = float(data[2])
                elif (tag == "HOST_GALAXY_GALID") or (tag == "HOSTGAL_OBJID"):
                    self.hostGalaxyID = int(data[0])
                elif (tag == "HOST_GALAXY_PHOTO-Z") or (tag == "HOSTGAL_PHOTOZ"):
                    self.zPhotHost = float(data[0])
                    self.zPhotHostErr = float(data[2])
                elif tag == "MJD_MAX_FLUX-CCF":
                    self.ccMjdMaxFlux = float(data[0])

        for b in self.lcsDict.keys():
            self.lcsDict[b].set_badCurve()

    def __cmp__(self, other):
        return 2*(self.zPhotHost - other.zPhotHost > 0) - 1

    def set_shifted_mjd(self, supernovaFitObj, destretchFlag=True):
        """Apply `set_shifted_mjd_2` from Supernova class to all the bands.
        Depending on `destretchFlag` can or cannot perform correction
        for time delay.
        """
        for b in self.lcsDict.keys():
            if destretchFlag:
                self.lcsDict[b].set_shifted_mjd_2(
                    supernovaFitObj.lcsDict[b],
                    (self.zSpec if self.zSpec else self.zPhotHost),
                    supernovaFitObj.ccMjdMaxFlux
                    )
            else:
                self.lcsDict[b].set_shifted_mjd(
                    supernovaFitObj.lcsDict[b])

    def calc_dereddened_flux(self, R):
