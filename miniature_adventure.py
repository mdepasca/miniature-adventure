import utilities as util
import argparse
# Parsing input from command line
parser = argparse.AgrumentParser(
	description = "SN lightcurve fitter and classifier.",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
	"-f", "--fitting", dest="fitting",
	type=bool, action="store_true",
	help="Fit lightcurves")

parser.add_argument(
	"-z", "--zero-point", dest=zeroPoint,
	type=bool, action="store_false",
	help="Set zero point in time. For each object it is the time, in MJD, \
	of maximum observed in r-band (tailored on simulated data from SNANA).")

parser.add_argument(
	"-d", "--distance-metric", dest="distMetric",
	type=bool, action="store_false",
	help="Calculate distance between lightcurves in same band. It is use to \
	build a diffusion map (see Coifman \& Lafon (2006) and Lafon \& Lee \
	(2006)).")