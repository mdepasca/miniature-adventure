import numpy as np
import gc
import argparse

if __name__ == "__main__":
	"""
	The routine works by reading in a file containig in each line the 
	time stamp of the sub-matrix file (WARNING: they have to be in the 
	correct order, 1st line -> sub-matrix 00; 2nd line -> sub-matrix 01 
	etc.).
	The time stamp is inserted in the file name. THIS COULD BETTER BE A LIST OF FILE NAMES!!
	"""


	parser = argparse.ArgumentParser(
        description = "SN lightcurve fitter and classifier: Bulding distance matrix from files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		'--path', '-p', dest='path',
		help='Path where to find files.')

	parser.add_argument(
		'--time-stamps', '-t', dest='timeStamps',
		default='',
		help='ASCII file containing time stamps of sub-matrix files. One per line in the correct order.')

	args = parser.parse_args()

if __name__ == "__main__":
	if args.timeStamps == '':
		raise SystemExit

	inFile = file(args.timeStamps, 'r')
	timeStamps = inFile.readlines()
	inFile.close()

	# path = 'results/SIMGEN_PUBLIC_FIT/RBF_test-length/distance_matrix/'
	fileName = 'dist_matrix_Sum_mothra_{:<14.3f}.txt'

	print 'mat00 ...'
	mat00 = np.loadtxt(args.path+fileName.format(float(timeStamps[0])))

	print 'mat01 ...'
	mat01 = np.loadtxt(args.path+fileName.format(float(timeStamps[1])))

	print 'mat02 ...'
	mat02 = np.loadtxt(args.path+fileName.format(float(timeStamps[2])))

	print 'mat03 ...'
	mat03 = np.loadtxt(args.path+fileName.format(float(timeStamps[3])))

	#---> hstacking ....
	print 'hstacking mat0 ...'
	mat0 = np.hstack((mat00, mat01, mat02, mat03))

	del mat00
	gc.collect()

	# ----- line 1 -----
	print 'mat10 ...'
	mat10 = np.transpose(mat01)
	
	del mat01
	gc.collect()

	print 'mat11 ...'
	mat11 = np.loadtxt(args.path+fileName.format(float(timeStamps[4])))

	print 'mat12 ...'
	mat12 = np.loadtxt(args.path+fileName.format(float(timeStamps[5])))

	print 'mat13 ...'
	mat13 = np.loadtxt(args.path+fileName.format(float(timeStamps[6])))

	#---> hstacking ....
	print 'stacking mat1 ...'
	mat1 = np.hstack((mat10, mat11, mat12, mat13))

	del mat10, mat11, mat12
	gc.collect()

	# ----- line 2 -----
	print 'mat20 ...'
	mat20 = np.transpose(mat02)

	del mat02

	print 'mat21 ...'
	mat21 = np.transpose(mat12)

	del mat12
	fc.collect()

	print 'mat22 ...'
	mat22 = np.loadtxt(args.path+fileName.format(float(timeStamps[7])))

	print 'mat23 ...'
	mat23 = np.loadtxt(args.path+fileName.format(float(timeStamps[8])))

	#---> stacking ....

	print 'stacking mat2 ...'
	mat2 = np.hstack((mat20, mat21, mat22, mat23))

	del mat20, mat21, mat22
	gc.collect()

	# ----- line 3 -----

	print 'mat30 ...'
	mat30 = np.transpose(mat03)

	print 'mat31 ...'
	mat31 = np.transpose(mat13)

	print 'mat32 ...'
	mat32 = np.transpose(mat23)

	del mat03, mat13, mat23
	
	print 'mat33 ...'
	mat33 = np.loadtxt(args.path+fileName.format(float(timeStamps[9])))

	#---> stacking ....
	print 'stacking mat3 ...'
	mat3 = np.hstack((mat30, mat31, mat32, mat33))

	del mat30, mat31, mat32, mat33
	gc.collect()

	#---> stacking ....
	print 'stacking mat ...'
	mat = np.vstack((mat0, mat1, mat2, mat3))

	del mat0, mat1, mat2, mat3
	gc.collect()

	print 'saveing mat ...'
	np.savetxt(args.path+'dist_matrix_Sum.txt', mat)

	del mat
	gc.collect()
	print 'END'
	#END
