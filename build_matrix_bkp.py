import numpy as np
import gc

path = 'results/SIMGEN_PUBLIC_FIT/RBF_test-length/distance_matrix/'
fileName = 'dist_matrix_Sum_mothra_{:<5.3f}.txt'

print 'mat00 ...'
mat00 = np.loadtxt(path+'dist_matrix_Sum_mothra_1424728136.080.txt')

print 'mat01 ...'
mat01 = np.loadtxt(path+'dist_matrix_Sum_mothra_1424761931.836.txt')

print 'mat02 ...'
mat02 = np.loadtxt(path+'dist_matrix_Sum_mothra_1424836678.486.txt')

print 'mat03 ...'
mat03 = np.loadtxt(path+'dist_matrix_Sum_mothra_1426935010.715.txt')#dist_matrix_Sum_mothra_1424835453.155.txt')

#---> hstacking ....
print 'hstacking mat0 ...'
mat0 = np.hstack((mat00, mat01, mat02, mat03))

del mat00, mat01, mat02, mat03
gc.collect()

print 'mat10 ...'
mat10 = np.transpose(mat01)
# ----- line 1 -----

print 'mat11 ...'
mat11 = np.loadtxt(path+'dist_matrix_Sum_mothra_1424976703.216.txt')

print 'mat12 ...'
mat12 = np.loadtxt(path+'dist_matrix_Sum_mothra_1425008330.511.txt')

print 'mat13 ...'
mat13 = np.loadtxt(path+'dist_matrix_Sum_mothra_1426934839.968.txt')#dist_matrix_Sum_mothra_1425030478.363.txt')

#---> hstacking ....
print 'stacking mat1 ...'
mat1 = np.hstack((mat10, mat11, mat12, mat13))

del mat10, mat11, mat12, mat13
gc.collect()

# ----- line 2 -----
print 'mat20 ...'
mat20 = np.transpose(mat02)

print 'mat21 ...'
mat21 = np.transpose(mat12)

print 'mat22 ...'
mat22 = np.loadtxt(path+'dist_matrix_Sum_mothra_1424934904.523.txt')

print 'mat23 ...'
mat23 = np.loadtxt(path+'dist_matrix_Sum_mothra_1426935152.075.txt')#dist_matrix_Sum_mothra_1425005559.641.txt')

#---> stacking ....

print 'stacking mat2 ...'
mat2 = np.hstack((mat20, mat21, mat22, mat23))

del mat20, mat21, mat22, mat23
gc.collect()

# ----- line 3 -----

print 'mat30 ...'
mat30 = np.transpose(mat03)

print 'mat31 ...'
mat31 = np.transpose(mat13)

print 'mat32 ...'
mat32 = np.transpose(mat23)

print 'mat33 ...'
mat33 = np.loadtxt(path+'dist_matrix_Sum_mothra_1424972850.416.txt')

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
np.savetxt(path+'dist_matrix_Sum.txt', mat)

del mat
gc.collect()
print 'END'
#END
