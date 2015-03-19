import numpy as np

print 'mat00 ...'
mat00 = np.loadtxt('products/distance_matrix/Oct28/dist_matrix_Sum_pc017261.ads.eso.org_1414549460.544.txt')

print 'mat01 ...'
mat01 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414592263.858.txt')

print 'mat02 ...'
mat02 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414591966.764.txt')

print 'mat03 ...'
mat03 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414590774.128.txt')

#---> hstacking ....
print 'hstacking mat0 ...'
mat0 = np.hstack((mat00, mat01, mat02, mat03))

# ----- line 1 -----
print 'mat10 ...'
mat10 = np.transpose(mat01)

print 'mat11 ...'
mat11 = np.loadtxt('products/distance_matrix/Oct28/dist_matrix_Sum_pc017261.ads.eso.org_1414547556.833.txt')

print 'mat12 ...'
mat12 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414590827.558.txt')

print 'mat13_00 ...'
mat13_00 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414612211.194.txt')

print 'mat13_01 ...'
mat13_01 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414612669.116.txt')

print 'mat13_10 ...'
mat13_10 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414612742.041.txt')

print 'mat13_11 ...'
mat13_11 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414612592.755.txt')

#---> stacking ....
print 'stacking mat13_0 ...'
mat13_0 = np.hstack((mat13_00, mat13_01))

print 'stacking mat13_1 ...'
mat13_1 = np.hstack((mat13_10, mat13_11))

print 'stacking mat13 ...'
mat13 = np.vstack((mat13_0, mat13_1))

#---> hstacking ....
print 'stacking mat1 ...'
mat1 = np.hstack((mat10, mat11, mat12, mat13))

# ----- line 2 -----
print 'mat20 ...'
mat20 = np.transpose(mat02)

print 'mat21 ...'
mat21 = np.transpose(mat12)

print 'mat22 ...'
mat22 = np.loadtxt('products/distance_matrix/Oct28/dist_matrix_Sum_pc017261.ads.eso.org_1414547515.958.txt')

print 'mat23_00 ...'
mat23_00 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414612861.709.txt')

print 'mat23_01 ...'
mat23_01 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414613331.182.txt')

print 'mat23_10 ...'
mat23_10 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414614055.800.txt')

print 'mat23_11 ...'
mat23_11 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_mothra_1414613381.914.txt')

#---> stacking ....
print 'stacking mat23_0 ...'
mat23_0 = np.hstack((mat23_00, mat23_01))

print 'stacking mat23_1 ...'
mat23_1 = np.hstack((mat23_10, mat23_11))

print 'stacking mat23 ...'
mat23 = np.vstack((mat23_0, mat23_1))

#---> stacking ....

print 'stacking mat2 ...'
mat2 = np.hstack((mat20, mat21, mat22, mat23))


# ----- line 3 -----

print 'mat30 ...'
mat30 = np.transpose(mat03)

print 'mat31 ...'
mat31 = np.transpose(mat13)

print 'mat32 ...'
mat32 = np.transpose(mat23)

print 'mat33_00 ...'
mat33_00 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_pc017261.ads.eso.org_1414587633.166.txt')

print 'mat33_01 ...'
mat33_01 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_pc017261.ads.eso.org_1414596438.697.txt')

print 'mat33_10 ...'
mat33_10 = np.transpose(mat33_01)

print 'mat33_11 ...'
mat33_11 = np.loadtxt('products/distance_matrix/dist_matrix_Sum_pc017261.ads.eso.org_1414587956.411.txt')

#---> stacking ....
print 'stacking mat33_0 ...'
mat33_0 = np.hstack((mat33_00, mat33_01))

print 'stacking mat33_1 ...'
mat33_1 = np.hstack((mat33_10, mat33_11))

print 'stacking mat33 ...'
mat33 = np.vstack((mat33_0, mat33_1))

#---> stacking ....
print 'stacking mat3 ...'
mat3 = np.hstack((mat30, mat31, mat32, mat33))

#---> stacking ....
print 'stacking mat ...'
mat = np.vstack((mat0, mat1, mat2, mat3))

print 'saveing mat ...'
np.savetxt('products/distance_matrix/dist_matrix_Sum.txt', mat)

print 'END'
#END
