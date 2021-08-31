# This script checks the validity of reflectivity calculations, using the ORSO validation tests
# see: https://github.com/reflectivity/analysis/tree/master/validation/test/unpolarised
from anaklasis import ref
import numpy as np
import os



# ORSO validation tests
print(' ')

# test0
data=np.loadtxt('ORSO_validation_test0.dat', usecols = (0,1),comments = "#")
layers=np.loadtxt('test0.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = np.zeros_like(Q)
tmp = ref.Reflectivity(Q, dq, system, 0, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=8e-5)
  print("ORSO validation test #0 passed!")
except:
  print("ORSO validation test #0 not passed!")
  exit()

# test1
data=np.loadtxt('ORSO_validation_test1.dat', usecols = (0,1),comments = "#")
layers=np.loadtxt('test1.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = np.zeros_like(Q)
tmp = ref.Reflectivity(Q, dq, system, 0, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=8e-5)
  print("ORSO validation test #1 passed!")
except:
  print("ORSO validation test #1 not passed!")
  exit()

# test2
data=np.loadtxt('ORSO_validation_test2.dat', usecols = (0,1),comments = "#")
layers=np.loadtxt('test2.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = np.zeros_like(Q)
tmp = ref.Reflectivity(Q, dq, system, 0, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=8e-5)
  print("ORSO validation test #2 passed!")
except:
  print("ORSO validation test #2 not passed!")
  exit()

# test3
data=np.loadtxt('ORSO_validation_test3.dat', usecols = (0,1),comments = "#")
layers=np.loadtxt('test3.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = np.zeros_like(Q)
tmp = ref.Reflectivity(Q, dq, system, 0, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=8e-5)
  print("ORSO validation test #3 passed!")
except:
  print("ORSO validation test #3 not passed!")
  exit()


# test4
data=np.loadtxt('ORSO_validation_test4.dat', usecols = (0,1,2,3),comments = "#")
layers=np.loadtxt('test0.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = data[:, 3] * 2 * np.sqrt(2 * np.log(2))
tmp = ref.Reflectivity(Q, dq, system, -1, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=0.03)
  print("ORSO validation test #4 passed!")
except:
  print("ORSO validation test #4 not passed!")
  exit()


# test5
data=np.loadtxt('ORSO_validation_test5.dat', usecols = (0,1,2,3),comments = "#")
layers=np.loadtxt('test1.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = data[:, 3] * 2 * np.sqrt(2 * np.log(2))
tmp = ref.Reflectivity(Q, dq, system, -1, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=0.03)
  print("ORSO validation test #5 passed!")
except:
  print("ORSO validation test #5 not passed!")
  exit()

# test6
data=np.loadtxt('ORSO_validation_test6.dat', usecols = (0,1),comments = "#")
layers=np.loadtxt('test6.layers', usecols = (0,1,2,3),comments = "#")
Q = data[:,0]
patches=[1.0]
model=[]
for i in range(layers.shape[0]-1):
       model.append([layers[i,1]*1e-6,layers[i,2]*1e-6,layers[i,0],layers[i+1,3],0,'layer'])
model.append([layers[layers.shape[0]-1,1]*1e-6,layers[layers.shape[0]-1,2]*1e-6,layers[layers.shape[0]-1,0],0,0,'layer'])

system=[model]
dq = np.zeros_like(Q)
tmp = ref.Reflectivity(Q, dq, system, 0, 0, 1, patches, 1) 
try:
  test=np.testing.assert_allclose(tmp[:, 1], data[:, 1], rtol=8e-5)
  print("ORSO validation test #6 passed!")
except:
  print("ORSO validation test #6 not passed!")
  exit()

print(' ')
print('More details on ORSO validation tests can be found here: https://github.com/reflectivity/analysis/tree/master/validation/test/unpolarised' )
print(' ')

