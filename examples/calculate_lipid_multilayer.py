# perform neutron reflectivity calculation for a lipid multilayer
from anaklasis import anaklasis

project='lipid_multilayer'

patches = [1.0] # single laterally uniform layer

system=[] # Create empty system list

# Append model to system list
# model numbering starts from 0 so this will be system[0]
system.append([
	 # Re_sld Im-sld thk rough solv description
	[  2.07e-6, 0.0, 0,  'p0', 0.0, 'Si'],
	[  3.50e-6, 0.0, 10, 'p0', 0.1, 'SiO2'],
	[  6.35e-6, 0.0, 5,  'p0', 1.0, 'D2O thin layer'],
	[  6.35e-6, 0.0, 0,  0.0,  1.0, 'D2O bulk'],
	])

# "for loop" for inseting bilayers between layer #2 and layer #3
# note that layer numbering starts from 0
for i in range(10): # i runs from 0 to 9
                           # Re_sld Im-sld thk rough solv  description
	system[0].insert(3+4*i,[ 1.7e-6,  0.0, 10, 'p0', 0.3, 'inner_heads'+str(i+1)])
	system[0].insert(4+4*i,[ -0.5e-6, 0.0, 35, 'p0', 0.0, 'tails'+str(i+1)])
	system[0].insert(5+4*i,[ 1.7e-6,  0.0, 10, 'p0', 0.3, 'outer_heads'+str(i+1)])
	system[0].insert(6+4*i,[ 6.35e-6, 0.0, 15, 'p0', 1.0, 'D2O gap'+str(i+1)])

global_param = [
	['p0', 2.5, 'roughness'],
	]

resolution=[0.012]
background = [5e-7]
scale = [1.0]
qmax = [0.3]

# Before calling the calculate function, one may use print(system[0]) in order
# to check model definition
print(system[0])

res = anaklasis.calculate(project,resolution,patches,system,global_param,background,scale,qmax, True)

# Note that the tuple (res) contains the results in the form of three lists
# res[0] list containg reflectivity (Q, R(Q), R(Q)Q^4)
# res[1][0] containing the sld profile (z, sld)
# res[2][0] containing the solvent profile (z, solv)