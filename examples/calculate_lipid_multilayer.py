# perform neutron reflectivity calculation for a lipid multilayer
from anaklasis import ref

project='lipid_multilayer'

# Create model list where 
# we first include only the fronting Si medium
# and the SiO2 and water thin layer on Si
model=[
	 # Re_sld Im-sld thk rough solv description
	[  2.07e-6, 0.0, 0,  'p0', 0.0, 'Si'],
	[  3.50e-6, 0.0, 10, 'p0', 0.1, 'SiO2'],
	[  6.35e-6, 0.0, 5,  'p0', 1.0, 'D2O thin layer'],
	]

# "for loop" for appending bilayers on top of Si/SiO2/thin water layers
for i in range(10): # i runs from 0 to 9
                 # Re_sld Im-sld thk rough solv  description
	model.append([ 1.7e-6,  0.0, 10, 'p0', 0.3, 'inner_heads'+str(i+1)])
	model.append([ -0.5e-6, 0.0, 35, 'p0', 0.0, 'tails'+str(i+1)])
	model.append([ 1.7e-6,  0.0, 10, 'p0', 0.3, 'outer_heads'+str(i+1)])
	model.append([ 6.35e-6, 0.0, 15, 'p0', 1.0, 'D2O gap'+str(i+1)])

# Finally add the D2O semi-infinite backcking
model.append([  6.35e-6, 0.0, 0,  0.0,  1.0, 'D2O bulk'])

patches = [1.0] # single laterally uniform layer
system=[model] # we have a single model

global_param = [
	['p0', 2.5, 'roughness'],
	]

resolution=[0.012]
background = [5e-7]
scale = [1.0]
qmax = [0.3]

# Before calling the calculate function, one may use print(model) in order
# to check model definition
print(model)

res = ref.calculate(project,resolution,patches,system,global_param,background,scale,qmax, True)

# Note that the dictionary (res) contains the results
# res[("reflectivity")] Numpy array containing reflectivity ([Q, R(Q), R(Q)Q^4])
# res[("profile")] Numpy array containing the sld profile ([z, sld])
# res[("solvent")] Numpy array containing the solvent profile ([z, solv])