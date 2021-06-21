# This is a template that helps in writing anaklasis scripts
# for defining interfacial models and for performing reflectivity calculations

from anaklasis import ref

project='my_project' #Here put the name of your project

# Define the model. Note that layer elements may be numerical
# values or symbolic expressions (string) involving parameters
model=[
	#  Re_sld  Im_sld   thk rough solv description
	[ 0.00e-6, 0.00e-6,  0 , 0.0, 0.0, 'fronting'],
	[ 3.45e-6, 0.00e-6,'p0',  0.0, 0.0, 'layer'],
	[ 2.07e-6, 0.00e-6,  0 , 0.0, 0.0, 'backing'],
	]

# You need to modify the next two lines only if you want
# to define a mixed area system with multiple models
patches=[1.0]
system=[model]

# Here you can set the model parameters from p0 up to p39
# together with their values and descriptions
global_param = [
	['p0', 50, 'layer_thickness'],
	]

resolution=[0.05] # Here set the dQ/Q instrumental resolution
background = [1.0e-9] # Here set the background
scale = [1.0] # Here set the scale
qmax = [0.7] # Here set the maximum Q value

res = ref.calculate(project, resolution, 
	patches, system, global_param, 
	background, scale, qmax, plot=True)

# Note that the dictionary (res) contains the results
# res[("reflectivity")] Numpy array containing reflectivity ([Q, R(Q), R(Q)Q^4])
# res[("profile")] Numpy array containing the sld profile ([z, sld])
# res[("solvent")] Numpy array containing the solvent profile ([z, solv])
