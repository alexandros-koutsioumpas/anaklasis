# This is a template that helps in writing anaklasis scripts
# for defining interfacial models and for performing comparions of
# theoretical reflectivity curves with actual data

from anaklasis import ref

project='my_project' #Here put the name of your project

input_file = 'my_data.dat' # enter the filename of input data
units = ['A'] # set units of the input file, Angstrom (A) or nanometers (nm)

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
				  # if your data file contains a 4-th column
				  # with dQ/Q, you may set resolution to -1
				  # for pointwise resolution
background = [1.0e-7] # Here set the background
scale = [1.0] # Here set the scale
qmax = [0.25] # Here set the maximum Q value

# Note that we set experror to True or False depending on the
# presence of a 3rd column in the input file containting dQ errors
res = ref.compare(project, input_file, units, resolution, 
	patches, system, global_param,background, scale, qmax, 
	experror=True, plot=True)

# Note that the dictionary (res) contains the results
# res[("reflectivity")] Numpy array containing reflectivity ([Q, R(Q), R(Q)Q^4])
# res[("profile")] Numpy array containing the sld profile ([z, sld])
# res[("solvent")] Numpy array containing the solvent profile ([z, solv])
# res["chi_square"] float number corresponding to the chi squared
