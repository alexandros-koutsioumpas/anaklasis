# This is a template that helps in writing anaklasis scripts
# for performing refinements of reflectivity data (single curve)
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

# Here you can define the model parameters from p0 up to p39
# together with the their bounds. Parameters can be of type
# uniform of normal.

  # param  min  max  description  type='uniform'
  # param  mean  sd  description  type='normal'
global_param = [
	['p0', 10,  50,  'layer_thickness', 'uniform'],
	]

# No multi-parameters since we have a single input curve
multi_param = []

# constraints are symbolic inequality expressions between
# defined parameters.
constraints = [] # no constraints

resolution=[0.05] # Here set the dQ/Q instrumental resolution
				  # if your data file contains a 4-th column
				  # with dQ/Q, you may set resolution to -1
				  # for pointwise resolution
		      # min max     type 
background = [[0.0,1.0e-6,'uniform']] # Here set the background bounds.

         # mean sd   type
scale = [[0.8,1.2,'uniform']] # Here set the scale bounds.

# Here we set the fit details (Figure of merit, fit weight, and method)
fit_mode=0 # using Chi squared with errors figure of merit
fit_weight=[1] 
method = 'simple' # May be simple, mcmc or bootstrap

# Note that we set experror to True or False depending on the
# presence of a 3rd column in the input file containing dQ errors
res = ref.fit(project, in_file, units, fit_mode,
fit_weight,method,resolution,patches, system,
global_param,multi_param, constraints,
background,scale,experror=True, plot=True,fast=True)

# Note that the dictionary (res) contains the results
# res[("reflectivity")] Numpy array containing reflectivity ([Q, R(Q), R(Q)Q^4])
# res[("profile")] Numpy array containing the sld profile ([z, sld])
# res[("solvent")] Numpy array containing the solvent profile ([z, solv])
# res["chi_square"] float number corresponding to the chi squared
# res["background"] fitted background list containing two floats [mean value, standard deviation (if available)]
# res["scale"] fitted scale list containing two floats [mean value, standard deviation (if available)]
# res[("global_parameters","pi")] fitted global parameter list containing two floats [mean value, standard deviation (if available)]

