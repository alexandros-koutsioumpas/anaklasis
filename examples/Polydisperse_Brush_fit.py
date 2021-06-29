# Fit of neutron reflectivity data from a polydisperse polymer brush
# this is essentially an extension of the Brush_fit example where we
# also include normally distributed brush length polydispersity in the model.
from anaklasis import ref

project='Polydisperse_Brush_fit'
in_file=['DA34424a.dat']
units=['A'] # Angstrom units
fit_mode=0 # using FOM1
fit_weight=[1] 
method = 'mcmc' # Markov Chain Monte Carlo


# We create model list and first add the semi-
# infinite fronting and the thin water layer
model=[
	# Re_sld Im_sld thk rough solv description
	[  4.32e-6, 0.0, 0, 5, 0.0, 'Quartz'],
	[ -0.56e-6, 0.0, 'p4', 3, 0.0, 'water']
	]

#Brush layer solvent volume fraction expression
#(n-1.5) x slice thickness gives the midpoint
# of each slice. The first slice in layer 2
# since the water layer is in position 1.

# Approximate solvent volume fraction integral using a sum.
# Refer to the related Jupyter-notebook for details
expr="""Sum(((1/(p2*(2*pi)**0.5))*exp(-0.5*(ii/10)**2))*
Min(1-p0+(p0/((p1+(ii/10)*p2)**p3))*((p1/10)*(n-1.5))**p3,1)
*(p2/10),(ii,-40,40))"""

# we add 20 slices for the brush layer with
# thickness brush_length/10 and roughness
# brush_length/20 to smooth the profile
for i in range(20): 
 	model.append([ 1.41e-6, 0.0e-6, 'p1/10', 'p1/20', expr, 'NP_layer'])

# Finally we add the semi-infinite backing (toluene)
model.append([5.86e-6, 0.0,  0, 0, 1.0, 'd-toluene'])

system=[model] # add single model in system list
patches=[1.0] # that completely covers the surface

global_param = [
    # param  min  max  description  type='uniform'
    # param  mean  sd  description  type='normal'
	['p0', 0.02, 0.20, 'phi0','uniform'],
	['p1', 100, 700, 'brush_length','uniform'],
	['p2', 0, 100, 'brush_length_sd','uniform'],
	['p3', 1, 4, 'exponent','uniform'],
	['p4', 5, 1, 'water thickness','normal'],
	]

multi_param = [] # no multi-parameters
constraints = ['4*p2<p1'] # Keep brush length bigger than 4 times the brush length standard deviation
resolution=[0.06] # dQ/Q=6%
background = [[0.0,0.0,'uniform']] #background corrected data
scale = [[1.0,1.0,'uniform']] #data correctly scaled

res = ref.fit(project, in_file, units, fit_mode, fit_weight,method,resolution,patches, system,
global_param,multi_param, constraints, background,scale,experror=True, plot=True,fast=True)

