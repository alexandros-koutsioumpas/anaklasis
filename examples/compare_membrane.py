from anaklasis import anaklasis

project='membrane_ref_data_comparison'

input_file = 'membrane.dat' # input curve
units = ['A'] # Q units in Angstrom

patches=[1.0] # single patch 100% covergae

# create system list and append the fronting/backing layers
# and 5 layers SiO2/water/heads/tails/heads
system = []
system.append([
	# Re_sld  Im_sld thk  rough  solv  description 
	[  2.07e-6, 0.0, 0,    1.09, 0.0,  'Si'],
	[  3.5e-6,  0.0, 12.0, 3.50, 0.0,  'SiO2'],
	[  6.15e-6, 0.0, 3.6,  'p0', 1.0,  'D2O'],
	[  1.7e-6,  0.0, 10.5, 'p0', 0.30, 'heads'],
	[ -0.4e-6,  0.0, 30.9, 'p0', 0.03, 'tails'],
	[  1.7e-6,  0.0, 6.3,  'p0', 0.20, 'heads'],
	[  6.15e-6, 0.0,  0,   0.0,  1.0,  'D2O'],
	])

global_param = [
	['p0', 3.6, 'roughness'],
	]

resolution=[-1] # pointwise resolution
background = [5.3e-7] # instrumental background
scale = [1.05] # small scale correction
qmax = [0.25] 

res = anaklasis.compare(project, input_file, units, resolution, 
	patches, system, global_param,background, scale, qmax, 
	experror=True, plot=True)