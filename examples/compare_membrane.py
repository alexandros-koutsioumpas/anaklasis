from anaklasis import anaklasis

project='membrane_ref_data_comparison'

input_file = ['membrane.dat'] # input curve
units = ['A'] # Q units in Angstrom

patches=[1.0] # single patch 100% covergae

# create system list and append the fronting/backing layers
# and 5 layers SiO2/water/heads/tails/heads
system = []
system.append([
	# Re_sld  Im_sld thk  rough  solv  description 
	[  2.07e-6, 0.0, 0,    0.89, 0.0,  'Si'],
	[  3.5e-6,  0.0, 12.6, 3.82, 0.0,  'SiO2'],
	[  6.15e-6, 0.0, 4.4,  'p0', 1.0,  'D2O'],
	[  1.7e-6,  0.0, 10.6, 'p0', 0.24, 'heads'],
	[ -0.4e-6,  0.0, 29.2, 'p0', 0.03, 'tails'],
	[  1.7e-6,  0.0, 9.0,  'p0', 0.47, 'heads'],
	[  6.15e-6, 0.0,  0,   0.0,  1.0,  'D2O'],
	])

global_param = [
	['p0', 2.67, 'roughness'],
	]

resolution=[0.1] # dQ/Q=10%
background = [1.9e-7] # instrumental background
scale = [0.978] # small scale correction
qmax = [0.25] 

res = anaklasis.compare(project, input_file, units, resolution, 
	patches, system, global_param,background, scale, qmax, 
	experror=True, plot=True)
