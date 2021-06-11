# fit of 7-contrast neutron reflectivity data from a lipid monolayer at the air/water interface 
from anaklasis import ref

project='lipid_monolayer'

# import seven experimental curves
in_file=['hd2o30.dat','d13d2o30.dat','d13acmw30.dat','d70d2o30.dat','d70acmw30.dat',
		'd83d2o30.dat','d83acmw30.dat']

# all units in Angstrom
units=['A','A','A','A','A','A','A']

fit_mode=0 # use linear Figure of Merit

# all curves equal fit weight
fit_weight=[1,1,1,1,1,1,1]

method = 'mcmc' # perform mcmc

# dQ/Q = 10% for all input curves
resolution=[0.1,0.1,0.1,0.1,0.1,0.1,0.1]

model= [
	# Re_sld Im_sld thk rough solv description
	[  0.0e-6, 0.0, 0, 'p0', 0.0, 'Air'],
	[ 'm2*p4/p5', 0.0, 'p3', 'p0', '1.0-p4', 'tails'],
	[ 'm1*((p3*p2*p4)/(p5*p1))/p2', 0.0, 'p1', 'p0', '1.0-(p3*p2*p4)/(p5*p1)', 'heads'],
	[  'm0', 0.0,  0, 0.0, 1.0, 'bulk'],
	]

patches=[1.0] # single uniform layer
system=[model]

# Definition of Global parameters
param = [
    # param  min  max  description
	['p0', 3.0, 5.0, 'roughness','uniform'],
	['p1', 8.0, 16.0, 'heads thickness','uniform'],
	['p2', 300, 380.0, 'head volume','uniform'],
	['p3', 10.0, 26.0, 'tail thickness','uniform'],
	['p4', 1.0, 1.0, 'tail vol fraction','uniform'],
	['p5',800, 1000, 'tail volume','uniform'],
	]

# Definition of multi parameters
# solvent sld, lipid head and tail scattering lengths
multi_param = [
	# param  min  max  min max ... description
	['m0', 6.35e-6, 6.35e-6, 6.35e-6, 6.35e-6, -0.00e-6, 0.00e-6,6.35e-6, 6.35e-6,
	 -0.00e-6, 0.00e-6,6.35e-6, 6.35e-6, -0.00e-6, 0.00e-6, 'solvent_sld','uniform'],
	['m1', 6.01e-4, 6.01e-4, 19.54e-4, 19.54e-4, 19.54e-4, 19.54e-4, 11.21e-4, 11.21e-4,
	 11.21e-4, 11.21e-4, 24.75e-4, 24.75e-4, 24.75e-4, 24.75e-4 ,'b_heads','uniform'],
	['m2', -3.58e-4, -3.58e-4, -3.58e-4, -3.58e-4, -3.58e-4, -3.58e-4, 69.32e-4, 69.32e-4,
	 69.32e-4, 69.32e-4, 69.32e-4, 69.32e-4, 69.32e-4, 69.32e-4 ,'b_tails','uniform'],
	]

# Solvent volume fraction should stay always
# larger than zero
constraints = [
	'(p3*p2*p4)/(p5*p1)<1',
	]

# background as free parameter
background = [
	[0.0e-7,4.0e-6,'uniform'],
	[0.0e-7,4.0e-6,'uniform'],
	[0.0e-7,4.0e-6,'uniform'],
	[0.0e-7,4.0e-6,'uniform'],
	[0.0e-7,4.0e-6,'uniform'],
	[0.0e-7,4.0e-6,'uniform'],
	[0.0e-7,4.0e-6,'uniform'],
	]

# curve scale as a free parameter
scale = [
	[0.10,0.25,'uniform'],
	[0.10,0.25,'uniform'],
	[0.10,0.25,'uniform'],
	[0.10,0.25,'uniform'],
	[0.10,0.25,'uniform'],
	[0.10,0.25,'uniform'],
	[0.10,0.25,'uniform'],
	]

res = ref.fit(project, in_file, units, fit_mode, fit_weight, method, resolution,
			 patches, system, param, multi_param, constraints,background,scale,
			 experror=True,plot=True,fast=True)
