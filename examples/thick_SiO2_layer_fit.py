# Example of neutron reflectivity dataset refinement from a
# thick SiO2 layer at the air/Si interface
# Measurments acquired at D17 instrument (ILL, Grenoble)
from anaklasis import ref

project='SiO_fit'
in_file=['D17_SiO.dat']
units=['A']

fit_mode=0 # 0 is for linear, 1 is for log
fit_weight=[1]
method = 'mcmc' # Markov Chain Monte Carlo Sampling

resolution=[-1] # pointwise resolution

model = [
	# Re_sld Im_sld thk rough solv description
	[  0.0, 0.0,  0, 'p0', 0.0, 'Air'],
	[ 'p1', 0.0, 'p2', 'p3', 0.0, 'SiOx'],
	[  2.07e-6, 0.0, 0, 0.0, 0.0, 'Si'],
	]

system=[model]
patches=[1.0]

global_param = [
    # param  min  max  description
    ['p0', 0, 20, 'air/SiOx_roughness','uniform'],
    ['p1', 3.3e-6, 3.7e-6, 'SiOx_sld','uniform'],
    ['p2', 0, 2000, 'SiOx_thickness','uniform'],
    ['p3', 0, 30, 'SiOx/Si_roughness','uniform'],
	]

multi_param = []
constraints = []

background = [
	[0.0e-11,1.0e-5,'uniform'],
	]

scale = [
	[0.8,1.1,'uniform'],
	]

res = ref.fit(project, in_file, units, fit_mode, fit_weight,method,resolution,patches, system,
global_param,multi_param, constraints, background,scale,experror=True, plot=True,fast=True)