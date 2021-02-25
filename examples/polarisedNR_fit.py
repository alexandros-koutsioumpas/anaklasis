from anaklasis import anaklasis

project='PNR_FeNi'
in_file=[]
in_file.append('FeNi_Au_ud.dat')
in_file.append('FeNi_Au_dd.dat')

units=['A','A'] # all input Q in inverse Angstrom

fit_mode=0 # 0 is linear, 1 is log
fit_weight=[1.0,1.0] # equal fit weiht for both curves

method = 'simple' # can be 'simple', or 'mcmc' or 'bootstrap'

patches=[1.0]

# Here we define the model Si/SiO2/FeNi/Au/D2O
# note the expression for the sld of the FeNi layer
# nuclear sld plus or minus magnetic sld
system = [[
	# Re_sld Im_sld thk rough solv description
	[  2.07e-6, 0.0, 0, 2.0, 0.0, 'Si'],
	[ 3.5e-6, 0.0, 12, 2.0, 0.0, 'SiO2'],
	[ 'p0+m0*p1', 0.0, 'p2', 'p4', 0.0, 'FeNi'],
	[ 4.66e-6, 0.0, 'p3', 'p5', 0.0, 'Au'],
	[  'p6', 0.0,  0, 0.0, 1.0, 'D2O'],
	]]
	
model_param = [
    # param  min  max  description type
	['p0', 7.57e-6, 9.57e-6, 'nuclear_sld','uniform'],
	['p1', 1.5e-6, 5e-6, 'magnetic_sld','uniform'],
	['p2',100,700,'FeNi_thickness','uniform'],
	['p3',50,150,'Au_thickness','uniform'],
	['p4',0,20,'FeNi_roughness','uniform'],
	['p5',0,20,'Au_roughness','uniform'],
	['p6',6.0e-6,6.4e-6,'solvent_sld','uniform'],
	]

# we define a multiparameter (m0) that assumes ge value -1 or +1
# for subtracting or adding the magnetic sld to the total
# sld depending on the polarsation of the incoming beam
multi_param = [
	# param  min  max  min max ... description
	['m0',  -1.0, -1.0, 1.0, 1.0, 'up/down','uniform']
	]

constraints = [] # no constraints
resolution=[0.1,0.1] # dQ/Q=10%

# background corrected data so evrything set to zero
background = [
	[0.0,0.0,'uniform'],
	[0.0,0.0,'uniform'],
	]
# scaling of the curve left slightly free
scale = [
	[0.9,1.1,'uniform'],
	[0.9,1.1,'uniform'],
	]

res = anaklasis.fit(project, in_file, units, fit_mode, fit_weight, method, resolution, patches, system, 
	model_param, multi_param, constraints,background,scale,experror=False,plot=True,fast=False)