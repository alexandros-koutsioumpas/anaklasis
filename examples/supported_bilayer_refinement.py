# 3-contrast neutron reflectometry fit, for a supported lipid bilayer at the Si/water interface
from anaklasis import ref

project='supported_bilayer_refinement'
in_file=['bilayer_D2O.dat','bilayer_SMW.dat','bilayer_H2O.dat']

units=['A','A','A'] # all units in Angstrom
fit_mode=0 # 0 is for linear, 1 is for log
fit_weight=[1,1,1] # equal fit weight for all three curves

method = 'mcmc' # perform MCMC

model = [
	# Re_sld Im_sld thk rough solv description
	[  2.07e-6, 0.0, 0, 2.0, 0.0, 'Si'],
	[ 3.5e-6, 0.0, 10, 'p0', 0.0, 'SiO2'],
	[ 'm0', 0.0, 'p1', 'p0', 1.0, 'water'],
	[ 'p5/(p9*p3)', 0.0, 'p3', 'p0', '1.0-p7/(p9*p3)', 'inner_heads'],
	[ 'p6/(p9*p2)', 0.0, 'p2', 'p0', '1.0-p8/(p9*p2)', 'inner_tails'],
	[ 'p6/(p9*p2)', 0.0, 'p2', 'p0', '1.0-p8/(p9*p2)', 'outer_tails'],
	[ 'p5/(p9*p4)', 0.0, 'p4', 'p0', '1.0-p7/(p9*p4)', 'outer_heads'],
	[  'm0', 0.0,  0, 0.0, 1.0, 'bulk'],
	]

patches=[1.0]
system=[model]

param = [
    # param  min  max  description, for type 'uniform'
    # param  mean sd   description  for type 'normal'
	['p0', 2.0, 0.5, 'roughness','normal'],
	['p1', 3.0, 15.0, 'water_d','uniform'],
	['p2', 7.0, 20.0, 'tail_d','uniform'],
	['p3', 5.0, 15.0, 'inner_head_d','uniform'],
	['p4', 5.0, 15.0, 'outer_head_d','uniform'],
	['p5', 6.01e-4, 6.01e-4, 'b_heads','uniform'],
	['p6', -2.92e-4, -2.92e-4, 'b_tails','uniform'],
	['p7', 319, 319, 'V_heads','uniform'],
	['p8', 782, 782, 'V_tails','uniform'],
	['p9', 40, 70, 'area_lipid','uniform'],
	]

# Note the multiparameter m0 represents the solvent sld that assumes a different 
# set of bounds for each input curve.
multi_param = [
	# param  min  max  min max min max description type
	['m0', 6.15e-6, 6.40e-6, 1.80e-6, 2.30e-6, -0.56e-6, 0.0e-6, 'solv_sld','uniform']
	]

constraints = [
	'1.0-p8/(p9*p2)>0', # With these three constaints we avoid that the solvent volume
	'1.0-p7/(p9*p4)>0', # fraction might become negative for the inner_head, inner_tail
	'1.0-p7/(p9*p3)>0', # outer_tail ot outer_head layer.
	]

resolution=[-1,-1,-1] # pointwise resolution
background = [
	[0.0,1.0e-5,'uniform'], # background is left free for all three curves
	[0.0,1.0e-5,'uniform'],
	[0.0,1.0e-5,'uniform'],
	]
scale = [
	[1.0,0.2,'normal'], # curve scale is also a free parmater for all curves
	[1.0,0.2,'normal'],
	[1.0,0.2,'normal'],
	]

res = ref.fit(project,in_file,units,fit_mode,fit_weight,method,resolution,patches,system,param
	,multi_param,constraints,background,scale,experror=True,plot=True,fast=True)
