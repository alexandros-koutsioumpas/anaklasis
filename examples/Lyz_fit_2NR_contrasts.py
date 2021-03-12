# 2-contrast neutron reflectometry fit, from an adsorbed Lysozyme layer at the Si/water interface
from anaklasis import anaklasis

project='Lyz_fit_2_contrasts'
in_file=['Sample6_Si_Lyz_D2O_refl_master.dat','Sample6_Si_Lyz_H2O_refl_master.dat']

units=['A','A'] # all units in Angstrom
fit_mode=0 # 0 is for linear, 1 is for log
fit_weight=[1,1] # equal fit weight for all three curves

method = 'mcmc' # perform MCMC

patches=[1.0]
system = [[
	# Re_sld Im_sld thk rough solv description
	[  2.07e-6, 0.0, 0, 2.0, 0.0, 'Si'],
	[ 3.5e-6, 0.0, 10, 'p0', 0.0, 'SiO2'],
	[ 'm1', 0.0, 'p1', 'p2', 'p3', 'protein'],
	[  'm0', 0.0,  0, 0.0, 1.0, 'bulk'],
	]]

param = [
    # param  min  max  description, for type 'uniform'
    # param  mean sd   description  for type 'normal'
	['p0', 3.0, 3.0, 'SiO2_roughness','uniform'],
	['p1', 25.0, 40.0, 'protein_thk','uniform'],
	['p2', 0, 10.0, 'protein_rough','uniform'],
	['p3', 0.0, 1.0, 'protein_solv','uniform'],
	]

# Note the multiparameter m0 represents the solvent sld that assumes a different 
# set of bounds for each input curve. We leave the sld a bit free to account
# for potential imperefect solvent exchange during the experiment.
# Note that multiparameter m1 represents the sld of the protein. For the two curves
# the sld is different due to the exchange of labile hydrogens
multi_param = [
	# param  min  max  min max  description type
	['m0', 6.0e-6, 6.40e-6, -0.56e-6, 0.5e-6, 'solv_sld','uniform'],
	['m1', 3.45e-6, 3.45e-6,1.98e-6, 1.98e-6, 'solv_sld','uniform']
	]

constraints = [
	]

resolution=[-1,-1] # pointwise resolution
background = [
	[0.0,1.0e-5,'uniform'], # background is left free for all curves
	[0.0,1.0e-5,'uniform'],
	]
scale = [
	[1.0,0.2,'normal'], # curve scale is also a free parameter for all curves
	[1.0,0.2,'normal'],
	]

res = anaklasis.fit(project,in_file,units,fit_mode,fit_weight,method,resolution,patches,system,param
	,multi_param,constraints,background,scale,experror=True,plot=True,fast=True)
