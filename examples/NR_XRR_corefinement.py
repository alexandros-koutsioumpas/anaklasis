# Co-refienement of XRR and polarised NR data from a MgO/Fe/Pt system
from anaklasis import ref

project='XRR_fit'
in_file=[]
in_file.append('Fe_Pt-XRR.dat') # XRR curve
in_file.append('Fe_Pt_UP.dat')  # NR curve spin up
in_file.append('Fe_Pt_DN.dat')  # NR curve spin down

units=['A','A','A'] # all data in Angstrom units

fit_mode=0 # linear FOM
fit_weight=[1,1,1] # equal fit weight for all curves

method = 'simple' # no MCMC or Bootstrap

model = [
	# Re_sld Im_sld thk rough solv description
	[  0.0e-6, 0.0, 0, 'p1', 0.0, 'air'],
	[ 'm0', 'm1', 'p2', 'p3', 0.0, 'Pt'],
	[ 'm2+m4*p0', 'm3', 'p4', 'p5', 0.0, 'Fe'],
	[ 'm5', 'm6',  0, 0.0, 0.0, 'MgO'],
	]

patches=[1.0]
system=[model]

param = [
    # param  min  max  description type
    ['p0', 4.5e-6, 4.5e-6, 'Fe- mag sld','uniform'],
    ['p1', 0, 6, 'air/Pt roughness','uniform'],
    ['p2', 20, 30, 'Pt thickness','uniform'],
    ['p3', 0, 6, 'Pt/Fe roughness','uniform'],
    ['p4', 200, 250, 'Fe thickness','uniform'],
    ['p5', 0, 6, 'Fe/MgO roughness','uniform'],
	]

multi_param = [
	# param  min_XRR  max_XRR  min_NR_up max_NR_up min_NR_down max_NR_down  description type
	['m0', 13.7e-5, 13.7e-5, 6.35e-6, 6.35e-6, 6.35e-6, 6.35e-6, 'Pt-Re(sld)','uniform'],
	['m1', 1.35e-5, 1.35e-5, 0.0, 0.0, 0.0, 0.0, 'Pt-Im(sld)','uniform'],
	['m2', 5.94e-5, 5.94e-5, 8.024e-6, 8.024e-6, 8.024e-6, 8.024e-6, 'Fe-nuclear Re(sld)','uniform'], 
	['m3', 7.69e-6, 7.69e-6, 0.0, 0.0, 0.0, 0.0, 'Fe-nuclear Im(sld)','uniform'],
	['m4', 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 'mag_sld contribution','uniform'],
	['m5', 3.06e-5, 3.06e-5, 6.01e-6, 6.01e-6, 6.01e-6, 6.01e-6, 'MgO-Re(sld)','uniform'],  
	['m6', 3.26e-7, 3.26e-7, 0.0, 0.0, 0.0, 0.0, 'MgO-Im(sld)','uniform'],
	]

constraints = [
	]

resolution=[0.0, 0.12, 0.12] # dQ/Q for XRR and the two NR curves

background = [
	[0.0e-11,1.0e-5,'uniform'],
	[0.0e-11,1.0e-5,'uniform'],
	[0.0e-11,1.0e-5,'uniform'],
	]

scale = [
	[1.5,2.5,'uniform'],
	[1.00,1.00,'uniform'],
	[1.00,1.00,'uniform'],
	]

res = ref.fit(project, in_file, units, fit_mode,fit_weight, method, resolution, patches, system, param, multi_param, constraints,background,scale,experror=False, plot=True,fast=False)

