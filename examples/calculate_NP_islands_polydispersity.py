# reflectivity calculations for polydisperse nanoparticle islands
from anaklasis import ref

project='nanoparticle_islands'

# model0 Si/D2O interface
model0=[
	#  Re_sld   Im_sld thk rough solv decription
	[ 2.07e-6,  0.0e-6, 0, 2.0,  0.0, 'Si'],
	[ 6.35e-6,  0.0e-6, 0, 0.0,  1.0, 'D2O'],
	]

# model1 Si/nanoparticles/D2O
# First we declare just Si semi-infinite fronting
model1=[
	#  Re_sld   Im_sld thk rough solv decription
	[ 2.07e-6,  0.0e-6, 0,  2.0, 0.0, 'Si'],
	]

# Approximate solvent volume fraction integral using a sum.
expr="""Sum(((1/(p2*(2*pi)**0.5))*exp(-0.5*(((p1+(ii/10)*p2)-p1)/p2)**2))*
Min(1-(4*p0/(p1+(ii/10)*p2)**2)*((p1+(ii/10)*p2)*(n-0.5)*(p1/100)-((n-0.5)*(p1/100))**2),1)
*(p2/10),(ii,-30,30))"""

# Define solvent volume fraction for each NP slice. Note that n is the layer number. 
# The product of (n-0.5) x slice thickness gives us the middle z point of each slice.

# append NP layer sliced in 170 slabs.
for i in range(170): 
 	model1.append([ 1.41e-6, 0.0e-6, 'p1/100', 0.0, expr, 'NP_layer'])

# Finally append D2O semi-infinite backing in model1
model1.append([ 6.35e-6,  0.0e-6, 0,  0.0, 1.0, 'D2O'])

# Define different patch coverages that should add up to unity
patches=[0.3,0.7]
system=[model0,model1] # Note we have two models

global_param = [
	['p0', 0.91, 'packing_constant'],
	['p1', 150, 'nanoparticle_diameter'],
	['p2', 40, 'nanoparticle_diameter_sd'],
	]

resolution=[0.05]
background = [1.0e-7]
scale = [1.0]
qmax = [0.3]

res = ref.calculate(project, resolution, patches, 
	system, global_param, 
	background,scale, qmax, plot=True)

# Note that the dictionary (res) contains the results
# res[("reflectivity")] Numpy array containing reflectivity ([Q, R(Q), R(Q)Q^4])
# res[("profile","model0")] Numpy array containing the sld profile ([z, sld]) of the Si/D2O model
# res[("solvent","model0")] Numpy array containing the solvent profile ([z, solv]) of the Si/D2O model
# res[("profile","model1")] Numpy array containing the sld profile ([z, sld]) of the Si/NP/D2O model
# res[("solvent","model1")] Numpy array containing the solvent profile ([z, solv]) of the Si/NP/D2O model