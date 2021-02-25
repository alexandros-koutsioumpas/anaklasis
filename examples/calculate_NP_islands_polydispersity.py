# reflectivity calculations for polydisperse nanoparticle islands
from anaklasis import anaklasis

project='nanoparticle_islands'

# Define different patch coverages that should add up to unity
patches=[0.3,0.7]

# Create empty system list where we will append models
system=[]

# Append first model Si/D2O interface, this will be system[0]
system.append([
	#  Re_sld   Im_sld thk rough solv decription
	[ 2.07e-6,  0.0e-6, 0, 2.0,  0.0, 'Si'],
	[ 6.35e-6,  0.0e-6, 0, 0.0,  1.0, 'D2O'],
	])

# Append second model Si/nanoparticles/D2O, this will be system[1]
# First declare just Si and D2O
system.append([
	#  Re_sld   Im_sld thk rough solv decription
	[ 2.07e-6,  0.0e-6, 0,  2.0, 0.0, 'Si'],
	[ 6.35e-6,  0.0e-6, 0,  0.0, 1.0, 'D2O'],
	])

# Approximate solvent volume fraction integral using a sum.
expr="""Sum(((1/(p2*(2*pi)**0.5))*exp(-0.5*(((p1+(ii/10)*p2)-p1)/p2)**2))
Min(1-(p0-(4*p0/(p1+(ii/10)*p2)**2)*((n-0.5)*(p1/100)-(p1+(ii/10)*p2)/2)**2),1)
*(p2/10),(ii,-30,30))"""

# Define solvent volume fraction for each NP slice. Note that n is the layer number. 
# The product of (n-0.5) x slice thickness gives us the middle z point of each slice.

# insert between Si & D2O, NP layer sliced in 170 slabs.
for i in range(170): 
	# we insert in system[1] because index counting # in Python starts from 0
 	system[1].insert(1,[ 1.41e-6, 0.0e-6, 'p1/100', 0.0, expr, 'NP_layer'])

global_param = [
	['p0', 0.91, 'packing_constant'],
	['p1', 150, 'nanoparticle_diameter'],
	['p2', 40, 'nanoparticle_diameter_sd'],
	]

resolution=[0.05]
background = [1.0e-7]
scale = [1.0]
qmax = [0.3]

res = anaklasis.calculate(project, resolution, patches, 
	system, global_param, 
	background,scale, qmax, plot=True)

# Note that the tuple (res) contains the results in the form of three lists
# res[0] list containg reflectivity (Q, R(Q), R(Q)Q^4)
# res[1][0] list containing the sld profile (z, sld) of the Si/D2O model
# res[1][1] list containing the sld profile (z, sld) of the Si/NP/D2O model
# res[2][0] list containing the solvent profile (z, solv) of the Si/D2O model
# res[2][1] list containing the solvent profile (z, solv) of the Si/NP/D2O model