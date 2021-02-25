# perform neutron reflectivity calculations for a bimodal polymer brush
from anaklasis import anaklasis

project='Bimodal_Brush'

patches = [1.0] # single laterally uniform layer

# Create empty system list and Append model
# model numbering starts from 0 so this will be system[0]
system=[]
system.append([
	 # Re_sld Im-sld thk rough solv description
	[  2.07e-6, 0.0,  0,  3,  0.0, 'Si'],
	[  3.5e-6,  0.0, 10,  3,  0.0, 'SiO2'],
	[  'p0',    0.0,  0,  0,  1.0, 'D2O'],
	])

# We define the expression for the solvent volume fraction within the brush layer
expr='1-Max(p1-(p1/(p2**p5))*((p4/50.0)*(n-1.5))**p5,0)-p3+(p3/(p4**p6))*((p4/50.0)*(n-1.5))**p6'

# we insert 50 "brush" slabs after the SiO2 layer.
                                      # Re_sld Im-sld  thk       rough      solv  description
for i in range(50): system[0].insert(2,[ 'p7', 0.0, 'p4/50.0', 'p4/100.0', expr, 'brush_layer'])

global_param = [
    # name value  description
	['p0', 6.35e-6, 'solvent sld'],
	['p1', 0.4, 'brush1_phi0'],
	['p2', 400, 'brush1_length'],
	['p3', 0.2, 'brush2_phi0'],
	['p4', 900, 'brush2_length'],
 	['p5', 2, 'exponent_brush1'],
  	['p6', 2, 'exponent_brush2'],
  	['p7', 0.7e-6, 'polymer_sld'],
	]

resolution=[0.07] # dQ/Q=7%
background = [0.0] # no background
scale = [1.0]
qmax = [0.1]


res = anaklasis.calculate(project,resolution,patches,system,global_param,background,scale,qmax,True)

# Note that the tuple (res) contains the results in the form of three lists
# res[0] list containg reflectivity (Q, R(Q), R(Q)Q^4)
# res[1][0] containing the sld profile (z, sld)
# res[2][0] containing the solvent profile (z, solv)