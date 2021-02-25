from anaklasis import anaklasis

project='XRR_supported_DPPC_membrane'

# We have a single uniform layer with full coverage
patches=[1.0]

# Create empty system list where we then append
# the model.
system=[]

system.append([
	#  Re_sld  Im_sld   thk rough solv description
	[ 9.41e-6, 0.00e-6,  0 , 3.0, 1.0, 'H2O'],
	[ 14.4e-6, 0.00e-6, 10,  3.0, 0.4, 'outer_heads'],
	[ 8.40e-6, 0.00e-6, 35 , 3.0, 0.0, 'tails'],
	[ 14.4e-6, 0.00e-6, 10,  3.0, 0.3, 'inner_heads'],
	[ 9.41e-6, 0.00e-6,  5 , 3.0, 1.0, 'thin_H2O_layer'],
	[ 19.5e-6, 0.00e-6, 10 , 1.0, 0.0, 'SiO2'],
	[ 18.5e-6, 0.00e-6,  0 , 0.0, 0.0, 'Si'],
	])

global_param = []

resolution=[0.001]
background = [1.0e-9]
scale = [1.0]
qmax = [0.7]

res = anaklasis.calculate(project, resolution, 
	patches, system, global_param, 
	background, scale, qmax, plot=True)

# Note that the tuple (res) contains the results in the form of three lists
# res[0] list containing reflectivity (Q, R(Q), R(Q)Q^4)
# res[1][0] containing the sld profile (z, sld)
# res[2][0] containing the solvent profile (z, solv)
