from anaklasis import ref

project='2layers'

# Create single model(patch)
model=[
	#  Re_sld  Im_sld   thk rough solv description
	[ 0.00e-5, 0.00e-7,  0 , 3.0, 0.0, 'air'],
	[ 12.4e-5, 1.28e-5, 50,  3.0, 0.0, 'Au'],
	[ 5.94e-5, 7.69e-6, 40 , 3.0, 0.0, 'Fe'],
	[ 2.00e-5, 4.58e-7,  0 , 0.0, 0.0, 'Si'],
	]

# We have a single uniform layer with full coverage
patches=[1.0]
system=[model]

global_param = []

resolution=[0.001]
background = [1.0e-9]
scale = [1.0]
qmax = [0.7]

res = ref.calculate(project, resolution, 
	patches, system, global_param, 
	background, scale, qmax, plot=True)

# Note that the dictionary (res) contains the results
# res[("reflectivity")] Numpy array containing reflectivity ([Q, R(Q), R(Q)Q^4])
# res[("profile")] Numpy array containing the sld profile ([z, sld])
# res[("solvent")] Numpy array containing the solvent profile ([z, solv])
