'''
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
import os
import numpy as np
import corner
import pprint
import sympy
import matplotlib.pyplot as plt
import numdifftools as nd
import emcee
import sys
from tqdm import tqdm
import os.path
from os import path
import shutil
import time
from multiprocessing import Pool
import multiprocessing
from scipy.optimize import differential_evolution
from scipy.optimize import minimize
import scipy
import warnings
from scipy import special
from fortran_ref import f_realref
# if os.name == 'nt':
# 	spath=np.__file__
# 	extra_dll_dir = os.path.join(spath[:-17], 'fortran_ref\.libs')
# 	if os.path.isdir(extra_dll_dir): 
# 		os.add_dll_directory(extra_dll_dir)


def calc_profile(z, LayerMatrix, Nlayers):
	#z, distance from the surface along the z axis
	#LayerMatrix, matrix containing the real (1st column) and imaginary (2nd column) scattering
	# length density of each layer. The 3rd column corresponds to the thickness (A) of
	# each layer and the 4th column to the roughness (A). 5th column layer hydration
	# row 1 and Nlayers+1 concern the fronting and backing material respectively
	# Nlayers, Number of layers in the model
	#from scipy import special
	#import numpy as np

	rho=LayerMatrix[0][0]
	for i in range(2,Nlayers+3):
		zi=0
		for j in range(1,i):
			zi=zi+LayerMatrix[j-1][2]

		if LayerMatrix[Nlayers+1][4] == 1.0 and LayerMatrix[0][4] == 0.0:
			lmi1=LayerMatrix[i-1][0]*(1.0-LayerMatrix[i-1][4])+LayerMatrix[Nlayers+1][0]*LayerMatrix[i-1][4]
			lmi2=LayerMatrix[i-2][0]*(1.0-LayerMatrix[i-2][4])+LayerMatrix[Nlayers+1][0]*LayerMatrix[i-2][4]
		if LayerMatrix[Nlayers+1][4] == 0.0 and LayerMatrix[0][4] == 1.0:
			lmi1=LayerMatrix[i-1][0]*(1.0-LayerMatrix[i-1][4])+LayerMatrix[0][0]*LayerMatrix[i-1][4]
			lmi2=LayerMatrix[i-2][0]*(1.0-LayerMatrix[i-2][4])+LayerMatrix[0][0]*LayerMatrix[i-2][4]
		if LayerMatrix[Nlayers+1][4] == 0.0 and LayerMatrix[0][4] == 0.0:
			lmi1=LayerMatrix[i-1][0]
			lmi2=LayerMatrix[i-2][0]

		if LayerMatrix[i-2][3] != 0:
			#rho=rho+((LayerMatrix[i-1][0]-LayerMatrix[i-2][0])/2.0)*(1.0+special.erf((z-zi)/(np.sqrt(2.0)*LayerMatrix[i-2][3])))
			rho=rho+((lmi1-lmi2)/2.0)*(1.0+special.erf((z-zi)/(np.sqrt(2.0)*LayerMatrix[i-2][3])))

		else:
			if z >= zi:
				#rho=rho+((LayerMatrix[i-1][0]-LayerMatrix[i-2][0])/2.0)*(2.0)
				rho=rho+((lmi1-lmi2)/2.0)*(2.0)
			else:
				rho=rho+0.0	

	return rho

def calc_solvent_penetration(z, LayerMatrix, Nlayers):
	#z, distance from the surface along the z axis
	#LayerMatrix, matrix containing the real (1st column) and imaginary (2nd column) scattering
	# length density of each layer. The 3rd column corresponds to the thickness (A) of
	# each layer and the 4th column to the roughness (A). 5th column layer hydration
	# row 1 and Nlayers+1 concern the fronting and backing material respectively
	# Nlayers, Number of layers in the model
	#from scipy import special
	#import numpy as np

	if LayerMatrix[Nlayers+1][4] == 1.0 and LayerMatrix[0][4] == 0.0: rho=0.0
	if LayerMatrix[Nlayers+1][4] == 0.0 and LayerMatrix[0][4] == 1.0: rho=1.0
	if LayerMatrix[Nlayers+1][4] == 0.0 and LayerMatrix[0][4] == 0.0: return 0.0

	for i in range(2,Nlayers+3):
		zi=0
		for j in range(1,i):
			zi=zi+LayerMatrix[j-1][2]

		if LayerMatrix[Nlayers+1][4] == 1.0 and LayerMatrix[0][4] == 0.0:
			lmi1=LayerMatrix[i-1][4]
			lmi2=LayerMatrix[i-2][4]
		if LayerMatrix[Nlayers+1][4] == 0.0 and LayerMatrix[0][4] == 1.0:
			lmi1=LayerMatrix[i-1][4]
			lmi2=LayerMatrix[i-2][4]

		if LayerMatrix[i-2][3] != 0:
			#rho=rho+((LayerMatrix[i-1][0]-LayerMatrix[i-2][0])/2.0)*(1.0+special.erf((z-zi)/(np.sqrt(2.0)*LayerMatrix[i-2][3])))
			rho=rho+((lmi1-lmi2)/2.0)*(1.0+special.erf((z-zi)/(np.sqrt(2.0)*LayerMatrix[i-2][3])))

		else:
			if z > zi:
				#rho=rho+((LayerMatrix[i-1][0]-LayerMatrix[i-2][0])/2.0)*(2.0)
				rho=rho+((lmi1-lmi2)/2.0)*(2.0)
			else:
				rho=rho+0.0	
	return rho

def profile(LayerMatrix, npoints):
	# Given a layer model (LayerMatrix), calculate sld profile with npoints between -10A -> Dmax+10A
	# and return an array with z, sld(z) format
	#import numpy as np
	global_total_d=0.0
	#for i in range(1,np.size(LayerMatrix,0)-1):
	for j in range(len(LayerMatrix)):
		total_d=0.0
		for i in range(1,np.size(LayerMatrix[j],0)):
			total_d=total_d+LayerMatrix[j][i][2]
		if total_d > global_total_d: global_total_d=total_d

	#z_bin = np.linspace(-10, total_d+10, npoints+1)
	#Profile = np.empty([np.size(LayerMatrix[j],0),npoints, 2])
	Profile = np.empty([np.size(LayerMatrix,0),npoints, 2])
	for j in range(len(LayerMatrix)):
		z_bin = np.linspace(-4*LayerMatrix[j][0][3]-5.0, global_total_d+4*LayerMatrix[j][np.size(LayerMatrix[j],0)-2][3]+5.0, npoints+1)
		for i in range(len(z_bin)-1):
			Profile[j][i][0]=z_bin[i]
			Profile[j][i][1]=calc_profile(z_bin[i],LayerMatrix[j], np.size(LayerMatrix[j],0)-2)

	return Profile

def Reflectivity(q_bin,res_bin, LayerMatrix, resolution,bkg,scale,patches):
	# Given a Q bin 1-d array, layer model (LayerMatrix) and instrumental resolution  calculate the specular reflectivity
	# and return an array with Q, R(Q), R(Q)Q^4 
	#import numpy as np
	Refl=np.empty([len(q_bin), 3])
	for i in range(len(q_bin)):
		Refl[i][0]=q_bin[i]
		#This is for f2py
		Refl[i][1]=0.0
		Refl[i][2]=0.0
		for j in range(len(LayerMatrix)):
			flayers = np.array(LayerMatrix[j])
			if resolution == -1:
				Refl[i][1]=Refl[i][1]+patches[j]*scale*f_realref(0,q_bin[i],res_bin[i]/q_bin[i],flayers[:,0:5],np.size(LayerMatrix[j],0)-2)+bkg
			else:
				Refl[i][1]=Refl[i][1]+patches[j]*scale*f_realref(0,q_bin[i],res_bin[i],flayers[:,0:5],np.size(LayerMatrix[j],0)-2)+bkg
			#The following line for pure Python calculation
			#Refl[i][1]=real_refl(q_bin[i], LayerMatrix, np.size(LayerMatrix,0)-2, res_bin[i])
		Refl[i][2]=Refl[i][1]*q_bin[i]**4

	return Refl

def solvent_penetration(LayerMatrix, npoints):
	# Given a layer model (LayerMatrix), calculate sld profile with npoints between -10A -> Dmax+10A
	# and return an array with z, sovent penetration from 0 to 1
	#import numpy as np
	global_total_d=0.0
	#for i in range(1,np.size(LayerMatrix,0)-1):
	for j in range(len(LayerMatrix)):
		total_d=0.0
		for i in range(1,np.size(LayerMatrix[j],0)):
			total_d=total_d+LayerMatrix[j][i][2]
		if total_d > global_total_d: global_total_d=total_d

	#z_bin = np.linspace(-10, total_d+10, npoints+1)
	#Profile = np.empty([np.size(LayerMatrix[j],0),npoints, 2])
	Profile = np.empty([np.size(LayerMatrix,0),npoints, 2])
	for j in range(len(LayerMatrix)):
		z_bin = np.linspace(-4*LayerMatrix[j][0][3]-5.0, global_total_d+4*LayerMatrix[j][np.size(LayerMatrix[j],0)-2][3]+5.0, npoints+1)
		for i in range(len(z_bin)-1):
			Profile[j][i][0]=z_bin[i]
			Profile[j][i][1]=calc_solvent_penetration(z_bin[i],LayerMatrix[j], np.size(LayerMatrix[j],0)-2)

	return Profile

def chi_square(data, LayerMatrix, resolution, bkg, scale, patches):
	# Given a set of data, a layer model (LayerMatrix) and instrumental resolution
	# calculate the discrapancy (chi square) between the model and experimental data.
	chi=0
	Nexp=np.size(data,0)
	for i in range(0, Nexp):
		if data[i][2] == 0:
			if resolution != -1:
				for j in range(len(LayerMatrix)):
					flayers = np.array(LayerMatrix[j])
					chi=chi+patches[j]*(1.0/float(Nexp))*((data[i][1]-scale*f_realref(1,data[i][0],resolution,flayers[:,0:5],np.size(LayerMatrix[j],0)-2)-bkg))**2
			else:
				for j in range(len(LayerMatrix)):
					flayers = np.array(LayerMatrix[j])
					chi=chi+patches[j]*(1.0/float(Nexp))*((data[i][1]-scale*f_realref(1,data[i][0],data[i][3]/data[i][0],flayers[:,0:5],np.size(LayerMatrix[j],0)-2)-bkg))**2
		else:
			if resolution != -1:
				for j in range(len(LayerMatrix)):
					flayers = np.array(LayerMatrix[j])
					chi=chi+patches[j]*(1.0/float(Nexp))*((data[i][1]-scale*f_realref(1,data[i][0],resolution,flayers[:,0:5],np.size(LayerMatrix[j],0)-2)-bkg)/(data[i][2]))**2
			else:
				for j in range(len(LayerMatrix)):
					flayers = np.array(LayerMatrix[j])
					chi=chi+patches[j]*(1.0/float(Nexp))*((data[i][1]-scale*f_realref(1,data[i][0],data[i][3]/data[i][0],flayers[:,0:5],np.size(LayerMatrix[j],0)-2)-bkg)/(data[i][2]))**2
	return chi


def fig_of_merit_sym(x, *argv):
	#
	data=argv[0]
	resolution=argv[1]
	Nlayers=argv[2]
	m_param=argv[3]
	layers_max=argv[4]
	model_param=argv[5]
	model_constraints=argv[6]
	fit_mode=argv[7]
	num_curves=argv[8]
	c_param=argv[9]
	multi_param=argv[10]
	f_layer_fun=argv[11]
	f_left_fun=argv[12]
	f_right_fun=argv[13]
	center_fun=argv[14]
	experror=argv[15]
	fit_weight=argv[16]
	patches=argv[17]

	#print(x)
	fom=0.0
	for curve in range(0,num_curves):

		vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,m_param):
			vp[i]=x[sum(Nlayers)*5+i]

		vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,c_param):
			vm[i]=x[sum(Nlayers)*5+m_param+2*num_curves+num_curves*i+curve]

		#Check constraints
		for i in range(0,len(model_constraints)):
			if '>' in str(model_constraints[i]) or '<' in str(model_constraints[i]):
				if center_fun[i] == '>':
					left=np.float((f_left_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					right=np.float((f_right_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))

					if left < right:
						fom=np.inf
						#print(f'FOM = {fom}\r', end="       ")
						return np.float('inf')

				if center_fun[i] == '<':
					left=np.float((f_left_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					right=np.float((f_right_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					if left > right:
						fom=np.inf
						#print(f'FOM = {fom}\r', end="       ")
						return np.float('inf')

		layers=[]
		for k in range(len(layers_max)):
			sub_layers=[]
			for i in range(0,Nlayers[k]):
				line=[]
				for j in range(0,5):
					if isinstance(layers_max[k][i][j], str):
						line.append(np.float((f_layer_fun(k,i,j,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],i,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))
					else:
						if k > 0:
							line.append(x[5*sum(Nlayers[0:k])+i*5+j])
						else:
							line.append(x[i*5+j])

				sub_layers.append(line)
			layers.append(sub_layers)

		bkg=x[sum(Nlayers)*5+m_param+curve]
		scl=x[sum(Nlayers)*5+m_param+num_curves+curve]

		h1=0.0
		h2=0.0
		Nexp=np.size(data[curve],0)
		for i in range(0, Nexp):

			if fit_mode == 0:
				# pure chi square
				# The folloiwng line for pure Python calculation
				#h3=real_refl_fast(data[i][0], layers, np.size(layers,0)-2, resolution[curve])
				#The folloiwing two lines for f2py
				
				if resolution[curve] != -1:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i,0],resolution[curve],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				else:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i,0],data[curve][i,3]/data[curve][i,0],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				if experror == False:
					h1=h1+(1.0/float(Nexp))*((data[curve][i,1]-h3)/(h3))**2
				else:
					h1=h1+(1.0/float(Nexp))*((data[curve][i,1]-h3)/data[curve][i,2])**2
				#h2=h2+1.0

			if fit_mode == 1:
				# pure chi square
				# The folloiwng line for pure Python calculation
				#h3=real_refl_fast(data[i][0], layers, np.size(layers,0)-2, resolution[curve])
				#The folloiwing two lines for f2py
				
				if resolution[curve] != -1:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i,0],resolution[curve],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				else:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i,0],data[curve][i,3]/data[curve][i,0],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				if experror == True:
					h1=h1+(1.0/float(Nexp))*((np.log10(data[curve][i,1])-np.log10(h3))/np.absolute((data[curve][i,2])/(np.log(10.0)*data[curve][i,1])))**2
				else:
					h1=h1+(1.0/float(Nexp))*((np.log10(data[curve][i,1])-np.log10(h3)))**2


		fom=fom+fit_weight[curve]*h1


	#print(f'FOM = {fom}\r', end="       ")

	return fom

def reduced_chi_sym(x, *argv):
	#
	data=argv[0]
	resolution=argv[1]
	Nlayers=argv[2]
	m_param=argv[3]
	layers_max=argv[4]
	model_param=argv[5]
	fit_mode=argv[6]
	num_curves=argv[7]
	c_param=argv[8]
	multi_param=argv[9]
	f_layer_fun=argv[10]
	fit_weight=argv[11]
	patches=argv[12]
	free_param=argv[13]


	fom=0.0

	for curve in range(0,num_curves):

		vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,m_param):
			vp[i]=x[sum(Nlayers)*5+i]

		vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,c_param):
			vm[i]=x[sum(Nlayers)*5+m_param+2*num_curves+num_curves*i+curve]

		layers=[]
		for k in range(len(layers_max)):
			sub_layers=[]
			for i in range(0,Nlayers[k]):
				line=[]
				for j in range(0,5):
					if isinstance(layers_max[k][i][j], str):
						line.append(np.float((f_layer_fun(k,i,j,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],i,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))
					else:
						if k > 0:
							line.append(x[5*sum(Nlayers[0:k])+i*5+j])
						else:
							line.append(x[i*5+j])

				sub_layers.append(line)
			layers.append(sub_layers)

		bkg=x[sum(Nlayers)*5+m_param+curve]
		scl=x[sum(Nlayers)*5+m_param+num_curves+curve]

		h1=0.0
		h2=0.0
		Nexp=np.size(data[curve],0)-free_param
		sum_weight=sum(fit_weight)
		for i in range(0, Nexp):

			if fit_mode == 0:
				# pure chi square
				# The folloiwng line for pure Python calculation
				#h3=real_refl_fast(data[i][0], layers, np.size(layers,0)-2, resolution[curve])
				#The folloiwing two lines for f2py
				
				if resolution[curve] != -1:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i,0],resolution[curve],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				else:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i,0],data[curve][i,3]/data[curve][i,0],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				#h1=h1+(1.0/float(Nexp))*((data[i][1]-h3)/(h3))**2
				h1=h1+(1.0/float(Nexp))*((data[curve][i,1]-h3)/data[curve][i,2])**2
				h2=h2+1.0
			if fit_mode == 1:
				# pure chi square
				# The folloiwng line for pure Python calculation
				#h3=real_refl_fast(data[i][0], layers, np.size(layers,0)-2, resolution[curve])
				#The folloiwing two lines for f2py

				if resolution[curve] != -1:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i][0],resolution[curve],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				else:
					h3=0.0
					for k in range(len(layers_max)):
						flayers = np.array(layers[k])
						h3=h3+patches[k]*scl*f_realref(1,data[curve][i][0],data[curve][i][3]/data[curve][i][0],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
				h1=h1+(1.0/float(Nexp))*((np.log10(data[curve][i,1])-np.log10(h3))/np.absolute((data[curve][i,2])/(np.log(10.0)*data[curve][i,1])))**2
				h2=h2+1.0

		fom=fom+h1*fit_weight[curve]/sum_weight

	return fom


def log_likelihood(theta, *argv):
	data=argv[0]
	resolution=argv[1]
	Nlayers=argv[2]
	m_param=argv[3]
	layers_max=argv[4]
	layers_min=argv[5]
	model_param=argv[6]
	background=argv[7]
	fit_mode=argv[8]
	num_curves=argv[9]
	c_param=argv[10]
	multi_param=argv[11]
	scale=argv[12]
	f_layer_fun=argv[13]
	fit_weight=argv[14]
	patches=argv[15]


	logl=0.0
	mk=0
	for curve in range(0,num_curves):

		k=0
		for l in range(len(layers_min)):
			for i in range(0,Nlayers[l]):
				for j in range(0,5):
					if layers_min[l][i][j] != layers_max[l][i][j]:
						k=k+1


		vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		l=0
		for i in range(0,m_param):
			if model_param[i][4] == 'uniform':
				if model_param[i][1] != model_param[i][2]:
					vp[i]=theta[k+l]
					l=l+1
				else:
					vp[i]=model_param[i][1]
			if model_param[i][4] == 'normal':
				if model_param[i][2] != 0.0:
					vp[i]=theta[k+l]
					l=l+1
				else:
					vp[i]=model_param[i][1]

		k=k+l

		if background[curve][2] == 'uniform':
			if background[curve][0] != background[curve][1]:
				bkg=theta[k+curve]
			else:
				bkg=background[curve][0]
		if background[curve][2] == 'normal':
			if background[curve][1] != 0.0:
				bkg=theta[k+curve]
			else:
				bkg=background[curve][0]

		l=0
		for i in range(0,num_curves):
			if background[i][2] == 'uniform':
				if background[i][0] != background[i][1]:
					l=l+1
			if background[i][2] == 'normal':
				if background[i][1] != 0.0:
					l=l+1
		k=k+l

		if scale[curve][2] == 'uniform':
			if scale[curve][0] != scale[curve][1]:
				scl=theta[k+curve]
			else:
				scl=scale[curve][0]
		if scale[curve][2] == 'normal':
			if scale[curve][1] != 0.0:
				scl=theta[k+curve]
			else:
				scl=scale[curve][0]

		l=0
		for i in range(0,num_curves):
			if scale[i][2] == 'uniform':
				if scale[i][0] != scale[i][1]:
					l=l+1
			if scale[i][2] == 'normal':
				if scale[i][1] != 0.0:
					l=l+1
		k=k+l

		vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		l=0

		for i in range(0,c_param):
			if multi_param[i][-1] == 'uniform':
				if multi_param[i][2*curve+1] != multi_param[i][2*curve+2]:
					vm[i]=theta[k+mk+l]
					l=l+1
				else:
					vm[i]=multi_param[i][2*curve+1]
			if multi_param[i][-1] == 'normal':
				if multi_param[i][2*curve+2] != 0.0:
					vm[i]=theta[k+mk+l]
					l=l+1
				else:
					vm[i]=multi_param[i][2*curve+1]
		mk=mk+l

		layers=[]
		k=0
		for l in range(len(layers_min)):
			sub_layers=[]
			for i in range(0,Nlayers[l]):
				line=[]
				for j in range(0,5):
					if isinstance(layers_max[l][i][j], str):
						line.append(np.float((f_layer_fun(l,i,j,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],i,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))
					else:
						if layers_min[l][i][j] != layers_max[l][i][j]:
							line.append(theta[k])
							k=k+1
						else:
							line.append(layers_max[l][i][j])

				sub_layers.append(line)
			layers.append(sub_layers)


		x=data[curve][:,0]
		y=data[curve][:,1]
		yerr=data[curve][:,2]
		Nexp=np.size(data[curve],0)
		model=np.zeros(Nexp)
		
		for i in range(0, Nexp):
			if resolution[curve] != -1:
				model[i]=0.0
				for k in range(len(layers_max)):
					flayers = np.array(layers[k])
					model[i] = model[i]+patches[k]*scl*f_realref(1,data[curve][i][0],resolution[curve],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
			else:
				model[i]=0.0
				for k in range(len(layers_max)):
					flayers = np.array(layers[k])
					model[i] = model[i]+patches[k]*scl*f_realref(1,data[curve][i][0],data[curve][i][3]/data[curve][i][0],flayers[:,0:5],np.size(layers[k],0)-2)+bkg
		if fit_mode == 0:
			sigma2 = yerr ** 2
			logl = logl + fit_weight[curve]*(-0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(2*np.pi*sigma2)))
		if fit_mode == 1:
			sigma2 = np.absolute(yerr/(y*np.log(10))) ** 2
			logl = logl + fit_weight[curve]*(-0.5 * np.sum((np.log10(y) - np.log10(model)) ** 2 / sigma2 + np.log(2*np.pi*sigma2)))

	return logl

def log_prior(theta, *argv):
	Nlayers=argv[0]
	m_param=argv[1]
	layers_max=argv[2]
	layers_min=argv[3]
	model_param=argv[4]	
	model_constraints=argv[5]
	background=argv[6]
	num_curves=argv[7]
	c_param=argv[8]
	multi_param=argv[9]
	scale=argv[10]
	f_left_fun=argv[11]
	f_right_fun=argv[12]
	center_fun=argv[13]
 	# check bound
	cpass=0
	l=0
	for k in range(len(layers_min)):
		for i in range(0,Nlayers[k]):
			for j in range(0,5): 
				if layers_min[k][i][j] != layers_max[k][i][j]:
					if np.float(layers_min[k][i][j]) < theta[l] < np.float(layers_max[k][i][j]):
			  			cpass=cpass+1
					l=l+1

	for i in range(0,m_param):
		if model_param[i][4] == 'uniform':
			if model_param[i][1] != model_param[i][2]:
				if np.float(model_param[i][1]) < theta[l] < np.float(model_param[i][2]):
			  		cpass=cpass+1
				l=l+1
		if model_param[i][4] == 'normal':
			if model_param[i][2] != 0.0:
				if np.float(model_param[i][1]-3*model_param[i][2]) < theta[l] < np.float(model_param[i][1]+3*model_param[i][2]):
			  		cpass=cpass+1
				l=l+1

	for i in range(0,num_curves):
		if background[i][2] == 'uniform':
			if background[i][0] != background[i][1]:
				if np.float(background[i][0]) < theta[l] < np.float(background[i][1]):
					cpass=cpass+1
				l=l+1
		if background[i][2] == 'normal':
			if background[i][1] != 0.0:
				if np.float(background[i][0]-3*background[i][1]) < theta[l] < np.float(background[i][0]+3*background[i][1]):
					cpass=cpass+1
				l=l+1

	for i in range(0,num_curves):
		if scale[i][2] == 'uniform':
			if scale[i][0] != scale[i][1]:
				if np.float(scale[i][0]) < theta[l] < np.float(scale[i][1]):
					cpass=cpass+1
				l=l+1
		if scale[i][2] == 'normal':
			if scale[i][1] != 0.0:
				if np.float(scale[i][0]-3*scale[i][1]) < theta[l] < np.float(scale[i][0]+3*scale[i][1]):
					cpass=cpass+1
				l=l+1

	for i in range(0,c_param):
		for j in range(0,num_curves):
			if multi_param[i][-1] == 'uniform':
				if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
					if np.float(multi_param[i][2*j+1]) < theta[l] < np.float(multi_param[i][2*j+2]):
						cpass = cpass + 1
					l=l+1
			if multi_param[i][-1] == 'normal':
				if multi_param[i][2*j+2] != 0.0:
					if np.float(multi_param[i][2*j+1]-3*multi_param[i][2*j+2]) < theta[l] < np.float(multi_param[i][2*j+1]+3*multi_param[i][2*j+2]):
						cpass = cpass + 1
					l=l+1


	if cpass != l:
		return -np.inf


	mk=0
	for curve in range(0,num_curves):
		k=0
		for l in range(len(layers_min)):
			for i in range(0,Nlayers[l]):
				for j in range(0,5):
					if layers_min[l][i][j] != layers_max[l][i][j]:
						k=k+1

		vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		l=0
		for i in range(0,m_param):
			if model_param[i][4] == 'uniform':
				if model_param[i][1] != model_param[i][2]:
					vp[i]=theta[k+l]
					l=l+1
				else:
					vp[i]=model_param[i][1]
			if model_param[i][4] == 'normal':
				if model_param[i][2] != 0.0:
					vp[i]=theta[k+l]
					l=l+1
				else:
					vp[i]=model_param[i][1]

		k=k+l
		for i in range(0,num_curves):
			if background[i][2] == 'uniform':
				if background[i][0] != background[i][1]:
					k=k+1
			if background[i][2] == 'normal':
				if background[i][1] != 0.0:
					k=k+1

		for i in range(0,num_curves):
			if scale[i][2] == 'uniform':
				if scale[i][0] != scale[i][1]:
					k=k+1
			if scale[i][2] == 'normal':
				if scale[i][1] != 0.0:
					k=k+1

		vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		l=0
		for i in range(0,c_param):
			if multi_param[i][-1] == 'uniform':
				if multi_param[i][2*curve+1] != multi_param[i][2*curve+2]:
					vm[i]=theta[k+mk+l]
					l=l+1
				else:
					vm[i]=multi_param[i][2*curve+1]
			if multi_param[i][-1] == 'normal':
				if multi_param[i][2*curve+2] != 0.0:
					vm[i]=theta[k+mk+l]
					l=l+1
				else:
					vm[i]=multi_param[i][2*curve+1]
		mk=mk+l
		#Check constraints
		for i in range(0,len(model_constraints)):
			if '>' in str(model_constraints[i]) or '<' in str(model_constraints[i]):
				if center_fun[i] == '>':
					left=np.float((f_left_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					right=np.float((f_right_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					if left < right:
						return -np.inf


				if center_fun[i] == '<':
					left=np.float((f_left_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					right=np.float((f_right_fun(i,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39])))
					if left > right:
						return -np.inf

	
	sum_prior=0.0

	l=0

	for i in range(0,m_param):
		if model_param[i][4] == 'uniform':
			if model_param[i][1] != model_param[i][2]:
				l=l+1
		if model_param[i][4] == 'normal':
			if model_param[i][2] != 0.0:
				sum_prior=sum_prior+np.log(1.0/(np.sqrt(2*np.pi)*model_param[i][2]))-0.5*(theta[l]-model_param[i][1])**2/model_param[i][2]**2
				l=l+1

	for i in range(0,num_curves):
		if background[i][2] == 'uniform':
			if background[i][0] != background[i][1]:
				l=l+1
		if background[i][2] == 'normal':
			if background[i][1] != 0.0:
				sum_prior=sum_prior+np.log(1.0/(np.sqrt(2*np.pi)*background[i][1]))-0.5*(theta[l]-background[i][0])**2/background[i][1]**2
				l=l+1

	for i in range(0,num_curves):
		if scale[i][2] == 'uniform':
			if scale[i][0] != scale[i][1]:
				l=l+1
		if scale[i][2] == 'normal':
			if scale[i][1] != 0.0:
				sum_prior=sum_prior+np.log(1.0/(np.sqrt(2*np.pi)*scale[i][1]))-0.5*(theta[l]-scale[i][0])**2/scale[i][1]**2
				l=l+1

	for i in range(0,c_param):
		for j in range(0,num_curves):
			if multi_param[i][-1] == 'uniform':
				if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
					l=l+1
			if multi_param[i][-1] == 'normal':
				if multi_param[i][2*j+2] != 0.0:
					sum_prior=sum_prior+np.log(1.0/(np.sqrt(2*np.pi)*multi_param[i][2*j+2]))-0.5*(theta[l]-multi_param[i][2*j+1])**2/multi_param[i][2*j+2]**2
					l=l+1

	return sum_prior


def log_probability(theta, *argv):
	data=argv[0]
	resolution=argv[1]
	Nlayers=argv[2]
	m_param=argv[3]
	layers_max=argv[4]
	layers_min=argv[5]
	model_param=argv[6]
	model_constraints=argv[7]
	background=argv[8]
	fit_mode=argv[9]
	num_curves=argv[10]
	c_param=argv[11]
	multi_param=argv[12]
	scale=argv[13]
	f_layer_fun=argv[14]
	f_left_fun=argv[15]
	f_right_fun=argv[16]
	center_fun=argv[17]
	fit_weight=argv[18]
	patches=argv[19]

	lp = log_prior(theta,Nlayers,m_param,layers_max,layers_min,model_param,model_constraints,background,num_curves,c_param,multi_param,scale,f_left_fun,f_right_fun,center_fun)
	if not np.isfinite(lp):
		return -np.inf
	#print(lp)
	return lp + log_likelihood(theta, data, resolution, Nlayers, m_param, layers_max, layers_min, model_param, background, fit_mode,num_curves,c_param,multi_param,scale,f_layer_fun,fit_weight,patches)



def f_layer_fun(k,i,j,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39):
	return np.float((layer_fun[k][i][j](p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39)))

def f_left_fun(i,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39):
	return np.float((left_fun[i](p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39)))

def f_right_fun(i,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39):
	return np.float((right_fun[i](p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39)))

def fit(project, in_file, units, fit_mode,fit_weight, method, resolution, patches, system, model_param, multi_param, model_constraints,background,scale,experror,plot,fast):
	"""
	--------------------------------------------------------------------
	Program ANAKLASIS - Fit Module for X-ray/Neutron reflection datasets
	version 1.2, February 2021
	developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ
	for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de
	--------------------------------------------------------------------
	This function performs fits of x-ray and neutron reflection data

	:param project: string containing the name of the project. All output files
	 are saved in a directory with the same name. If project name is 'none'
	 no output files are written on disk. 

	:param in_file: list of M string elements. 

	 Syntax:

	 in_file = ['file1','file2',...,'fileM']

	 Each element is the path to an ASCII file containing reflectivity data. 
	 The file structure should be in column format ( Q, Reflectivity, error in 
	 reflectivity (dR), error in Q (dQ)). If a third column (dR) is not present
	 in the data, set experror=False (see below). If a 4th column (dQ) is 
	 present you may take pointwise resolution into account by setting resolution
	 equal to -1 (see below). Additional columns will be ignored. 
	 Lines beginning with '#' can be present and are considered as comments. 
	 Reflectivity should be footprint corrected.

	:param units: list of M string elements that can be either 'A' inverse
	 Angstrom or 'nm' inverse nanometers, describing the units of momentum
	 transfer (Q) in each of the corresponding input files.  

	 Syntax:

	 units = ['units1','units2',...,'unitsM']

	:fit_mode: integer that can take values 0 (linear figure of merit) or 
	 1 (log10 figure of merit)

	:method: string that can be either 'simple', either 'mcmc' or 'bootstrap'
	With 'simple' method, parameter uncertainty calculation is performed by
	calculation of the Hessian Matrix. With 'mcmc' method a Markov Chain Monte
	Carlo run is performed and together with parameter uncertainty estimation
	a corner plot with parameter correlations is exported. With 'bootstrap' 
	also parameter correlation is exported however note that this methods is
	far more computationally intensive than 'mcmc'.

	:resolution: list of M float elements. Each element is the dQ/Q resolution of
	the corresponding input data in param in_file

	Syntax:

		resolution = [res1, res2, ... , resM]

	 for example res1 describes the resolution of the data present in 'file1'.

	Note that if a 4th dQ column is present in the in_file, you may set res = -1
	so that pointwise dQ/Q resolution is taken into account.

	:patches: list of surface coverage of each defined model. In case of a single
	 defined model (most usual case) the definition has the following syntax

	Syntax:
		patches = [1.0]

	in case of K defined models

	Syntax:
	 	patches = [coverage_1, coverage_2 ... coverage_K]

	where the sum of all coverages should add up to 1.

	:system: List of lists containing defined models.
	 Each model is represented as a list of N+1 lists(lines) that contain 6 elements.


	Syntax:
		system = [[
			[  Re_sld0, Im_sld0, d0, sigma0, solv0, 'layer0'],
			[  Re_sld1, Im_sld1, d1, sigma1, solv1, 'layer1'],
			[  Re_sld2, Im_sld2, d2, sigma2, solv2, 'layer2'],
			.
			.
			.
			[  Re_sldN, Im_sldN, dN, sigmaN, solvN, 'layerN'],
			]]

	 If more than one models needs to be defined, more lists can be appended as:

	 	system.append([
			[  Re_sld0, Im_sld0, d0, sigma0, solv0, 'layer0'],
			[  Re_sld1, Im_sld1, d1, sigma1, solv1, 'layer1'],
			[  Re_sld2, Im_sld2, d2, sigma2, solv2, 'layer2'],
			.
			.
			.
			[  Re_sldN, Im_sldN, dN, sigmaN, solvN, 'layerN'],
	 		])

	 Each line (6 element list) represents a layer, from layer 0 (semi-infinite
	 fronting) to layer N (semi-infinite backing). The elements of the list 
	 correspond to Real sld (in 10^-6 A^-2) , Imaginary sld in 10^-6 A^-2,
	 thickness (in Angstrom), roughness (in Angstrom), solvent volumer fraction
	 (0 to 1) and layer description (string) respectively. All elements (except
	 description) can be numerical values or SymPy expressions (string) involving
	 global and multi parameters. Additionally in the SymPy expressions the 
	 integer (n) can be used, that represents the number of the layer from 0 
	 to N, and/or the summation integers ii,jj,kk, and/or the variables x,y,z
	 than may be used in SymPy integrals or derivatives.

	 When solv0 = 0 and solvN = 1 (fronting and backing solvent volume fraction)
	 then the solvent volume fraction parameter assumes that the backing layer
	 represents a semi-infinite liquid medium and that the liquid may
	 penetrate layers 1 to N-1 (usual in measurements at the solid/liquid or
	 air/liquid interface). 

	 When solv0 = 1 and solvN = 0 (fronting and backing solvent volume fraction)
	 then the solvent volume fraction parameter assumes that the fronting layer
	 represents a semi-infinite liquid medium and that the liquid may
	 penetrate layers 1 to N-1 (usual in measurements at the solid/liquid or
	 air/liquid interface). 

	 When solv0 = 0 and solvN = 0 (fronting and backing solvent volume fraction)
	 all solv values should be set zero. Any non zero value is ignored. 

	 Note that sigma_i represents the roughness between layer_i and layer_(i+1) 

	:model_param: Global parameter list of X 5-element lists.

	Syntax:
		param = [
			['p0', p0_A, p0_B, 'p0_description','p0_type'],
			['p1', p1_A, p1_B, 'p1_description','p1_type'],
			.
			.
			.
			['pX', pX_A, pX_B, 'pX_name''pX_type'],
			]

	 Up to X=20 global parameters can be defined. The names should be strings of
	 the form 'p0','p1' ... 'p19' respectively. The two last elements of each 
	 global parameter are strings defining the description and type of the 
	 parameter. The type can be either 'uniform' or 'normal'. In case of uniform
	 type, The two values p_A/P_B represent the min/max bound of a 
	 uniform distribution. In case of normal type the values p_A/p_B represent
	 the mean value and standard deviation (sd) of a normal distribution.
	 If min is equal to max value or standard deviation is zero, then the global 
	 parameter is considered as fixed, otherwise is left free to vary within the 
	 min/max bounds or within mean - 3*sd/mean + 3*sd respectively.

	:multi_param: Multi parameter list of Y 2+2*M-element lists, where M is the
	 number of input curves.

	Syntax:
		param = [
			['m0', m0_A1, m0_B1, ..., m0_AM, m0_BM, 'm0_description','m0_type'],
			['m1', m1_A1, m1_B1, ..., m1_AM, m1_BM, 'm1_description','m1_type'],
			.
			.
			.
			['mY', mY_A1, mY_B1, ..., mY_AM, mY_BM, 'mY_description','mY_type'],
			]


	 Up to Y=20 multi parameters can be defined. The names should be strings of
	 the form 'm0','m1' ... 'm19' respectively. The two last elements of each 
	 multi parameter are strings defining the description and type of the 
	 parameter. The type can be either 'uniform' or 'normal'. In case of uniform
	 type, the two values m_A/m_B represent the min/max bounds of a 
	 uniform distribution. In case of normal type the values m_A/m_B represent
	 the mean value and standard deviation (sd) of a normal distribution.
	 If min is equal to max value or standard deviation is zero, then the multi 
	 parameter is considered as fixed, otherwise is left free to vary within the 
	 min/max bounds or within mean - 3*sd/mean + 3*sd respectively.

	 The difference with global parameters is than the m_A/M_B value of a multi
	 parameter is different for each input curve, thus M (m_A/M_B) value pairs
	 need to be defined, where M is the number of input curves.

	:model_constraints: list of strings containing inequality expressions.

	Syntax:

		model_constraints = [
			'expression1',
			'expression2'.
			.
			.
			.
			'expressionZ',]

	 Any sympy valid expression (string) can be used involving defined global and 
	 multi parameters. An example expression is the following '2*p1/p0>1'.

	:background: list of M lists that contain two numbers and type (normal or uniform)
	 For uniform type these two numbers represent the lower and upper bound of the 
	 instrumental background of the  corresponding data. If the upper and lower value
	 is the same, the parameter is considered fixed during the fit, otherwise it is 
	 considered as a free to vary within the min/max bounds. In case were the
	 type is normal these two numbers represent the mean and standard deviation of the 
	 instrumental background of the  corresponding data. If the standard deviation is
	 zero, the parameter is considered fixed during the fit, otherwise it is 
	 considered as a free to vary within mean-3*sd/mean+3*sd.

	Syntax: 

		background = [
			[bkg_A1,bkg_B1,'type'],
			[bkg_A2,bkg_B2,'type'],
			.
			.
			.
			.
			[bkg_AM,bkg_BM],'type'],
			]	

	 for example [bkg_min1,bkg_max1] describe the lower and upper bound of the
	 background of the data present in 'file1'.

	 Note: Theoretical reflectivity is calculated as

	 	R = scale * R(Q) + background

	:scale: list of M lists that contain two numbers and type (normal or uniform)
	 For uniform type these two numbers represent the lower and upper bound of the 
	 scale of the corresponding data. If the upper and lower value is the same, 
	 the parameter is considered fixed during the fit, otherwise it is 
	 considered as a free to vary within the min/max bounds. In case were the
	 type is normal these two numbers represent the mean and standard deviation of the 
	 scale of the  corresponding data. If the standard deviation is
	 zero, the parameter is considered fixed during the fit, otherwise it is 
	 considered as a free to vary within mean-3*sd/mean+3*sd.

	Syntax: 

		scale = [
			[scale_A1,scale_B1,'type'],
			[scale_A2,scale_B2,'type'],
			.
			.
			.
			.
			[scale_AM,scale_BM,'type'],]

	 for example [scale_min1,scale_max1] describe the lower and upper bound of
	 the scale of the data present in 'file1'.

	 Note: Theoretical reflectivity is calculated as

	 	R = scale * R(Q) + background

	:param experror: Boolean. Set as True, if all input files contain a third 
	column with Reflectivity uncertainty. Set as False, if at least one of the
	input files contains no Reflectivity uncertainty. Also, if all data contain
	a third column but for some reason Reflectivity uncertainty is considered
	as non-reliable, you may set the parameter to False thus ignoring errors.

	:param plot: Boolean. If True, an interactive plot is displayed at the end
	of the calculations.

	:param fast: Boolean. If True, a small differential evolution population size
	is used. If False a large differential evolution population size is used.

	:return: a tuple containing 8 lists: 

	return[0][i][j][Q,R,RxQ^4] <- reflectivity list
	return[1][i][j][z,sld,sd sld*] <- sld list
	return[2][i][j][z,solv, sd solv*] <- solvent volume fraction list
	return[3][k][mean, sd*] <- global_parameter fit list
	return[4][l][i][mean, sd*] <- multi_parameter fit list
	return[5][i][mean, sd*] <- background fit list
	return[6][i][mean, sd*] <- scale fit list
	return[7][i][value] <- chi squared list 

	* if available 

	where i the curve number starting from 0
	      j the model number starting from 0
	      k the global_parameter number starting from 0
	      l the multi_parameter number starting from 0

	 In case of error a string that describes the error that occurred, is returned.

	 Output files are written in the created directory 'project' and include
	 a log file, reflectivity curves in R vs Q and R vs Q^4, solvent volume 
	 fraction and scattering length density profiles. Also corresponding PDF
	 figures are saved together with the ASCII data files.
	"""

	#Increase recursion depth for Sympy
	#sys.setrecursionlimit(100000)

	# Multitasking only on macOS and Linux
	if os.name == 'posix':
		if multiprocessing.get_start_method(allow_none=True) != 'fork':
			multiprocessing.set_start_method('fork') # This is needed for Pyhton versions above 3.7!!!!!
		mp=-1
	else:
		mp=1

	warnings.filterwarnings("ignore", category=RuntimeWarning)  # This to filter RuntimeWarning especially from numdifftools
	num_curves=np.size(in_file,0)


	if method == 'mcmc':
		mcmc=True
	else:
		mcmc=False

	if method == 'bootstrap':
		bootstraping = True
	else:
		bootstraping = False


	layers_max = system.copy() # new version 1.1, stop using min/max list input but keeping functionality in the code
	layers_min = system.copy()

	print('--------------------------------------------------------------------')
	print('Program ANAKLASIS - Fit Module for X-ray/Neutron reflection datasets')
	print('version 1.2, February 2021')
	print('developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ')
	print('for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de')
	print('--------------------------------------------------------------------')

	print_buffer=''
	print('Project name: '+project)
	if experror == True:
		if fit_mode == 0:
			print('Using chi squared with errors figure of merit (FOM)')
			print_buffer=print_buffer+'Using chi squared with errors figure of merit (FOM)\n'
		else:
			print('Using log10 with errors figure of merit (FOM)')
			print_buffer=print_buffer+'Using log10 with errors figure of merit (FOM)\n'
	else:
		if fit_mode == 0:
			print('Using chi / R squared figure of merit (FOM)')
			print_buffer=print_buffer+'Using chi / R squared figure of merit (FOM)\n'
		else:
			print('Using log10 figure of merit (FOM)')
			print_buffer=print_buffer+'Using log10 figure of merit (FOM)\n'		

	if fast == True:
		print('Fast fit mode. Using small differential evolution population size')
		print_buffer=print_buffer+'Fast fit mode. Using small differential evolution population size\n'
	else:
		print('Slow fit mode. Using large differential evolution population size')
		print_buffer=print_buffer+'Slow fit mode. Using large differential evolution population size\n'

	if experror == False and method == 'mcmc':
		method = 'simple'
		print('Cannot perform Markov Chain Monte Carlo without provided exprimental errors (dR)')
		print_buffer=print_buffer+'Cannot perform Markov Chain Monte Carlo without provided exprimental errors (dR)\n'

	if experror == False and method == 'bootstrap':
		method = 'simple'
		print('Cannot perform Bootstrap without provided exprimental errors (dR)')
		print_buffer=print_buffer+'Cannot perform Bootstrap without provided exprimental errors (dR)\n'

	if experror == False:
		print('No parameter uncertainty estimation will be performed')
		print_buffer=print_buffer+'No parameter uncertainty estimation will be performed\n'

	if method == 'mcmc':
		print('Parameter uncertainity calculation using Markov Chain Monte Carlo')
		print_buffer=print_buffer+'Parameter uncertainity calculation using Markov Chain Monte Carlo\n'

	if method == 'bootstrap':
		print('Parameter uncertainity calculation using Bootstrap method')
		print_buffer=print_buffer+'Parameter uncertainity calculation using Bootstrap method\n'

	if method != 'mcmc' and method !='bootstrap' and experror == True:
		print('Parameter uncertainity calculation by Hessian matrix estimation')
		print_buffer=print_buffer+'Parameter uncertainity calculation by Hessian matrix estimation\n'

	if plot == True:
		print('A plot summarizing results will pop-up after the end of the calculation')
		print_buffer=print_buffer+'A plot summarizing results will pop-up after the end of the calculation\n'
	else:
		print('No plot will pop-up after the end of the calculation')
		print_buffer=print_buffer+'No plot will pop-up after the end of the calculation\n'

	print('\n')

	if project != 'none':
		project.strip('/')
		folder="project-"+project
		if not os.path.exists(folder):
	   		os.makedirs(folder)
		else:
			print('Directory already exists.. overwriting data..')
			print_buffer=print_buffer+'Directory already exists.. overwriting data..\n'

	for i in range(0,num_curves):	
		if os.path.isfile(in_file[i]):
			if project != 'none': shutil.copy(in_file[i], folder+'/input_curve'+str(i+1)+'.dat')
		else:
			print('error: input file '+str(in_file[i])+' does not exit in current directory..')
			return 'error: input file '+str(in_file[i])+' does not exit in current directory..'

	#exit_code = process.wait()	

	# Load data and sort according to asceding Q
	data=[]
	for i in range(0,num_curves):
		if experror == True:
			if resolution[i] == -1:
				data.append(np.loadtxt(in_file[i], usecols = (0,1,2,3),comments = "#"))
			else:
				data.append(np.loadtxt(in_file[i], usecols = (0,1,2),comments = "#"))
				zerocol = np.zeros((np.size(data[i],0),1))
				data[i] = np.append(data[i], zerocol, axis=1)
				data[i][:,3]=resolution[i]
		else:
			if resolution[i] == -1:
				# data.append(np.loadtxt(in_file[i], usecols = (0,1),comments = "#"))
				# rescol = np.loadtxt(in_file[i], usecols = (3),comments = "#")
				# zerocol = np.zeros((np.size(data[i],0),1))
				# data[i] = np.append(data[i], zerocol, axis=1)
				# data[i] = np.append(data[i], rescol, axis=1)
				data.append(np.loadtxt(in_file[i], usecols = (0,1,3),comments = "#"))
				zerocol = np.zeros((np.size(data[i],0),1))
				data[i] = np.insert(data[i], 2,zerocol, axis=1)
			else:
				data.append(np.loadtxt(in_file[i], usecols = (0,1),comments = "#"))
				zerocol = np.zeros((np.size(data[i],0),1))
				data[i] = np.append(data[i], zerocol, axis=1)
				data[i] = np.append(data[i], zerocol, axis=1)	
				data[i][:,3]=resolution[i]
		
		#data.append(np.loadtxt(in_file[i], usecols = (0,1,2),comments = "#"))
		#print(data)


	#Check data
	for i in range(0,num_curves):
		for j in range(np.size(data[i],0)):
			if data[i][j,0] <= 0:
				print('Invalid wave-vector tranfer (Q) value present in input data!')
				return 'Invalid wave-vector tranfer (Q) value present in input data!'

			if experror == True and data[i][j,2] <= 0:
				print('Invalid experimental error (dR) value present in input data!')
				return 'Invalid experimental error (dR) value present in input data!'

			if resolution[i] == -1 and data[i][j,3] <0:
				print('Invalid (dQ) value present in input data!')
				return 'Invalid (dQ) value present in input data'
	# Check units
	if np.size(units,0) != num_curves:
		print('Incosistent number of entries for wavevector tranfer unints!')
		return 'Incosistent number of entries for wavevector tranfer unints!'

	for i in range(num_curves):
		if units[i] != 'A' and units[i] != 'a' and units[i] != 'nm' and units[i] != 'NM':
			print('Q unints can only be in inverse Angstroms (A) or nanometers (nm)!')
			return 'Q unints can only be in inverse Angstroms (A) or nanometers (nm)!'


		#Delete data with negative reflectivity
		row_to_delete=[]
		for j in range(np.size(data[i],0)):
			if data[i][j,1] <= 0: #or data[i][j,0] > qmax:
				row_to_delete.append(j)

		data[i] = np.delete(data[i], row_to_delete, axis=0)		

	for i in range(0,num_curves):
		data[i] = data[i][data[i][:,0].argsort()]

	for i in range(0,num_curves):
		if units[i] == 'nm' or units[i] == 'NM':
			for j in range(np.size(data[i],0)):
				data[i][j,0]=data[i][j,0]/10.0
				data[i][j,3]=data[i][j,3]/10.0


	#Check of defined model

	if np.size(resolution,0) != num_curves:
		print('Incosistent number of entries for instrumental resolution!')
		return 'Incosistent number of entries for instrumental resolution!'

	if np.size(background,0) != num_curves:
		print('Incosistent number of entries for instrumental background!')
		return 'Incosistent number of entries for instrumental background!'

	for i in range(np.size(background,0)):
		if np.size(background[i]) != 3:
			print('Defined background for each data curve needs exactly two entries (min,max,uniform type) or (mean,sd,normal type)!')
			return 'Defined background for each data curve needs exactly two entries (min,max,uniform type) or (mean,sd,normal type)!'

	if np.size(scale,0) != num_curves:
		print('Incosistent number of entries for reflectivity scaling!')
		return 'Incosistent number of entries for reflectivity scaling!'	

	for i in range(np.size(scale,0)):
		if np.size(scale[i]) != 3:
			print('Defined scale for each data curve needs exactly two entries (min,max,uniform type) or (mean,sd,normal type)!')
			return 'Defined scale for each data curve needs exactly two entries (min,max,uniform type) or (mean,sd,normal type)!'

	if fit_mode !=0 and fit_mode !=1:
		print('Fit mode can only be equal to 0 (linear) or 1 (log10)!')
		return 'Fit mode can only be equal to 0 (linear) or 1 (log10)!'	

	if len(layers_min) != len(layers_max):
		print('Defined min/max models do not have the same number!')
		return 'Defined min/max models do not have the same number!'			

	for j in range(0,len(layers_min)):
		if np.size(layers_min[j],0) != np.size(layers_max[j],0):
			print('Defined min/max models do not have the same number of layers!')
			return 'Defined min/max models do not have the same number of layers!'		

	for j in range(0,len(layers_min)):
		for i in range(np.size(layers_min[j],0)):
			if np.size(layers_min[j][i]) != 6:
				print('Defined min layer model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1))
				print('correct syntax is: [ real sld, imaginary sld, thickness, roughness, solvent volume fraction, name],')
				return 'Defined min layer model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1)				
			if not isinstance(layers_min[j][i][5], str):
				print('Name entry should be a string, min model layer #'+str(i)+', patch #'+str(j+1))
				return 'Name entry should be a string, min model layer #'+str(i)+', patch #'+str(j+1)	

	for j in range(0,len(layers_max)):
		for i in range(np.size(layers_max[j],0)):
			if np.size(layers_max[j][i]) != 6:
				print('Defined max layer model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1))
				print('correct syntax is: [ real sld, imaginary sld, thickness, roughness, solvent volume fraction, name],')
				return 'Defined max layer model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1)				
			if not isinstance(layers_max[j][i][5], str):
				print('Name entry should be a string, max model layer #'+str(i)+', patch #'+str(j+1))
				return 'Name entry should be a string, max model layer #'+str(i)+', patch #'+str(j+1)		

	for k in range(0,len(layers_max)):
		for i in range(np.size(layers_min[k],0)):
			if not isinstance(layers_min[k][i][4], str):
				if layers_min[k][i][4] < 0 or layers_min[k][i][4] > 1:
					print('Invalid solvent volume fraction in layer #'+str(i)+', patch #'+str(k+1))
					print('it should be between 0 and 1!')
					return 'Invalid solvent volume fraction in layer #'+str(i)+', patch #'+str(k+1)
			for j in range(6):
				if isinstance(layers_min[k][i][j], str):
					if layers_min[k][i][j] != layers_max[k][i][j]:
						print('Parametric and name entries in the model should be identical!')
						print('Check layer #'+str(i)+', entry number '+str(j+1)+', patch #'+str(k+1))
						return 'Parametric and name entries in the model should be identical!'
				if not isinstance(layers_min[k][i][j], str):
					if layers_min[k][i][j] > layers_max[k][i][j]:	
						print('Min value cannot be larger than max value!')
						print('Check layer #'+str(i)+', entry number '+str(j+1)+', patch #'+str(k+1))
						return 'Min value cannot be larger than max value!'


	if np.size(model_param,0) > 40:
		print('maximum number of model parameters is equal to 40')
		return 'maximum number of model parameters is equal to 40'

	for i in range(np.size(model_param,0)):
		if model_param[i][0] != 'p'+str(i):
			print('parameter #'+str(i)+' should be named p'+str(i))
			return 'parameter #'+str(i)+' should be named p'+str(i)
		if not isinstance(model_param[i][3], str):
			print('description of parameter #'+str(i)+' should be a string!')
			return 'description of parameter #'+str(i)+' should be a string!'	
		if isinstance(model_param[i][1], str) and model_param[i][-1] == 'uniform':	
			print('min value of parameter #'+str(i)+' should be a number!')
			return 'min value of parameter #'+str(i)+' should be a number!'	
		if isinstance(model_param[i][2], str) and model_param[i][-1] == 'uniform':	
			print('max value of parameter #'+str(i)+' should be a number!')
			return 'max value of parameter #'+str(i)+' should be a number!'	
		if model_param[i][1] > model_param[i][2] and model_param[i][-1] == 'uniform':	
			print('min value of parameter #'+str(i)+' should be larger than max value!')
			return 'min value of parameter #'+str(i)+' should be larger than max value!'
		if isinstance(model_param[i][1], str) and model_param[i][-1] == 'normal':	
			print('mean value of parameter #'+str(i)+' should be a number!')
			return 'mean value of parameter #'+str(i)+' should be a number!'	
		if isinstance(model_param[i][2], str) and model_param[i][-1] == 'normal':	
			print('sd value of parameter #'+str(i)+' should be a number!')
			return 'sd value of parameter #'+str(i)+' should be a number!'		
		if np.size(model_param[i]) != 5:
			print('Number of entries for p'+str(i)+' is wrong!')
			print('correct syntax is: [parameter name, min or mean, max or sd, parameter description, type],')
			return 'Number of entries for p'+str(i)+' is wrong!'

	if np.size(multi_param,0) > 40:
		print('maximum number of multi parameters is equal to 40')
		return 'maximum number of multi parameters is equal to 40'
	for i in range(np.size(multi_param,0)):
		if multi_param[i][0] != 'm'+str(i):
			print('multi parameter #'+str(i)+' should be named m'+str(i))
			return 'multi parameter #'+str(i)+' should be named m'+str(i)
		if np.size(multi_param[i]) != 2*num_curves+3:
			print('Number of entries for m'+str(i)+' is wrong!')
			print('correct syntax is: [parameter name, min1 or mean1, max1 or sd1, ..., minN or meanN, maxN or sdN, parameter description, type uniform or normal],')
			print('where N is the number of input curves.')
			return 'Number of entries for m'+str(i)+' is wrong!'			
		if not isinstance(multi_param[i][-2], str):
			print('description of multi parameter #'+str(i)+' should be a string!')
			return 'description of multi parameter #'+str(i)+' should be a string!'	
		for j in range(1,2*num_curves):
			if isinstance(multi_param[i][j], str):
				print('multi parameter #'+str(i)+' min/max or mean/sd values should be numbers')
				return 'multi parameter #'+str(i)+' min/max or mean/sd values should be numbers'
		for j in range(1, 2*num_curves,2):
			if multi_param[i][j] > multi_param[i][j+1] and multi_param[i][-1] == 'uniform':
				print('multi parameter #'+str(i)+', min value #'+str(j)+' cannot be larger than max value #'+str(j+1))
				return 'multi parameter #'+str(i)+', min value #'+str(j)+' cannot be larger than max value #'+str(j+1)									

	for i in range(np.size(model_constraints)):
		if not '>' in model_constraints[i] and not '<' in model_constraints[i]:
			print('Constraint expression #'+str(i+1)+' is not an inequality!')
			return 'Constraint expression #'+str(i+1)+' is not an inequality!'

		if '=' in model_constraints[i]:
			print('Constraint expression #'+str(i+1)+' is not a pure inequality!')
			return 'Constraint expression #'+str(i+1)+' is not a pure inequality!'

	bt_iter=1000


	if project != 'none': os.chdir(folder)
	start_time = time.time()

	for i in range(0,num_curves):
		print("file#"+str(i+1)+" experimental points: "+str(np.size(data[i],0)))
		print('Q units in inverse '+units[i])
		if resolution[i] != -1:
			print('dQ/Q = ',resolution[i])
		else:
			print('dQ/Q pointwise')
		print("fit weight: ",fit_weight[i])
	np.set_printoptions(precision=3)


	# These values of the semi-infinite fronting and backing have no physical meaning 
	for i in range(len(layers_min)):
		layers_min[i][0][2]=0.0
		layers_min[i][-1][2]=0.0
		layers_min[i][-1][3]=0.0

		layers_max[i][0][2]=0.0
		layers_max[i][-1][2]=0.0
		layers_max[i][-1][3]=0.0

	free_param=0
	list_free_param=''
	bounds =[]
	corner_bounds=[]
	corner_labels=[]
	global layer_fun
	layer_fun=[]
	ii,jj,kk=sympy.symbols('ii jj kk', integer=True)
	x,y=sympy.symbols('x y')
	p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n=sympy.symbols('p0 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 n')
	m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39=sympy.symbols('m0 m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 m20 m21 m22 m23 m24 m25 m26 m27 m28 m29 m30 m31 m32 m33 m34 m35 m36 m37 m38 m39')


	for k in range(len(layers_min)):
		sub_layer_fun=[]
		for i in range(0,np.size(layers_min[k],0)):
			line=[]
			for j in range(0,5):
				if isinstance(layers_min[k][i][j], str):
					bounds.append((0,0))
					expr = layers_max[k][i][j]
					line.append((sympy.lambdify([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39],expr,'numpy')))
				else:
					bounds.append((layers_min[k][i][j],layers_max[k][i][j]))
					expr = '0'
					line.append((sympy.lambdify([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39],expr,'numpy')))
				if layers_min[k][i][j] != layers_max[k][i][j]:
					free_param=free_param+1
					if j == 0: 
						if len(layers_min) == 1:
							corner_labels.append(r"$\Re sld$-"+layers_max[k][i][5])
							list_free_param=list_free_param+'Re(sld) '+layers_max[k][i][5]+','
						else:
							corner_labels.append(r"$\Re sld$-"+layers_max[k][i][5]+'_patch #'+str(k+1))
							list_free_param=list_free_param+'Re(sld) '+layers_max[k][i][5]+'_patch #'+str(k+1)+','
					if j == 1:
						if len(layers_min) == 1: 
							corner_labels.append(r"$\Im sld-$ "+layers_max[k][i][5])
							list_free_param=list_free_param+'Im(sld) '+layers_max[k][i][5]+','							
						else:
							corner_labels.append(r"$\Im sld-$ "+layers_max[k][i][5]+'_patch #'+str(k+1))
							list_free_param=list_free_param+'Im(sld) '+layers_max[k][i][5]+'_patch #'+str(k+1)+','
					if j == 2:
						if len(layers_min) == 1:  
							corner_labels.append(r"$d-$ "+layers_max[k][i][5])
							list_free_param=list_free_param+'thickness '+layers_max[k][i][5]+','
						else:
							corner_labels.append(r"$d-$ "+layers_max[k][i][5]+'_patch #'+str(k+1))
							list_free_param=list_free_param+'thickness '+layers_max[k][i][5]+'_patch #'+str(k+1)+','							
					if j == 3:
						if len(layers_min) == 1:   
							corner_labels.append(r"$\sigma-$ "+layers_max[k][i][5]+"/"+layers_max[k][i+1][5])
							list_free_param=list_free_param+'roughness '+layers_max[k][i][5]+"/"+layers_max[k][i+1][5]+','
						else:
							corner_labels.append(r"$\sigma-$ "+layers_max[k][i][5]+"/"+layers_max[k][i+1][5]+'_patch #'+str(k+1))
							list_free_param=list_free_param+'roughness '+layers_max[k][i][5]+"/"+layers_max[k][i+1][5]+'_patch #'+str(k+1)+','
					if j == 4: 
						if len(layers_min) == 1:   
							corner_labels.append(r"$\phi-$ "+layers_max[k][i][5])
							list_free_param=list_free_param+'solvent volume fraction '+layers_max[k][i][5]+','
						else:
							corner_labels.append(r"$\phi-$ "+layers_max[k][i][5]+'_patch #'+str(k+1))
							list_free_param=list_free_param+'solvent volume fraction '+layers_max[k][i][5]+'_patch #'+str(k+1)+','
					corner_bounds.append((layers_min[k][i][j],layers_max[k][i][j]))
			sub_layer_fun.append(line)
		layer_fun.append(sub_layer_fun)

	# This function definition needs to be here in case this function in the main program
	# def f_layer_fun(i,j,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,n,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9):
	# 	return np.float((layer_fun[i][j](p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,n,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9)))

	global left_fun
	global right_fun
	left_fun=[]
	right_fun=[]
	center_fun=[]
	for i in range(0,len(model_constraints)):
		if '>' in str(model_constraints[i]) or '<' in str(model_constraints[i]):
			if '>' in str(model_constraints[i]):
				left_expr = model_constraints[i].split('>')[0]
				right_expr = model_constraints[i].split('>')[1]
				left_fun.append((sympy.lambdify([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39],left_expr,'numpy')))
				right_fun.append((sympy.lambdify([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39],right_expr,'numpy')))
				center_fun.append('>')



			if '<' in str(model_constraints[i]):
				left_expr = model_constraints[i].split('<')[0]
				right_expr = model_constraints[i].split('<')[1]
				left_fun.append((sympy.lambdify([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39],left_expr,'numpy')))
				right_fun.append((sympy.lambdify([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39],right_expr,'numpy')))
				center_fun.append('<')

# Needs to be here in case this function acts as the main program
# def f_left_fun(i,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9):
# 	return np.float((left_fun[i](p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9)))

# def f_right_fun(i,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9):
# 	return np.float((right_fun[i](p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9)))

	m_param=0
	for i in range(0,np.size(model_param,0)):
		m_param=m_param+1
		if model_param[i][4] != 'uniform' and model_param[i][4] != 'normal':
			print('Type of global parameter p'+str(i)+' can only be uniform or normal!')
			return 'Type of global parameter p'+str(i)+' can only be uniform or normal!'		 
		if model_param[i][4] == 'uniform':
			bounds.append((model_param[i][1],model_param[i][2]))
			if model_param[i][1] != model_param[i][2]:
				free_param=free_param+1
				#corner_labels.append(r"param-"+model_param[i][3])
				corner_labels.append(r""+model_param[i][3])
				corner_bounds.append((model_param[i][1],model_param[i][2]))
				list_free_param=list_free_param+'p'+str(i)+' ('+model_param[i][3]+'),'
		if model_param[i][4] == 'normal':
			bounds.append((model_param[i][1]-3*model_param[i][2],model_param[i][1]+3*model_param[i][2]))
			if model_param[i][2] != 0:
				free_param=free_param+1
				#corner_labels.append(r"param-"+model_param[i][3])
				corner_labels.append(r""+model_param[i][3])
				corner_bounds.append((model_param[i][1]-3*model_param[i][2],model_param[i][1]+3*model_param[i][2]))
				list_free_param=list_free_param+'p'+str(i)+' ('+model_param[i][3]+'),'

	for i in range(0,num_curves):
		if background[i][2] != 'uniform' and background[i][2] != 'normal':
			print('Type of background for curve #'+str(i+1)+' can only be uniform or normal!')
			return 'Type of background for curve #'+str(i+1)+' can only be uniform or normal!'
		if background[i][2] == 'uniform':	
			bounds.append((background[i][0],background[i][1]))
			if background[i][0] != background[i][1]:
				free_param=free_param+1
				corner_labels.append(r"bkg_curve"+str(i+1))
				corner_bounds.append((background[i][0],background[i][1]))
				list_free_param=list_free_param+'bkg curve#'+str(i+1)+','
		if background[i][2] == 'normal':	
			bounds.append((background[i][0]-3*background[i][1],background[i][0]+3*background[i][1]))
			if background[i][1] != 0:
				free_param=free_param+1
				corner_labels.append(r"bkg_curve"+str(i+1))
				corner_bounds.append((background[i][0]-3*background[i][1],background[i][0]+3*background[i][1]))
				list_free_param=list_free_param+'bkg curve#'+str(i+1)+','

	for i in range(0,num_curves):
		if scale[i][2] != 'uniform' and scale[i][2] != 'normal':
			print('Type of scale for curve #'+str(i+1)+' can only be uniform or normal!')
			return 'Type of scale for curve #'+str(i+1)+' can only be uniform or normal!'
		if scale[i][2] == 'uniform':
			bounds.append((scale[i][0],scale[i][1]))
			if scale[i][0] != scale[i][1]:
				free_param=free_param+1
				corner_labels.append(r"scale_curve"+str(i+1))
				corner_bounds.append((scale[i][0],scale[i][1]))
				list_free_param=list_free_param+'scale curve#'+str(i+1)+','
		if scale[i][2] == 'normal':
			bounds.append((scale[i][0]-3*scale[i][1],scale[i][0]+3*scale[i][1]))
			if scale[i][1] != 0:
				free_param=free_param+1
				corner_labels.append(r"scale_curve"+str(i+1))
				corner_bounds.append((scale[i][0]-3*scale[i][1],scale[i][0]+3*scale[i][1]))
				list_free_param=list_free_param+'scale curve#'+str(i+1)+','				

	c_param=0
	for i in range(0,np.size(multi_param,0)):
		c_param=c_param+1
		if multi_param[i][-1] != 'uniform' and multi_param[i][-1] != 'normal':
			print('Type of multi parameter m'+str(i)+' can only be uniform or normal!')
			return 'Type of multi parameter m'+str(i)+' can only be uniform or normal!'	
		if multi_param[i][-1] == 'uniform':
			for j in range(0,num_curves):
				bounds.append((multi_param[i][2*j+1],multi_param[i][2*j+2]))
				if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
					free_param = free_param + 1
					#corner_labels.append(r"param-"+multi_param[i][2*num_curves+1]+"_curve"+str(j+1))
					corner_labels.append(r""+multi_param[i][2*num_curves+1]+"_curve"+str(j+1))
					corner_bounds.append((multi_param[i][2*j+1],multi_param[i][2*j+2]))	
					list_free_param=list_free_param+'m'+str(i)+' curve#'+str(j+1)+' ('+multi_param[i][-2]+'),'			
		if multi_param[i][-1] == 'normal':
			for j in range(0,num_curves):
				bounds.append((multi_param[i][2*j+1]-3*multi_param[i][2*j+2],multi_param[i][2*j+1]+3*multi_param[i][2*j+2]))
				if multi_param[i][2*j+2] != 0:
					free_param = free_param + 1
					#corner_labels.append(r"param-"+multi_param[i][2*num_curves+1]+"_curve"+str(j+1))
					corner_labels.append(r""+multi_param[i][2*num_curves+1]+"_curve"+str(j+1))
					corner_bounds.append((multi_param[i][2*j+1]-3*multi_param[i][2*j+2],multi_param[i][2*j+1]+3*multi_param[i][2*j+2]))	
					list_free_param=list_free_param+'m'+str(i)+' curve#'+str(j+1)+' ('+multi_param[i][-2]+'),'	

	list_free_param=list_free_param+'\n'
	#population_size=int(15*free_param/(5*np.size(layers_min,0)))
	mnlayers=[]
	for i in range(len(layers_min)):
		mnlayers.append(np.size(layers_min[i],0))
	population_size=int(15*(free_param)/(m_param+c_param+(5*sum(mnlayers))))
	if population_size == 0:
		population_size = 1

	if fast == True:
		population_size = 1*population_size
	else:
		population_size = 5*population_size

	print('free parameters = ',free_param)
	print('\n')
	print('list of free parameters: '+list_free_param)
	#print('pop_size = ',population_size)


	print('Running differential evolution minimization...')
	mnlayers=[]
	for i in range(len(layers_min)):
		mnlayers.append(np.size(layers_min[i],0))
	results = differential_evolution(fig_of_merit_sym, bounds, (data, resolution, mnlayers, m_param,layers_max,model_param,model_constraints,fit_mode,num_curves,c_param,multi_param,f_layer_fun,f_left_fun,f_right_fun,center_fun,experror,fit_weight,patches),tol=0.1,popsize=population_size,maxiter=200,polish=True,updating='deferred',workers=mp)
	de_results = results.x.copy()

	print('/n')
	print('Success: '+str(results.success))
	print('Number of evaluation: '+str(results.nfev))
	print('Number of iterations: '+str(results.nit))
	print('FOM: '+str(results.fun))


	#if not bootstraping and not mcmc and experror==True:
	if not mcmc and experror==True:
		print('Estimating Hessian matrix...')
		mnlayers=[]
		for i in range(len(layers_min)):
			mnlayers.append(np.size(layers_min[i],0))
		df=nd.Hessdiag(reduced_chi_sym,num_steps=50)(de_results,data,resolution,mnlayers,m_param,layers_max,model_param,fit_mode,num_curves,c_param,multi_param,f_layer_fun,fit_weight,patches,free_param)

		df_counter=0
		df_layers=[]
		df_vm=[]
		df_vp=[]
		df_bkg=[]
		df_scl=[]

		df_init=[]

		for k in range(len(layers_min)):
			sub_df_layers=[]
			for i in range(0,np.size(layers_min[k],0)):
				line=[]
				for j in range(0,5):
					if layers_min[k][i][j] != layers_max[k][i][j]:
						line.append(np.sqrt(2.0/df[df_counter]))
						df_init.append(np.sqrt(2.0/df[df_counter]))
					else:
						line.append(np.float(0.0))
						df_init.append(np.float(0.0))
					df_counter = df_counter + 1
				sub_df_layers.append(line)
			df_layers.append(sub_df_layers)

		for i in range(0,m_param):
			if model_param[i][1] != model_param[i][2]:
		 		df_vp.append(np.sqrt(2.0/df[df_counter]))
		 		df_init.append(np.sqrt(2.0/df[df_counter]))
			else:
		 		df_vp.append(np.float(0.0))
		 		df_init.append(np.float(0.0))
			df_counter = df_counter + 1

		for i in range(0,num_curves):
			if background[i][0] != background[i][1]:
				df_bkg.append(np.sqrt(2.0/df[df_counter]))
				df_init.append(np.sqrt(2.0/df[df_counter]))
			else:
				df_bkg.append(np.float(0.0))
				df_init.append(np.float(0.0))
			df_counter = df_counter + 1

		for i in range(0,num_curves):
			if scale[i][0] != scale[i][1]:
				df_scl.append(np.sqrt(2.0/df[df_counter]))
				df_init.append(np.sqrt(2.0/df[df_counter]))
			else:
				df_scl.append(np.float(0.0))
				df_init.append(np.float(0.0))
			df_counter = df_counter + 1

		for i in range(0,np.size(multi_param,0)):
			for j in range(0,num_curves):
				if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
					df_vm.append(np.sqrt(2.0/df[df_counter]))
					df_init.append(np.sqrt(2.0/df[df_counter]))
				else:
					df_vm.append(np.float(0.0))
					df_init.append(np.float(0.0))
				df_counter = df_counter + 1



		boot_bounds=[]
		for i in range((5*sum(mnlayers)+m_param+c_param*num_curves+2*num_curves)):
			if not np.isnan(df_init[i]):
				if df_init[j] != 0.0:
					lower=de_results[i]-3*df_init[i]
					if lower < bound[i][0]: lower = bound[i][0]
					upper=de_results[i]+3*df_init[i]
					if upper > bound[i][1]: upper = bound[i][1]
					boot_bounds.append((lower,upper))
				else:
					boot_bounds.append(bounds[i])
			else:
				boot_bounds.append(bounds[i])


	q_bin=[]
	res_bin=[]
	Refl2=[]
	Refl=[]
	Profile=[]
	Solvent=[]
	chi_s=[]

	for curve in range(0,num_curves):

		vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,m_param):
			vp[i]=results.x[sum(mnlayers)*5+i]

		vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,c_param):
			vm[i]=results.x[sum(mnlayers)*5+m_param+2*num_curves+num_curves*i+curve]


		layers=[]
		for k in range(len(layers_min)):
			sub_layers=[]
			for i in range(0,np.size(layers_min[k],0)):
				line=[]
				for j in range(0,5):
					if isinstance(layers_max[k][i][j], str):
						line.append(np.float((f_layer_fun(k,i,j,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],i,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))

					else:
						if k > 0:
							line.append(results.x[5*sum(mnlayers[0:k])+i*5+j])
						else:
							line.append(results.x[i*5+j])
				line.append(layers_max[k][i][5])
				sub_layers.append(line)
			layers.append(sub_layers)

		bkg=results.x[sum(mnlayers)*5+m_param+curve]
		scl=results.x[sum(mnlayers)*5+m_param+num_curves+curve]


		q_bin.append(np.linspace(np.min(data[curve][:,0]), np.max(data[curve][:,0]), 1001))
		res_bin.append(np.zeros((np.size(q_bin[curve],0),1)))
		for i in range(np.size(res_bin[curve],0)):
			if resolution[curve] != -1:
				res_bin[curve][i] = resolution[curve]
			else:
				res_bin[curve][i] = np.interp(q_bin[curve][i], data[curve][:,0],data[curve][:,3])


		Refl2.append(Reflectivity(data[curve][:,0], data[curve][:,3], layers, resolution[curve],bkg,scl,patches)) # needed for bootstrapping


		Refl.append(Reflectivity(q_bin[curve], res_bin[curve], layers, resolution[curve],bkg,scl,patches))

		Profile.append(profile(layers, 1000))

		Solvent.append(solvent_penetration(layers, 1000))

		chi_s.append(chi_square(data[curve], layers, resolution[curve], bkg, scl,patches))




	if mcmc:
		#Create sol tuple containging the fitted solution if the parameters
		print('Running Markov Chain Monte Carlo sampling...')
		hsol=[]
		min_hsol=[]
		max_hsol=[]
		count=0
		for k in range(len(layers_min)):
			for i in range(0,np.size(layers_min[k],0)):
				for j in range(0,5):
					if layers_min[k][i][j] != layers_max[k][i][j]:
						if k > 0:
							hsol.append(results.x[5*sum(mnlayers[0:k])+i*5+j])
						else:
							hsol.append(results.x[i*5+j])
						min_hsol.append(bounds[count][0])
						max_hsol.append(bounds[count][1])
					count=count+1

		for i in range(0,m_param):
			if model_param[i][4] == 'uniform':
				if model_param[i][1] != model_param[i][2]:
					hsol.append(results.x[sum(mnlayers)*5+i])
					min_hsol.append(bounds[count][0])
					max_hsol.append(bounds[count][1])
			if model_param[i][4] == 'normal':
				if model_param[i][2] != 0.0:
					hsol.append(results.x[sum(mnlayers)*5+i])
					min_hsol.append(bounds[count][0])
					max_hsol.append(bounds[count][1])
			count=count+1

		for i in range(0,num_curves):
			if background[i][2] == 'uniform':
				if background[i][0] != background[i][1]:
					hsol.append(results.x[sum(mnlayers)*5+m_param+i])
					min_hsol.append(bounds[count][0])
					max_hsol.append(bounds[count][1])
			if background[i][2] == 'normal':
				if background[i][1] != 0.0:
					hsol.append(results.x[sum(mnlayers)*5+m_param+i])
					min_hsol.append(bounds[count][0])
					max_hsol.append(bounds[count][1])
			count=count+1

		for i in range(0,num_curves):
			if scale[i][2] == 'uniform':
				if scale[i][0] != scale[i][1]:
					hsol.append(results.x[sum(mnlayers)*5+m_param+num_curves+i])
					min_hsol.append(bounds[count][0])
					max_hsol.append(bounds[count][1])
			if scale[i][2] == 'normal':
				if scale[i][1] != 0.0:
					hsol.append(results.x[sum(mnlayers)*5+m_param+num_curves+i])
					min_hsol.append(bounds[count][0])
					max_hsol.append(bounds[count][1])
			count=count+1

		for i in range(0,np.size(multi_param,0)):
			if multi_param[i][-1] == 'uniform':
				for j in range(0,num_curves):
					if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
						hsol.append(results.x[sum(mnlayers)*5+m_param+2*num_curves+c_param*i+j])	
						min_hsol.append(bounds[count][0])
						max_hsol.append(bounds[count][1])
					count=count+1
			if multi_param[i][-1] == 'normal':
				for j in range(0,num_curves):
					if multi_param[i][2*j+2] != 0.0:
						hsol.append(results.x[sum(mnlayers)*5+m_param+2*num_curves+c_param*i+j])	
						min_hsol.append(bounds[count][0])
						max_hsol.append(bounds[count][1])
					count=count+1

		sol=tuple(hsol)


		print('Initializing MCMC walkers..')
		num_walkers=5*free_param
		nsol=np.zeros((num_walkers, free_param))
		for i in range(0,num_walkers):
			while True:
				for j in range(0,free_param):
					while True:
						if hsol[j] != 0.0: 
							nsol[i][j]=np.float(hsol[j])*np.float(np.random.randn(1))
							#nsol[i][j]=np.float(hsol[j])+(np.float(max_hsol[j]-min_hsol[j])/2.0)*np.float(np.random.randn(1))
						else:
							nsol[i][j]=(np.float(max_hsol[j]-min_hsol[j])/2.0)*np.float(np.random.randn(1))


						if min_hsol[j] <= hsol[j]+1e-2*nsol[i][j] <=  max_hsol[j]: 
							break

				if log_prior(hsol[:]+1e-2*nsol[i][:],mnlayers,m_param,layers_max,layers_min,model_param,model_constraints,background,num_curves,c_param,multi_param,scale,f_left_fun,f_right_fun,center_fun) != -np.inf:
					break

		#pos = sol + 1e-4 * np.random.randn(32, free_param)
		pos = sol + 1e-2 * nsol
		nwalkers, ndim = pos.shape

		print('MCMC test run..')
		os.environ["OMP_NUM_THREADS"] = "1"
		mcmc_iter=500
		if mp == -1:
			with Pool() as pool:
				sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,pool=pool, args=(data,resolution,mnlayers,m_param,layers_max,layers_min,model_param,model_constraints,background,fit_mode,num_curves,c_param,multi_param,scale,f_layer_fun,f_left_fun,f_right_fun,center_fun,fit_weight,patches))
				sampler.run_mcmc(pos, mcmc_iter, progress=True)
			#pool.join()

		#	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data,resolution,np.size(layers_min,0),m_param,layers_max,layers_min,model_param,model_constraints,background,fit_mode,num_curves,c_param,multi_param))
		#	sampler.run_mcmc(pos, mcmc_iter, progress=True);
				tau = sampler.get_autocorr_time(quiet='True')
				#print('MCMC autocorellation time:')
				#print(tau)
				maxtau=max(tau)
				mcmc_iter=60*int(maxtau)
				#time.sleep(5)
				print('MCMC production run..')
				#sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,pool=pool, args=(data,resolution,np.size(layers_min,0),m_param,layers_max,layers_min,model_param,model_constraints,background,fit_mode,num_curves,c_param,multi_param))
				sampler.run_mcmc(None, mcmc_iter, progress=True)
		else:
			sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data,resolution,mnlayers,m_param,layers_max,layers_min,model_param,model_constraints,background,fit_mode,num_curves,c_param,multi_param,scale,f_layer_fun,f_left_fun,f_right_fun,center_fun,fit_weight,patches))
			sampler.run_mcmc(pos, mcmc_iter, progress=True)


	#	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data,resolution,np.size(layers_min,0),m_param,layers_max,layers_min,model_param,model_constraints,background,fit_mode,num_curves,c_param,multi_param))
	#	sampler.run_mcmc(pos, mcmc_iter, progress=True);
			tau = sampler.get_autocorr_time(quiet='True')
			#print('MCMC autocorellation time:')
			#print(tau)
			maxtau=max(tau)
			mcmc_iter=60*int(maxtau)
			#time.sleep(5)
			print('MCMC production run..')
			#sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,pool=pool, args=(data,resolution,np.size(layers_min,0),m_param,layers_max,layers_min,model_param,model_constraints,background,fit_mode,num_curves,c_param,multi_param))
			sampler.run_mcmc(None, mcmc_iter, progress=True)
		#flat_samples = sampler.get_chain(discard=500, thin=15, flat=True)
		flat_samples = sampler.get_chain(discard=int(8*maxtau), thin=int(maxtau/2), flat=True)


		samples=flat_samples.copy()

		count=0
		for i in range(0,np.size(flat_samples,0),int(maxtau/2)):
			count=count+1

		r_flat_samples=np.zeros((count,free_param))

		count=0
		for i in range(0,np.size(flat_samples,0),int(maxtau/2)):
			r_flat_samples[count,:]=flat_samples[i,:]
			count=count+1


		print('Calculating mean curves...')
		samples_Refl=np.empty([num_curves,1001, np.size(r_flat_samples,0)+1])
		samples_Refl2=[]
		for i in range(0,num_curves):
			samples_Refl[i,:,0]=q_bin[i][:]
			samples_Refl2.append(np.empty([np.size(data[i],0),np.size(r_flat_samples,0)+1]))

		samples_Profile=np.empty([num_curves,len(layers_min),1000, 2*np.size(r_flat_samples,0)+1])
		samples_Solvent=np.empty([num_curves,len(layers_min),1000, 2*np.size(r_flat_samples,0)+1])


		samples[:,:]=flat_samples[:,:]

		Profile_mean=[]
		Solvent_mean=[]

		with tqdm(total=(num_curves)*(np.size(r_flat_samples,0)), file=sys.stdout) as pbar:
			#Nlayers=np.size(layers_min,0)
			for m in range(0,np.size(r_flat_samples,0)):
				mk=0
				for curve in range(0,num_curves):

					k=0
					for l in range(len(layers_min)):
						for i in range(0,mnlayers[l]):
							for j in range(0,5):
								if layers_min[l][i][j] != layers_max[l][i][j]:
									k=k+1


					vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
					l=0
					for i in range(0,m_param):
						if model_param[i][4] == 'uniform':
							if model_param[i][1] != model_param[i][2]:
								vp[i]=r_flat_samples[m,k+l]
								l=l+1
							else:
								vp[i]=model_param[i][1]
						if model_param[i][4] == 'normal':
							if model_param[i][2] != 0:
								vp[i]=r_flat_samples[m,k+l]
								l=l+1
							else:
								vp[i]=model_param[i][1]

					k=k+l

					if background[curve][2] == 'uniform':
						if background[curve][0] != background[curve][1]:
							bkg=r_flat_samples[m,k+curve]
						else:
							bkg=background[curve][0]
					if background[curve][2] == 'normal':
						if background[curve][1] != 0:
							bkg=r_flat_samples[m,k+curve]
						else:
							bkg=background[curve][0]

					for i in range(0,num_curves):
						if background[curve][2] == 'uniform':
							if background[i][0] != background[i][1]:
								k=k+1
						if background[curve][2] == 'normal':
							if background[i][1] != 0:
								k=k+1

					if scale[i][2] == 'uniform':
						if scale[curve][0] != scale[curve][1]:
							scl=r_flat_samples[m,k+curve]
						else:
							scl=scale[curve][0]

					if scale[i][2] == 'normal':
						if scale[curve][1] != 0:
							scl=r_flat_samples[m,k+curve]
						else:
							scl=scale[curve][0]

					if scale[i][2] == 'uniform':
						for i in range(0,num_curves):
							if scale[i][0] != scale[i][1]:
								k=k+1

					if scale[i][2] == 'normal':
						for i in range(0,num_curves):
							if scale[i][1] != 0:
								k=k+1


					vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
					l=0
					for i in range(0,c_param):
						if multi_param[i][-1] == 'uniform':
							if multi_param[i][2*curve+1] != multi_param[i][2*curve+2]:
								vm[i]=r_flat_samples[m,k+mk+l]
								l=l+1
							else:
								vm[i]=multi_param[i][2*curve+1]
						if multi_param[i][-1] == 'normal':
							if multi_param[i][2*curve+2] != 0:
								vm[i]=r_flat_samples[m,k+mk+l]
								l=l+1
							else:
								vm[i]=multi_param[i][2*curve+1]
					mk=mk+l


					layers=[]
					k=0
					for l in range(len(layers_min)):
						sub_layers=[]
						for i in range(0,mnlayers[l]):
							line=[]
							for j in range(0,5):
								if isinstance(layers_max[l][i][j], str):
									line.append(np.float((f_layer_fun(l,i,j,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],i,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))

								else:
									if layers_min[l][i][j] != layers_max[l][i][j]:
										line.append(r_flat_samples[m,k])
										k=k+1
									else:
										line.append(layers_max[l][i][j])
							line.append(layers_max[l][i][5])
							sub_layers.append(line)
						layers.append(sub_layers)


					cRefl=Reflectivity(q_bin[curve],res_bin[curve], layers, resolution[curve], bkg, scl,patches)
					samples_Refl[curve,:,m+1]=cRefl[:,1]
					samples_Refl2[curve][:,0]=data[curve][:,0]
					samples_Refl2[curve][:,m+1]=Reflectivity(data[curve][:,0],data[curve][:,3], layers, resolution[curve], bkg, scl,patches)[:,1]
					cProfile=profile(layers, 1000)
					for j in range(len(layers_min)):
						samples_Profile[curve,j,:,2*m]=cProfile[j][:,0]
						samples_Profile[curve,j,:,2*m+1]=cProfile[j][:,1]
					cSolvent=solvent_penetration(layers, 1000)
					for j in range(len(layers_min)):
						samples_Solvent[curve,j,:,2*m]=cSolvent[j][:,0]
						samples_Solvent[curve,j,:,2*m+1]=cSolvent[j][:,1]
					pbar.update(1)

		print('Plotting curves... please wait...')
		for curve in range(0,num_curves):

			Refl[curve][:,1]=0.0
			Refl2[curve][:,1]=0.0
			for k in range(0,np.size(r_flat_samples,0)):
				Refl[curve][:,1]=Refl[curve][:,1]+samples_Refl[curve,:,k+1]/np.float(np.size(r_flat_samples,0))
				Refl2[curve][:,1]=Refl2[curve][:,1]+samples_Refl2[curve][:,k+1]/np.float(np.size(r_flat_samples,0))


			chi_s[curve]=0.0
			Nexp=np.size(data[curve],0)
			for i in range(0, Nexp):
				chi_s[curve]=chi_s[curve]+(1.0/float(Nexp))*((data[curve][i][1]-Refl2[curve][i,1])/(data[curve][i][2]))**2

			#Calculate mean profiles
			z_max=0
			z_min=0
			for k in range(len(layers_min)):
				for i in range(0,np.size(r_flat_samples,0)):
					for j in range(0,1000):
						if samples_Profile[curve,k,j,2*i] > z_max: z_max=samples_Profile[curve,k,j,2*i]
						if samples_Profile[curve,k,j,2*i] < z_min: z_min=samples_Profile[curve,k,j,2*i]		

			z_bin = np.linspace(z_min, z_max, 1000)
			Profile_mean.append(np.zeros([len(layers_min),len(z_bin), 3]))
			Solvent_mean.append(np.zeros([len(layers_min),len(z_bin), 3]))
			for k in range(len(layers_min)):
				Profile_mean[curve][k][:,0]=z_bin[:]
				Solvent_mean[curve][k][:,0]=z_bin[:]


			for i in range(0,np.size(r_flat_samples,0)):
				for k in range(len(layers_min)):
					for j in range(0,len(z_bin)):
						Profile_mean[curve][k][j,1]=Profile_mean[curve][k][j,1]+np.interp(Profile_mean[curve][k][j,0], samples_Profile[curve,k,:,2*i],samples_Profile[curve,k,:,2*i+1])/float(np.size(r_flat_samples,0))
						Solvent_mean[curve][k][j,1]=Solvent_mean[curve][k][j,1]+np.interp(Solvent_mean[curve][k][j,0], samples_Solvent[curve,k,:,2*i],samples_Solvent[curve,k,:,2*i+1])/float(np.size(r_flat_samples,0))

			for i in range(0,np.size(r_flat_samples,0)):
				for k in range(len(layers_min)):
					for j in range(0,len(z_bin)):
						Profile_mean[curve][k][j,2]=Profile_mean[curve][k][j,2]+((Profile_mean[curve][k][j,1]-np.interp(Profile_mean[curve][k][j,0], samples_Profile[curve,k,:,2*i],samples_Profile[curve,k,:,2*i+1]))**2)
						Solvent_mean[curve][k][j,2]=Solvent_mean[curve][k][j,2]+((Solvent_mean[curve][k][j,1]-np.interp(Solvent_mean[curve][k][j,0], samples_Solvent[curve,k,:,2*i],samples_Solvent[curve,k,:,2*i+1]))**2)

			for j in range(0,len(z_bin)):
				for k in range(len(layers_min)):
					Profile_mean[curve][k][j,2]=np.sqrt(Profile_mean[curve][k][j,2]/float(np.size(r_flat_samples,0)-1))
					Solvent_mean[curve][k][j,2]=np.sqrt(Solvent_mean[curve][k][j,2]/float(np.size(r_flat_samples,0)-1))

		samples_truths=[]
		samples_truths_plus=[]
		samples_truths_minus=[]


		sd_vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		count=0
		for k in range(len(layers_min)):
			for i in range(0,np.size(layers_min[k],0)):
				for j in range(0,5):
					if layers_min[k][i][j] != layers_max[k][i][j]:
						layers[k][i][j]=np.mean(r_flat_samples[:,count])
						samples_truths.append(np.mean(r_flat_samples[:,count]))
						samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
						samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))

						count=count+1

		for j in range(0,m_param):
			if model_param[j][4] == 'uniform' and model_param[j][1] != model_param[j][2]:
				vp[j]=np.mean(r_flat_samples[:,count])
				sd_vp[j]=np.std(r_flat_samples[:,count])
				samples_truths.append(np.mean(r_flat_samples[:,count]))
				samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
				samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))

				count=count+1
			if model_param[j][4] == 'normal' and model_param[j][2] != 0:
				vp[j]=np.mean(r_flat_samples[:,count])
				sd_vp[j]=np.std(r_flat_samples[:,count])
				samples_truths.append(np.mean(r_flat_samples[:,count]))
				samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
				samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))

				count=count+1

		mbkg=np.zeros(num_curves)
		sd_bkg=np.zeros(num_curves)

		for i in range(0,num_curves):
			if background[i][2] == 'uniform': 
				if background[i][0] != background[i][1]:
					mbkg[i]=np.mean(r_flat_samples[:,count])
					sd_bkg[i]=np.std(r_flat_samples[:,count])
					samples_truths.append(np.mean(r_flat_samples[:,count]))
					samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
					samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))
					count=count+1
				else:
					mbkg[i]=background[i][0]
					sd_bkg[i]=0.0
			if background[i][2] == 'normal': 
				if background[i][1] != 0:
					mbkg[i]=np.mean(r_flat_samples[:,count])
					sd_bkg[i]=np.std(r_flat_samples[:,count])
					samples_truths.append(np.mean(r_flat_samples[:,count]))
					samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
					samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))
					count=count+1
				else:
					mbkg[i]=background[i][0]
					sd_bkg[i]=0.0

		mscl=np.zeros(num_curves)
		sd_scl=np.zeros(num_curves)

		for i in range(0,num_curves):
			if scale[i][2] == 'uniform': 
				if scale[i][0] != scale[i][1]:
					mscl[i]=np.mean(r_flat_samples[:,count])
					sd_scl[i]=np.std(r_flat_samples[:,count])
					samples_truths.append(np.mean(r_flat_samples[:,count]))
					samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
					samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))
					count=count+1
				else:
					mscl[i]=scale[i][0]
					sd_scl[i]=0.0
			if scale[i][2] == 'normal': 
				if scale[i][1] != 0:
					mscl[i]=np.mean(r_flat_samples[:,count])
					sd_scl[i]=np.std(r_flat_samples[:,count])
					samples_truths.append(np.mean(r_flat_samples[:,count]))
					samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
					samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))
					count=count+1
				else:
					mscl[i]=scale[i][0]
					sd_scl[i]=0.0

		m_vm=[]
		sd_vm=[]
		for i in range(0,num_curves):
			m_vm.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
			sd_vm.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		for i in range(0,c_param):
			for j in range(0,num_curves):
				if multi_param[i][-1] == 'uniform': 
					if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
						m_vm[j][i]=np.mean(r_flat_samples[:,count])
						sd_vm[j][i]=np.std(r_flat_samples[:,count])
						samples_truths.append(np.mean(r_flat_samples[:,count]))
						samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
						samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))	

						count=count+1
					else:
						m_vm[j][i]=multi_param[i][2*j+1]
						sd_vm[j][i]=0.0

				if multi_param[i][-1] == 'normal':
					if multi_param[i][2*j+2] != 0:
						m_vm[j][i]=np.mean(r_flat_samples[:,count])
						sd_vm[j][i]=np.std(r_flat_samples[:,count])
						samples_truths.append(np.mean(r_flat_samples[:,count]))
						samples_truths_plus.append(np.mean(r_flat_samples[:,count])+np.std(r_flat_samples[:,count]))
						samples_truths_minus.append(np.mean(r_flat_samples[:,count])-np.std(r_flat_samples[:,count]))	

						count=count+1
					else:
						m_vm[j][i]=multi_param[i][2*j+1]
						sd_vm[j][i]=0.0

	if bootstraping:
		samples_Refl=np.empty([num_curves,1001, bt_iter+1])
		samples_Refl2=[]
		for i in range(0,num_curves):
			samples_Refl[i,:,0]=q_bin[i][:]
			samples_Refl2.append(np.empty([np.size(data[i],0),bt_iter+1]))

		samples_Profile=np.empty([num_curves,len(layers_min),1000, 2*bt_iter+1])
		samples_Solvent=np.empty([num_curves,len(layers_min),1000, 2*bt_iter+1])
		samples=np.empty([bt_iter,free_param])
		samples_sigma=np.empty([bt_iter+1])

		Profile_mean=[]
		Solvent_mean=[]

		for i in range(0,bt_iter):
			print('Bootstrap #',i+1)
			databoot=[]
			for curve in range(0,num_curves):
				databoot.append(np.zeros([np.size(data[curve],0),4]))

			for curve in range(0,num_curves):
				for j in range(0,np.size(data[curve],0)):
					databoot[curve][j,0]=data[curve][j,0]
					databoot[curve][j,1]=np.random.normal(loc=data[curve][j,1], scale=data[curve][j,2]) #data[j][1]+np.random.normal(scale=data[j][2])  #Refl2[j][1]+np.random.normal(scale=data[j][2])  
					databoot[curve][j,2]=data[curve][j,2]
					databoot[curve][j,3]=data[curve][j,3]

			population_size=int(3*(free_param)/(m_param+c_param+(5*sum(mnlayers))))
			if population_size == 0:
				population_size = 1	  
			results = differential_evolution(fig_of_merit_sym, boot_bounds, (databoot, resolution, mnlayers, m_param,layers_max,model_param,model_constraints,fit_mode,num_curves,c_param,multi_param,f_layer_fun,f_left_fun,f_right_fun,center_fun,experror,fit_weight,patches),tol=0.1,popsize=population_size,maxiter=100,polish=True,updating='deferred',workers=mp)


			print('/n')
			print('Success: '+str(results.success))
			print('Number of evaluation: '+str(results.nfev))
			print('Number of iterations: '+str(results.nit))
			print('FOM: '+str(results.fun))

			count=0
			for l in range(len(layers_min)):
				for j in range(0,np.size(layers_min[l],0)):
					for k in range(0,5):
						if layers_min[l][j][k] != layers_max[l][j][k]:
							if l>0:
								samples[i][count] = results.x[5*sum(mnlayers[0:l])+j*5+k]
							else:
								samples[i][count] = results.x[j*5+k]
							count=count+1

			for j in range(0,m_param):
				if model_param[j][4] == 'uniform' and model_param[j][1] != model_param[j][2]:
					samples[i][count] = results.x[(sum(mnlayers)-0)*5+j]
					count=count+1
				if model_param[j][4] == 'normal' and model_param[j][2] != 0:
					samples[i][count] = results.x[(sum(mnlayers)-0)*5+j]
					count=count+1

			for curve in range(0,num_curves):
				if background[curve][2] == 'uniform' and background[curve][0] != background[curve][1]:
					samples[i][count] = results.x[(sum(mnlayers)-0)*5+m_param+curve]
					count=count+1
				if background[curve][2] == 'normal' and background[curve][1] != 0:
					samples[i][count] = results.x[(sum(mnlayers)-0)*5+m_param+curve]
					count=count+1

			for curve in range(0,num_curves):
				if scale[curve][2] == 'uniform' and scale[curve][0] != scale[curve][1]:
					samples[i][count] = results.x[(sum(mnlayers)-0)*5+m_param+num_curves+curve]
					count=count+1
				if scale[curve][2] == 'normal' and scale[curve][1] != 0:
					samples[i][count] = results.x[(sum(mnlayers)-0)*5+m_param+num_curves+curve]
					count=count+1

			for curve in range(0,num_curves):
				for j in range(0,c_param):
					if multi_param[j][-1] == 'uniform' and multi_param[j][2*curve+1] != multi_param[j][2*curve+2]:
						samples[i][count] = results.x[sum(mnlayers)*5+m_param+2*num_curves+c_param*curve+j]
						count=count+1
					if multi_param[j][-1] == 'normal' and multi_param[j][2*curve+2] != 0:
						samples[i][count] = results.x[sum(mnlayers)*5+m_param+2*num_curves+c_param*curve+j]
						count=count+1

			for curve in range(0,num_curves):

				vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
				for j in range(0,m_param):
					vp[j]=results.x[sum(mnlayers)*5+j]

				vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
				for j in range(0,c_param):
					vm[j]=results.x[sum(mnlayers)*5+m_param+2*num_curves+c_param*curve+j]


				layers=[]
				for l in range(len(layers_min)):
					sub_layers=[]
					for j in range(0,np.size(layers_min[l],0)):
						line=[]
						for k in range(0,5):
							if isinstance(layers_max[l][j][k], str):
								line.append(np.float((f_layer_fun(l,j,k,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],j,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))
							else:
								if l>0:
									line.append(results.x[5*sum(mnlayers[0:l])+j*5+k])
								else:
									line.append(results.x[j*5+k])
						line.append(layers_max[l][j][5])
						sub_layers.append(line)
					layers.append(sub_layers)


				bkg=results.x[sum(mnlayers)*5+m_param+curve]
				scl=results.x[sum(mnlayers)*5+m_param+num_curves+curve]


				cRefl=Reflectivity(q_bin[curve],res_bin[curve], layers, resolution[curve], bkg,scl,patches)
				samples_Refl[curve,:,i+1]=cRefl[:,1]
				samples_Refl2[curve][:,0]=data[curve][:,0]
				samples_Refl2[curve][:,i+1]=Reflectivity(data[curve][:,0],data[curve][:,3], layers, resolution[curve], bkg, scl, patches)[:,1]
				cProfile=profile(layers, 1000)
				for j in range(len(layers_min)):
					samples_Profile[curve,j,:,2*i]=cProfile[j][:,0]
					samples_Profile[curve,j,:,2*i+1]=cProfile[j][:,1]
				cSolvent=solvent_penetration(layers, 1000)
				for j in range(len(layers_min)):
					samples_Solvent[curve,j,:,2*i]=cSolvent[j][:,0]
					samples_Solvent[curve,j,:,2*i+1]=cSolvent[j][:,1]



		print('Plotting curves... please wait...')

		for curve in range(0,num_curves):

			Refl[curve][:,1]=0.0
			Refl2[curve][:,1]=0.0
			for k in range(0,bt_iter):
				Refl[curve][:,1]=Refl[curve][:,1]+samples_Refl[curve,:,k+1]/np.float(bt_iter)
				Refl2[curve][:,1]=Refl2[curve][:,1]+samples_Refl2[curve][:,k+1]/np.float(bt_iter)


			chi_s[curve]=0.0
			Nexp=np.size(data[curve],0)
			for i in range(0, Nexp):
				chi_s[curve]=chi_s[curve]+(1.0/float(Nexp))*((data[curve][i][1]-Refl2[curve][i,1])/(data[curve][i][2]))**2

			#Calculate mean profiles
			z_max=0
			z_min=0
			for k in range(len(layers_min)):
				for i in range(0,bt_iter):
					for j in range(0,1000):
						if samples_Profile[curve,k,j,2*i] > z_max: z_max=samples_Profile[curve,k,j,2*i]
						if samples_Profile[curve,k,j,2*i] < z_min: z_min=samples_Profile[curve,k,j,2*i]		

			z_bin = np.linspace(z_min, z_max, 1000)
			Profile_mean.append(np.zeros([len(layers_min),len(z_bin), 3]))
			Solvent_mean.append(np.zeros([len(layers_min),len(z_bin), 3]))
			for k in range(len(layers_min)):
				Profile_mean[curve][k][:,0]=z_bin[:]
				Solvent_mean[curve][k][:,0]=z_bin[:]


			for i in range(0,bt_iter):
				for k in range(len(layers_min)):
					for j in range(0,len(z_bin)):
						Profile_mean[curve][k][j,1]=Profile_mean[curve][k][j,1]+np.interp(Profile_mean[curve][k][j,0], samples_Profile[curve,k,:,2*i],samples_Profile[curve,k,:,2*i+1])/float(bt_iter)
						Solvent_mean[curve][k][j,1]=Solvent_mean[curve][k][j,1]+np.interp(Solvent_mean[curve][k][j,0], samples_Solvent[curve,k,:,2*i],samples_Solvent[curve,k,:,2*i+1])/float(bt_iter)

			for i in range(0,bt_iter):
				for k in range(len(layers_min)):
					for j in range(0,len(z_bin)):
						Profile_mean[curve][k][j,2]=Profile_mean[curve][k][j,2]+((Profile_mean[curve][k][j,1]-np.interp(Profile_mean[curve][k][j,0], samples_Profile[curve,k,:,2*i],samples_Profile[curve,k,:,2*i+1]))**2)
						Solvent_mean[curve][k][j,2]=Solvent_mean[curve][k][j,2]+((Solvent_mean[curve][k][j,1]-np.interp(Solvent_mean[curve][k][j,0], samples_Solvent[curve,k,:,2*i],samples_Solvent[curve,k,:,2*i+1]))**2)

			for j in range(0,len(z_bin)):
				for k in range(len(layers_min)):
					Profile_mean[curve][k][j,2]=np.sqrt(Profile_mean[curve][k][j,2]/float(bt_iter-1))
					Solvent_mean[curve][k][j,2]=np.sqrt(Solvent_mean[curve][k][j,2]/float(bt_iter-1))

		samples_truths=[]
		samples_truths_plus=[]
		samples_truths_minus=[]


		sd_vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		count=0
		for k in range(len(layers_min)):
			for i in range(0,np.size(layers_min[k],0)):
				for j in range(0,5):
					if layers_min[k][i][j] != layers_max[k][i][j]:
						layers[k][i][j]=np.mean(samples[:,count])
						samples_truths.append(np.mean(samples[:,count]))
						samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
						samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))

						count=count+1

		for j in range(0,m_param):
			if model_param[j][4] == 'uniform' and model_param[j][1] != model_param[j][2]:
				vp[j]=np.mean(samples[:,count])
				sd_vp[j]=np.std(samples[:,count])
				samples_truths.append(np.mean(samples[:,count]))
				samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
				samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))

				count=count+1
			if model_param[j][4] == 'normal' and model_param[j][2] != 0:
				vp[j]=np.mean(samples[:,count])
				sd_vp[j]=np.std(samples[:,count])
				samples_truths.append(np.mean(samples[:,count]))
				samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
				samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))

				count=count+1

		mbkg=np.zeros(num_curves)
		sd_bkg=np.zeros(num_curves)

		for i in range(0,num_curves):
			if background[i][2] == 'uniform':
				if background[i][0] != background[i][1]:
					mbkg[i]=np.mean(samples[:,count])
					sd_bkg[i]=np.std(samples[:,count])
					samples_truths.append(np.mean(samples[:,count]))
					samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
					samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))
					count=count+1
				else:
					mbkg[i]=background[i][0]
					sd_bkg[i]=0.0
			if background[i][2] == 'normal':
				if background[i][1] != 0:
					mbkg[i]=np.mean(samples[:,count])
					sd_bkg[i]=np.std(samples[:,count])
					samples_truths.append(np.mean(samples[:,count]))
					samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
					samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))
					count=count+1
				else:
					mbkg[i]=background[i][0]
					sd_bkg[i]=0.0

		mscl=np.zeros(num_curves)
		sd_scl=np.zeros(num_curves)

		for i in range(0,num_curves):
			if scale[i][2] == 'uniform':
				if scale[i][0] != scale[i][1]:
					mscl[i]=np.mean(samples[:,count])
					sd_scl[i]=np.std(samples[:,count])
					samples_truths.append(np.mean(samples[:,count]))
					samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
					samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))
					count=count+1
				else:
					mscl[i]=scale[i][0]
					sd_scl[i]=0.0
			if scale[i][2] == 'normal':
				if scale[i][1] != 0:
					mscl[i]=np.mean(samples[:,count])
					sd_scl[i]=np.std(samples[:,count])
					samples_truths.append(np.mean(samples[:,count]))
					samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
					samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))
					count=count+1
				else:
					mscl[i]=scale[i][0]
					sd_scl[i]=0.0

		m_vm=[]
		sd_vm=[]
		for i in range(0,num_curves):
			m_vm.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
			sd_vm.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		for i in range(0,c_param):
			for j in range(0,num_curves):
				if multi_param[i][-1] == 'uniform':
					if multi_param[i][2*j+1] != multi_param[i][2*j+2]:
						m_vm[j][i]=np.mean(samples[:,count])
						sd_vm[j][i]=np.std(samples[:,count])
						samples_truths.append(np.mean(samples[:,count]))
						samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
						samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))	

						count=count+1
					else:
						m_vm[j][i]=multi_param[i][2*j+1]
						sd_vm[j][i]=0.0
				if multi_param[i][-1] == 'normal':
					if multi_param[i][2*j+2] != 0:
						m_vm[j][i]=np.mean(samples[:,count])
						sd_vm[j][i]=np.std(samples[:,count])
						samples_truths.append(np.mean(samples[:,count]))
						samples_truths_plus.append(np.mean(samples[:,count])+np.std(samples[:,count]))
						samples_truths_minus.append(np.mean(samples[:,count])-np.std(samples[:,count]))	

						count=count+1
					else:
						m_vm[j][i]=multi_param[i][2*j+1]
						sd_vm[j][i]=0.0

	elapsed_time = time.time() - start_time


	# Export profiles
	if project != 'none':
		for curve in range(0,num_curves):
			np.savetxt(project+"_reflectivity_curve"+str(curve+1)+".dat", Refl2[curve], fmt='%1.4e', header="Q (A^-1),  R, RQ^4 (A^-4)")
			if mcmc or bootstraping:
				for i in range(len(layers_min)):
					if len(layers_min) == 1:
						np.savetxt(project+"_sld_profile_curve"+str(curve+1)+".dat", Profile_mean[curve][i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2), sld standard deviation")
					else:
						np.savetxt(project+"_sld_profile_curve"+str(curve+1)+"_patch#"+str(i+1)+".dat", Profile_mean[curve][i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2), sld standard deviation")

			else:
				for i in range(len(layers_min)):
					if len(layers_min) == 1:
						np.savetxt(project+"_sld_profile_curve"+str(curve+1)+".dat", Profile[curve][i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2)")
					else:
						np.savetxt(project+"_sld_profile_curve"+str(curve+1)+"_patch#"+str(i+1)+".dat", Profile[curve][i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2)")
		if mcmc or bootstraping:
			for i in range(len(layers_min)):
				if len(layers_min) == 1:
					np.savetxt(project+"_solvent_profile.dat", Solvent_mean[curve][i], fmt='%1.4e',header="z (A),  solvent volume fraction, solvent volume fraction standard deviation")
				else:
					np.savetxt(project+"_solvent_profile"+"_patch#"+str(i+1)+".dat", Solvent_mean[curve][i], fmt='%1.4e',header="z (A),  solvent volume fraction, solvent volume fraction standard deviation")
		else:
			for i in range(len(layers_min)):
				if len(layers_min) == 1:
					np.savetxt(project+"_solvent_profile.dat", Solvent[curve][i], fmt='%1.4e', header="z (A),  solvent volume fraction")
				else:
					np.savetxt(project+"_solvent_profile"+"_patch#"+str(i+1)+".dat", Solvent[curve][i], fmt='%1.4e', header="z (A),  solvent volume fraction")

		#print(Refl)
		pcolors=['black','blue','green','cyan','magenta','yellow', 'red']
		scolors=['.k','.b','.g','.c','.m','.y','.r']
		plinestyle=['-','--',':','-.','.','--',':']
		#plt.figure(figsize=(9, 7))
		plt.figure(figsize=(7, 5.5))
		grid = plt.GridSpec(6, 2)
		ax1=plt.subplot(grid[0:3, 0]) 
		plt.xlim([0.0,data[0][-1,0]])
		if bootstraping:
			for curve in range(0,num_curves):
				for i in range(0,bt_iter):
					plt.plot(samples_Refl[curve,:,0],samples_Refl[curve,:,i+1]/(10.0**(2*curve)),color=pcolors[curve],alpha=0.05)
		if mcmc:
			for curve in range(0,num_curves):
				for i in range(0,np.size(r_flat_samples,0)):
					plt.plot(samples_Refl[curve,:,0],samples_Refl[curve,:,i+1]/(10.0**(2*curve)),color=pcolors[curve],alpha=0.05)
		for curve in range(0,num_curves):
			plt.errorbar(data[curve][:,0], data[curve][:,1]/(10.0**(2*curve)), data[curve][:,2]/(10.0**(2*curve)),fmt=scolors[curve],ecolor=pcolors[curve],zorder=0,alpha=0.7)
			plt.plot(Refl[curve][:,0],Refl[curve][:,1]/(10.0**(2*curve)),color='red')
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)$')
		#plt.title(r'$R(Q)$')
		plt.grid(True)

		#plt.subplot(1,2,2)
		plt.subplot(grid[0:, 1])
		plt.xlim([0.0,data[0][-1,0]])
		if bootstraping:
			for curve in range(0,num_curves):
				for i in range(0,bt_iter):
					plt.plot(samples_Refl[curve,:,0],(samples_Refl[curve,:,i+1]/(10.0**(2*curve)))*samples_Refl[curve,:,0]**4,color=pcolors[curve],alpha=0.05)
		if mcmc:
			for curve in range(0,num_curves):
				for i in range(0,np.size(r_flat_samples,0)):
					plt.plot(samples_Refl[curve,:,0],(samples_Refl[curve,:,i+1]/(10.0**(2*curve)))*samples_Refl[curve,:,0]**4,color=pcolors[curve],alpha=0.05)
		for curve in range(0,num_curves):
			plt.errorbar(data[curve][:,0], (data[curve][:,1]/(10.0**(2*curve)))*data[curve][:,0]**4, (data[curve][:,2]/(10.0**(2*curve)))*data[curve][:,0]**4,fmt=scolors[curve],ecolor=pcolors[curve],zorder=0,label=r'$\chi_{'+str(curve+1)+'}^2 = $'+'{:06.3e}'.format(chi_s[curve]),alpha=0.7)
			plt.plot(Refl[curve][:,0],(Refl[curve][:,1]/(10.0**(2*curve)))*Refl[curve][:,0]**4,color='red')
		plt.legend(framealpha=0.4)
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)Q^4$')
		plt.grid(True)


		ax3=plt.subplot(grid[3, 0]) 
		plt.xlim([0.0,data[0][-1,0]])
		plt.axhline(y=0, lw=1,ls='dashed', color='k')
		for curve in range(0,num_curves):
			if experror == True:
				plt.scatter(data[curve][:,0],(-Refl2[curve][:,1]+data[curve][:,1])/(data[curve][:,2]),2,color=pcolors[curve])
			else:
				plt.scatter(data[curve][:,0],(-Refl2[curve][:,1]+data[curve][:,1])/(Refl2[curve][:,1]),2,color=pcolors[curve])
		#plt.yscale('log')
		if experror == True:
			plt.ylabel(r'$\Delta / \sigma$')
		else:
			plt.ylabel(r'$\Delta / R$')
		plt.xlabel(r'$Q(\AA^{-1})$')

		#plt.subplot(3,2,5)
		plt.subplot(grid[4:, 0])
		if bootstraping:
			for i in range(len(layers_min)):
				for curve in range(0,num_curves):
					plt.plot(Profile_mean[curve][i][:,0],Profile_mean[curve][i][:,1]*1e6,color=pcolors[curve],linestyle=plinestyle[i])
					plt.fill_between(Profile_mean[curve][i][:,0],(Profile_mean[curve][i][:,1]-Profile_mean[curve][i][:,2])*1e6,(Profile_mean[curve][i][:,1]+Profile_mean[curve][i][:,2])*1e6,color=pcolors[curve],alpha=0.3)
		if mcmc:
			for i in range(len(layers_min)):
				for curve in range(0,num_curves):
					plt.plot(Profile_mean[curve][i][:,0],Profile_mean[curve][i][:,1]*1e6,color=pcolors[curve],linestyle=plinestyle[i])
					plt.fill_between(Profile_mean[curve][i][:,0],(Profile_mean[curve][i][:,1]-Profile_mean[curve][i][:,2])*1e6,(Profile_mean[curve][i][:,1]+Profile_mean[curve][i][:,2])*1e6,color=pcolors[curve],alpha=0.3)
		if not bootstraping and not mcmc:
			for curve in range(0,num_curves):
				for i in range(len(layers_min)):
					plt.plot(Profile[curve][i][:,0],Profile[curve][i][:,1]*1e6,color=pcolors[curve],linestyle=plinestyle[i])

		plt.xlabel(r'$z(\AA)$')
		plt.ylabel(r'$sld(10^{-6}\AA^{-2})$',color='black')

		no_solvent=True
		for i in range(len(layers_min)):
			if system[i][-1][4] != 0.0 or system[i][0][4] != 0.0:
				no_solvent=False

		if no_solvent == False:
			plt.twinx()
			if bootstraping:
				for i in range(len(layers_min)):
					for curve in range(0,num_curves):
						plt.plot(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1],color='orange',linestyle=plinestyle[i])
						plt.fill_between(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1]-Solvent_mean[curve][i][:,2],Solvent_mean[curve][i][:,1]+Solvent_mean[curve][i][:,2],color='orange',alpha=0.3)
			if mcmc:
				for i in range(len(layers_min)):
					for curve in range(0,num_curves):
						plt.plot(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1],color='orange',linestyle=plinestyle[i])
						plt.fill_between(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1]-Solvent_mean[curve][i][:,2],Solvent_mean[curve][i][:,1]+Solvent_mean[curve][i][:,2],color='orange',alpha=0.3)
			if not bootstraping and not mcmc:
				for curve in range(0,num_curves):
					for i in range(len(layers_min)):
						plt.plot(Solvent[curve][i][:,0],Solvent[curve][i][:,1],color='orange',linestyle=plinestyle[i])
			plt.ylabel('solvent',color='orange')

		ax1.get_shared_x_axes().join(ax1, ax3)
		plt.tight_layout()
		plt.draw() 
		plt.savefig(project+'_fit_summary.pdf')
		if plot is True: plt.show()
		plt.close()  

		#individual plots
		# Reflectivity
		plt.figure()
		grid = plt.GridSpec(4, 1)
		ax1=plt.subplot(grid[0:3, 0])
		plt.xlim([0.0,data[0][-1,0]])
		if bootstraping:
			for curve in range(0,num_curves):
				for i in range(0,bt_iter):
					plt.plot(samples_Refl[curve,:,0],samples_Refl[curve,:,i+1]/(10.0**(2*curve)),color=pcolors[curve],alpha=0.05)
		if mcmc:
			for curve in range(0,num_curves):
				for i in range(0,np.size(r_flat_samples,0)):
					plt.plot(samples_Refl[curve,:,0],samples_Refl[curve,:,i+1]/(10.0**(2*curve)),color=pcolors[curve],alpha=0.05)
		for curve in range(0,num_curves):
			plt.errorbar(data[curve][:,0], data[curve][:,1]/(10.0**(2*curve)), data[curve][:,2]/(10.0**(2*curve)),fmt=scolors[curve],ecolor=pcolors[curve],zorder=0,label=r'$\chi_{'+str(curve+1)+'}^2 = $'+'{:06.3e}'.format(chi_s[curve]),alpha=0.7)
			plt.plot(Refl[curve][:,0],Refl[curve][:,1]/(10.0**(2*curve)),color='red')
		plt.legend(framealpha=0.4)
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)$')
		plt.grid(True)

		ax2=plt.subplot(grid[3, 0])
		plt.xlim([0.0,data[0][-1,0]])
		plt.axhline(y=0, lw=1,ls='dashed', color='k')
		for curve in range(0,num_curves):
			if experror == True:
				plt.scatter(data[curve][:,0],(-Refl2[curve][:,1]+data[curve][:,1])/(data[curve][:,2]),2,color=pcolors[curve])
			else:
				plt.scatter(data[curve][:,0],(-Refl2[curve][:,1]+data[curve][:,1])/(Refl2[curve][:,1]),2,color=pcolors[curve])

		if experror == True:
			plt.ylabel(r'$\Delta / \sigma$')
		else:
			plt.ylabel(r'$\Delta / R$')
		ax1.get_shared_x_axes().join(ax1, ax2)
		plt.tight_layout()
		plt.draw() 
		plt.savefig(project+'_fit_reflectivity.pdf')

		#Reflecticity x Q^4
		plt.figure()
		plt.xlim([0.0,data[0][-1,0]])
		if bootstraping:
			for curve in range(0,num_curves):
				for i in range(0,bt_iter):
					plt.plot(samples_Refl[curve,:,0],(samples_Refl[curve,:,i+1]/(10.0**(2*curve)))*samples_Refl[curve,:,0]**4,color=pcolors[curve],alpha=0.05)
		if mcmc:
			for curve in range(0,num_curves):
				for i in range(0,np.size(r_flat_samples,0)):
					plt.plot(samples_Refl[curve,:,0],(samples_Refl[curve,:,i+1]/(10.0**(2*curve)))*samples_Refl[curve,:,0]**4,color=pcolors[curve],alpha=0.05)
		for curve in range(0,num_curves):
			plt.errorbar(data[curve][:,0], (data[curve][:,1]/(10.0**(2*curve)))*data[curve][:,0]**4, (data[curve][:,2]/(10.0**(2*curve)))*data[curve][:,0]**4,fmt=scolors[curve],ecolor=pcolors[curve],zorder=0,label=r'$\chi_{'+str(curve+1)+'}^2 = $'+'{:06.3e}'.format(chi_s[curve]),alpha=0.7)
			plt.plot(Refl[curve][:,0],(Refl[curve][:,1]/(10.0**(2*curve)))*Refl[curve][:,0]**4,color='red')
		plt.legend(framealpha=0.4)
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)Q^4$')
		plt.grid(True)
		plt.draw() 
		plt.savefig(project+'_fit_RQ^4.pdf')

		#sld profiles
		plt.figure()
		if bootstraping:
			for i in range(len(layers_min)):
				for curve in range(0,num_curves):
					plt.plot(Profile_mean[curve][i][:,0],Profile_mean[curve][i][:,1]*1e6,color=pcolors[curve],linestyle=plinestyle[i])
					plt.fill_between(Profile_mean[curve][i][:,0],(Profile_mean[curve][i][:,1]-Profile_mean[curve][i][:,2])*1e6,(Profile_mean[curve][i][:,1]+Profile_mean[curve][i][:,2])*1e6,color=pcolors[curve],alpha=0.3)
		if mcmc:
			for i in range(len(layers_min)):
				for curve in range(0,num_curves):
					plt.plot(Profile_mean[curve][i][:,0],Profile_mean[curve][i][:,1]*1e6,color=pcolors[curve],linestyle=plinestyle[i])
					plt.fill_between(Profile_mean[curve][i][:,0],(Profile_mean[curve][i][:,1]-Profile_mean[curve][i][:,2])*1e6,(Profile_mean[curve][i][:,1]+Profile_mean[curve][i][:,2])*1e6,color=pcolors[curve],alpha=0.3)
		if not bootstraping and not mcmc:
			for curve in range(0,num_curves):
				for i in range(len(layers_min)):
					plt.plot(Profile[curve][i][:,0],Profile[curve][i][:,1]*1e6,color=pcolors[curve],linestyle=plinestyle[i])

		plt.xlabel(r'$z(\AA)$')
		plt.ylabel(r'$sld(10^{-6}\AA^{-2})$',color='black')
		plt.draw() 
		plt.savefig(project+'_fit_sld_profile.pdf')

		#solvent profile
		plt.figure()
		if bootstraping:
			for i in range(len(layers_min)):
				for curve in range(0,num_curves):
					plt.plot(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1],color='black',linestyle=plinestyle[i])
					plt.fill_between(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1]-Solvent_mean[curve][i][:,2],Solvent_mean[curve][i][:,1]+Solvent_mean[curve][i][:,2],color='black',alpha=0.3)
		if mcmc:
			for i in range(len(layers_min)):
				for curve in range(0,num_curves):
					plt.plot(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1],color='black',linestyle=plinestyle[i])
					plt.fill_between(Solvent_mean[curve][i][:,0],Solvent_mean[curve][i][:,1]-Solvent_mean[curve][i][:,2],Solvent_mean[curve][i][:,1]+Solvent_mean[curve][i][:,2],color='black',alpha=0.3)
		if not bootstraping and not mcmc:
			for curve in range(0,num_curves):
				for i in range(len(layers_min)):
					plt.plot(Solvent[curve][i][:,0],Solvent[curve][i][:,1],color='black',linestyle=plinestyle[i])
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel('solvent',color='black')
		plt.draw() 
		plt.savefig(project+'_fit_solvent_profile.pdf')

		#corner plot
		if bootstraping:
			#corner_figure = corner.corner(samples,labels=corner_labels,quantiles=[0.16, 0.5, 0.84],show_titles=True,truths=samples_truths, title_kwargs={"fontsize": 12})
			#plt.figure(figsize=(7, 7))
			corner_figure = corner.corner(samples,labels=corner_labels,show_titles=False,truths=samples_truths, title_kwargs={"fontsize": 18},label_kwargs={"fontsize": 18})
			corner_figure.savefig(project+'_corner_plot.pdf')
			#corner_figure.show()

		if mcmc:
			#plt.figure(figsize=(7, 7))
			plt.tight_layout()
			corner_figure = corner.corner(flat_samples,labels=corner_labels,show_titles=False,truths=samples_truths,title_kwargs={"fontsize": 18},label_kwargs={"fontsize": 18})
			corner_figure.savefig(project+'_corner_plot.pdf')
			#corner_figure.show()	

		f = open(project+"_final_parameters.log", "w")
		f.write('--------------------------------------------------------------------\n')
		f.write('Program ANAKLASIS - Fit Module for X-ray/Neutron reflection datasets\n')
		f.write('version 1.2, February 2021\n')
		f.write('developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ\n')
		f.write('for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de\n')
		f.write('--------------------------------------------------------------------\n')


		f.write('Project name: '+project+'\n')
		f.write(print_buffer)

		for i in range(0,num_curves):
			f.write('input file #'+str(i+1)+': '+in_file[i]+'\n')


		f.write('free parameters = '+str(free_param))
		f.write('\n')
		f.write('list of free parameters: '+list_free_param+'\n')

		f.write('--------------------------------------------------------------------\n')
		f.write('Final model parameters\n')
		f.write('--------------------------------------------------------------------\n')
		for curve in range(0,num_curves):
			ii,jj,kk=sympy.symbols('ii jj kk', integer=True)
			x,y=sympy.symbols('x y')
			p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n=sympy.symbols('p0 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 n')
			m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39=sympy.symbols('m0 m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 m20 m21 m22 m23 m24 m25 m26 m27 m28 m29 m30 m31 m32 m33 m34 m35 m36 m37 m38 m39')

			vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
			for i in range(0,m_param):
				vp[i]=results.x[sum(mnlayers)*5+i]

			vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
			for i in range(0,c_param):
				#vm[i]=results.x[np.size(layers_max,0)*5+m_param+2*num_curves+c_param*curve+i]
				vm[i]=results.x[sum(mnlayers)*5+m_param+2*num_curves+num_curves*i+curve]


			layers=[]
			for k in range(len(layers_min)):
				sub_layers=[]
				for i in range(0,np.size(layers_min[k],0)):
					line=[]
					for j in range(0,5):
						if isinstance(layers_max[k][i][j], str):
							line.append(np.float((f_layer_fun(k,i,j,vp[0],vp[1],vp[2],vp[3],vp[4],vp[5],vp[6],vp[7],vp[8],vp[9],vp[10],vp[11],vp[12],vp[13],vp[14],vp[15],vp[16],vp[17],vp[18],vp[19],vp[20],vp[21],vp[22],vp[23],vp[24],vp[25],vp[26],vp[27],vp[28],vp[29],vp[30],vp[31],vp[32],vp[33],vp[34],vp[35],vp[36],vp[37],vp[38],vp[39],i,vm[0],vm[1],vm[2],vm[3],vm[4],vm[5],vm[6],vm[7],vm[8],vm[9],vm[10],vm[11],vm[12],vm[13],vm[14],vm[15],vm[16],vm[17],vm[18],vm[19],vm[20],vm[21],vm[22],vm[23],vm[24],vm[25],vm[26],vm[27],vm[28],vm[29],vm[30],vm[31],vm[32],vm[33],vm[34],vm[35],vm[36],vm[37],vm[38],vm[39]))))
						else:
							if k > 0:
								line.append(results.x[5*sum(mnlayers[0:k])+i*5+j])
							else:
								line.append(results.x[i*5+j])
					line.append(layers_max[k][i][5])
					sub_layers.append(line)
				layers.append(sub_layers)

			bkg=results.x[sum(mnlayers)*5+m_param+curve]
			scl=results.x[sum(mnlayers)*5+m_param+num_curves+curve]

			f.write('\n')
			f.write('Curve #'+str(curve+1))
			f.write('\n')
			if resolution[curve] !=-1:
				f.write('Instrumental Resolution, dQ/Q: '+str(resolution[curve])+'\n')
			else:
				f.write('Instrumental Resolution, dQ/Q pointwise\n')
			f.write('fit weigth: '+str(fit_weight[curve])+'\n')
			f.write('--------------------------------------------------------------------\n')
			count=0
			for k in range(len(layers_min)):
				if len(layers_min) > 1:
					f.write('Patch #'+str(k+1)+', coverage: '+str(patches[k])+'\n')
					f.write('\n')
				for i in range(0,np.size(layers_min[k],0)):
					for j in range(0,5):
						if layers_min[k][i][j] == layers_max[k][i][j] and not isinstance(layers_min[k][i][j], str):
							if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (fixed)\n')
							if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (fixed)\n')
							if j == 2 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (fixed)\n')
							if j == 3 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (fixed)\n')
							if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (fixed)\n')
						elif layers_min[k][i][j] == layers_max[k][i][j] and isinstance(layers_min[k][i][j], str):
							if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (parametric), -> '+layers_min[k][i][0]+'\n')
							if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (parametric), -> '+layers_min[k][i][1]+'\n')
							if j == 2 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (parametric), -> '+layers_min[k][i][2]+'\n')
							if j == 3 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (parametric), -> '+layers_min[k][i][3]+'\n')
							if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (parametric), -> '+layers_min[k][i][4]+'\n')
						else:
							if bootstraping or mcmc:
								if j == 0: 
									f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {samples[:,count].mean()/1e-6:.2f} \u00B1 {np.std(samples[:,count])/1e-6:.2f}, bounds: {float(layers_min[k][i][0]/1e-6):.2f} -> {float(layers_max[k][i][0]/1e-6):.2f}\n')
									count=count+1
								if j == 1: 
									f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {samples[:,count].mean()/1e-6:.2f} \u00B1 {np.std(samples[:,count])/1e-6:.2f}, bounds: {float(layers_min[k][i][1]/1e-6):.2f} -> {float(layers_max[k][i][1]/1e-6):.2f}\n')
									count=count+1
								if j == 2 and i+1 != np.size(layers_min[k],0): 
									f.write(layers[k][i][5]+f' thickness (A)              = {samples[:,count].mean():.2f} \u00B1 {np.std(samples[:,count]):.2f}, bounds: {float(layers_min[k][i][2]):.2f} -> {float(layers_max[k][i][2]):.2f}\n')
									count=count+1
								if j == 3 and i+1 != np.size(layers_min[k],0): 
									f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {samples[:,count].mean():.2f} \u00B1 {np.std(samples[:,count]):.2f}, bounds: {float(layers_min[k][i][3]):.2f} -> {float(layers_max[k][i][3]):.2f}\n')
									count=count+1
								if j == 4: 
									f.write(layers[k][i][5]+f' solvent volume fraction    = {samples[:,count].mean():.2f} \u00B1 {np.std(samples[:,count]):.2f}, bounds: {float(layers_min[k][i][4]):.2f} -> {float(layers_max[k][i][4]):.2f}\n')
									count=count+1
							else:
								if experror == True:
									if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} \u00B1 {float(df_layers[k][i][0]/1e-6):.2f}, bounds: {(layers_min[k][i][0]/1e-6):.2f} -> {(layers_max[k][i][0]/1e-6):.2f}\n')
									if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} \u00B1 {float(df_layers[k][i][1]/1e-6):.2f}, bounds: {(layers_min[k][i][1]/1e-6):.2f} -> {(layers_max[k][i][1]/1e-6):.2f}\n')
									if j == 2 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} \u00B1 {float(df_layers[k][i][2]):.2f}, bounds: {float(layers_min[k][i][2]):.2f} -> {float(layers_max[k][i][2]):.2f}\n')
									if j == 3 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} \u00B1 {float(df_layers[k][i][3]):.2f}, bounds: {float(layers_min[k][i][3]):.2f} -> {float(layers_max[k][i][3]):.2f}\n')
									if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} \u00B1 {float(df_layers[k][i][4]):.2f}, bounds: {float(layers_min[k][i][4]):.2f} -> {float(layers_max[k][i][4]):.2f}\n')
								else:
									if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f}, bounds: {(layers_min[k][i][0]/1e-6):.2f} -> {(layers_max[k][i][0]/1e-6):.2f}\n')
									if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f}, bounds: {(layers_min[k][i][1]/1e-6):.2f} -> {(layers_max[k][i][1]/1e-6):.2f}\n')
									if j == 2 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f}, bounds: {float(layers_min[k][i][2]):.2f} -> {float(layers_max[k][i][2]):.2f}\n')
									if j == 3 and i+1 != np.size(layers_min[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f}, bounds: {float(layers_min[k][i][3]):.2f} -> {float(layers_max[k][i][3]):.2f}\n')
									if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f}, bounds: {float(layers_min[k][i][4]):.2f} -> {float(layers_max[k][i][4]):.2f}\n')	
					f.write('--------------------------------------------------------------------\n')
			f.write(' \n')

			if bootstraping or mcmc:
				if background[curve][2] == 'uniform': 
					f.write('Background: '+f' = {mbkg[curve]:.2e} \u00B1 {float(sd_bkg[curve]):.2e}, bounds: '+ str(float(background[curve][0])) +' -> '+str(float(background[curve][1]))+' (uniform)'+'\n')
				if background[curve][2] == 'normal': 
					f.write('Background: '+f' = {mbkg[curve]:.2e} \u00B1 {float(sd_bkg[curve]):.2e}, bounds: '+ str(float(background[curve][0])) +' \u00B1 '+str(float(background[curve][1]))+' (normal)'+'\n')
			else:
				if experror == True:
					if background[curve][2] == 'uniform':
						f.write('Background: '+f' = {bkg:.2e} \u00B1 {float(df_bkg[curve]):.2e}, bounds: '+ str(float(background[curve][0])) +' -> '+str(float(background[curve][1]))+' (uniform)'+'\n')
					if background[curve][2] == 'normal':
						f.write('Background: '+f' = {bkg:.2e} \u00B1 {float(df_bkg[curve]):.2e}, bounds: '+ str(float(background[curve][0])) +' \u00B1 '+str(float(background[curve][1]))+' (normal)'+'\n')
				else:
					if background[curve][2] == 'uniform':
						f.write('Background: '+f' = {bkg:.2e}, bounds: '+ str(float(background[curve][0])) +' -> '+str(float(background[curve][1]))+' (uniform)'+'\n')
					if background[curve][2] == 'normal':
						f.write('Background: '+f' = {bkg:.2e}, bounds: '+ str(float(background[curve][0])) +' \u00B1 '+str(float(background[curve][1]))+' (normal)'+'\n')

			if bootstraping or mcmc:
				if scale[curve][2] == 'uniform':
					f.write('Scale: '+f' = {mscl[curve]:.2e} \u00B1 {float(sd_scl[curve]):.2e}, bounds: '+ str(float(scale[curve][0])) +' -> '+str(float(scale[curve][1]))+' (uniform)'+'\n')
				if scale[curve][2] == 'normal':
					f.write('Scale: '+f' = {mscl[curve]:.2e} \u00B1 {float(sd_scl[curve]):.2e}, bounds: '+ str(float(scale[curve][0])) +' \u00B1 '+str(float(scale[curve][1]))+' (normal)'+'\n')
			else:
				if experror == True:
					if scale[curve][2] == 'uniform':
						f.write('Scale: '+f' = {scl:.2e} \u00B1 {float(df_scl[curve]):.2e}, bounds: '+ str(float(scale[curve][0])) +' -> '+str(float(scale[curve][1]))+' (uniform)'+'\n')
					if scale[curve][2] == 'normal':
						f.write('Scale: '+f' = {scl:.2e} \u00B1 {float(df_scl[curve]):.2e}, bounds: '+ str(float(scale[curve][0])) +' \u00B1 '+str(float(scale[curve][1]))+' (normal)'+'\n')
				else:
					if scale[curve][2] == 'uniform':
						f.write('Scale: '+f' = {scl:.2e}, bounds: '+ str(float(scale[curve][0])) +' -> '+str(float(scale[curve][1]))+' (uniform)'+'\n')
					if scale[curve][2] == 'normal':
						f.write('Scale: '+f' = {scl:.2e}, bounds: '+ str(float(scale[curve][0])) +' \u00B1 '+str(float(scale[curve][1]))+' (normal)'+'\n')

			#f.write('Multi Parameters:\n')
			for i in range(0,c_param):
				if bootstraping or mcmc:
					if multi_param[i][-1] == 'uniform':
						f.write(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {m_vm[curve][i]:.2e} \u00B1 {float(sd_vm[curve][i]):.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' -> '+str(float(multi_param[i][2*curve+2]))+' (uniform)'+'\n')
					if multi_param[i][-1] == 'normal':
						f.write(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {m_vm[curve][i]:.2e} \u00B1 {float(sd_vm[curve][i]):.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' \u00B1 '+str(float(multi_param[i][2*curve+2]))+' (normal)'+'\n')
				else:
					if experror == True:
						if multi_param[i][-1] == 'uniform':
							f.write(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e} \u00B1 {float(df_vm[i*num_curves+curve]):.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' -> '+str(float(multi_param[i][2*curve+2]))+' (uniform)'+'\n')
						if multi_param[i][-1] == 'normal':
							f.write(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e} \u00B1 {float(df_vm[i*num_curves+curve]):.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' \u00B1 '+str(float(multi_param[i][2*curve+2]))+' (normal)'+'\n')
					else:
						if multi_param[i][-1] == 'uniform':
							f.write(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' -> '+str(float(multi_param[i][2*curve+2]))+' (uniform)'+'\n')
						if multi_param[i][-1] == 'normal':
							f.write(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' \u00B1 '+str(float(multi_param[i][2*curve+2]))+' (normal)'+'\n')
		f.write('\n')
		f.write('Parameters:\n')
		for i in range(0,m_param):
			if bootstraping or mcmc:
				if model_param[i][4] == 'uniform':
					f.write(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {sd_vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' -> '+str(float(model_param[i][2]))+' (uniform)'+'\n')
				if model_param[i][4] == 'normal':
					f.write(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {sd_vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' \u00B1 '+str(float(model_param[i][2]))+' (normal)'+'\n')
			else:
				if experror == True:
					if model_param[i][4] == 'uniform':
						f.write(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {float(df_vp[i]):.2e}, bounds: '+ str(float(model_param[i][1])) +' -> '+str(float(model_param[i][2]))+' (uniform)'+'\n')
					if model_param[i][4] == 'normal':
						f.write(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {float(df_vp[i]):.2e}, bounds: '+ str(float(model_param[i][1])) +' \u00B1 '+str(float(model_param[i][2]))+' (normal)'+'\n')
				else:
					if model_param[i][4] == 'uniform':
						f.write(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' -> '+str(float(model_param[i][2]))+' (uniform)'+'\n')
					if model_param[i][4] == 'normal':
						f.write(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' \u00B1 '+str(float(model_param[i][2]))+' (normal)'+'\n')


		for curve in range(0,num_curves):
			f.write(' ')
			f.write('curve#'+str(curve+1)+f' chi^2 = {chi_s[curve]:.2e}\n')


	#export final model results

	fitted_global=[]
	fitted_multi=[]
	fitted_bkg=[]
	fitted_scale=[]
	print('--------------------------------------------------------------------')
	print('Final model parameters')
	print('--------------------------------------------------------------------')
	for curve in range(0,num_curves):
		ii,jj,kk=sympy.symbols('ii jj kk', integer=True)
		x,y=sympy.symbols('x y')

		p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n=sympy.symbols('p0 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 n')
		m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39=sympy.symbols('m0 m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 m20 m21 m22 m23 m24 m25 m26 m27 m28 m29 m30 m31 m32 m33 m34 m35 m36 m37 m38 m39')

		vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,m_param):
			vp[i]=results.x[sum(mnlayers)*5+i]

		vm=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		for i in range(0,c_param):
			#vm[i]=results.x[np.size(layers_max,0)*5+m_param+2*num_curves+c_param*curve+i]
			vm[i]=results.x[sum(mnlayers)*5+m_param+2*num_curves+num_curves*i+curve]


		layers=[]
		for k in range(len(layers_min)):
			sub_layers=[]
			for i in range(0,np.size(layers_min[k],0)):
				line=[]
				for j in range(0,5):
					if isinstance(layers_max[k][i][j], str):
						expr = sympy.sympify(layers_max[k][i][j])
						line.append(np.float(expr.subs([(n,i),(p0,vp[0]),(p1,vp[1]),(p2,vp[2]),(p3,vp[3]),(p4,vp[4]),(p5,vp[5]),(p6,vp[6]),(p7,vp[7]),(p8,vp[8]),(p9,vp[9]),(p10,vp[10]),(p11,vp[11]),(p12,vp[12]),(p13,vp[13]),(p14,vp[14]),(p15,vp[15]),(p16,vp[16]),(p17,vp[17]),(p18,vp[18]),(p19,vp[19]),(m0,vm[0]),(m1,vm[1]),(m2,vm[2]),(m3,vm[3]),(m4,vm[4]),(m5,vm[5]),(m6,vm[6]),(m7,vm[7]),(m8,vm[8]),(m9,vm[9]),(m10,vm[10]),(m11,vm[11]),(m12,vm[12]),(m13,vm[13]),(m14,vm[14]),(m15,vm[15]),(m16,vm[16]),(m17,vm[17]),(m18,vm[18]),(m19,vm[19])])))
					else:
						if k > 0:
							line.append(results.x[5*sum(mnlayers[0:k])+i*5+j])
						else:
							line.append(results.x[i*5+j])
				line.append(layers_max[k][i][5])
				sub_layers.append(line)
			layers.append(sub_layers)

		bkg=results.x[sum(mnlayers)*5+m_param+curve]
		scl=results.x[sum(mnlayers)*5+m_param+num_curves+curve]


		print('\n')
		print('Curve #'+str(curve+1))
		print('\n')
		if resolution[curve] != -1:
			print('Instrumental Resolution, dQ/Q: '+str(resolution[curve]))
		else:
			print('Instrumental Resolution, dQ/Q: pointwise')
		print('fit weigth: '+str(fit_weight[curve]))
		print('--------------------------------------------------------------------')
		count=0
		for k in range(len(layers_min)):
			if len(layers_min) > 1:
				print('Patch #'+str(k+1)+', coverage: '+str(patches[k]))
				print(' ')
			for i in range(0,np.size(layers_min[k],0)):
				for j in range(0,5):
					if layers_min[k][i][j] == layers_max[k][i][j] and not isinstance(layers_min[k][i][j], str):
						if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (fixed)')
						if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (fixed)')
						if j == 2 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (fixed)')
						if j == 3 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (fixed)')
						if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (fixed)')
					elif layers_min[k][i][j] == layers_max[k][i][j] and isinstance(layers_min[k][i][j], str):
						if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (parametric), -> '+layers_min[k][i][0])
						if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (parametric), -> '+layers_min[k][i][1])
						if j == 2 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (parametric), -> '+layers_min[k][i][2])
						if j == 3 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (parametric), -> '+layers_min[k][i][3])
						if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (parametric), -> '+layers_min[k][i][4])
					else:
						if bootstraping or mcmc:
							if j == 0: 
								print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {samples[:,count].mean()/1e-6:.2f} \u00B1 {np.std(samples[:,count])/1e-6:.2f}, bounds: {float(layers_min[k][i][0]/1e-6):.2f} -> {float(layers_max[k][i][0]/1e-6):.2f}')
								count=count+1
							if j == 1: 
								print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {samples[:,count].mean()/1e-6:.2f} \u00B1 {np.std(samples[:,count])/1e-6:.2f}, bounds: {float(layers_min[k][i][1]/1e-6):.2f} -> {float(layers_max[k][i][1]/1e-6):.2f}')
								count=count+1
							if j == 2 and i+1 != np.size(layers_min[k],0): 
								print(layers[k][i][5]+f' thickness (A)              = {samples[:,count].mean():.2f} \u00B1 {np.std(samples[:,count]):.2f}, bounds: {float(layers_min[k][i][2]):.2f} -> {float(layers_max[k][i][2]):.2f}')
								count=count+1
							if j == 3 and i+1 != np.size(layers_min[k],0): 
								print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {samples[:,count].mean():.2f} \u00B1 {np.std(samples[:,count]):.2f}, bounds: {float(layers_min[k][i][3]):.2f} -> {float(layers_max[k][i][3]):.2f}')
								count=count+1
							if j == 4: 
								print(layers[k][i][5]+f' solvent volume fraction    = {samples[:,count].mean():.2f} \u00B1 {np.std(samples[:,count]):.2f}, bounds: {float(layers_min[k][i][4]):.2f} -> {float(layers_max[k][i][4]):.2f}')
								count=count+1
						else:
							if experror == True:
								if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} \u00B1 {(df_layers[k][i][0]/1e-6):.2f}, bounds: {float(layers_min[k][i][0]/1e-6):.2f} -> {float(layers_max[k][i][0]/1e-6):.2f}')
								if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} \u00B1 {(df_layers[k][i][1]/1e-6):.2f}, bounds: {float(layers_min[k][i][1]/1e-6):.2f} -> {float(layers_max[k][i][1]/1e-6):.2f}')
								if j == 2 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} \u00B1 {(df_layers[k][i][2]):.2f}, bounds: {float(layers_min[k][i][2]):.2f} -> {float(layers_max[k][i][2]):.2f}')
								if j == 3 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} \u00B1 {(df_layers[k][i][3]):.2f}, bounds: {float(layers_min[k][i][3]):.2f} -> {float(layers_max[k][i][3]):.2f}')
								if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} \u00B1 {(df_layers[k][i][4]):.2f}, bounds: {float(layers_min[k][i][4]):.2f} -> {float(layers_max[k][i][4]):.2f}')
							else:
								if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f}, bounds: {float(layers_min[k][i][0]/1e-6):.2f} -> {float(layers_max[k][i][0]/1e-6):.2f}')
								if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f}, bounds: {float(layers_min[k][i][1]/1e-6):.2f} -> {float(layers_max[k][i][1]/1e-6):.2f}')
								if j == 2 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f}, bounds: {float(layers_min[k][i][2]):.2f} -> {float(layers_max[k][i][2]):.2f}')
								if j == 3 and i+1 != np.size(layers_min[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f}, bounds: {float(layers_min[k][i][3]):.2f} -> {float(layers_max[k][i][3]):.2f}')
								if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f}, bounds: {float(layers_min[k][i][4]):.2f} -> {float(layers_max[k][i][4]):.2f}')
				print('--------------------------------------------------------------------')
		print(' \n')

		if bootstraping or mcmc:
			if background[curve][2] == 'uniform': 
				print('Background: '+f' = {mbkg[curve]:.2e} \u00B1 {sd_bkg[curve]:.2e}, bounds: '+ str(float(background[curve][0])) +' -> '+str(float(background[curve][1]))+' (uniform)')
			if background[curve][2] == 'normal': 
				print('Background: '+f' = {mbkg[curve]:.2e} \u00B1 {sd_bkg[curve]:.2e}, bounds: '+ str(float(background[curve][0])) +' \u00B1 '+str(float(background[curve][1]))+' (normal)')
			fitted_bkg.append([mbkg[curve],sd_bkg[curve]])
		else:
			if experror == True:
				if background[curve][2] == 'uniform': 
					print('Background: '+f' = {bkg:.2e} \u00B1 {df_bkg[curve]:.2e}, bounds: '+ str(float(background[curve][0])) +' -> '+str(float(background[curve][1]))+' (uniform)')
				if background[curve][2] == 'normal': 
					print('Background: '+f' = {bkg:.2e} \u00B1 {df_bkg[curve]:.2e}, bounds: '+ str(float(background[curve][0])) +' \u00B1 '+str(float(background[curve][1]))+' (normal)')
				fitted_bkg.append([bkg,df_bkg[curve]])
			else:
				if background[curve][2] == 'uniform': 
					print('Background: '+f' = {bkg:.2e}, bounds: '+ str(float(background[curve][0])) +' -> '+str(float(background[curve][1]))+' (uniform)')	
				if background[curve][2] == 'normal': 
					print('Background: '+f' = {bkg:.2e}, bounds: '+ str(float(background[curve][0])) +' \u00B1 '+str(float(background[curve][1]))+' (normal)')
				fitted_bkg.append([bkg,'not available'])				

		if bootstraping or mcmc:
			if scale[curve][2] == 'uniform':
				print('Scale: '+f' = {mscl[curve]:.2e} \u00B1 {sd_scl[curve]:.2e}, bounds: '+ str(float(scale[curve][0])) +' -> '+str(float(scale[curve][1]))+' (uniform)')
			if scale[curve][2] == 'normal':
				print('Scale: '+f' = {mscl[curve]:.2e} \u00B1 {sd_scl[curve]:.2e}, bounds: '+ str(float(scale[curve][0])) +' \u00B1 '+str(float(scale[curve][1]))+' (normal)')
			fitted_scale.append([mscl[curve],sd_scl[curve]])
		else:
			if experror == True:
				if scale[curve][2] == 'uniform':
					print('Scale: '+f' = {scl:.2e} \u00B1 {df_scl[curve]:.2e}, bounds: '+ str(float(scale[curve][0])) +' -> '+str(float(scale[curve][1]))+' (uniform)')
				if scale[curve][2] == 'normal':
					print('Scale: '+f' = {scl:.2e} \u00B1 {df_scl[curve]:.2e}, bounds: '+ str(float(scale[curve][0])) +' \u00B1 '+str(float(scale[curve][1]))+' (normal)')
				fitted_scale.append([scl,df_scl[curve]])
			else:
				if scale[curve][2] == 'uniform':
					print('Scale: '+f' = {scl:.2e}, bounds: '+ str(float(scale[curve][0])) +' -> '+str(float(scale[curve][1]))+' (uniform)')
				if scale[curve][2] == 'normal':
					print('Scale: '+f' = {scl:.2e}, bounds: '+ str(float(scale[curve][0])) +' \u00B1 '+str(float(scale[curve][1]))+' (normal)')
				fitted_scale.append([scl,'not available'])

		#print('Multi Parameters:\n')
		h_list=[]
		for i in range(0,c_param):
			if bootstraping or mcmc:
				if multi_param[i][-1] == 'uniform':
					print(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {m_vm[curve][i]:.2e} \u00B1 {float(sd_vm[curve][i]):.2e}, bounds: '+ str(float(multi_param[i][1+2*curve])) +' -> '+str(float(multi_param[i][2+2*curve]))+' (uniform)')
				if multi_param[i][-1] == 'normal':
					print(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {m_vm[curve][i]:.2e} \u00B1 {float(sd_vm[curve][i]):.2e}, bounds: '+ str(float(multi_param[i][1+2*curve])) +' \u00B1 '+str(float(multi_param[i][2+2*curve]))+' (normal)')
				h_list.append([m_vm[curve][i],sd_vm[curve][i]])
			else:
				if experror == True:
					if multi_param[i][-1] == 'uniform':
						print(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e} \u00B1 {float(df_vm[i*num_curves+curve]):.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' -> '+str(float(multi_param[i][2*curve+2]))+' (uniform)')
					if multi_param[i][-1] == 'normal':
						print(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e} \u00B1 {float(df_vm[i*num_curves+curve]):.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' \u00B1 '+str(float(multi_param[i][2*curve+2]))+' (normal)')
					h_list.append([vm[i],df_vm[i*num_curves+curve]])
				else:
					if multi_param[i][-1] == 'uniform':
						print(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' -> '+str(float(multi_param[i][2*curve+2]))+' (uniform)')
					if multi_param[i][-1] == 'normal':
						print(str(multi_param[i][1+2*num_curves])+': '+'m'+str(i)+f' = {vm[i]:.2e}, bounds: '+ str(float(multi_param[i][2*curve+1])) +' \u00B1 '+str(float(multi_param[i][2*curve+2]))+' (normal)')
					h_list.append([vm[i],'not available'])
		fitted_multi.append(h_list)

	print('\n')
	print('Parameters:\n')
	for i in range(0,m_param):
		if bootstraping or mcmc:
			if model_param[i][4] == 'uniform':
				print(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {sd_vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' -> '+str(float(model_param[i][2]))+' (uniform)')
			if model_param[i][4] == 'normal':
				print(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {sd_vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' \u00B1 '+str(float(model_param[i][2]))+' (normal)')
			fitted_global.append([vp[i],sd_vp[i]])
		else:
			if experror == True:
				if model_param[i][4] == 'uniform':
					print(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {float(df_vp[i]):.2e}, bounds: '+ str(float(model_param[i][1])) +' -> '+str(float(model_param[i][2]))+' (uniform)')
				if model_param[i][4] == 'normal':
					print(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e} \u00B1 {float(df_vp[i]):.2e}, bounds: '+ str(float(model_param[i][1])) +' \u00B1 '+str(float(model_param[i][2]))+' (normal)')
				fitted_global.append([vp[i],df_vp[i]])
			else:
				if model_param[i][4] == 'uniform':
					print(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' -> '+str(float(model_param[i][2]))+' (uniform)')
				if model_param[i][4] == 'normal':
					print(str(model_param[i][3])+': '+'p'+str(i)+f' = {vp[i]:.2e}, bounds: '+ str(float(model_param[i][1])) +' \u00B1 '+str(float(model_param[i][2]))+' (normal)')
				fitted_global.append([vp[i],'not available'])

	print('\n')
	for curve in range(0,num_curves):
		print('curve#'+str(curve+1)+f' chi^2 = {chi_s[curve]:.2e}')
	print('\n')


	if experror==True:
		print('\n')
		print('Note: reported chi^2 is given by 1/N x sum[((R-R_exp)/sigma_exp)^2]')
		print('      where: N is the number of experimental points, sigma_exp the experimental uncertainty')
		print('      R and R_exp the theoretical and experimental reflectivity respectively.')
		if project != 'none':
			f.write('\n')
			f.write('Note: reported chi^2 is given by 1/N x sum[((R-R_exp)/sigma_exp)^2]\n')
			f.write('      where: N is the number of experimental points, sigma_exp the experimental uncertainty\n')
			f.write('      R and R_exp the theoretical and experimental reflectivity respectively.\n')
	else:
		print('\n')
		print('Note: reported chi^2 is given by 1/N x sum[(R-R_exp)^2]')
		print('      where: N is the number of experimental points,')
		print('      R and R_exp the theoretical and experimental reflectivity respectively.')
		if project != 'none':
			f.write('\n')
			f.write('Note: reported chi^2 is given by 1/N x sum[(R-R_exp)^2]\n')
			f.write('      where: N is the number of experimental points,\n')
			f.write('      R and R_exp the theoretical and experimental reflectivity respectively.\n')

	print('\n')
	print("Total calculation time (sec): ", np.int(elapsed_time))
	print(' ')

	if project != 'none':
		f.write('\n')
		f.write("Total calculation time (sec): "+str(np.int(elapsed_time)))
		f.write('\n')


	print('Library versions used for the calculations:')
	print('numpy: '+np.__version__)
	print('scipy: '+scipy.__version__)
	print('numdifftools: '+nd.__version__)
	print('sympy: '+sympy.__version__)
	print('emcee: '+emcee.__version__)

	if project != 'none':
		f.write('Library versions used for the calculations:\n')
		f.write('numpy: '+np.__version__+'\n')
		f.write('scipy: '+scipy.__version__+'\n')
		f.write('numdifftools: '+nd.__version__+'\n')
		f.write('sympy: '+sympy.__version__+'\n')
		f.write('emcee: '+emcee.__version__+'\n')

		f.close()
		os.chdir('..')

	if bootstraping or mcmc: 
		return Refl2,Profile_mean,Solvent_mean,fitted_global,fitted_multi,fitted_bkg,fitted_scale,chi_s
	else:
		return Refl2,Profile,Solvent,fitted_global,fitted_multi,fitted_bkg,fitted_scale,chi_s

def calculate(project,resolution, patches, system, system_param, background, scale, qmax, plot):
	"""
	--------------------------------------------------------------------
	Program ANAKLASIS - Calculation Module for X-ray/Neutron reflection
	version 1.2, February 2021
	developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ
	for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de
	--------------------------------------------------------------------
	This function performs x-ray and neutron reflectivity curve calculations 

	:param project: string containing the name of the project. All output files
	 are saved in a directory with the same name. If project name is 'none'
	 no output files re written on disk. 

	:resolution: list of single float. The element is the dQ/Q resolution

	Syntax:

		resolution = [res1]

	:patches: list of surface coverage of each defined model. In case of a single
	 defined model (most usual case) the definition has the following syntax

	Syntax:
		patches = [1.0]

	in case of K defined models

	Syntax:
	 	patches = [coverage_1, coverage_2 ... coverage_K]

	where the sum of all coverages should add up to 1.

	:system: List of lists containing defined models.
	 Each model is represented as a list of N+1 lists(lines) that contain 6 elements.

	Syntax:
		system = [[
			[  Re_sld0, Im_sld0, d0, sigma0, solv0, 'layer0'],
			[  Re_sld1, Im_sld1, d1, sigma1, solv1, 'layer1'],
			[  Re_sld2, Im_sld2, d2, sigma2, solv2, 'layer2'],
			.
			.
			.
			[  Re_sldN, Im_sldN, dN, sigmaN, solvN, 'layerN'],
			]]

	 If more than one models needs to be defined, more lists can be appended as:

	 	system.append([
			[  Re_sld0, Im_sld0, d0, sigma0, solv0, 'layer0'],
			[  Re_sld1, Im_sld1, d1, sigma1, solv1, 'layer1'],
			[  Re_sld2, Im_sld2, d2, sigma2, solv2, 'layer2'],
			.
			.
			.
			[  Re_sldN, Im_sldN, dN, sigmaN, solvN, 'layerN'],
	 		])

	 Each line (6 element list) represents a layer, from layer 0 (semi-infinite
	 fronting) to layer N (semi-infinite backing). The elements of the list 
	 correspond to Real sld (in 10^-6 A^-2) , Imaginary sld in 10^-6 A^-2,
	 thickness (in Angstrom), roughness (in Angstrom), solvent volumer fraction
	 (0 to 1) and layer description (string) respectively. All elements (except
	 description) can be numerical values or sympy expressions (string) involving
	 global and multi parameters. Additionally in the sympy expressions the 
	 integer (n) can be used, that represents the number of the layer from 0 
	 to N, and/or the summation integers ii,jj,kk, and/or the variables x,y,z
	 than may be used in SymPy integrals or derivatives.


	 When solv0 = 0 and solvN = 1 (fronting and backing solvent volume fraction)
	 then the solvent volume fraction parameter assumes that the backing layer
	 represents a semi-infinite liquid medium and that the liquid may
	 penetrate layers 1 to N-1 (usual in measurements at the solid/liquid or
	 air/liquid interface). 

	 When solv0 = 1 and solvN = 0 (fronting and backing solvent volume fraction)
	 then the solvent volume fraction parameter assumes that the fronting layer
	 represents a semi-infinite liquid medium and that the liquid may
	 penetrate layers 1 to N-1 (usual in measurements at the solid/liquid or
	 air/liquid interface). 

	 When solv0 = 0 and solvN = 0 (fronting and backing solvent volume fraction)
	 all solv values should be set zero. Any non zero value is ignored. 

	 Note that sigma_i represents the roughness between layer_i and layer_(i+1) 

	:system_param: Global parameter list of X 3-element lists.

	Syntax:
		system_param = [
			['p0', p0_value, 'p0_description'],
			['p1', p1_value, 'p1_description'],
			.
			.
			.
			['pX', pX_valuex, 'pX_name'],
			]

	 Up to X=40 global parameters can be defined. The name should be strings of
	 the form 'p0','p1' ... 'p39' respectively. The last element of each global
	 parameter is also a string (description). The middle elements are floats
	 corresponding to the value of the parameter.

	:background: list with single numerical element that corresponds to the 
	background.

	Syntax: 

		background = [bkg_value]	

	 Note: Theoretical reflectivity is calculated as

	 	R = scale * R(Q) + background

	:scale: list with single numerical element that corresponds to the 
	scaling of the reflectivity curve.

	Syntax: 

		scale = [scale_value]

	 Note: Theoretical reflectivity is calculated as

	 	R = scale * R(Q) + background

	:param plot: Boolean. If True, an interactive plot is displayed at the end
	of the calculations. 

	:qmax: list with single numerical element that corresponds to the max
	Q in inverse Angstrom.

	Syntax:

		qmax = [Qmax_value]

	:return: a tuple containing 3 lists: 

	return[0][i][Q,R,RxQ^4] <- reflectivity list
	return[1][i][z,sld] <- sld list
	return[2][i][z,solv] <- solvent volume fraction list

	where i is the model number starting from 0

	 In case of error a string that describes the error that occurred, is returned.

	 Output files are written in the created directory 'project' and include
	 a log file, reflectivity curves in R vs Q and R vs Q^4, solvent volume 
	 fraction and scattering length density profiles. Also corresponding PDF
	 figures are saved together with the ASCII data files.
	"""

	#Increase recursion depth for Sympy
	#sys.setrecursionlimit(100000)

	print('--------------------------------------------------------------------')
	print('Program ANAKLASIS - Calculation Module for X-ray/Neutron reflection ')
	print('version 1.2, February 2021')
	print('developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ')
	print('for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de')
	print('--------------------------------------------------------------------')


	num_curves=1
	#Check of defined model
	if np.size(resolution,0) != num_curves:
		print('Incosistent number of entries for instrumental resolution!')
		return 'Incosistent number of entries for instrumental resolution!'

	if np.size(background,0) != num_curves:
		print('Incosistent number of entries for instrumental background!')
		return 'Incosistent number of entries for instrumental background!'

	for i in range(np.size(background,0)):
		if np.size(background[i]) != 1:
			print('Defined background needs exactly one entry!')
			return 'Defined background needs exactly one entry!'

	if np.size(scale,0) != num_curves:
		print('Incosistent number of entries for reflectivity scaling!')
		return 'Incosistent number of entries for reflectivity scaling!'	

	for i in range(np.size(scale,0)):
		if np.size(scale[i]) != 1:
			print('Defined scale needs exactly one entry!')
			return 'Defined scale needs exactly one entry!'
	
	if len(system) != len(patches):
		print('Coverage fractions in patches should be the same as the defined number of systems!')
		return 'Coverage fractions in patches should be the same as the defined number of systems!'

	total_coverage=0.0
	for i in range(np.size(patches,0)):
		total_coverage=total_coverage+patches[i]

	if total_coverage != 1.0:
		print('Total coverage from all patches should be equal to 1!')
		return 'Total coverage from all patches should be equal to 1!'		

	for j in range(0,len(system)):
		for i in range(np.size(system[j],0)):
			if not isinstance(system[j][i][4], str):
				if system[j][i][4] < 0 or system[j][i][4] > 1:
					print('Invalid solvent volume fraction in layer #'+str(i)+', patch #'+str(j+1))
					print('it should be between 0 and 1!')
					return 'Invalid solvent volume fraction in layer #'+str(i)+', patch #'+str(j+1)
			if np.size(system[j][i]) != 6:
				print('Defined model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1))
				print('correct syntax is: [ real sld, imaginary sld, thickness, roughness, solvent volume fraction, name],')
				return 'Defined model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1)				
			if not isinstance(system[j][i][5], str):
				print('Name entry should be a string, model layer #'+str(i)+', patch #'+str(j+1))
				return 'Name entry should be a string, model layer #'+str(i)+', patch #'+str(j+1)		


	if np.size(system_param,0) > 40:
		print('maximum number of model parameters is equal to 20')
		return 'maximum number of model parameters is equal to 20'

	for i in range(np.size(system_param,0)):
		if system_param[i][0] != 'p'+str(i):
			print('parameter #'+str(i)+' should be named p'+str(i))
			return 'parameter #'+str(i)+' should be named p'+str(i)
		if not isinstance(system_param[i][2], str):
			print('description of parameter #'+str(i)+' should be a string!')
			return 'description of parameter #'+str(i)+' should be a string!'	
		if isinstance(system_param[i][1], str):	
			print('value of parameter #'+str(i)+' should be a number!')
			return 'value of parameter #'+str(i)+' should be a number!'	
		if np.size(system_param[i]) != 3:
			print('Number of entries for p'+str(i)+' is wrong!')
			print('correct syntax is: [parameter name, value, parameter description],')
			return 'Number of entries for p'+str(i)+' is wrong!'

	# These values of the semi-infinite fronting and backing have no physical meaning 
	for i in range(len(system)):
		system[i][0][2]=0.0
		system[i][-1][2]=0.0
		system[i][-1][3]=0.0

	if project != 'none':
		project.strip('/')
		folder="project-"+project
		if not os.path.exists(folder):
	   		os.makedirs(folder)
		else:
			print('Directory already exists.. overwriting data..')
		os.chdir(folder)

		f = open(project+"_calculation_parameters.log", "w")
		f.write('--------------------------------------------------------------------\n')
		f.write('Program ANAKLASIS - Calculation Module for X-ray/Neutron reflection \n')
		f.write('version 1.2, February 2021\n')
		f.write('developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ\n')
		f.write('for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de\n')
		f.write('--------------------------------------------------------------------\n')

	ii,jj,kk=sympy.symbols('ii jj kk', integer=True)
	x,y=sympy.symbols('x y')
	p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n=sympy.symbols('p0 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 n')

	vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	for i in range(0,np.size(system_param,0)):
		vp[i]=system_param[i][1]


	layers=[]
	for k in range(len(system)):
		sub_layers=[]
		for i in range(0,np.size(system[k],0)):
			line=[]
			for j in range(0,5):
				if isinstance(system[k][i][j], str):
					expr = sympy.sympify(system[k][i][j])
					line.append(np.float(expr.subs([(n,i),(p0,vp[0]),(p1,vp[1]),(p2,vp[2]),(p3,vp[3]),(p4,vp[4]),(p5,vp[5]),(p6,vp[6]),(p7,vp[7]),(p8,vp[8]),(p9,vp[9]),(p10,vp[10]),(p11,vp[11]),(p12,vp[12]),(p13,vp[13]),(p14,vp[14]),(p15,vp[15]),(p16,vp[16]),(p17,vp[17]),(p18,vp[18]),(p19,vp[19]),(p20,vp[20]),(p21,vp[21]),(p22,vp[22]),(p23,vp[23]),(p24,vp[24]),(p25,vp[25]),(p26,vp[26]),(p27,vp[27]),(p28,vp[28]),(p29,vp[29]),(p30,vp[30]),(p31,vp[31]),(p32,vp[32]),(p33,vp[33]),(p34,vp[34]),(p35,vp[35]),(p36,vp[36]),(p37,vp[37]),(p38,vp[38]),(p39,vp[39])])))
				else:
					line.append(system[k][i][j])
			line.append(system[k][i][5])
			sub_layers.append(line)
		layers.append(sub_layers)


	#print('dQ/Q = ',resolution)
	#np.set_printoptions(threshold=sys.maxsize)
	#pprint.pprint(layers)

	q_bin = np.linspace(0.001, np.float(qmax[0]), 1001)
	res_bin = np.zeros(np.size(q_bin,0))
	for i in range(np.size(res_bin,0)):
		res_bin[i] = resolution[0]


	Refl=Reflectivity(q_bin, res_bin, layers, resolution[0], np.float(background[0]),scale[0],patches)

	Profile=profile(layers, 1000)

	Solvent=solvent_penetration(layers, 1000)

	if project != 'none':
		np.savetxt(project+"_reflectivity_curve.dat", Refl, fmt='%1.4e', header="Q (A^-1),  R, RQ^4 (A^-4)")
		for i in range(len(system)):
			if len(system) > 1:
				np.savetxt(project+"_sld_profile_patch#"+str(i+1)+".dat", Profile[i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2)")
				np.savetxt(project+"_solvent_profile_patch#"+str(i+1)+".dat", Solvent[i], fmt='%1.4e',header="z (A),  solvent volume fraction")
			else:
				np.savetxt(project+"_sld_profile.dat", Profile[i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2)")
				np.savetxt(project+"_solvent_profile.dat", Solvent[i], fmt='%1.4e',header="z (A),  solvent volume fraction")

	print('\n')
	print('Instrumental Resolution, dQ/Q: '+str(resolution[0]))
	print('--------------------------------------------------------------------')
	for k in range(0,len(system)):
		if len(system) > 1:
			print('Patch #'+str(k+1)+', coverage: '+str(patches[k]))
			print('')
		for i in range(0,np.size(system[k],0)):
			for j in range(0,5):
				if not isinstance(system[k][i][j], str):
					if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} ')
					if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} ')
					if j == 2 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} ')
					if j == 3 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} ')
					if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} ')
				else:
					if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (parametric), -> '+system[k][i][0])
					if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (parametric), -> '+system[k][i][1])
					if j == 2 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (parametric), -> '+system[k][i][2])
					if j == 3 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (parametric), -> '+system[k][i][3])
					if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (parametric), -> '+system[k][i][4])
			print('--------------------------------------------------------------------')
	print(' \n')
	print('Background: '+f' = {background[0]:.2e}')
	print('Scale: '+f' = {scale[0]:.2e}')
	print('\n')
	print('Parameters:\n')
	for i in range(0,np.size(system_param,0)):
		print(str(system_param[i][2])+': '+'p'+str(i)+f' = {vp[i]:.2e}')

	if project != 'none':
		f.write('\n')
		f.write('Instrumental Resolution, dQ/Q: '+str(resolution[0])+'\n')
		f.write('--------------------------------------------------------------------\n')
		for k in range(0,len(system)):
			if len(system) > 1:
				f.write('Patch #'+str(k+1)+', coverage: '+str(patches[k])+'\n')
				f.write('\n')
			for i in range(0,np.size(system[k],0)):
				for j in range(0,5):
					if not isinstance(system[k][i][j], str):
						if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} \n')
						if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} \n')
						if j == 2 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} \n')
						if j == 3 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} \n')
						if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} \n')
					else:
						if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (parametric), -> '+system[k][i][0]+'\n')
						if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (parametric), -> '+system[k][i][1]+'\n')
						if j == 2 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (parametric), -> '+system[k][i][2]+'\n')
						if j == 3 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (parametric), -> '+system[k][i][3]+'\n')
						if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (parametric), -> '+system[k][i][4]+'\n')
				f.write('--------------------------------------------------------------------\n')
		f.write(' \n')
		f.write('Background: '+f' = {background[0]:.2e}'+'\n')
		f.write('Scale: '+f' = {scale[0]:.2e}'+'\n')
		f.write('\n')
		f.write('Parameters:\n')
		for i in range(0,np.size(system_param,0)):
			f.write(str(system_param[i][2])+': '+'p'+str(i)+f' = {vp[i]:.2e}\n')

		plinestyle=['-','--',':','-.','.','--',':']
		#print(Refl)
		#plt.figure(figsize=(9, 7))
		plt.subplot(2,2,1)
		plt.plot(Refl[:,0],Refl[:,1],color='red')
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)$')
		#plt.title(r'$R(Q)$')
		plt.grid(True)

		#ax2=plt.add_axes((.1,.1,.8,.2))
		#ax2.plot(Refl[:,0],1.0,color='red')

		plt.subplot(2,2,2)
		plt.plot(Refl[:,0],Refl[:,2],color='red')
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)Q^4$')
		#plt.title(r'$R(Q)Q^4$')
		plt.grid(True)
		#plt.savefig("test.png")

		plt.subplot(2,2,3)
		for i in range(len(system)):
			plt.plot(Profile[i][:,0],Profile[i][:,1]*1e6,color='green',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel(r'$sld(10^{-6}\AA^{-2})$')

		plt.subplot(2,2,4)
		for i in range(len(system)):
			plt.plot(Solvent[i][:,0],Solvent[i][:,1],color='blue',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel('solvent')

		plt.tight_layout()
		plt.draw() 
		plt.savefig(project+'_summary.pdf')
		#plt.savefig(project+'_summary.ps')
		if plot is True: plt.show()

		plt.figure()
		plt.plot(Refl[:,0],Refl[:,1],color='red')
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)$')
		#plt.title(r'$R(Q)$')
		plt.grid(True)
		plt.draw() 
		plt.savefig(project+'_reflectivity.pdf')	

		plt.figure()
		plt.plot(Refl[:,0],Refl[:,2],color='red')
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)Q^4$')
		#plt.title(r'$R(Q)Q^4$')
		plt.grid(True)
		plt.draw() 
		plt.savefig(project+'_reflectivity_RQ^4.pdf')

		plt.figure()
		for i in range(len(system)):
			plt.plot(Profile[i][:,0],Profile[i][:,1]*1e6,color='green',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel(r'$sld(10^{-6}\AA^{-2})$')
		plt.grid(True)
		plt.draw() 
		plt.savefig(project+'_sld_profile.pdf')

		plt.figure()
		for i in range(len(system)):
			plt.plot(Solvent[i][:,0],Solvent[i][:,1],color='blue',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel('solvent')
		plt.draw() 
		plt.savefig(project+'_solvent_profile.pdf')

		plt.close()  

		os.chdir('..') 


	print('Library versions used for the calculations:')
	print('numpy: '+np.__version__)
	print('scipy: '+scipy.__version__)
	print('numdifftools: '+nd.__version__)
	print('sympy: '+sympy.__version__)

	if project != 'none':
		f.write('Library versions used for the calculations:\n')
		f.write('numpy: '+np.__version__+'\n')
		f.write('scipy: '+scipy.__version__+'\n')
		f.write('numdifftools: '+nd.__version__+'\n')
		f.write('sympy: '+sympy.__version__+'\n')

	return Refl,Profile,Solvent

def compare(project, in_file, units, resolution, patches, system, system_param, background, scale, qmax,  experror, plot):
	"""
	--------------------------------------------------------------------
	Program ANAKLASIS - Comparison Module for X-ray/Neutron reflection
	version 1.2, February 2021
	developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ
	for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de
	--------------------------------------------------------------------
	This function performs comparison of x-ray and neutron reflection data
	with theoretical reflectivity from a defined model.

	:param project: string containing the name of the project. All output files
	 are saved in a directory with the same name. If project name is 'none'
	 no output files re written on disk. 

	:param in_file: list of single string element. 

	 Syntax:

	 in_file = ['file']

	 which is the path to an ASCII file containing the reflectivity data. 
	 The file structure should be in column format ( Q, Reflectivity, error in 
	 reflectivity (dR), error in Q (dQ)). If a third column (dR) is not present
	 in the data, set experror=False (see below). If a 4th column (dQ) is 
	 present you may take pointwise resolution into account by setting resolution
	 equal to -1 (see below). Additional columns will be ignored. 
	 Lines beginning with '#' can be present and are considered as comments. 
	 Reflectivity should be footprint corrected.

	:param units: list of a single string element that can be either 'A' inverse
	 Angstrom or 'nm' inverse nanometers, describing the units of momentum
	 transfer (Q) in the input file.  

	 Syntax:

	 units = ['units1']

	:resolution: list of single float element corresponding to the dQ/Q resolution 
	of the corresponding input data

	Syntax:

		resolution = [res_value]

	Note that if a 4th dQ column is present in the in_file, you may set res_value
	 = -1, so that pointwise dQ/Q resolution is taken into account.

	:patches: list of surface coverage of each defined model. In case of a single
	 defined model (most usual case) the definition has the following syntax

	Syntax:
		patches = [1.0]

	in case of K defined models

	Syntax:
	 	patches = [coverage_1, coverage_2 ... coverage_K]

	where the sum of all coverages should add up to 1.

	:system: List of lists containing defined models.
	 Each model is represented as a list of N+1 lists that contain 6 elements.


	Syntax:
		system = [[
			[  Re_sld0, Im_sld0, d0, sigma0, solv0, 'layer0'],
			[  Re_sld1, Im_sld1, d1, sigma1, solv1, 'layer1'],
			[  Re_sld2, Im_sld2, d2, sigma2, solv2, 'layer2'],
			.
			.
			.
			[  Re_sldN, Im_sldN, dN, sigmaN, solvN, 'layerN'],
			]]

	 If more than one models needs to be defined, more lists can be appended as:

	 	system.append([
			[  Re_sld0, Im_sld0, d0, sigma0, solv0, 'layer0'],
			[  Re_sld1, Im_sld1, d1, sigma1, solv1, 'layer1'],
			[  Re_sld2, Im_sld2, d2, sigma2, solv2, 'layer2'],
			.
			.
			.
			[  Re_sldN, Im_sldN, dN, sigmaN, solvN, 'layerN'],
	 		])

	 Each line (6 element list) represents a layer, from layer 0 (semi-infinite
	 fronting) to layer N (semi-infinite backing). The elements of the list 
	 correspond to Real sld (in 10^-6 A^-2) , Imaginary sld in 10^-6 A^-2,
	 thickness (in Angstrom), roughness (in Angstrom), solvent volumer fraction
	 (0 to 1) and layer description (string) respectively. All elements (except
	 description) can be numerical values or sympy expressions (string) involving
	 global and multi parameters. Additionally in the sympy expressions the 
	 integer (n) can be used, that represents the number of the layer from 0 
	 to N, and/or the summation integers ii,jj,kk, and/or the variables x,y,z
	 than may be used in SymPy integrals or derivatives.


	 When solv0 = 0 and solvN = 1 (fronting and backing solvent volume fraction)
	 then the solvent volume fraction parameter assumes that the backing layer
	 represents a semi-infinite liquid medium and that the liquid may
	 penetrate layers 1 to N-1 (usual in measurements at the solid/liquid or
	 air/liquid interface). 

	 When solv0 = 1 and solvN = 0 (fronting and backing solvent volume fraction)
	 then the solvent volume fraction parameter assumes that the fronting layer
	 represents a semi-infinite liquid medium and that the liquid may
	 penetrate layers 1 to N-1 (usual in measurements at the solid/liquid or
	 air/liquid interface). 

	 When solv0 = 0 and solvN = 0 (fronting and backing solvent volume fraction)
	 all solv values should be set zero. Any non zero value is ignored. 

	 Note that sigma_i represents the roughness between layer_i and layer_(i+1) 

	:system_param: Global parameter list of X 3-element lists.

	Syntax:
		system_param = [
			['p0', p0_value, 'p0_description'],
			['p1', p1_value, 'p1_description'],
			.
			.
			.
			['pX', pX_valuex, 'pX_name'],
			]

	 Up to X=40 global parameters can be defined. The name should be strings of
	 the form 'p0','p1' ... 'p39' respectively. The last element of each global
	 parameter is also a string (description). The middle elements are floats
	 corresponding to the value of the parameter.

	:background: list with single numerical element that corresponds to the 
	background.

	Syntax: 

		background = [bkg_value]	

	 Note: Theoretical reflectivity is calculated as

	 	R = scale * R(Q) + background

	:scale: list with single numerical element that corresponds to the 
	scaling of the reflectivity curve.

	Syntax: 

		scale = [scale_value]

	 Note: Theoretical reflectivity is calculated as

	 	R = scale * R(Q) + background

	:param plot: Boolean. If True, an interactive plot is displayed at the end
	of the calculations.

	:return: a tuple containing 4 lists: 

	return[0][i][Q,R,RxQ^4] <- reflectivity list
	return[1][i][z,sld] <- sld list
	return[2][i][z,solv] <- solvent volume fraction list
	return[3][value] <- chi squared list 

	where i the model number starting from 0

	In case of error a string that describes the error that occurred, is returned.

	:param experror: Boolean. Set as True, if all input files contain a third 
	column with Reflectivity uncertainty. Set as False, if at least one of the
	input files contains no Reflectivity uncertainty. Also, if all data contain
	a third column but for some reason Reflectivity uncertainty is considered
	as non-reliable, you may set the parameter to False thus ignoring errors.

	:param plot: Boolean. If True, an interactive plot is displayed at the end
	of the calculations.

	:param fast: Boolean. If True, a small differential evolution population size
	is used. If False a large differential evolution population size is used.

	:return: Reflectivity list, Profile List, Solvent volume fraction list and
	chi squared list. Reflectivity list contains (Q,R,RQ^4) data. Profile and
	Solvent list have dimensions model X 2 (z, Profile).
	In case of error a string that describes the error that occured. is returned. 

	 Output files are written in the created directory 'project' and include
	 a log file, reflectivity curves in R vs Q and R vs Q^4, solvent volume 
	 fraction and scattering length density profiles. Also corresponding PDF
	 figures are saved together with the ASCII data files.
	"""

	#Increase recursion depth for Sympy
	#sys.setrecursionlimit(100000)

	print('--------------------------------------------------------------------')
	print('Program ANAKLASIS - Comparison Module for X-ray/Neutron reflection ')
	print('version 1.2, February 2021')
	print('developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ')
	print('for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de')
	print('--------------------------------------------------------------------')

	num_curves=1
	#Check of defined model
	if np.size(units,0) != num_curves:
		print('Incosistent number of entries for momentum transfer units!')
		return 'Incosistent number of entries for for momentum transfer units!'

	if units[0] != 'A' and units[0] != 'a' and units[0] != 'nm' and units[0] != 'NM':
		print('Q unints can only be in inverse Angstroms (A) or nanometers (nm)!')
		return 'Q unints can only be in inverse Angstroms (A) or nanometers (nm)!'

	if np.size(resolution,0) != num_curves:
		print('Incosistent number of entries for instrumental resolution!')
		return 'Incosistent number of entries for instrumental resolution!'

	if np.size(background,0) != num_curves:
		print('Incosistent number of entries for instrumental background!')
		return 'Incosistent number of entries for instrumental background!'

	for i in range(np.size(background,0)):
		if np.size(background[i]) != 1:
			print('Defined background needs exactly one entry!')
			return 'Defined background needs exactly one entry!'

	if np.size(scale,0) != num_curves:
		print('Incosistent number of entries for reflectivity scaling!')
		return 'Incosistent number of entries for reflectivity scaling!'	

	for i in range(np.size(scale,0)):
		if np.size(scale[i]) != 1:
			print('Defined scale needs exactly one entry!')
			return 'Defined scale needs exactly one entry!'

	if len(system) != len(patches):
		print('Coverage fractions in patches should be the same as the defined number of systems!')
		return 'Coverage fractions in patches should be the same as the defined number of systems!'

	total_coverage=0.0
	for i in range(np.size(patches,0)):
		total_coverage=total_coverage+patches[i]

	if total_coverage != 1.0:
		print('Total coverage from all patches should be equal to 1!')
		return 'Total coverage from all patches should be equal to 1!'

	for j in range(0,len(system)):
		for i in range(len(system[j])):
			if not isinstance(system[j][i][4], str):
				if system[j][i][4] < 0 or system[j][i][4] > 1:
					print('Invalid solvent volume fraction in layer #'+str(i)+', patch #'+str(j+1))
					print('it should be between 0 and 1!')
					return 'Invalid solvent volume fraction in layer #'+str(i)+', patch #'+str(j+1)
			if len(system[j][i]) != 6:
				print('Defined model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1))
				print('correct syntax is: [ real sld, imaginary sld, thickness, roughness, solvent volume fraction, name],')
				return 'Defined model has an invalid number of entries, layer #'+str(i)+', patch #'+str(j+1)				
			if not isinstance(system[j][i][5], str):
				print('Name entry should be a string, model layer #'+str(i)+', patch #'+str(j+1))
				return 'Name entry should be a string, model layer #'+str(i)+', patch #'+str(j+1)		


	if np.size(system_param,0) > 40:
		print('maximum number of model parameters is equal to 40')
		return 'maximum number of model parameters is equal to 40'

	for i in range(np.size(system_param,0)):
		if system_param[i][0] != 'p'+str(i):
			print('parameter #'+str(i)+' should be named p'+str(i))
			return 'parameter #'+str(i)+' should be named p'+str(i)
		if not isinstance(system_param[i][2], str):
			print('description of parameter #'+str(i)+' should be a string!')
			return 'description of parameter #'+str(i)+' should be a string!'	
		if isinstance(system_param[i][1], str):	
			print('value of parameter #'+str(i)+' should be a number!')
			return 'value of parameter #'+str(i)+' should be a number!'	
		if np.size(system_param[i]) != 3:
			print('Number of entries for p'+str(i)+' is wrong!')
			print('correct syntax is: [parameter name, value, parameter description],')
			return 'Number of entries for p'+str(i)+' is wrong!'


	# These values of the semi-infinite fronting and backing have no physical meaning 
	for i in range(len(system)):
		system[i][0][2]=0.0
		system[i][-1][2]=0.0
		system[i][-1][3]=0.0

	if project != 'none':
		project.strip('/')
		folder="project-"+project
		if not os.path.exists(folder):
	   		os.makedirs(folder)
		else:
			print('Directory already exists.. overwriting data..')

	if experror == True:
		if resolution[0] == -1:
			data = np.loadtxt(in_file, usecols = (0,1,2,3),comments = "#")
		else:
			data = np.loadtxt(in_file, usecols = (0,1,2),comments = "#")
			zerocol = np.zeros((np.size(data,0),1))
			data = np.append(data, zerocol, axis=1)
			data[:,3]=resolution[0]
	else:
		if resolution[0] == -1:
			data = np.loadtxt(in_file, usecols = (0,1,3),comments = "#")
			zerocol = np.zeros((np.size(data,0),1))
			data = np.insert(data, 2,zerocol, axis=1)
		else:
			data = np.loadtxt(in_file, usecols = (0,1),comments = "#")
			zerocol = np.zeros((np.size(data,0),1))
			data = np.append(data, zerocol, axis=1)
			data = np.append(data, zerocol, axis=1)
			data[:,3]=resolution[0]		

	#Check data
	for i in range(np.size(data,0)):
		if data[i][0] <= 0:
			print('Invalid wave-vector tranfer (Q) value present in input data!')
			return 'Invalid wave-vector tranfer (Q) value present in input data!'

		if experror == True and data[i][2] <= 0:
			print('Invalid experimental error (dR) value present in input data!')
			return 'Invalid experimental error (dR) value present in input data!'

		if resolution[0] == -1 and data[j][3] <0:
			print('Invalid (dQ) value present in input data!')
			return 'Invalid (dQ) value present in input data'	

	#Delete data with negative reflectivity
	row_to_delete=[]
	for i in range(np.size(data,0)):
		if data[i][1] <= 0 or data[i][0] > qmax:
			row_to_delete.append(i)

	data = np.delete(data, row_to_delete, axis=0)

	if units[0] == 'nm' or units[0] == 'NM':
		for j in range(np.size(data,0)):
			data[j,0]=data[j,0]/10.0
			data[j,3]=data[j,3]/10.0

	if os.path.isfile(in_file):
		if project != 'none': shutil.copy(in_file, folder+'/input_curve'+'.dat')
	else:
		print('error: input file '+str(in_file)+' does not exit in current directory..')
		return 'error: input file '+str(in_file)+' does not exit in current directory..'

	if project != 'none':
		os.chdir(folder)
		f = open(project+"_comparison_parameters.log", "w")
		f.write('--------------------------------------------------------------------\n')
		f.write('Program ANAKLASIS - Comparison Module for X-ray/Neutron reflection \n')
		f.write('version 1.2, February 2021\n')
		f.write('developed by Dr. Alexandros Koutsioumpas. JCNS @ MLZ\n')
		f.write('for bugs and requests contact: a.koutsioumpas[at]fz-juelich.de\n')
		f.write('--------------------------------------------------------------------\n')
	
	ii=sympy.symbols('ii', integer=True)
	x=sympy.symbols('x')
	p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,n=sympy.symbols('p0 p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 n')

	vp=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	for i in range(0,np.size(system_param,0)):
		vp[i]=system_param[i][1]


	layers=[]
	for k in range(len(system)):
		sub_layers=[]
		for i in range(0,np.size(system[k],0)):
			line=[]
			for j in range(0,5):
				if isinstance(system[k][i][j], str):
					expr = sympy.sympify(system[k][i][j])
					line.append(np.float(expr.subs([(n,i),(p0,vp[0]),(p1,vp[1]),(p2,vp[2]),(p3,vp[3]),(p4,vp[4]),(p5,vp[5]),(p6,vp[6]),(p7,vp[7]),(p8,vp[8]),(p9,vp[9]),(p10,vp[10]),(p11,vp[11]),(p12,vp[12]),(p13,vp[13]),(p14,vp[14]),(p15,vp[15]),(p16,vp[16]),(p17,vp[17]),(p18,vp[18]),(p19,vp[19]),(p20,vp[20]),(p21,vp[21]),(p22,vp[22]),(p23,vp[23]),(p24,vp[24]),(p25,vp[25]),(p26,vp[26]),(p27,vp[27]),(p28,vp[28]),(p29,vp[29]),(p30,vp[30]),(p31,vp[31]),(p32,vp[32]),(p33,vp[33]),(p34,vp[34]),(p35,vp[35]),(p36,vp[36]),(p37,vp[37]),(p38,vp[38]),(p39,vp[39])])))
				else:
					line.append(system[k][i][j])
			line.append(system[k][i][5])
			sub_layers.append(line)
		layers.append(sub_layers)

	#print('dQ/Q = ',resolution[0])
	#np.set_printoptions(threshold=sys.maxsize)
	#pprint.pprint(layers)



	q_bin = np.linspace(data[0][0], qmax[0], 1001)
	res_bin = np.zeros(np.size(q_bin,0))
	for i in range(np.size(res_bin,0)):
		if resolution[0] != -1:
			res_bin[i] = resolution[0]
		else:
			res_bin[i] = np.interp(q_bin[i], data[:,0], data[:,3])

	#print(res_bin)

	Refl=Reflectivity(q_bin, res_bin, layers, resolution[0], np.float(background[0]),scale[0],patches)
	Refl2=Reflectivity(data[:,0], data[:,3], layers, resolution[0], np.float(background[0]),scale[0],patches)

	Profile=profile(layers, 1000)

	Solvent=solvent_penetration(layers, 1000)

	chi_s=chi_square(data, layers, resolution[0], np.float(background[0]),scale[0],patches)


	print('\n')
	if resolution[0] != -1:
		print('Instrumental Resolution, dQ/Q: '+str(resolution[0]))
	else:
		print('Instrumental Resolution, dQ/Q: pointwise')
	print('--------------------------------------------------------------------')
	for k in range(0,len(system)):
		if len(system) > 1:
			print('Patch #'+str(k+1)+', coverage: '+str(patches[k]))
			print('')
		for i in range(0,np.size(system[k],0)):
			for j in range(0,5):
				if not isinstance(system[k][i][j], str):
					if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} ')
					if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} ')
					if j == 2 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} ')
					if j == 3 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} ')
					if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} ')
				else:
					if j == 0: print(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (parametric), -> '+system[k][i][0])
					if j == 1: print(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (parametric), -> '+system[k][i][1])
					if j == 2 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (parametric), -> '+system[k][i][2])
					if j == 3 and i+1 != np.size(layers[k],0): print(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (parametric), -> '+system[k][i][3])
					if j == 4: print(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (parametric), -> '+system[k][i][4])
			print('--------------------------------------------------------------------')
	print(' \n')
	print('Background: '+f' = {background[0]:.2e}')
	print('Scale: '+f' = {scale[0]:.2e}')
	print('\n')
	print('Parameters:\n')
	for i in range(0,np.size(system_param,0)):
		print(str(system_param[i][2])+': '+'p'+str(i)+f' = {vp[i]:.2e}')


	if project != 'none':
		f.write('\n')
		if resolution[0] != -1:
			f.write('Instrumental Resolution, dQ/Q: '+str(resolution[0])+'\n')
		else:
			f.write('Instrumental Resolution, dQ/Q: pointwise'+'\n')
		f.write('--------------------------------------------------------------------\n')
		for k in range(0,len(system)):
			if len(system) > 1:
				f.write('Patch #'+str(k+1)+', coverage: '+str(patches[k])+'\n')
				f.write('\n')
			for i in range(0,np.size(system[k],0)):
				for j in range(0,5):
					if not isinstance(system[k][i][j], str):
						if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} \n')
						if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} \n')
						if j == 2 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} \n')
						if j == 3 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} \n')
						if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} \n')
					else:
						if j == 0: f.write(layers[k][i][5]+f' real sld (10^-6 A^-2)      = {(layers[k][i][0]/1e-6):.2f} (parametric), -> '+system[k][i][0]+'\n')
						if j == 1: f.write(layers[k][i][5]+f' imaginary sld (10^-6 A^-2) = {(layers[k][i][1]/1e-6):.2f} (parametric), -> '+system[k][i][1]+'\n')
						if j == 2 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+f' thickness (A)              = {layers[k][i][2]:.2f} (parametric), -> '+system[k][i][2]+'\n')
						if j == 3 and i+1 != np.size(layers[k],0): f.write(layers[k][i][5]+'/'+layers[k][i+1][5]+f' roughness (A) = {layers[k][i][3]:.2f} (parametric), -> '+system[k][i][3]+'\n')
						if j == 4: f.write(layers[k][i][5]+f' solvent volume fraction    = {layers[k][i][4]:.2f} (parametric), -> '+system[k][i][4]+'\n')
				f.write('--------------------------------------------------------------------\n')
		f.write(' \n')
		f.write('Background: '+f' = {background[0]:.2e}'+'\n')
		f.write('Scale: '+f' = {scale[0]:.2e}'+'\n')
		f.write('\n')
		f.write('Parameters:\n')
		for i in range(0,np.size(system_param,0)):
			f.write(str(system_param[i][2])+': '+'p'+str(i)+f' = {vp[i]:.2e}\n')

	print('chi^2 = ',chi_s)
	if project != 'none': f.write('chi^2 = '+str(chi_s)+'\n')

	if experror==True:
		print('\n')
		print('Note: reported chi^2 is given by 1/N x sum[((R-R_exp)/sigma_exp)^2]')
		print('      where: N is the number of experimental points, sigma_exp the experimental uncertainty')
		print('      R and R_exp the theoretical and experimental reflectivity respectively.')
		if project != 'none':
			f.write('\n')
			f.write('Note: reported chi^2 is given by 1/N x sum[((R-R_exp)/sigma_exp)^2]\n')
			f.write('      where: N is the number of experimental points, sigma_exp the experimental uncertainty\n')
			f.write('      R and R_exp the theoretical and experimental reflectivity respectively.\n')
	else:
		print('\n')
		print('Note: reported chi^2 is given by 1/N x sum[(R-R_exp)^2]')
		print('      where: N is the number of experimental points,')
		print('      R and R_exp the theoretical and experimental reflectivity respectively.')
		if project != 'none':
			f.write('\n')
			f.write('Note: reported chi^2 is given by 1/N x sum[(R-R_exp)^2]\n')
			f.write('      where: N is the number of experimental points,\n')
			f.write('      R and R_exp the theoretical and experimental reflectivity respectively.\n')

	if project != 'none':
		np.savetxt(project+"_reflectivity_curve.dat", Refl, fmt='%1.4e', header="Q (A^-1),  R, RQ^4 (A^-4)")
		for i in range(len(system)):
			if len(system) > 1:
				np.savetxt(project+"_sld_profile_patch#"+str(i+1)+".dat", Profile[i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2)")
				np.savetxt(project+"_solvent_profile_patch#"+str(i+1)+".dat", Solvent[i], fmt='%1.4e',header="z (A),  solvent volume fraction")
			else:
				np.savetxt(project+"_sld_profile.dat", Profile[i], fmt='%1.4e', header="z (A),  sld (10^-6 A^-2)")
				np.savetxt(project+"_solvent_profile.dat", Solvent[i], fmt='%1.4e',header="z (A),  solvent volume fraction")

	if project != 'none':
		plinestyle=['-','--',':','-.','.','--',':']
		#print(Refl)
		plt.figure(figsize=(9, 7))
		grid = plt.GridSpec(6, 2)
		ax1=plt.subplot(grid[0:3, 0])
		plt.xlim([data[0,0]-data[0,0]/2.0,data[-1,0]])
		plt.errorbar(data[:,0], data[:,1], data[:,2],fmt='.k',ecolor='grey',zorder=0)
		plt.plot(Refl[:,0],Refl[:,1],color='red',label=r'$\chi^2 = $'+'{:06.3e}'.format(chi_s))
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)$')
		#plt.title(r'$R(Q)$')
		plt.grid(True)

		plt.subplot(grid[0:, 1])
		plt.xlim([data[0,0]-data[0,0]/2.0,data[-1,0]])
		plt.errorbar(data[:,0], data[:,1]*data[:,0]**4, data[:,2]*data[:,0]**4,fmt='.k',ecolor='grey',zorder=0,label=r'$\chi^2 = $'+'{:06.3e}'.format(chi_s))
		plt.plot(Refl[:,0],Refl[:,2],color='red')
		plt.legend()
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)Q^4$')
		#plt.title(r'$R(Q)Q^4$')
		plt.grid(True)
		#plt.savefig("test.png")

		ax3=plt.subplot(grid[3, 0])
		plt.xlim([data[0,0]-data[0,0]/2.0,data[-1,0]])
		plt.axhline(y=0, lw=1,ls='dashed', color='k')
		if experror == True:
			plt.scatter(data[:,0],(Refl2[:,1]-data[:,1])/(data[:,2]),2,color='red')
		else:
			plt.scatter(data[:,0],(Refl2[:,1]-data[:,1])/Refl2[:,1],2,color='red')
		#plt.yscale('log')
		if experror == True:
			plt.ylabel(r'$\Delta / \sigma$')
		else:
			plt.ylabel(r'$\Delta / R$')
		plt.xlabel(r'$Q(\AA^{-1})$')


		plt.subplot(grid[4:, 0])
		for i in range(len(system)):
			plt.plot(Profile[i][:,0],Profile[i][:,1]*1e6,color='green',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel(r'$sld(10^{-6}\AA^{-2})$',color='green')
		plt.twinx()
		for i in range(len(system)):
			plt.plot(Solvent[i][:,0],Solvent[i][:,1],color='blue',linestyle=plinestyle[i])
		plt.ylabel('solvent',color='blue')

		#plt.subplot(3,2,7)
		#plt.plot(Solvent[:,0],Solvent[:,1],color='green')
		#plt.xlabel(r'$z(\AA)$')
		#plt.ylabel('solvent')
		ax1.get_shared_x_axes().join(ax1, ax3)
		plt.tight_layout()
		plt.draw() 
		plt.savefig(project+'_summary.pdf')
		if plot is True: plt.show()

		plt.figure()
		grid = plt.GridSpec(4, 1)
		ax1=plt.subplot(grid[0:3, 0])
		plt.xlim([data[0,0]-data[0,0]/2.0,data[-1,0]])
		plt.errorbar(data[:,0], data[:,1], data[:,2],fmt='.k',ecolor='grey',zorder=0)
		plt.plot(Refl[:,0],Refl[:,1],color='red',label=r'$\chi^2 = $'+'{:06.3e}'.format(chi_s))
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)$')
		plt.legend()
		#plt.title(r'$R(Q)$')
		plt.grid(True)
		ax3=plt.subplot(grid[3, 0])
		plt.xlim([data[0,0]-data[0,0]/2.0,data[-1,0]])
		plt.axhline(y=0, lw=1,ls='dashed', color='k')
		if experror == True:
			plt.scatter(data[:,0],(Refl2[:,1]-data[:,1])/(data[:,2]),2,color='red')
		else:
			plt.scatter(data[:,0],(Refl2[:,1]-data[:,1])/Refl2[:,1],2,color='red')
		#plt.yscale('log')
		if experror == True:
			plt.ylabel(r'$\Delta / \sigma$')
		else:
			plt.ylabel(r'$\Delta / R$')
		plt.xlabel(r'$Q(\AA^{-1})$')	
		ax1.get_shared_x_axes().join(ax1, ax3)
		plt.tight_layout()
		plt.draw() 
		plt.savefig(project+'_reflectivity.pdf')

		plt.figure()
		plt.xlim([data[0,0]-data[0,0]/2.0,data[-1,0]])
		plt.errorbar(data[:,0], data[:,1]*data[:,0]**4, data[:,2]*data[:,0]**4,fmt='.k',ecolor='grey',zorder=0,label=r'$\chi^2 = $'+'{:06.3e}'.format(chi_s))
		plt.plot(Refl[:,0],Refl[:,2],color='red')
		plt.legend()
		plt.yscale('log')
		plt.xlabel(r'$Q(\AA^{-1})$')
		plt.ylabel(r'$R(Q)Q^4$')
		#plt.title(r'$R(Q)Q^4$')
		plt.grid(True)
		plt.draw() 
		plt.savefig(project+'_reflectivity_RQ^4.pdf')

		plt.figure()
		for i in range(len(system)):
			plt.plot(Profile[i][:,0],Profile[i][:,1]*1e6,color='green',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel(r'$sld(10^{-6}\AA^{-2})$',color='green')
		plt.savefig(project+'_sld_profile.pdf')

		plt.figure()
		for i in range(len(system)):
			plt.plot(Solvent[i][:,0],Solvent[i][:,1],color='blue',linestyle=plinestyle[i])
		#plt.yscale('log')
		plt.xlabel(r'$z(\AA)$')
		plt.ylabel('solvent',color='blue')
		plt.savefig(project+'_solvent_profile.pdf')

		plt.close()  


	print('Library versions used for the calculations:')
	print('numpy: '+np.__version__)
	print('scipy: '+scipy.__version__)
	print('numdifftools: '+nd.__version__)
	print('sympy: '+sympy.__version__)

	if project != 'none':
		f.write('Library versions used for the calculations:\n')
		f.write('numpy: '+np.__version__+'\n')
		f.write('scipy: '+scipy.__version__+'\n')
		f.write('numdifftools: '+nd.__version__+'\n')
		f.write('sympy: '+sympy.__version__+'\n')

	return chi_s,Refl,Profile,Solvent
