"""

PACKAGE WITH MISCELLANEOUS FUNCTIONS.

"""

import datetime
import time 	
import numpy as np


def round_time(t=None, round_to=60):

	"""
	Rounds off datetime.datetime.now() to nearest second.
	
	Parameters
	----------
	t : 		Current date and time.
	round_to : 	int, float.
				Nearest second dt is to be rounded to.

	Example
	-------
	print("%s: in my_function." %str(round_time(datetime.datetime.now(),round_to=60)
	"""
	if t == None: 
		t = datetime.datetime.now()
	seconds = (t - t.min).seconds
	rounding = (seconds+round_to/2) // round_to * round_to
	return t + datetime.timedelta(0,rounding-seconds,-t.microsecond)


def index_of_value(a,value):
	"""
	Get value of index that is closest to value in the array/list a.
	"""

	return min(range(len(a)),key=lambda i: abs(a[i] - value))


def write_var(par,var,var_extra_title_str):
	
	"""Writing log to file outfile"""
	
	filename_ = var_extra_title_str + par.filename 

	np.save(filename_,var)	# file extension automatically .npy
	print("%s: saved datalog file %s" %(round_time(datetime.datetime.now(),round_to=60),(filename_ + ".npy")))


def hopping_lw(site0, site1):
	"""
	Hopping linewidth between site0 and site1. 

	Uses the fact that to site naming convention, sites above/below each other and next to each other must have the y coupling. Else, they have the x coupling.
	"""
	
	return 0.04 if np.abs(site0-site1)==1 else 0.1


def check_PHS(sys_):
	par.Ez = 0.025
	s = kwant.solvers.default.smatrix(sys_,energy=0.105,args=[par])
	s_ee = s.submatrix(0,0)#(0,0),(0,0))
	s_hh = s.submatrix(1,0)#(0,1),(0,1))
	print('s_ee: \n', np.round(s_ee, 3))	# 3: to 3 decimals
	print('s_hh: \n', np.round(s_hh[::-1, ::-1], 3))
	print('s_ee - s_hh^*: \n',np.round(s_ee - s_hh[::-1, ::-1].conj(), 3), '\n')


def set_filename_skeleton(par, ppar):
	"""
	Function setting a general filename to be appended to the name of the specific quantities being saved.

	Parameters
	----------
	par : 	Simplenamespace obect.
			Contains values of the physical parameters of the system.

	ppar : 	Simplenamespace object.
			Contains program specific parameters that may be useful to have in filename.
			 - Example: "Linear" scale being used in plotting conductances or "seis" being the colormap used in the plot.
			 - Example: resolutions being used for the variables biasenergies, Ez and mu; k being the number of lowest eigenvalues and states the diagonalization is performed for; n0 and n0' being the indices of the eigenstates and energies that have the lowest absolute value, which is used to calculate the Cooper charge and E1 and Em1 corresponding to these specified states.

	"""


	if par.pincher == True:
		if ppar.doubleres == False:
			ppar.filename = str.format("Nx%g_Ny%g_LM%g_LL%g_LR%g_LY%g_ax%g_ay%g_mu%g_Gamma%g_pL%g_pR%g_k%g_biE%g_%g_Ez%g_%g_mu%g_%g_biEres%g_Ezres%g_mures%g_pinchcoef%g_tsys%g_" %(par.Nx,par.Ny,par.LM,par.LL,par.LR,par.LY,par.ax,par.ay,par.mu,par.Gamma,par.phiL,par.phiR,par.k,par.biasenergies[0],par.biasenergies[-1],par.Ez_values[0],par.Ez_values[-1],par.mu_values[0],par.mu_values[-1],par.biasEres,par.Ezres,par.mures,par.t_p_coeff,par.ty))	# name is for both fig and datalog files and specifies only the parameters such that can be used for any general file if needed; adding file extension later when saving the respective files
		if ppar.doubleres == True:
			# make equivalent par.filename_2 with the res_2 values. This is used when saving datalog files for the dependent variables that are a function of the independent variables whose resolutions are given by the second resolution.
			ppar.filename_2 = str.format("Nx%g_Ny%g_LM%g_LL%g_LR%g_LY%g_ax%g_ay%g_mu%g_Gamma%g_pL%g_pR%g_k%g_biE%g_%g_Ez%g_%g_mu%g_%g_biEres%g_Ezres%g_mures%g_pinchcoef%g_tsys%g_" %(par.Nx,par.Ny,par.LM,par.LL,par.LR,par.LY,par.ax,par.ay,par.mu,par.Gamma,par.phiL,par.phiR,par.k,par.biasenergies_2[0],par.biasenergies_2[-1],par.Ez_values_2[0],par.Ez_values_2[-1],par.mu_values_2[0],par.mu_values_2[-1],par.biasEres_2,par.Ezres_2,par.mures_2,par.t_p_coeff,par.ty))
	else:
		if ppar.doubleres == False:
			ppar.filename = str.format("Nx%g_Ny%g_LM%g_LL%g_LR%g_LY%g_ax%g_ay%g_mu%g_Gamma%g_pL%g_pR%g_k%g_biE%g_%g_Ez%g_%g_mu%g_%g_biEres%g_Ezres%g_mures%g_ty%g_" %(par.Nx,par.Ny,par.LM,par.LL,par.LR,par.LY,par.ax,par.ay,par.mu,par.Gamma,par.phiL,par.phiR,par.k,par.biasenergies[0],par.biasenergies[-1],par.Ez_values[0],par.Ez_values[-1],par.mu_values[0],par.mu_values[-1],par.biasEres,par.Ezres,par.mures,par.ty))
		if ppar.doubleres == True:
			ppar.filename = str.format("Nx%g_Ny%g_LM%g_LL%g_LR%g_LY%g_ax%g_ay%g_mu%g_Gamma%g_pL%g_pR%g_k%g_biE%g_%g_Ez%g_%g_mu%g_%g_biEres%g_Ezres%g_mures%g_ty%g_" %(par.Nx,par.Ny,par.LM,par.LL,par.LR,par.LY,par.ax,par.ay,par.mu,par.Gamma,par.phiL,par.phiR,par.k,par.biasenergies_2[0],par.biasenergies_2[-1],par.Ez_values_2[0],par.Ez_values_2[-1],par.mu_values_2[0],par.mu_values_2[-1],par.biasEres_2,par.Ezres_2,par.mures_2,par.ty))

	if ppar.make_N_wire:
		ppar.filename = "SM_wire_" + ppar.filename 

	if ppar.SClead:
		if par.tx_SC_lead and par.tx_SC_pincher:
			ppar.filename = "txySClead_%.1f_%.1f_txySCp_%.1f_%.1f_tNlead_%.1fmeV_gSgN_%g_" %(par.tx_SC_lead,par.ty_SC_lead,par.tx_SC_pincher,par.ty_SC_pincher,par.t_N_lead,par.g_S_over_g_N) + ppar.filename
		else:
			ValueError("To make ppar.filename skeleton, need both SC lead and pincher t-coefficients when ppar.SClead = True.")
	if ppar.make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight:
		ppar.filename = "ToCaliSys_"+ ppar.filename 	# the system to be calibrated, N-N(S)
	elif ppar.make_1D_NLeft_1D_S_Heff_No_NRight:
		ppar.filename = "CaliSys_" + ppar.filename 		# the system to calibrate against, N-S (Heff in S region)
	else:
		pass
	# if len(par.Ez_values) == 1:
	# 	ppar.filename = ppar.filename + "Ezconst_%g_" %par.Ez
	# 	if ppar.doubleres == True:
	# 		ppar.filename_2 = ppar.filename_2 + "Ezconst_%g_" %par.Ez_2

	# elif len(par.mu_values) == 1:
	# 	ppar.filename = ppar.filename + "muconst_%g_" %par.mu
	# 	if ppar.doubleres == True:
	# 		ppar.filename_2 = ppar.filename_2 + "muconst_%g_" %par.mu_2
	
	if ppar.doubleres == False:
		ppar.filename_11 = "G_11_%s_"%(ppar.var_name) + ppar.filename+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_12 = "G_12_%s_"%(ppar.var_name) + ppar.filename+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_11_S = "G_11_S_%s_"%(ppar.var_name) + ppar.filename+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_11_A = "G_11_A_%s_"%(ppar.var_name) + ppar.filename+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_12_S = "G_12_S_%s_"%(ppar.var_name) + ppar.filename+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_12_A = "G_12_A_%s_"%(ppar.var_name) + ppar.filename+"%s_%g"%(ppar.par_const_name,ppar.par_const)
	
	elif ppar.doubleres == True:
		ppar.filename_11 = "G_11_%s_"%(ppar.var_name) + ppar.filename_2+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_12 = "G_12_%s_"%(ppar.var_name) + ppar.filename_2+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_11_S = "G_11_S_%s_"%(ppar.var_name) + ppar.filename_2+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_11_A = "G_11_A_%s_"%(ppar.var_name) + ppar.filename_2+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_12_S = "G_12_S_%s_"%(ppar.var_name) + ppar.filename_2+"%s_%g"%(ppar.par_const_name,ppar.par_const)
		ppar.filename_12_A = "G_12_A_%s_"%(ppar.var_name) + ppar.filename_2+"%s_%g"%(ppar.par_const_name,ppar.par_const)


def addto_filename(ppar, *args):
	"""

	Parameters
	----------
	ppar: 				Simplenamespace Object ppar that contains the current ppar.filename in use.

	args : 				dict of strings.
						each string corresponds to what one wants to add to filename. These need to be one of the following:
						 - "n": the index/indices of the eigenstate for which one is considering.

	Notes
	-----
	 - only call this function once per parameter value you want to add to the filename. Else, will have duplicates of same parameter value in the name.
	"""

	if "n" in args:
		n = ppar.n
		ppar.filename = ppar.filename + "n0%g_n0'%g_"%(n[0],n[1])
		if ppar.doubleres == True:
			ppar.filename_2 = ppar.filename_2 + "n0%g_n0'%g_"%(n[0],n[1])
