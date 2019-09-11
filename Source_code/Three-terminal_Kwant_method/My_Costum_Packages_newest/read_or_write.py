"""

PACKAGE FOR READING AND WRITING DATA FILES

"""

from My_Costum_Packages_newest import misc, calc, plot
import scipy.sparse.linalg
import kwant
import datetime
import numpy as np


def write_var(ppar,var,var_extra_title_str):
	
	"""Writing log to file outfile"""
	
	filename_ = var_extra_title_str + ppar.filename 

	np.save(filename_,var)	# file extension automatically .npy
	# print("%s: -saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),(filename_ + ".npy")))

def write_var_2(ppar,var,var_extra_title_str):
	
	"""Writing log to file outfile with filename skeleton ppar.filename_2"""
	
	filename_ = var_extra_title_str + ppar.filename_2

	np.save(filename_,var)


def read_or_write_G_11_12_S_A_var(sys,par,ppar,generating_G_11_12_S_A_var,biasenergies_asymm,var_name,filename_11,filename_12,filename_11_S,filename_11_A,filename_12_S,filename_12_A):

	"""	Read or write conductances for a general variable named var_name.
		var_name: string of variable name. The variable has to be exactly the name of a parameter used when constructing the Hamiltonian in the Builder object sys.
		Suggested filename format: G_11_var_%g_%g_"%(par.mu_values[0],par.mu_values[-1]) + ppar.filename
	"""

	import time
	# from My_Costum_Packages import calc
	# import kwant

	print("%s: in read_or_write_G_11_12_S_A_var, about to calc/get conductance(s) depending on %s" %(misc.round_time(datetime.datetime.now(),round_to=60),var_name))

	exec("par.var_values = par.%s_values" %var_name)

	G_11_var = np.zeros((len(par.var_values),len(par.biasenergies)))
	G_12_var = np.zeros((len(par.var_values),len(par.biasenergies)))
	G_11_A_var = np.zeros((len(par.var_values),len(par.biasenergies)))
	G_11_S_var = np.zeros((len(par.var_values),len(par.biasenergies)))
	G_12_A_var = np.zeros((len(par.var_values),len(par.biasenergies)))
	G_12_S_var = np.zeros((len(par.var_values),len(par.biasenergies)))

	if generating_G_11_12_S_A_var[0] == True:
		print(" - Generating conductances")
		if biasenergies_asymm == True:
			print(" - User input biasenergies_asymm=True.")
			if np.round(par.biasenergies[0],10) == -np.round(par.biasenergies[-1],10):
				print(" - Check: This is the case for first and last element.")
			else:
				Inputerror(" - First element of p.biasenergies was not the negative of the last value. Ensure p.biasenergies has an odd number of entries and is antisymmetric about zero.")
		start_time = time.time()
		if ppar.oneNLead == True:	# one lead means G_12 not generated, so array axes won't match when executing the string in the 'else' block below.
			exec_str = "for par.%s, i in zip(par.%s_values,range(len(par.%s_values))):"%(var_name,var_name,var_name) + \
						"\n\telapsed_time = time.time() - start_time" + \
						"\n\t[G_11_var[i,:],G_12_var,G_11_S_var[i,:],G_11_A_var[i,:],G_12_S_var,G_12_A_var] = calc.calc_G_11_12_S_A(sys,par,ppar,i,var_name,elapsed_time,biasenergies_asymm=%s)" %(biasenergies_asymm)

		else:
			exec_str = "for par.%s, i in zip(par.%s_values,range(len(par.%s_values))):"%(var_name,var_name,var_name) + \
						"\n\telapsed_time = time.time() - start_time" + \
						"\n\t[G_11_var[i,:],G_12_var[i,:],G_11_S_var[i,:],G_11_A_var[i,:],G_12_S_var[i,:],G_12_A_var[i,:]] = calc.calc_G_11_12_S_A(sys,par,ppar,i,var_name,elapsed_time,biasenergies_asymm=%s)" %(biasenergies_asymm)

		exec(exec_str)
		np.save(filename_11,G_11_var)
		np.save(filename_12,G_12_var)
		np.save(filename_11_S,G_11_S_var)
		np.save(filename_11_A,G_11_A_var)
		np.save(filename_12_S,G_12_S_var)
		np.save(filename_12_A,G_12_A_var)

	if generating_G_11_12_S_A_var[0] == False:	
		print(" - Reading conductances")
		if filename_12 == None:
			G_11_var = np.load(filename_11 + '.npy' )
		else:
			G_11_var = np.load(filename_11 + '.npy' )
			G_12_var = np.load(filename_12 + '.npy' )
			G_11_S_var = np.load(filename_11_S + '.npy')
			G_11_A_var = np.load(filename_11_A + '.npy')
			G_12_S_var = np.load(filename_12_S + '.npy')
			G_12_A_var = np.load(filename_12_A + '.npy')

	return G_11_var,G_12_var,G_11_S_var,G_11_A_var,G_12_S_var,G_12_A_var




def read_or_write_energies_Ez(sys,par,ppar,*args):
	"""
	Reads from file or calculates eigenenergies and eigenvectors as a function of the parallel applied magnetic field par.Ez.
	Parameters
	----------
	args : 	if want to specify other filename than EvsEz_ + ppar.filename.

	Notes
	-----
	 - Using arr.argsort() to get indices that sorts the array arr. Here, using this on the eigenenergies in order to also be able to sort the eigenvectors, such that, later, one may access the k'th eigenvector.
	"""
	if ppar.generate_En_Ez == False:	# NB: need to have generated file with appropriate filename before setting use_energies_from_file = True.
		par.energies_Ez = np.load("EvsEz_" + ppar.filename + ".npy")	# .npy file faster to load with numpy
		par.Ez=par.Ez_values[0]	# just need this one defined in order for sys to work later
		if ppar.Eonly == True:
			pass
		else:
			par.eigenvectors_Ez = np.load("EigenvectorsvsEz_" + ppar.filename + ".npy")


	else:
		if par.printtype == True:	# if you want to print parts of the Hamiltonian. Only goes through the first iteration, when counter_Hprint = 0.
			Hprint(counter_Hprint)	
		"""Calculating energies for a range of applied fields, Ez"""
		par.energies_Ez=[]
		par.eigenvectors_Ez = []
		for par.Ez,i in zip(par.Ez_values,range(len(par.Ez_values))): # par is the class '__main__.SimpleNamespace'. par.Ez in Ez_values adds Ez=Ez_values[0] etc. to the par dictionary. In the loop, this happens one by one, so B is written over for every iteration of the loop.
			
			print(" - Ez-value number:", i)	

			H = sys.hamiltonian_submatrix(args=[par],sparse=True)
			H = H.tocsc()
			eigs = scipy.sparse.linalg.eigsh(H,k=par.k,sigma=0)

			sorting_indices = eigs[0].argsort()
			evals = eigs[0][sorting_indices]
			evecs = eigs[1][:,sorting_indices]
			
			par.energies_Ez.append(evals)
			par.eigenvectors_Ez.append(evecs)

		if ppar.doubleres== False:
			write_var(ppar,var=par.energies_Ez,var_extra_title_str="EvsEz_")
			print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("EvsEz_" + ppar.filename + ".npy")))
			write_var(ppar,var=par.eigenvectors_Ez,var_extra_title_str="EigenvectorsvsEz_")
			print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("EigenvectorsvsEz_" + ppar.filename + ".npy")))		
		else:
			write_var_2(ppar,var=par.energies_Ez,var_extra_title_str="EvsEz_")
			print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("EvsEz_" + ppar.filename + ".npy")))
			write_var_2(ppar,var=par.eigenvectors_Ez,var_extra_title_str="EigenvectorsvsEz_")
			print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("EigenvectorsvsEz_" + ppar.filename + ".npy")))	


def read_or_write_energies_mu(sys,par,ppar):
	"""
	Reads from file or calculates eigenenergies and eigenvectors as a function of only mu_values, for some constant biasenergy and par.Ez.

	Notes
	-----
	 - Using arr.argsort() to get indices that sorts the array arr. Here, using this on the eigenenergies in order to also be able to sort the eigenvectors, such that, later, one may access the k'th eigenvector.
	 - If ppar.pincher_SC_lead is True: In this case, the total Kwant system is the normal metal plus the SC pincher layer that runs along the long axis of the normal metal. But we want the energy of the normal metal wire only. Therefore, when ppar.pincher_SC_lead is True, then only the subhamiltonian at the sites of the normal metal are used to compute the energies of the system.
	"""	
	print("%s: In read_or_write_energies_mu" %(misc.round_time(datetime.datetime.now(),round_to=60)))

	if ppar.generate_En_mu == False:

		if ppar.doubleres==False:
			par.energies_mu = np.load("Evsmu_" + ppar.filename + ".npy")
			print(" - loaded Evsmu datalog file")
			if ppar.Eonly == True:
				pass
			else:
				par.eigenvectors_mu = np.load("Eigenvectorsvsmu_" + ppar.filename + ".npy")
				print(" - loaded Eigenvectorsvsmu datalog file")
		else:
			par.energies_mu = np.load("Evsmu_" + ppar.filename_2 + ".npy")
			print(" - loaded Evsmu datalog file")
			if ppar.Eonly == True:
				pass
			else:
				par.eigenvectors_mu = np.load("Eigenvectorsvsmu_" + ppar.filename_2 + ".npy")
				print(" - loaded Eigenvectorsvsmu datalog file")

	else:
		calc.calc_energies_mu(sys,par,ppar)

	if ppar.doubleres == False:
		write_var(ppar,var=par.energies_mu,var_extra_title_str="Evsmu_")
		print("%s: - saved Evsmu datalog file" %(misc.round_time(datetime.datetime.now(),round_to=60)))
		write_var(ppar,var=par.eigenvectors_mu,var_extra_title_str="Eigenvectorsvsmu_")
		print("%s: - saved Eigenvectorsvsmu datalog file" %(misc.round_time(datetime.datetime.now(),round_to=60)))		
	else:
		write_var_2(ppar,var=par.energies_mu,var_extra_title_str="Evsmu_")
		print("%s: - saved Evsmu datalog file" %(misc.round_time(datetime.datetime.now(),round_to=60)))
		write_var_2(ppar,var=par.eigenvectors_mu,var_extra_title_str="Eigenvectorsvsmu_")
		print("%s: - saved Eigenvectorsvsmu datalog file" %(misc.round_time(datetime.datetime.now(),round_to=60)))	

def read_or_write_energies_biasenergies(sys,par,ppar):
	"""
	Reads from file or calculates eigenenergies and eigenvectors as a function of only biasenergies, for some constant par.mu and par.Ez.
	"""
	if ppar.use_energies_from_file_biasenergies == True:
		par.energies_biasenergies = np.load("Evsmu_" + ppar.filename + ".npy")
		par.eigenvectors_biasenergies = np.load("Eigenvectorsvsbiasenergies_" + ppar.filename + ".npy")
		par.Ez = par.Ez_values[0]

	else:
		par.energies_biasenergies = []
		par.eigenvectors_biasenergies = []#np.zeros((len(par.mu_values),par.k))
		for par.biasenergy, i in zip(par.biasenergies,range(len(par.biasenergies))):
			print(i,par.biasenergy)

			H = sys.hamiltonian_submatrix(args=[par],sparse=True)
			H = H.tocsc()
			Hpbias = H + np.array(par.biasenergy*np.ones(np.shape(H)))

			eigs = scipy.sparse.linalg.eigsh(Hpbias,k=par.k,sigma=0)

			sorting_indices = eigs[0].argsort()
			evals = eigs[0][sorting_indices]
			evecs = eigs[1][:,sorting_indices]
			
			par.energies_mu.append(evals)
			par.eigenvectors_mu.append(evecs)

		write_var(ppar,var=par.energies_biasenergies,var_extra_title_str="Evsbiasenergies_")
		print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("Evsbiasenergies_" + ppar.filename + ".npy")))
		write_var(ppar,var=par.eigenvectors_biasenergies,var_extra_title_str="Eigenvectorsvsbiasenergies_")
		print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("Eigenvectorsvsbiasenergies_" + ppar.filename + ".npy")))	


def read_or_write_energies_mu_bias(sys,par,ppar):

	"""
	Reads from file or calculated eigenenergies and eigenvectors as a function of mu and biasenergy for the system.

	Parameters
	----------
	sys : 	kwant.Builder.
	par : 	parameters of system.
	generate_En_mu_bias : 	Boolean (True/False)
	"""

	print(misc.round_time(datetime.datetime.now(),round_to=60),": In read_or_write_energies_mu_bias")
	if ppar.generate_En_mu_bias == False:
		par.energies_mu_bias = np.load("Evsmubias_" + ppar.filename + ".npy")
		par.eigenvectors_mu_bias = np.load("Eigenvectorsvsmubias_" + ppar.filename + ".npy")
		par.Ez = par.Ez_values[0]

	else:
		par.energies_mu_bias = []
		par.eigenvectors_mu_bias = np.zeros((len(par.mu_values),len(par.biasenergies),par.Ny*2*2,par.k),dtype=np.complex_)
		for par.mu, i in zip(par.mu_values,range(len(par.mu_values))):
			print(" - %.0f'th mu-value" %i)
			for par.biasenergy,j in zip(par.biasenergies,range(len(par.biasenergies))):
				H = sys.hamiltonian_submatrix(args=[par],sparse=True)
				H = H.tocsc()
				Hpbias = H + np.array(par.biasenergy*np.ones(np.shape(H)))

				"""Trying: just not sorting the energies:"""
				eigs = scipy.sparse.linalg.eigsh(Hpbias,k=par.k,sigma=0)
				# eigs = np.sort(np.array(eigs),axis=0).all()
				par.energies_mu_bias.append(np.sort(eigs[0]))
				par.eigenvectors_mu_bias[i,j,:,:] = eigs[1]

		write_var(ppar,var=par.energies_mu_bias,var_extra_title_str="Evsmubias_")
		print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("Evsmubias_" + ppar.filename + ".npy")))
		write_var(ppar,var=par.eigenvectors_mu_bias,var_extra_title_str="Eigenvectorsvsmubias_")
		print("%s: saved datalog file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),("Eigenvectorsvsmubias_" + ppar.filename + ".npy")))	


def read_or_write_g012a_mu(sys,par,ppar):
	if ppar.generating_g012a == None:
		# if not doing anything with g0_12_a
		pass
	else:
		if ppar.doubleres == True:
			mu_values_ = par.mu_values_2
			biasenergies_ = par.biasenergies_2
			filename_ = ppar.filename_2
			var_values_str_ = "mu_values_2"
		else:
			mu_values_ = par.mu_values
			biasenergies_ = par.biasenergies
			filename_ = ppar.filename
			var_values_str_ = "mu_values"

		if ppar.generating_g012a[0] == True:
			g0_12_a = calc.calc_g0_12_asymm(sys, par, ppar, lambda_=mu_values_, t_L=par.t_pincher, t_R=par.t_pincher, omega=biasenergies_, En=par.energies_mu[:,ppar.n[0]], evec_n=par.eigenvectors_mu[:,:,ppar.n[0]],filename=ppar.filename)
			if ppar.generating_g012a[1] == True and len(ppar.generating_g012a) == 2:
				# plotting conductance only
				plot.plot_map(np.transpose(g0_12_a),mu_values_,biasenergies_,"g0_12_asym_"+filename_+".pdf","$g_{12,a}$","$\mu\ [meV]$","$E_{bias}\ [meV]$","seismic")
			elif ppar.generating_g012a[1] == True and len(ppar.generating_g012a) == 3:
				# plotting conductance together with energies and Cooper charges of nth and n'th eigenstates
				legend = ["$E_1\, [meV]$","$E_{-1}\, [meV]$","$(u^2_1-v^2_1)$","$(u^2_{-1}-v^2_{-1})$"]

				plot.plot_G_ij_var(par,scaling="Linear",G_ij=np.transpose(g0_12_a),var_values_str=var_values_str_,filename="g012a_nus_vsEbimu_seis_"+filename_,figtitle="$g^0_{12,a}\ [e^2/h]$",xlabel="$\mu\ [meV]$",ylabel="$E_{bias}\ [meV]$",cm="seismic",u2mv2_factor=np.abs(np.max(par.E1_mu))/np.abs(np.max(par.u2mv2_E1_mu)),var_values_=mu_values_, E1_var=par.E1_mu, Em1_var=par.Em1_mu, u2mv2_E1_var=par.u2mv2_E1_mu, u2mv2_Em1_var=par.u2mv2_Em1_mu,legend=legend)

		elif ppar.generating_g012a[0] == False:
			g0_12_a = np.load("g0_12_asym_"+ppar.filename+".npy")
			if ppar.generating_g012a[1] == True:
				# if plotting
				plot.plot_map(np.transpose(g0_12_a),par.mu_values,par.biasenergies,"g0_12_asym_"+ppar.filename+".pdf","$g_{12,a}$","$\mu\ [meV]$","$E_{bias}\ [meV]$","seismic")


def read_or_write_g012a_Ez(sys,par,ppar):
	if ppar.generating_g012a == None:
		# if not doing anything with g0_12_a
		pass
	else:
		if ppar.generating_g012a[0] == True:
			g0_12_a = calc.calc_g0_12_asymm(sys, par, ppar, lambda_=par.Ez_values, t_L=par.t_pincher, t_R=par.t_pincher, omega=par.biasenergies, En=par.energies_Ez[:,ppar.n[0]], evec_n=par.eigenvectors_Ez[:,:,ppar.n[0]],filename=ppar.filename)
			if ppar.generating_g012a[1] == True and len(ppar.generating_g012a) == 2:
				# plotting conductance only
				plot.plot_map(np.transpose(g0_12_a),par.Ez_values,par.biasenergies,"g0_12_asym_"+ppar.filename+".pdf","$g_{12,a}$","$E_z\ [meV]$","$E_{bias}\ [meV]$","seismic")
			elif ppar.generating_g012a[1] == True and len(ppar.generating_g012a) == 3:
				# plotting conductance together with energies and Cooper charges of nth and n'th eigenstates
				legend = ["$E_1\, [meV]$","$E_{-1}\, [meV]$","$(u^2_1-v^2_1)$","$(u^2_{-1}-v^2_{-1})$"]

				plot.plot_G_ij_var(par,scaling="Linear",G_ij=np.transpose(g0_12_a),var_values_str="Ez_values",filename="g012a_nus_vsEbiEz_seis_"+ppar.filename,figtitle="$g^0_{12,a}\ [e^2/h]$",xlabel="$E_z\ [meV]$",ylabel="$E_{bias}$",cm="seismic",u2mv2_factor=np.abs(np.max(par.E1_Ez))/np.abs(np.max(par.u2mv2_E1_Ez)),var_values_=par.Ez_values, E1_var=par.E1_Ez, Em1_var=par.Em1_Ez, u2mv2_E1_var=par.u2mv2_E1_Ez, u2mv2_Em1_var=par.u2mv2_Em1_Ez,legend=legend)

		elif ppar.generating_g012a[0] == False:
			g0_12_a = np.load("g0_12_asym_"+ppar.filename+".npy")
			if ppar.generating_g012a[1] == True:
				# if plotting
				plot.plot_map(np.transpose(g0_12_a),par.Ez_values,par.biasenergies,"g0_12_asym_"+ppar.filename+".pdf","$g_{12,a}$","$Ez\ [meV]$","$E_{bias}\ [meV]$","seismic")


def read_or_write_u2mv2_u2_v2(par,ppar):
	if ppar.doubleres == True:
		filename_ = ppar.filename_2
	else:
		filename_ = ppar.filename

	if ppar.generate_E1m1_mu == True:
		if ppar.doubleres == False:
			par.u2mv2_1_mu, par.u2_1_mu, par.v2_1_mu, par.u2_mu_sites_1, par.v2_mu_sites_1 = calc.calc_u2mv2_mu(par.evec1_mu, par.mu_values, par.Ny, filename_ = "1_"+filename_+".npy", filename="u2mv2_E1_vs_mu_"+ppar.filename+".npy")
			par.u2mv2_m1_mu, par.u2_m1_mu, par.v2_m1_mu, par.u2_mu_sites_m1, par.v2_mu_sites_m1 = calc.calc_u2mv2_mu(par.evecm1_mu, par.mu_values, par.Ny, filename_ = "m1_"+filename_+".npy", filename="u2mv2_Em1_vs_mu_"+ppar.filename+".npy")

			np.save("u2_E1_vs_mu_" + ppar.filename + ".npy",par.u2_1_mu)
			np.save("v2_E1_vs_mu_" + ppar.filename + ".npy",par.v2_1_mu)

			np.save("u2_Em1_vs_mu_" + ppar.filename + ".npy",par.u2_m1_mu)
			np.save("v2_Em1_vs_mu_" + ppar.filename + ".npy",par.v2_m1_mu)
			print(" - u2mv2 etc. vs. mu files saved (same res as conductances)")
		else:
			par.u2mv2_1_mu, par.u2_1_mu, par.v2_1_mu, par.u2_mu_sites_1, par.v2_mu_sites_1  = calc.calc_u2mv2_mu(par.evec1_mu, par.mu_values_2, par.Ny, filename_ = "1_"+filename_+".npy", filename="u2mv2_E1_vs_mu_"+ppar.filename_2+".npy")
			par.u2mv2_m1_mu, par.u2_m1_mu, par.v2_m1_mu, par.u2_mu_sites_m1, par.v2_mu_sites_m1 = calc.calc_u2mv2_mu(par.evecm1_mu, par.mu_values_2, par.Ny, filename_ = "m1_"+filename_+".npy", filename="u2mv2_Em1_vs_mu_"+ppar.filename_2+".npy")

			np.save("u2_E1_vs_mu_" + ppar.filename_2 + ".npy",par.u2_1_mu)
			np.save("v2_E1_vs_mu_" + ppar.filename_2 + ".npy",par.v2_1_mu)

			np.save("u2_Em1_vs_mu_" + ppar.filename_2 + ".npy",par.u2_m1_mu)
			np.save("v2_Em1_vs_mu_" + ppar.filename_2 + ".npy",par.v2_m1_mu)
			print(" - u2mv2 etc. vs. mu files saved (different res from conductances)")

	elif ppar.generate_E1m1_mu == False:
		if ppar.doubleres == False:
			par.u2mv2_E1_mu = np.load("u2mv2_E1_vs_mu_"+ppar.filename+".npy")
			par.u2mv2_Em1_mu = np.load("u2mv2_Em1_vs_mu_"+ppar.filename+".npy")

			par.u2_E1_mu = np.load("u2_E1_vs_mu_" + ppar.filename + ".npy")
			par.v2_E1_mu = np.load("v2_E1_vs_mu_" + ppar.filename + ".npy")

			par.u2_Em1_mu = np.load("u2_Em1_vs_mu_" + ppar.filename + ".npy")
			par.v2_Em1_mu = np.load("v2_Em1_vs_mu_" + ppar.filename + ".npy")

			print(" - u2mv2 etc. vs. mu files loaded (same res as conductances)")
		else:
			par.u2mv2_E1_mu = np.load("u2mv2_E1_vs_mu_"+ppar.filename_2+".npy")
			par.u2mv2_Em1_mu = np.load("u2mv2_Em1_vs_mu_"+ppar.filename_2+".npy")

			par.u2_E1_mu = np.load("u2_E1_vs_mu_" + ppar.filename_2 + ".npy")
			par.v2_E1_mu = np.load("v2_E1_vs_mu_" + ppar.filename_2 + ".npy")

			par.u2_Em1_mu = np.load("u2_Em1_vs_mu_" + ppar.filename_2 + ".npy")
			par.v2_Em1_mu = np.load("v2_Em1_vs_mu_" + ppar.filename_2 + ".npy")

			print(" - u2mv2 etc. vs. mu files loaded (different res from conductances)")

		# Note: filename_ already defined depending on whether ppar.doubleres is True or False	
		par.u2_mu_sites_1 = np.load("u2_mu_sites_1_"+filename_+".npy")
		par.v2_mu_sites_1 = np.load("v2_mu_sites_1_"+filename_+".npy")

		par.u2_mu_sites_m1 = np.load("u2_mu_sites_m1_"+filename_+".npy")
		par.v2_mu_sites_m1 = np.load("v2_mu_sites_m1_"+filename_+".npy")

	elif ppar.generate_E1m1_mu == None:
		pass

	if ppar.generate_E1m1_Ez == True and hasattr(par,"evec1_Ez"):
		par.u2mv2_1_Ez, par.u2_1_Ez, par.v2_1_Ez = calc.calc_u2mv2_var(par.evec1_Ez, var_values=par.Ez_values, Ny=par.Ny, filename="E1_vs_Ez_"+ppar.filename+".npy", savefiles=True)
		par.u2mv2_m1_Ez, par.u2_m1_Ez, par.v2_m1_Ez = calc.calc_u2mv2_var(par.evecm1_Ez, var_values=par.Ez_values, Ny=par.Ny, filename="Em1_vs_Ez_"+ppar.filename+".npy", savefiles=True)

		if ppar.doubleres == False:
			np.save("u2_E1_vs_Ez_" + ppar.filename + ".npy",par.u2_1_Ez)
			np.save("v2_E1_vs_Ez_" + ppar.filename + ".npy",par.v2_1_Ez)

			np.save("u2_Em1_vs_Ez_" + ppar.filename + ".npy",par.u2_m1_Ez)
			np.save("v2_Em1_vs_Ez_" + ppar.filename + ".npy",par.v2_m1_Ez)
			print(" - u2mv2 etc. vs. Ez files saved (same res as conductances)")
		else:
			np.save("u2_E1_vs_Ez_" + ppar.filename_2 + ".npy",par.u2_1_Ez)
			np.save("v2_E1_vs_Ez_" + ppar.filename_2 + ".npy",par.v2_1_Ez)

			np.save("u2_Em1_vs_Ez_" + ppar.filename_2 + ".npy",par.u2_m1_Ez)
			np.save("v2_Em1_vs_Ez_" + ppar.filename_2 + ".npy",par.v2_m1_Ez)
			print(" - u2mv2 etc. vs. Ez files saved (different res from conductances)")

	elif ppar.generate_E1m1_Ez == False:
		if ppar.doubleres == False:
			par.u2mv2_E1_Ez = np.load("u2mv2_E1_vs_Ez_"+ppar.filename+".npy")
			par.u2mv2_Em1_Ez = np.load("u2mv2_Em1_vs_Ez_"+ppar.filename+".npy")

			par.u2_E1_Ez = np.load("u2_E1_vs_Ez_" + ppar.filename + ".npy")
			par.v2_E1_Ez = np.load("v2_E1_vs_Ez_" + ppar.filename + ".npy")

			par.u2_Em1_Ez = np.load("u2_Em1_vs_Ez_" + ppar.filename + ".npy")
			par.v2_Em1_Ez = np.load("v2_Em1_vs_Ez_" + ppar.filename + ".npy")
			print(" - u2mv2 etc. vs. Ez files loaded (same res as conductances)")
		else:
			par.u2mv2_E1_Ez = np.load("u2mv2_E1_vs_Ez_"+ppar.filename_2+".npy")
			par.u2mv2_Em1_Ez = np.load("u2mv2_Em1_vs_Ez_"+ppar.filename_2+".npy")

			par.u2_E1_Ez = np.load("u2_E1_vs_Ez_" + ppar.filename_2 + ".npy")
			par.v2_E1_Ez = np.load("v2_E1_vs_Ez_" + ppar.filename_2 + ".npy")

			par.u2_Em1_Ez = np.load("u2_Em1_vs_Ez_" + ppar.filename_2 + ".npy")
			par.v2_Em1_Ez = np.load("v2_Em1_vs_Ez_" + ppar.filename_2 + ".npy")
			print(" - u2mv2 etc. vs. Ez files loaded (different res from conductances)")


def read_or_write_E1m1(par,ppar):
	"""
	Function reading or writing values for E1 and Em1 for par and ppar for the desired variable(s), mu or Ez in this current implementation.

	Parameters:
	-----------
	ppar.generate_E1m1_mu: 	Boolean.

	ppar.generate_E1m1_Ez: 	Boolean.

	"""

	if ppar.doubleres == False:
		pass
	else:
		mu_values_ = par.mu_values_2
	print(" - shape of evec0prime = ", np.shape(par.eigenvectors_mu[:,:,ppar.n[1]]))
	print(" - shape of eigenvectors_mu = ", np.shape(par.eigenvectors_mu))
	print(" - shape of energies_mu = ", np.shape(par.energies_mu))

	if ppar.generate_E1m1_mu == True:

		par.E1_mu, par.Em1_mu, par.evec1_mu, par.evecm1_mu = calc.calc_E1_Em1_from_E0_E0prime_var(E0 = par.energies_mu[:,ppar.n[0]], E0prime = par.energies_mu[:,ppar.n[1]], evec0 = par.eigenvectors_mu[:,:,ppar.n[0]], evec0prime = par.eigenvectors_mu[:,:,ppar.n[1]], var_values = mu_values_,N=par.N,filename="vs_mu_"+ppar.filename,savefiles=False) # N: because 4 DOFs at each 'site'
		if ppar.doubleres == False:
			np.save("E1_vs_mu_" + ppar.filename + ".npy", par.E1_mu)
			np.save("Em1_vs_mu_" + ppar.filename + ".npy", par.Em1_mu)
			print(" - E1_vs_mu and Em1_vs_mu files saved")
		else:
			np.save("E1_vs_mu_" + ppar.filename_2 + ".npy", par.E1_mu)
			np.save("Em1_vs_mu_" + ppar.filename_2 + ".npy", par.Em1_mu)
			print(" - E1_vs_mu and Em1_vs_mu files saved")

	elif ppar.generate_E1m1_mu == False:
		if ppar.doubleres == False:
			par.E1_mu = np.load("E1_vs_mu_" + ppar.filename + ".npy")
			par.Em1_mu = np.load("Em1_vs_mu_" + ppar.filename + ".npy")
			print(" - E1_vs_mu and Em1_vs_Ez files loaded")
		else:
			par.E1_mu = np.load("E1_vs_mu_" + ppar.filename_2 + ".npy")
			par.Em1_mu = np.load("Em1_vs_mu_" + ppar.filename_2 + ".npy")
			print(" - E1_vs_mu and Em1_vs_Ez files loaded")
	elif ppar.generate_E1m1_mu == None:
		pass

	if hasattr(par, "energies_Ez") and ppar.generate_E1m1_Ez == True:
		# if ppar.generate_E1m1_Ez == True:
		############ VARIABLE EZ: #################
		"""
		Notes
		-----
		 - this runs only when energies and eigenvectors are read from file. Else, get error that need slices, not tuple to index the energy/eigenvector arrays.
		"""
		par.E1_Ez, par.Em1_Ez, par.evec1_Ez, par.evecm1_Ez = calc.calc_E1_Em1_from_E0_E0prime_var(E0 = par.energies_Ez[:][:,ppar.n[0]], E0prime = par.energies_Ez[:][:,ppar.n[1]], evec0 = par.eigenvectors_Ez[:][:,:,ppar.n[0]], evec0prime = par.eigenvectors_Ez[:][:,:,ppar.n[1]], var_values = par.Ez_values,N=4*par.Ny,filename="vs_Ez_"+ppar.filename,savefiles=True)
		print(" - E1_vs_Ez and Em1_vs_Ez files saved")

	elif ppar.generate_E1m1_Ez == False:
		par.E1_Ez = np.load("E1_vs_Ez_" + ppar.filename + ".npy")
		par.Em1_Ez = np.load("Em1_vs_Ez_" + ppar.filename + ".npy")
		print(" - E1_vs_Ez and Em1_vs_Ez files loaded")

	elif ppar.generate_E1m1_Ez is None:
		pass

	else:
		TypeError("ppar.generate_E1m1_Ez needs to be True (while energies_Ez is an attribute of par), False or None. Neither was the case.")
	
