import kwant, datetime, time
import numpy as np 
import matplotlib.pylab as plt
from My_Costum_Packages_newest import  set_params, calc, pincher_FWHM, read_or_write, misc, make, plot
from concurrent import futures  # Available from futures package.

###


"""Classes and functions:"""

class SimpleNamespace(object):
	"""Contains parameters of the problem"""
	def __init__(self,**kwargs):
		self.__dict__.update(kwargs)


###


def main_En_gG12A(ppar):

	"""
	Main function that focuses on calculating, saving/reading and plotting the following quantities for a system defined by the parameters in par and the program specific parameters from ppar:
	 - The par.k lowest |Energy eigenvalues| and the corresponding eigenvectors of the system. 
	 - The energy branches E_1 and E_m1 of the ppar.n'th states (ppar.n has to be an iterable of the lowest pair of energies that are symmetric around the lambda-axis, such that the energy branches can be correctly defined.)
	 - The Cooper charge u^2-v^2 of the ppar.n'th states (ppar.n is an iterable of integers or may in general be changed after determining E_1 and E_m1 to specify only one integer).
	 - The non-local asymmetric conductance G_12A of the tight-binding (TB) model from Kwant for the parameters in par.
	 - The non-local asymmetric conductance g_12A of the analytical model for the same parameters (parameters in par).

	Notes
	-----
	 - filename: REMEMBER to update the filename if new parameters are introduced or similarly! If not, writing files will overwrite old ones AND reading files, you might read the wrong file and there will be a lot of mess in the system.

	"""
	ppar.doubleres = False
	misc.set_filename_skeleton(par,ppar)

	"""Making system:"""
	sys = make.make_system(par,ppar)
	sys = sys.finalized()
	print("%s: sys finalized" %misc.round_time(datetime.datetime.now(),round_to=60))

	""" - Plotting system graph:"""
	import matplotlib.pylab as plt
	from kwant import wraparound

	plt.rc('text', usetex=True)
	font_size = 10
	plt.rcParams.update({'font.size': font_size})
	ax = kwant.plot(sys,site_symbol='o',hop_color=['orange'],site_edgecolor='red',site_lw=0.01, hop_lw=misc.hopping_lw)
	ax.savefig("N-SC_sys.pdf")


	"""	Calculating and loading, as a function of mu or Ez
		 - eigenenergies
		 - eigenvectors (if ppar.Eonly = False)
	"""
	ppar.Eonly = False
	if ppar.generate_En_mu != None:
		if ppar.generate_En_mu == False:
			read_or_write.read_or_write_energies_mu(sys,par,ppar)
		elif ppar.generate_En_mu == True:
			for ppar.generate_En_mu in [True,False]:
				read_or_write.read_or_write_energies_mu(sys,par,ppar)
		plt.figure()
		plt.tight_layout()
		plt.plot(par.mu_values,par.energies_mu)
		plt.xlabel("$\mu\ [meV]$")
		plt.ylabel("$E_n\ [meV]$")
		# plt.ylim([-5,5])
		# plt.savefig("energies_mu_ylimm55_%s.pdf"%ppar.filename)
		plt.savefig("energies_mu_%s.pdf"%ppar.filename)
	# adfadsf
	if ppar.generate_En_Ez != None:
		if ppar.generate_En_Ez == False:
			read_or_write.read_or_write_energies_Ez(sys,par,ppar)
		elif ppar.generate_En_Ez == True:
			for ppar.generate_En_Ez in [True,False]:
				read_or_write.read_or_write_energies_Ez(sys,par,ppar)

	# plt.plot(par.Ez_values,par.energies_Ez)


	"""	Calculating/loading 
		 - g0_12_a : analytical conductance inserted for TB-model parameters;
		 - E1, Em1, evec1, evecm1;
		 - u, v, u^2 - v^2 for the n'th pair of eigenvectors.
	"""
	# misc.addto_filename(ppar,"n")

	if ppar.var_name == "mu":
		read_or_write.read_or_write_g012a_mu(sys,par,ppar)
	elif ppar.var_name == "Ez":
		read_or_write.read_or_write_g012a_Ez(sys,par,ppar)

	# read_or_write.read_or_write_E1m1(par,ppar)

	# read_or_write.read_or_write_u2mv2_u2_v2(par,ppar)

	"""	Computing/reading 
		 - G_ij(_S/A) : 	conductances as fn of some variable for some parameter, versus biasenergy:"""

	exec("par.%s = %g" %(ppar.par_const_name,ppar.par_const))

	exec_str = "par.G_11_%s,par.G_12_%s,par.G_11_S_%s,par.G_11_A_%s,par.G_12_S_%s,par.G_12_A_%s = "%(ppar.var_name,ppar.var_name,ppar.var_name,ppar.var_name,ppar.var_name,ppar.var_name) +\
				"read_or_write.read_or_write_G_11_12_S_A_var(sys,par,ppar,ppar.generating_G_11_12_S_A_var,ppar.generating_G_11_12_S_A_var,ppar.var_name,ppar.filename_11,ppar.filename_12,ppar.filename_11_S,ppar.filename_11_A,ppar.filename_12_S,ppar.filename_12_A)"
	exec(exec_str)


	if ppar.generating_G_11_12_S_A_var[1] == True:
		plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_12_mu,var_values_str="mu_values",filename="G_12_vsEbias1mu_"+ppar.filename,figtitle="$G_{12}\ [e^2/h]$",xlabel="$\mu\ [meV]$",ylabel="$E_{bias}$",cm="seismic")
		plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_11_mu,var_values_str="mu_values",filename="G_11_vsEbias1mu_"+ppar.filename,figtitle="$G_{11}\ [e^2/h]$",xlabel="$\mu\ [meV]$",ylabel="$E_{bias}$",cm="seismic")#,u2mv2_factor=np.abs(np.max(par.E1_mu))/np.abs(np.max(par.u2mv2_E1_mu)),var_values_=par.mu_values, E1_var=par.E1_mu, Em1_var=par.Em1_mu, u2mv2_E1_var=par.u2mv2_E1_mu, u2mv2_Em1_var=par.u2mv2_Em1_mu)
		plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_12_A_mu,var_values_str="mu_values",filename="G_12_A_vsEbias1mu_"+ppar.filename,figtitle="$G_{12}^A\ [e^2/h]$",xlabel="$\mu\ [meV]$",ylabel="$E_{bias}$",cm="seismic")


###


"""Preparing to start program (main()):"""
ppar = SimpleNamespace(	make_N_wire=True, \
						make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight=False, \
						make_1D_NLeft_1D_S_Heff_No_NRight=False, \
						# make_1D_NLeft_1D_N_L_and_R_2D_S_2D_SMiddle = True , \
							SClead=True,pincher_SC_lead=False,SC_pincher_1D=False, \
							#
							generate_En_Ez=None, \
							generate_En_mu=True, \
							generate_E1m1_mu=None, \
							generate_E1m1_Ez=None, \
							generating_G_11_12_S_A_var=[True,True], \
							generating_g012a=[True,True], \
							#
							filename = "",one_D=True,counter=0)
# ppar = SimpleNamespace(SClead=False,generate_En_Ez=True,generate_En_mu=False,generate_E1m1_mu=False,generate_E1m1_Ez=True,generating_G_11_12_S_A_var=[True,True],generating_g012a=[True,True],filename = "",one_D=True,counter=0)



const, par, ppar = set_params.set_params(ppar)
print("Nx=%s, Ny=%.0f; ax=%g, ay=%g" %(par.Nx,par.Ny,par.ax,par.ay))
# print("tSClead=%s, muSClead=%s, tSCpincher=%s, muSCpincher=%s" %(par.tx_SC_lead,par.mu_SC_lead,par.tx_SC_pincher,par.mu_SC_pincher))
print("Top region: mu<%s meV" %(np.sqrt(par.Ez_values[0]**2-par.Gamma**2)))

t0 = time.time()
par.printtype = False 	# Not printing Hamiltonian/parts of it

"""Main program:"""
main_En_gG12A(ppar) # generate_En to be sent into main() such that not running diagonalization over again when not necessary.

"""After program:"""
print("time taken to run main(): ", time.time() - t0)
ppar.counter += 1
print("Counter = ", ppar.counter)
print("Program finished %s" %str(misc.round_time(datetime.datetime.now(),round_to=60)))
