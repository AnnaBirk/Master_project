import kwant
import datetime
import time
import numpy as np
import matplotlib.pylab as plt
from My_Costum_Packages_newest import set_params, calc, pincher_FWHM, read_or_write, misc, make, plot
from concurrent import futures  # Available from futures package.

###


"""Classes and functions:"""


class SimpleNamespace(object):
    """Contains parameters of the problem"""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


###


def main_En_gG12A(ppar):
	"""
	Main function that focuses on calculating, saving/reading and plotting the following quantities 
	for a system defined by the parameters in par and the program specific parameters from ppar:
	        - The par.k lowest |Energy eigenvalues| (absolute value) and the corresponding eigenvectors 
	                 of the system. 
	        - The eigenenergies E_1 and E_m1 of the ppar.n'th states (ppar.n has to be an iterable 
	                of the lowest pair of energies that are symmetric around the lambda(chemical potential or
	                magnetic field)-axis, such that the energies can be correctly defined.)
	        - The Cooper charge u^2-v^2 of the ppar.n'th states (ppar.n is an iterable of integers or 
	                may in general be changed after determining E_1 and E_m1 to specify only one integer).
	        - The non-local antisymmetric conductance G_12A of the tight-binding (TB) model from Kwant, 
	                for the parameters in par.
	        - The non-local antisymmetric conductance g_12A of the analytical model for the same parameters 
	                (parameters in par).	
	Notes
	-----
	 - Both par and ppar are SimpleNamespace objects.
	 - filename: REMEMBER to update the filename if new parameters are introduced or similarly. 
	                If not, writing files will overwrite old ones, and reading files, you might read the wrong file.
	"""	
	
	""" Creating a filename based on the rules specified in the following function : """
	misc.set_filename_skeleton(par, ppar)	
	""" 'Building' system, using the Kwant package : 
			The way the system is created is specified in ppar (see bottom of this file).
			These parameters in ppar are input in the make_system function. This function
			will make different types of systems depending on how ppar is specified. """
	sys = make.make_system(par, ppar)
	sys = sys.finalized()
	print("%s: sys finalized" % misc.round_time(
	    datetime.datetime.now(), round_to=60))
	par.N = 4*len(sys.sites)	    # par.N is used in calc_E1_Em1_from_E0_E0prime_var
								    # in order to set correct dimension of evec1_cut and evecm1_cut.	
	# - Plotting system graph :
	import matplotlib.pylab as plt
	from kwant import wraparound	
	ax = kwant.plot(sys,
	                site_symbol='o',
	                hop_color=['orange'],
	                site_edgecolor='red',
	                site_lw=0.01,
	                hop_lw=misc.hopping_lw)
	ax.savefig("N-SC_sys.pdf")	
	"""	Calculating and loading, as a function of mu or Ez
		 - eigenenergies
		 - eigenvectors (if ppar.Eonly = False)
	"""
	ppar.Eonly = False
	if ppar.generate_En_mu != None:		# reading or writing energy data vs. chem. pot.
	    if ppar.generate_En_mu == False:
	        read_or_write.read_or_write_energies_mu(sys, par, ppar)
	    elif ppar.generate_En_mu == True:
	        for ppar.generate_En_mu in [True, False]:
				# first time in loop: generating energy data
				# second time in loop: reading the before generated energy data
	            read_or_write.read_or_write_energies_mu(sys, par, ppar)
	        ppar.generate_En_mu = True	
	    if ppar.doubleres == False:
			# plots energies on one set of mu-axes
			# (there is no second sets of data for the energy. 
			# The user may use ppar.doubleres is True, if wanting
			# a second set of data that is zoomed in on one area in
			# mu-space)
	        plt.figure()
	        plt.tight_layout()
	        plt.plot(par.mu_values, par.energies_mu)
	        plt.xlabel("$\mu\ [meV]$")
	        plt.ylabel("$E_n\ [meV]$")
	        plt.savefig("energies_mu_%s.pdf" % ppar.filename)
	        plt.plot(par.mu_values_2, par.energies_mu)	
	    # Cuts and pastes from the energies that are only positive, 
		# with those that are only negative, for the pair of two 
		# lowest energy modes, which oscillate around zero bias :
	    read_or_write.read_or_write_E1m1(par, ppar)	
	if ppar.generate_En_Ez != None:		# reading or writing energy data vs. field
	    if ppar.generate_En_Ez == False:
	        read_or_write.read_or_write_energies_Ez(sys, par, ppar)
	    elif ppar.generate_En_Ez == True:
	        for ppar.generate_En_Ez in [True, False]:
				# first time in loop: generating data
				# second time in loop: reading already generated data.
	            read_or_write.read_or_write_energies_Ez(sys, par, ppar)	
	    # cuts and pastes from the energies as above :
	    read_or_write.read_or_write_E1m1(par, ppar)	
	"""	Calculating/loading 
		 - g0_12_a : 	analytical conductance at zero temperature (see arXiv:1905.05438),
		 				inserted for TB-model parameters.
		 - E1, Em1, evec1, evecm1 : 	two lowest eigenenergies, 'cut and pasted' such that
		 								they oscillate around zero bias in the topological region.
										The corresponding eigenvectors, similarly 'cut and pasted'.
		 - u, v, u^2 - v^2 : 			Coherence factors and BCS charge for the two lowest energy
		 								eigenvectors above.
	"""
	# misc.addto_filename(ppar,"n")		## Optional	
	if ppar.var_name == "mu":
	    read_or_write.read_or_write_g012a_mu(sys, par, ppar)
	elif ppar.var_name == "Ez":
	    read_or_write.read_or_write_g012a_Ez(sys, par, ppar)	
	read_or_write.read_or_write_u2mv2_u2_v2(par, ppar)
	if ppar.generate_E1m1_mu == True:
	    ppar.generate_E1m1_mu = False
	    read_or_write.read_or_write_u2mv2_u2_v2(par, ppar)	
	""" Plotting eigenenergies if they are read or generated, for mu 
			(can do similarly for field variable, Ez) : """	
	if ppar.generate_En_mu != None:
	    if ppar.doubleres == True:
	        plot.pplot_E_vs_var(par, ppar, par.energies_mu, par.mu_values_2, 
	                            figfilename="energies_mu_%s" % ppar.filename,
	                            xlim=[par.mu_values_2[0], par.mu_values_2[-1]],
	                            ylim=[np.min(par.energies_mu),
	                            np.max(par.energies_mu)], 
								title="", 
								xlabel="$\mu\ [meV]$", ylabel="$E_n\ [meV]$", 
								add_zoom=True, 
								ylim_zoom=[-par.bias_maxmin, par.bias_maxmin], 
								yticks=[-par.bias_maxmin, 0, par.bias_maxmin], 
								file_extension=".pdf", 
								u2mv2=par.u2mv2_E1_mu, 
								u2_mu=[], v2_mu=[]
								)	
	"""	Computing/reading 
		 - G_ij(_S/A) : 	
		 	local/nonlocal conductances as functions of some variable for some set of parameters, versus biasenergy 
				- all these variables are specified in the par SimpleNamespace object 
				- how the user chooses to calculate and represent the data is specified 
				in the ppar Simple namespace object : """	
	# Specifying the constant and the independent variable in this manner 
	# is general (the user may for instance change between mu and field being constant/variable) :
	exec("par.%s = %g" % (ppar.par_const_name, ppar.par_const))	
	# Generating or reading the local and nonlocal conductances :
	exec_str = \
		"par.G_11_%s,par.G_12_%s,par.G_11_S_%s,par.G_11_A_%s,par.G_12_S_%s,par.G_12_A_%s = "\
			 % (ppar.var_name, ppar.var_name, ppar.var_name, ppar.var_name, ppar.var_name, ppar.var_name) +\
			"read_or_write.read_or_write_G_11_12_S_A_var(	\
				sys,par,ppar,								\
				ppar.generating_G_11_12_S_A_var,			\
				ppar.generating_G_11_12_S_A_var,			\
				ppar.var_name,								\
				ppar.filename_11,							\
				ppar.filename_12,							\
				ppar.filename_11_S,							\
				ppar.filename_11_A,							\
				ppar.filename_12_S,							\
				ppar.filename_12_A)"
	exec(exec_str)	
	# Plotting contour plot of the differential conductances and
	# their symmetry-decomposed versions, in chem.pot.-bias space :
	if ppar.generating_G_11_12_S_A_var[1] == True:
		for G_, G_filename_, G_title_str in zip(\
			# nonlocal conductance, local conductance, local antisymmetric conductance, local symmetric conductance,
			# nonlocal antisymmetric conductance, nonlocal symmetric conductance, respectively.
			[par.G_12_mu, par.G_11_mu, par.G_11_A_mu, par.G_11_S_mu, par.G_12_A_mu, par.G_12_S_mu],\
				["G_12_vsEbias1mu_", "G_11_vsEbias1mu_", "G_11_A_vsEbias1mu_", "G_11_S_vsEbias1mu_", "G_12_A_vsEbias1mu_", "G_12_S_vsEbias1mu_"],\
					["$G_{LR}^0\ [e^2/h]$", "$G_{LL}^0\ [e^2/h]$", \
						"$G_{LL}^{0,asym}\ [e^2/h]$","$G_{LL}^{0,sym}\ [e^2/h]$",\
							"$G_{LR}^{0,asym}\ [e^2/h]$", "$G_{LR}^{0,sym}\ [e^2/h]$"]\
							):
			plot.plot_G_ij_var(	par, scaling="SymLogNorm", 
								G_ij=G_, 
								var_values_str="mu_values", 
								filename=G_filename_+ppar.filename,
								figtitle=G_title_str, 
								xlabel="$\mu\ [meV]$", 
								ylabel="$V\ [meV]$", 
								cm="seismic", 
								colorbar_ticks=[-1e-1, -1e-2, 0, 1e-2, 1e-1])  

	""" mu_cut_index is the index in the chemical potential variable (mu), 
	where the user chooses to make a line-cut. This cut is then plotted 
	as a function of the biasenergy. """	
	G_full_absmax = np.max(
		np.abs(par.G_11_mu[ppar.mu_cut_index, :], par.G_12_mu[ppar.mu_cut_index, :]))	
	plot.plot_G_1_G_2_varconst(par, G_1=par.G_11_mu, G_2=par.G_12_S_mu+par.G_12_A_mu, index_of_varconst=ppar.mu_cut_index, filename="G_11_12_mu_0_vs_bias"+ppar.filename, figtitle=" ",
								xlabel=r"$V\ [meV]$", ylabel=r"$G_{\alpha \beta}^0(V)\ [e^2/h]$", legend=[r"$G_{LL}^0$", r"$G_{LR}^0$"], ylim=[-1.1*G_full_absmax, 1.1*G_full_absmax], colors=['--k', "C0"])	
	plot.plot_G_1_G_2_varconst(par,
								G_1=par.G_12_S_mu, G_2=par.G_12_A_mu,
								G_3=par.G_11_S_mu, G_4=par.G_11_A_mu,
								index_of_varconst=ppar.mu_cut_index, filename="G_11_12_A_S_mu_0_vs_bias", figtitle=" ", xlabel=r"$V\ [meV]$", ylabel=r"$G_{\alpha \beta}^{0,sym/asym}(V)\ [e^2/h]$",
								legend=[r"$G_{LR}^{0,sym}$", r"$G_{LR}^{0,asym}$", r"$G_{LL}^{0,sym}$", r"$G_{LL}^{0,asym}$"], ylim=[-1.1*G_full_absmax, 1.1*G_full_absmax], colors=['C1', "C0", "k", "gray"]
								)	
	#############################
	# MAY ALSO BE USED : plotting map of normalized charge in real space. Please also consider NOTE [*1].
	# par.normalization_allsites = par.u2_mu_sites_1+par.v2_mu_sites_1
	# par.u2mv2_allsites = (par.u2_mu_sites_1-par.v2_mu_sites_1)/par.normalization_allsites	
	# plot.plot_map(par.u2mv2_allsites,par.mu_values_2,np.arange(0,par.Ny),"u2mv2_mu_i_"+ppar.filename_2+".pdf","$(u^2_1-v^2_1)/(u^2_1+v^2_1)$","$\mu\ [meV]$", "$site$","seismic",(4, 4))	
	# plot.plot_map(par.u2mv2_allsites[:,:int(round(par.Ny))],par.mu_values_2[:],np.arange(0,int(round(par.Ny))),"u2mv2_mu_i_"+ppar.filename_2+"_zoom_saturated.pdf","$(u^2_1-v^2_1)/(u^2_1+v^2_1)$","$\mu\ [meV]$", "$site$","seismic",figsize_inches = (4, 4),vmin_max=(-1,1))#,add_vlines=(0.707,0.711,0.720,0.724,0.768,0.772))
	#############################

###


"""	Preparing to start program (main()):
	ppar gives a set of instructions to set_params.set_params. These instructions include
	- Category (i) :  	what type of closed system / graph should be built using Kwant
	- Category (ii) : 	what type of leads should be attached to the closed system graph
	- Category (iii) : 	what data should be 
							generated (give the entry : True), 
							read (False), or 
							ignored (None).
						When given as a list, the first element refers to the data 
						(whether it is to be read, generated or ignored), while the second element 
						refers to whether or not the data should be plotted.
	See set_params.set_params for details.
"""

ppar = SimpleNamespace(
    # Category (i) :
    make_1DNLeft_1DN_2DSMiddle_No_NRight=False,\
    make_1D_Heff_LR_leads=True, \
    make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight=False, \
    make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight=False, \
    make_1D_NLeft_1D_S_Heff_No_NRight=False, \
    # Category (ii) :
    SClead=False, pincher_SC_lead=False, SC_pincher_1D=False, \
    # Category (iii) :
    generate_En_Ez=None, \
    generate_En_mu=None, \
    generate_E1m1_mu=None, \
    generate_E1m1_Ez=None, \
    generating_G_11_12_S_A_var=[False, True], \
    generating_g012a=[None, None], \
    # add something here if you want to add to the filename 
	# generated by misc.set_filename_skeleton(par,ppar) : 
    filename="",
    # category (i) (see set_params.set_params for details) :
    one_D=True,
    counter=0,				# a parameter that is relevant if you want to print blocks 
							# of the Hamiltonian (see Hprint.py)
	mu_cut_index = 667		# the index of the array containing the chemical potential.
							# It is the index for which plot_G_1_G_2_varconst in main()
							# plots a line-cut of the conductance signal.
)


###


""" 
Setting const (simple namespace object containing the physical constants 
for the problem), par, and ppar :
"""
const, par, ppar = set_params.set_params(ppar)
# printing the tight-binding parameters and the theoretically predicted 
# topological region for the corresponding infinite wire (see Oreg-Lutchyn model) : 
#  - Number of sites in the x and y directions; lattice constant in the x and y directions : 
print("Nx=%s, Ny=%.0f; ax=%g, ay=%g" % (par.Nx, par.Ny, par.ax, par.ay))
# The topological region when the chemical potential is the variable
# is for (chem. pot)^2 < (field)^2 - (induced superconducting gap)^2 :
print("Top region: mu<%s meV" % (np.sqrt(par.Ez_values[0]**2-par.Gamma**2)))
# The starting time, used to time how long it takes to finish the main program :
t0 = time.time()
par.printtype = False 	# Not printing blocks of the Hamiltonian (with Hprint.py) when False


"""Main program : """

main_En_gG12A(ppar)


"""After program : """

print("time taken to run main(): ", time.time() - t0)
ppar.counter += 1 		# relevant parameter if Hprint.py is called (if par.printtype is not False).
print("Counter = ", ppar.counter)
print("Program finished %s" %
      str(misc.round_time(datetime.datetime.now(), round_to=60)))