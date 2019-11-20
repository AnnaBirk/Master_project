"""

This script specifies program specific parameters /
what the program will do (ppar), 
physical parameters (par), and physical constants (const).
ppar, par, and const are SimpleNamespace objects.

"""

"""
							PARAMETERS
							----------
const : 		
	Simplenamespace object containing relevant physical constants, 
	such as the speed of light, the electron mass and Planch's reduced constant.

"""

"""
par : 	SimpleNamespace object containing relevant physical parameters.

Notes
-----
 - Parameters special for conductance calculation:
	Beware that the resolution of the biasenergies affects the 
	numerical value of the conductance! 
	If the resolution is too small, one may not resolve the conductance peaks. 
	The max conudctance value saturates such that max conductance = 4 e^2/h 
	for some (high) resolution.
"""

"""
ppar : 			
	Simplenamespace object containing "program specific parameters". 
	Sorted into three categories, (i), (ii), and (iii).

	(i) The closed system type / geometry :
	--------------------------------------
	- one_D : 	
		Boolean. 
		If True (False), a one(two)-dimensional system will be built 
		in make.make_system().

	- make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight : 		
		Boolean. 
		If defined as instance of ppar and given as True, make.make_system 
		will call 
			make.make_1D_N_2D_SC_system.
		The resultant system builder object is a one-dimensional semiconducting wire 
		in proximity to a 2D superconductor, both in a magnetic field 
		(with the corresponding Zeeman energy, p.Ez), and with Rashba spin-orbit coupling. 
		A superconducitng pincher is also constructed at the barrier between the 
		1D semiconductor and the 2D superconducting lead, given that
			pincher_SC_lead
		is True.

	(ii) The leads being attached to the closed system :
	---------------------------------------------------
	- pincher_SC_lead : 	
		Boolean. 
		If true, a pincher is added in make.make_1D_N_2D_SC_system. Construction of the 
		components of the graph that are at and connected to the interface sites (the
		sites in the superconductor that connect to the semiconductor, and implements the
		pincher), is made from onsite, hoppingx, and hoppingy_2D_superconductor_pincher. 
		They are the components of the Hamiltonian, as specified in make.make_system().

	(iii) How the program is run / what to calculate, read, and ignore :
	--------------------------------------------------------------------
	- generating_g012a : 	
		None (do nothing with g012a) or length-2 or length-3 list/tuple/array.
		If not equal to None, the ...
			... first element is a boolean: 
					True means calculate g012a, 
					False means its datalog file (which should have been generated 
					previously) is to be read and saved in the par.g0_12a variable;
			... second element is a boolean referring to whether or not g0_12a is 
			to be plotted agains biasenergies and the chosen independent variable 
			(which may have been given as a string with ppar.var_name).
			... (optional) third element is an arbitrary value. If the third element 
			is given, g012a is plotted, also with energies and Cooper charges on top of it.
			This representation may be useful in initial calculations, for checking
			that the calculation works. 
	# NAMES #
	 - ppar.var_name : 				str.
									Name of the variable that is varied (e.g. "mu" or "Ez")
	 - ppar.ppar.par_const_name : 	str.
									Name of the variable that is held constant (e.g. "Ez" 
									or "mu" if the ppar.var_name is "mu" or "Ez", respectively.)
	 - ppar.par_const : 			float.
									Numerical value of the constant.

"""




import numpy as np
class SimpleNamespace(object):
    """Contains parameters of the problem"""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def set_params(ppar):
    """Physical constants:"""
	# nm = nanometers
	# meV = milli electron volts
	# c = speed of light
	# J = Joules
	# s = seconds
	# m = meters
	# kg = kilograms
    const = SimpleNamespace(
		# Speed of light :
        c = 2.99792458e17,			# [c]=nm/s
		# Mass of an electron :
        m_e = 0.5109989461e9,		# [m_e]=meV/c^2
		# Reduced Planch's constant : 
        hbar = 6.582119514e-13,  	# [hbar]=meV*s
		# Reduced Planch's constant in SI units : 
        hbar_SI = 1.0547e-34,  		# [hbar_SI]=J*s
		# Speed of light in SI units :
        c_SI = 2.99792458e8,  		# [c_SI]=m/s
		# Mass of an electron in SI units : 
        m_e_SI = 9.10938356e-31  	# [m_e_SI]=kg
    )

	# If the closed system graph is a one-dimensional system : 
    if ppar.one_D == True:
		# The global parameter class : 
        par = SimpleNamespace(
			# No regions ('middle' or 'right') with
			# sites in the x-direction : 
			LM=0, 
			LR=0, 
			# Length of the system in the y-direction :
			LY=1500.,	# [LY]=nm
			# Default value for the chemical potential.
			# May be changed further down in this script : 
			mu=0, 		# [mu]=meV
			# alpha * hbar * c : 
			# [alpha]=nm*meV
			# [hbar]= meV*s
			# [c]=nm/s
			# [alpha*hbar*c]
			# 	=nm*meV*meV*s*nm/s
			# 	=(nm*meV)^2			 
			alphahbar=1.42e-4*const.c*const.hbar,
			Gamma=180e-3, # [Gamma]=meV
			phiL=0, # [phi]=1 (unitless phase)
			phiR=0, 
			k=30	# integer number of lowest energy eigenenergies
					# to numerically diagonalize for the closed
					# system.
			) 
		# Number of sites in the x-direction : 
        par.Nx = 1		# is 1 since 1D system. Unitless.
		# Number of sites in the y-direction : 
        par.Ny = 60		
		# Lattice constant in the y-direction
		# is calculated from the length and the
		# number of sites in said direction : 
        par.ay = float(par.LY)/par.Ny
		# Square lattice implementation : 
        par.ax = par.ay
		# The only region ('left') in the x-direction where
		# sites are occupied : 
        par.LL = par.ax
	# If not 1D, implementing the 2D system described in Hell : 
    else:
        par = SimpleNamespace(
			LM=250., 
			LL=1.e3, 
			LR=1.e3, 
			LY=4.e3, 
			mu=0, 
			alphahbar=1.42e-4*const.c*const.hbar, 
			Gamma=180e-3, 	# Induced gap from self-energy due to
							# the proximity effect (see Hell)
			phiL=0, 
			phiR=0, 
			k=30
			)
        par.Nx = Nx 		# IMHERE <<--- should be a problem. 
        par.Ny = int(round(float(par.Nx/1.625)))  # rounding off to nearest int
        par.ax = float(par.LM+par.LL+par.LR)/par.Nx
        par.ay = float(par.LY)/par.Ny

    par.GammaL = par.Gamma
    par.GammaR = par.Gamma 	# GammaL/R used in make.gamma_phi_fn(site,p)
    par.m_star = 0.023*const.m_e/(const.c)**2
    if ppar.make_1DNLeft_1DN_2DSMiddle_No_NRight or ppar.make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight:
        # !!!		# effective mass in N ~ 20 times effective mass in S, typically [Proximity effect exercise]
        par.m_N = 20.*par.m_star
    if ppar.make_1D_Heff_LR_leads:
        # effective mass is m star, no distinction betwenn N and S mass because now Heff model
        par.m_N = par.m_star
    # [t]=meV(=(meV*s)^2/(meV/(nm/s)^2*nm^2))
    par.tx = const.hbar**2/(2*par.m_N*par.ax**2)
    par.ty = const.hbar**2/(2*par.m_N*par.ay**2)
    par.hbar = const.hbar

    ### REAL SPACE DIMENSIONS ###################################################
    # Setting symmetric or asymmetric (par.right=[]) device;
    # Def. lattice numbers for the left, middle and right part. starts from 0.
    # Defined such that NL=NR, by making NM even/odd depending on Nx:
    #############################################################################
    par.left = np.arange(0, round(float(par.LL)/par.ax))
    par.middle = np.arange(round(float(par.LL)/par.ax),
                           round(float(par.LL+par.LM)/par.ax))
    # only symm device now. for both, do --> #for par.right in (np.arange(round(float(par.LL+par.LM)/par.ax),round(float(par.LL+par.LM+par.LR)/par.ax)),[]):	# NL,NM,NR:lattice numbers, starting from 0 in NL
    par.right = np.arange(round(float(par.LL+par.LM)/par.ax),
                          round(float(par.LL+par.LM+par.LR)/par.ax))

    #############################################################################
    ######### DOMAINS AND RESOLUTIONS OF INDEPENDENT VARIABLES AND PARAMETERS ###
    ######### *** ###############################################################
    par.biasEres = 1001  # 251####1001#301
    # vs. mu:
    par.Ezres = 1
    par.mures = 1001  # 251####1001#200
    par.Ez = 0.8  # 1.6####	# Ez is the constant, mu is the variable. Ez = 2*V_Z in master thesis
    # 0.74,0.76,par.mures)#54,par.mures)		# [mu] = meV
    par.mu_values = np.linspace(-0.6, 0.6, par.mures)
    par.mu = par.mu_values[0]
    par.Ez_values = np.array([par.Ez])
    ppar.var_name, ppar.par_const_name, ppar.par_const = "mu", "Ez", par.Ez

    # ###### vs. Ez:
    # par.Ezres = 100
    # par.mures = 1
    # par.mu = 0.	# Now, mu is the constant, not Ez
    # par.Ez_values = np.linspace(6.5,8.3,par.Ezres) 	# [Ez] = meV
    # par.Ez = par.Ez_values[0] 		# !!! just so that not error later due to par.Ez undefined
    # par.mu_values = np.array([par.mu])
    # ppar.var_name,ppar.par_const_name,ppar.par_const = "Ez","mu",par.mu
    # #############################################################################

    # Screen Gonly:
    # par.biasenergies = np.linspace(-2.,2.,par.biasEres) 	# [meV]
    # Screen "inside gap":
    par.bias_maxmin = 0.25
    par.biasenergies = np.linspace(-par.bias_maxmin,
                                   par.bias_maxmin, par.biasEres) 	# [meV]

    #### PINCHER AND GAMMAS #########
    # Notes
    # -----
    # - as long as Gamma_L/R is specified, dos is not used in calculating g012A
    #################################
    par.pincher = True
    par.t_p_coeff = 3
    par.t_pincher = par.t_p_coeff*par.ty

    """
	########## v16: pinchertype specified here, then implemented in make.make_1D_system. It can take the following values:
	 - "onsite"
	 - "hopping"
	 - "onsite and hopping"
	######################
	"""
    ppar.pinchertype = "onsite"

    # one-level approx model
    par.Gamma_L = np.pi*np.abs(par.ty)**2*np.ones(par.mures)  # = 3.37 meV
    par.Gamma_R = np.pi*np.abs(par.ty)**2*np.ones(par.mures)
    ######################

    #### DOUBLERES: 															 #########
    #### used when wanting to load conductance files with different resolutions. #########
    ######################################################################################
    ppar.doubleres = True
    if ppar.doubleres == True:
        ###### Vs. Ez : ################
        # par.mures_2 = 1
        # par.Ezres_2 = 200
        # par.biasEres_2 = 201
        # par.mu_2 = 0.
        # par.mu_values_2 = np.array([par.mu_2])
        # par.Ez_values_2 = np.linspace(0,10,par.Ezres_2)
        # par.biasenergies_2 = np.linspace(-0.5,0.5,par.biasEres_2)
        ###### Vs. mu : #################
        par.mures_2 = 1001  # 251#IMHERE##1001
        par.Ezres_2 = 1
        par.biasEres_2 = 1001  # 251#IMHERE##1001
        par.Ez_2 = 0.8  # 1.6#IMHERE##0.8 			# didived by 2 gives V_Z in thesis!!!
        par.mu_values_2 = np.linspace(-.6, .6, par.mures_2)
        par.Ez_values_2 = np.array([par.Ez_2])
        par.biasenergies_2 = np.linspace(-par.bias_maxmin,
                                         par.bias_maxmin, par.biasEres_2)

        par.Gamma_L = np.pi*np.abs(par.ty)**2 * \
            np.ones(par.mures_2)  # = 3.37 meV
        par.Gamma_R = np.pi*np.abs(par.ty)**2*np.ones(par.mures_2)

    """
	FOR Calculating/loading 
	 			E1, Em1, evec1, evecm1,
				u, v, u^2 - v^2 for the n'th pair of eigenvectors:

	Parameters
	----------
	n : 		index of two eigenstates that have same absolute value of their energy for every value of the variable.
	"""
    # n = [24,25] 	# given k = 50
    # n = [9,10]	# given k = 20
    # n = [18,19] 	# gien k = 38
    n = [14, 15]		# given k = 30 OK
    ppar.n = n

    ##################################################################################
    ###### IF ADDING SC LEAD TO 1D N WIRE                                  ##########
    ###### - set the separate semiconducting and superconducting parameters ##########
    ######
    # 2tx_SC + 2ty_SC + mu_SC = 0
    ######	mu_N = 0
    ######
    ######
    ######
    ######
    ######
    ######
    ##################################################################################
    ### calculating k_F, lambda_F - unnecessary ###
    # par.m_SI = 0.02*const.m_e_SI
    # par.alpha_SI = const.c_SI*1e-4
    # par.kF = 2*par.m_SI*par.alpha_SI/(const.hbar_SI)
    # par.lambda_F = np.pi*const.hbar_SI/(2*par.m_SI*par.alpha_SI)

    ######### par.t_N_lead matters to make pattern backgrounds white (in conductance) ######
    par.t_N_lead = 1.  # meV
    if ppar.make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight or ppar.make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight or ppar.make_1D_NLeft_1D_S_Heff_No_NRight or ppar.make_1DNLeft_1DN_2DSMiddle_No_NRight:

        # ppar.pincher_SC_lead = True

        par.t_N = par.ty 	# is used for semiconductor hamiltonian functions in make module

        par.t_SC_lead_coeff = 0.0001
        par.tx_SC_lead = par.t_SC_lead_coeff*par.ty
        par.ty_SC_lead = par.t_SC_lead_coeff*par.tx
        par.mu_SC_lead = -2.0*(par.tx_SC_lead+par.ty_SC_lead)
        par.Delta_SC_lead = 0.180*1000.

        ######### Defining Ez of S and S lead ################################################
        ######### Differs from that of N because effective g-factor is different 	##########
        par.g_S_over_g_N = -1./5
        par.Ez_2D_S_lead = par.g_S_over_g_N*par.Ez
        par.Ez_2D_S = par.g_S_over_g_N*par.Ez

        # that particular system is implemented as 1D N attached to 2D sc lead, so pincher t == SC lead t.
        if ppar.make_1DNLeft_1DN_2DSMiddle_No_NRight == False:
            par.ty_SC_pincher_coeff = 1.  # 1./0.023
            par.ty_SC_pincher = par.ty_SC_pincher_coeff*par.ty
            par.tx_SC_pincher = par.ty_SC_pincher_coeff*par.tx
            par.mu_SC_pincher = par.mu_SC_lead
            par.Delta_pincher = par.Delta_SC_lead

    if ppar.make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight or ppar.make_1D_NLeft_1D_S_Heff_No_NRight:
        # to be used in calc to only calculate G_11, and not G_12 as well.
        ppar.oneNLead = True
    else:
        ppar.oneNLead = False

    return const, par, ppar
    ###
