[![DOI](https://zenodo.org/badge/192598166.svg)](https://zenodo.org/badge/latestdoi/192598166)

# Python and MatLab scripts used in Master thesis
## Thesis title: Localand nonlocal differential conductance spectroscopy of Andreev bound states in three-terminal superconducting devices 


### Chapter 3
#### Figures 3.9 and 3.10: 
	Run the script 
		
		Source_code/Single_Andreev_level_analytical/ABCD_plot.py



### Chapter 6
#### Figure 6.6:
	Majorana wire eigenenergies plot : Run the script

		Examples/Majorana_wire.py

	Differential conductance plot : Run the script

		Examples/NS_boundary_dIdV_Zbarriers.py	


### Chapter 7
#### Figure 7.3:
To plot already generated data for the symmetric and the antisymmetric decives : 
		
		Run : 		Examples/2D_platform_plot.py
	 	
		Output :	Two energy vs. E_Z plots, one for the symmetric, another for the asymmetric device.
				Filenames:
					EvsEz_Sym_<date-of-today>_260_160_0.8_0_250_1000_1000_4000_8.65385_25_0_'alphabar'_0.18_0_0
					EvsEz_Asym_<date-of-today>_260_160_0.8_0_250_1000_1000_4000_8.65385_25_0_'alphabar'_0.18_0_0
		
		These two output plots are also given in the folder
				"Examples"
		with the filenames
				"EvsEz_Asym_2019-09-06_260_160_0.8_0_250_1000_1000_4000_8.65385_25_0_'alphabar'_0.18_0_0.pdf"
				"EvsEz_Sym_2019-09-06_260_160_0.8_0_250_1000_1000_4000_8.65385_25_0_'alphabar'_0.18_0_0.pdf"
				
		This data has the same parameters as used in the analysis in chapter 7. :
				chemical potential = 0 meV
				superconducting phases in left and right 2D layers = 0
				lattice spacing in x-direction, a_x = 9 nm
				lattice spacing in y-direction, a_y = 25 cm
				magnetic field, E_Z = range(0,0.8,80) (a higher resolution was used in figure 7.3)
				Length in middle = 250 nm
				Length on left / right sides = 1000 nm
				Length in y-direction = 4000 nm
				alpha_R (Rashba spin-orbit strength) 	= 0.28 eV Å
									= 1.42*10^(-4)*c*hbar
					where 	c = speed of light [nm/s]
						hbar = reduced Planck constant [meV*s]
				effective mass, m_star = 0.023*m_e
					where m_e = electron mass [meV/c^2]
		

To generate the data from scratch : 
	
		Run : 		Source_code/Hell_2D_platform_paper/2D_platform_diagonalization.py
		
		Output : 	Looping over the symmetric and asymmetric device, it will show and save a plot when the diagonalization of each one is done.

		To change the resolution/values of fields, E_Z, looped over:
				change line 192 in the function "main()".
		
		To change the system parameters : 
				change the "par" object. The essential parameters are given in
					line 226
				where "par" is initially defined.
				
				the effective mass of the superconductor, m_star, is defined in 
					line 238.

		To change the system aspect ratio : 
				change line 234. The number dividing "par.Nx" gives the aspect ratio
					aspect_ratio = par.Nx/par.Ny.
		
		
			


### Chapters 8 and 9
#### Equations (8.26)-(8.29) :
Run the Mathematica stript
	
	"Source_code/Single_Andreev_level_analytical/S-matrix_one-level_approximation_web.nb"

#### Figures : 
####		8.3 - Nonlocal antisymmetric differential conductance, N = 800 sites
####		8.4 - Nonlocal antisymmetric and symmetric conductances, traced along one of the two lowest eigenenergies from diagonalizing the closed system, N = 800 sites
####		9.2 - Probability density versus lambda = (chemical potential, or magnetic field), and site number (real space), N = 200 (length same as above, namely 1500 nm)
####		9.3 - Spin polarization versus lambda, and site number, N = 200 sites
#### 		9.4 - BCS charge versus lambda, and site number, N = 200 sites

The code used to generate the data was adapted from a package written by Esben Bork Hansen.
This whole package will not be uploaded here. In the future, Esben plans to make the package public.
Please see his PhD thesis for further details about the kinds of calculations he has performed with his code.

To generate these figures, run the script : 

	Source_code/Three-terminal_selfenergy_method/figures/thesis_plots_Anna/Three-terminal_selfenergy_method_plots.py

This script makes use of the following library and data : 
	
	The data used in the thesis are given in : 

			"Source_code/Three-terminal_selfenergy_method/figures/thesis_plots_Anna/data/wire_101"

	One library I wrote, which is used to process the raw data (differential conductances), and to represent the data as done in the figures at hand, is given in : 

			"Source_code/Three-terminal_selfenergy_method/thesis_plot_library/thesis_plot_library.py"

	Arrays containing the effective topological gap, plotted in figure 9.2, is given in :

			"Source_code/Three-terminal_selfenergy_method/Anna\ thesis\ plots/Dtop_mu_1.npz"
			"Source_code/Three-terminal_selfenergy_method/Anna\ thesis\ plots/Dtop_b_1.npz"

	versus the chemical potential and the magnetic field, respectively.



#### Figure : 
####            8.5 - Line-cut of the in-gap conductance calculated for the same system as in the main results, but now using Kwant. 

The discretization is chosen differently, partially due to computational cost, but also because the point in this figure was to demonstrate

	G_LR^asym(V) = G_LL^asym(V)

which should be independent of the validity of the tight-binding approximation.
Discretization parameters used :  

	N = 60 sites, Length=1500 nm, inter-site distance = 25 nm.

To plot the figure from already generated data, run the script :

	Source_code/Three-terminal_Kwant_method/expsys_v16_webv02.py

This script generates several figures:

	N-SC_sys.pdf			 
					 - a plot of the graph for the system
	G_11_vsEbias1mu_Nx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1__SymLogNorm.pdf
					 - local total conductance
	G_12_vsEbias1mu_Nx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1__SymLogNorm.pdf
					 - nonlocal total conductance
	G_11_S_vsEbias1mu_Nx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1__SymLogNorm.pdf
					 - symmetric component of local conductance
	G_11_A_vsEbias1mu_Nx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1__SymLogNorm.pdf
					 - antisymmetric conponent of local conductance
	G_12_S_vsEbias1mu_Nx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1__SymLogNorm.pdf
					 - symmetric component of nonlocal conductance
	G_12_A_vsEbias1mu_Nx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1__SymLogNorm.pdf
					 - antisymmetric component of nonlocal conductance
	G_11_12_mu_0_vs_biasNx1_Ny60_LM0_LL25_LR0_LY1500_ax25_ay25_mu-0.6_Gamma0.18_pL0_pR0_k30_biE-0.25_0.25_Ez0.8_0.8_mu-0.6_0.6_biEres1001_Ezres1_mures1001_pinchcoef3_tsys2.65042_tLeads1_.pdf
					 - line-cut in the total local and the total nonlocal differential conductances, at a chemical potential of 0.8 meV (the 667'th element of the 1001 long chemical potential-axis, in this particular example data). 
	G_11_12_A_S_mu_0_vs_bias.pdf
					 - line-cut of the symmetry-decomposed local and nonlocal differential conductances. Zooming in around zero bias in this plot, you get figure 8.5.


The corresponding arrays containing the data for the conductances :

	The names are identical, except for the file-type. Substitute ".pdf" for ".npy" in the above.
