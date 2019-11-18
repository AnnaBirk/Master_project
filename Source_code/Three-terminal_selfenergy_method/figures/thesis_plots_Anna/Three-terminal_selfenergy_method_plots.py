"""
In this file, we plot already generated data for self-energy model of three-terminal device,
using functions from the module thesis_plot_library.
"""
from thesis_plot_library.thesis_plot_library import \
                                plot_pcolormesh_4subplots_wzoom_a, \
                                calc_plot_wf2_vs_var_and_site, \
                                calc_dEdvar, \
                                calc_QLR, \
                                calc_E1_Em1_from_E0_E0prime_var, \
                                plot_E0_gS_gA_and_QR_dE0_gSovergA_4subplots, \
                                calc_Q_dE_mu_b_phase_shift, calc_wf_data, plot_Pl_and_cut_wf_spinpol_Q
import numpy as np
import sys
from pathlib import Path
sys.path.append(str(Path('.').absolute().parent))
sys.path.append('../../')


""" 
                PART 1 : ANTISYMMETRIC DIFFERENTIAL CONDUCTANCE

Load antisymmetric differential conductances: """
var_mu = "mu"
var_b = "b"
# filepath skeleton, path specified by folders whose names are specified by
# the system parameters : 
filepath_skeleton = "wire_101/sites=800/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"

# filename of differential conductance for chemical potential in 
# the domain from -600 to 600 micro electron volts : 
filename_mu = "diff_cond_12_of_mu_in_-600-600_t=47118.61755164704_B=400 (copy)"
filepath_mu = (filepath_skeleton + "/diff_cond_of_" + var_mu)
# loading differential conductance for chemical potential in said domain :
data_mu = np.load(
            "data/" + filepath_mu + "/data_" + filename_mu +  ".npz"
            )
# loading dfferential conductance for chemical potential in shorter domain,
# inside the topological region. The charge can be extracted from the oscillations
# in this domain :
filename_mu_zoom = "diff_cond_12_of_mu_in_-150-150_t=47118.61755164704_B=400 (copy)"
filepath_mu_zoom = (filepath_skeleton + "/diff_cond_of_" + var_mu)
data_mu_zoom = np.load(
            "data/" + filepath_mu_zoom + "/data_" + filename_mu_zoom +  ".npz"
            )
# loading differential conductance for magnetic field in the domain
# 0 to 800 micro electron volts ... :
filename_b = "diff_cond_12_of_b_in_0-800_t=47118.61755164704 (copy)"
filepath_b = (filepath_skeleton + "/diff_cond_of_" + var_b)
data_b = np.load(
            "data/" + filepath_b + "/data_" + filename_b +  ".npz"
            )
# ... and in the smaller domain, from 200 to 800 micro electron volts,
# where oscillations in the differential conductance are clear, and where
# the charge can be probed :
filename_b_zoom = "diff_cond_12_of_b_in_200-800_t=47118.61755164704 (copy)"
filepath_b_zoom = (filepath_skeleton + "/diff_cond_of_" + var_b)
data_b_zoom = np.load(
            "data/" + filepath_b_zoom + "/data_" + filename_b_zoom +  ".npz"
            )
# the data loaded above is in an .nzip folder containing three different
# arrays. These arrays can be accessed as you would a dictionary. The
# arrays have the following names :
#       "x_range" :     np.array, containing all values of the independent
#                       variable, here the chemical potential, or the magnetic
#                       field.
#       "e_range" :     np.array, containing all values of the energies the
#                       differential conductance is calculated for.
#       "diff_cond" :   np.array, containing all differential conductance
#                       data points, calculated in "x_range"-"e_range" "space".
#                       It has the shape
#                               ( np.shape(e_range), np.shape(x_range) )
#                       in the following :
mu_range, e_range_mu, diff_cond_mu = \
        data_mu["x_range"], data_mu["e_range"], data_mu["diff_cond"]
mu_range_zoom, e_range_mu_zoom, diff_cond_mu_zoom = \
        data_mu_zoom["x_range"], data_mu_zoom["e_range"], data_mu_zoom["diff_cond"]
b_range, e_range_b, diff_cond_b = \
        data_b["x_range"], data_b["e_range"], data_b["diff_cond"]
b_range_zoom, e_range_b_zoom, diff_cond_b_zoom = \
        data_b_zoom["x_range"], data_b_zoom["e_range"], data_b_zoom["diff_cond"]


""" Plot Asymmetric differential conductances : """
fig, ax, cbar = plot_pcolormesh_4subplots_wzoom_a(\
        # This function plots differential conductance in four panels :
        #  - upper left and right panels : 
        #       zoomed "out" differential conductance, in mu versus energy (left)
        #       and b versus energy (right) space
        #  - lower left and right panels :
        #       zoomed "in" differential conductance, --''--.
                        # 1st variable (left upper and lower panels):
                        mu_range, e_range_mu, diff_cond_mu, 
                        mu_range_zoom, e_range_mu_zoom, diff_cond_mu_zoom,\
                        # 2nd variable (right upper and lower panels):
                        b_range, e_range_b, diff_cond_b,\
                        b_range_zoom, e_range_b_zoom, diff_cond_b_zoom,\
                        # plotting specifications:
                        colormap_log_2_value =0.003,    # lower cut-off value for log-scale
                                                        # (see matplotlib documentation)
                        colormap_limits=[False,-0.011,0.011],   # If first element is False, use automatic limits set
                                                                # in matplotlib. If first element is True, use the
                                                                # lower and upper limit, specified by the second and
                                                                # third elements, respectively.
                        colormap_norm=[True, 0],        # Sets the midpoint in the colormap to be at the value
                                                        # zero. With a scale with different colors for values above
                                                        # and below zero, one can easily see the difference between
                                                        # positive and negative differential conductance in the plot.
                        colormap_name='seismic',        # colormapping between values of data
                                                        # and colors plotted.
                        colormap_log=[True, 0.0001])    # If first element is True, the colormap is scaled
                                                        # logarithmically. The specific scaling is then 
                                                        # specified by the second element. If the first element
                                                        # is false, the map is scaled linearly.

ft = "pdf"
fig.savefig("figures/" + filepath_skeleton + "/thesis_plots_Anna" + "/fig_" \
                                + filename_mu + "." + ft, format=ft, dpi=300)
fig.show()



"""     
        PART 2 : DIFFERENTIAL CONDUCTANCE SPECTROSCOPY

Loading data :
"""   
filepath_skeleton = "wire_101/sites=200/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"

filepath_b = (filepath_skeleton + "/Psi_components_vs_b")
filename_b = "Psi_=b_in_0-800"
data_b =np.load(
                "data/" + filepath_b + "/data_" + filename_b + ".npz"
                )

filepath_mu = (filepath_skeleton + "/Psi_components_vs_mu/")
filename_mu = "Psi_=mu_in_-600-600_B=400" # new trial domain larger
data_mu = np.load(
                 "data/" + filepath_mu + "/data_" + filename_mu + ".npz"
                 )

mu = data_mu["var_range"]
wf_mu = data_mu["psi_components"]
E0_mu = data_mu["energies"][:,0]

b = data_b["var_range"]
wf_b = data_b["psi_components"]
E0_b = data_b["energies"][:,0]

# Regions of well-behaved energy, found by inspection and hard-coded
# (level-crossings, with not-near zero energy modes, have not occured inside said regions) : 
b_domain = 500
b_domain_2 = 975
mu_domain_2 = 792 
mu_domain = 317

""" 
'Cutting and pasting' lowest-energy modes 
with the function calc_E1_Em1_from_E0_E0prime_var
(Notes
 -----
 Before : one of the energies is always positive, the other always negative. 
 After : two energies what moth oscillate, out of phase, around zero energy) 
"""
E0_mu, Em0_mu, wf_mu[:,:,:,0], wf_mu[:,:,:,1] = \
                                calc_E1_Em1_from_E0_E0prime_var(
                                        E0_mu[:],
                                        data_mu["energies"][:,1][:],mu[:], 
                                        evecs = [wf_mu[:,:,:,0], wf_mu[:,:,:,1]])

E0_b, Em0_b, wf_b[:,:,:,0], wf_b[:,:,:,1] = \
                                calc_E1_Em1_from_E0_E0prime_var(
                                        E0_b[:], 
                                        data_b["energies"][:,1][:],b[:], 
                                        evecs= [wf_b[:,:,:,0], wf_b[:,:,:,1]])

"""
Calculating the derivative of the lowest eigenenergies : 
"""
dEdmu, dEdb = calc_dEdvar(mu[:],b[:],E0_mu,E0_b)

"""
Calculating the end charges, q_L or q_R, end weights, n_L or n_R, 
and normalized end charges, Q_L or Q_R :
"""
q_E0_R_mu, n_E0_R_mu, q_E0_L_mu, n_E0_L_mu, Q_E0_R_mu, Q_E0_L_mu,\
q_Em0_R_mu, n_Em0_R_mu, q_Em0_L_mu, n_Em0_L_mu, Q_Em0_R_mu, Q_Em0_L_mu = calc_QLR(wf_mu[:])

q_E0_R_b, n_E0_R_b, q_E0_L_b, n_E0_L_b, Q_E0_R_b, Q_E0_L_b,\
q_Em0_R_b, n_Em0_R_b, q_Em0_L_b, n_Em0_L_b, Q_Em0_R_b, Q_Em0_L_b = calc_QLR(wf_b[:])


"""
Calculating the shift between integrated and local charges 
(used to estimate how localized the wavefunction may be, away from the ends of the wire.
Not plotted in thesis) : 
"""
QR_mu_i, dEdmu_i, ph_mu_i,\
QR_b_i, dEdb_i, ph_b_i = calc_Q_dE_mu_b_phase_shift(
                                mu[mu_domain:mu_domain_2], 
                                b[b_domain:b_domain_2], 
                                Q_E0_R_mu[mu_domain:mu_domain_2], 
                                dEdmu[mu_domain:mu_domain_2], 
                                Q_E0_R_b[b_domain:b_domain_2], 
                                dEdb[b_domain:b_domain_2])

"""
Loading the effective topological gap, 
calculated from an infinite Oreg-Lutchyn wire with the same parameters as in this modes : 
"""
Dtop_mu_data = np.load("Anna thesis plots/Dtop_mu_1.npz")
Dtop_b_data = np.load("Anna thesis plots/Dtop_b_1.npz")

"""
                PROBABILITY DENSITY/SPIN POLARIZATION/BCS CHARGE IN REAL SPACE
Plotting the 
        probability density, spin polarization, and the charge, as functions of 
                lambda = chemical potential, or magnetic field, on the x-axis, and
                sites (real space) on the y-axis. 
"""
fig, ax, fig_spin, fig_Q, wf2_P_l_A_L_mu_b, spin_P_l_A_L_mu_b, Q_P_l_A_L_mu_b \
                                = calc_plot_wf2_vs_var_and_site(   
                                        mu = 		data_mu["var_range"],
                                        wf_mu = 		wf_mu,
                                        b = 		data_b["var_range"], 
                                        wf_b = 		wf_b,
                                        phase_data=[    False,
                                                        (QR_mu_i+(dEdmu_i))/2., ph_mu_i,
                                                        (QR_b_i+(dEdb_i))/2., ph_b_i
                                                        ],
                                        Q_realspace=[	True,
                                                        dEdmu,
                                                        dEdb
                                                        ],
                                        Dtop2 = [	True,  
                                                        Dtop_mu_data["var_range"], 
                                                        Dtop_mu_data["D_top"] ,\
                                                        Dtop_b_data["var_range"], 
                                                        Dtop_b_data["D_top"]
                                                        ],
                                        level_cross_mu = [True,409,'gray'] )


"""             DIFFERENTIAL CONDUCTANCE SPECTROSCOPY
                        along the two lowest eigenenergies """
# loading data, similarly to earlier :                         
filepath_skeleton = "wire_101/sites=800/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"

filepath_b = (filepath_skeleton + "/Psi_components_vs_b")
filename_b = "Psi_=b_in_200-800_gS_and_gA_line_E0_Q__tunnel_t=47118.61755164704 (another copy)"
data_b =np.load(
                "data/" + filepath_b + "/data_" + filename_b + ".npz"
                )

filepath_mu = (filepath_skeleton + "/Psi_components_vs_mu")
filename_mu = "Psi_=mu_in_-150-150_B=400_gS_and_gA_line_E0_Q__tunnel_t=47118.61755164704 (another copy)"
data_mu = np.load(
                 "data/" + filepath_mu + "/data_" + filename_mu + ".npz"
                 )
# plotting traced symmetric (gS) and antisymmetric (gA) components of the nonlocal conductance (upper panel),
# the energy E0 along which the trace has been performed,
# the end charge at the right end, QR, 
# the derivative of the eigenenergy,
# and the charge probe, gS/gA : 
fig, ax, data_phaseshifts = plot_E0_gS_gA_and_QR_dE0_gSovergA_4subplots(   
                                            mu = data_mu["x_range"],
                                            E0_mu = data_mu["E0"],
                                            gS_mu = data_mu["G_LR_S"],
                                            gA_mu = data_mu["G_LR_A"],
                                            QR_mu = data_mu["QR_E0"],
                                            b = data_b["x_range"], 
                                            E0_b = data_b["E0"],
                                            gS_b = data_b["G_LR_S"],
                                            gA_b = data_b["G_LR_A"],
                                            QR_b = data_b["QR_E0"]
                                            )

""" 
ANNOTATING THE TWO TYPES OF ZEROS 
OBSERVED IN THE ANTISYMMETRIC CONDUCTANCE 
(I.E. THE MAXIMA AND MINIMA OF THE CHARGE, 
WHICH WE HAVE EASILY AVAILABLE FROM DATA_PHASESHIFTS) : 
"""
QR_mu_zeros_and_extrema = data_phaseshifts[2]
star_point = QR_mu_zeros_and_extrema[2]         # x-value of a 'dip' in gA
triangle_point = QR_mu_zeros_and_extrema[3]     # x-value of a zero energy crossing in gA

ax[0,0].annotate(r'$\star$',
            xy=(star_point,-0.2),
            xytext=(star_point,-2),
            arrowprops=dict(arrowstyle="->",relpos=(0, 0)),
            horizontalalignment='left',
            verticalalignment='bottom'
            )
ax[0,0].annotate(r'$\triangle\hspace{5pt}$',
            xy=(triangle_point,-0.2),
            xytext=(triangle_point,-2),
            arrowprops=dict(arrowstyle="->",relpos=(0, 0)),
            horizontalalignment='left',
            verticalalignment='bottom'
            )

filename_mu = "E0_gS_gA_and_QR_dE0_gSovergA_4subplots"
print(" - Filepath where conductance trace plot is saved : ", filepath_skeleton)
ft = "pdf"
fig.savefig("figures/" + filepath_skeleton + "/thesis_plots_Anna" + "/fig_" \
                                + filename_mu + "." + ft, format=ft, dpi=300)
fig.show()



""" NOT INCLUDED IN THESIS : 
plotting the phase function (P_l) which may be used to estimate the 'degree of localization' function (A_L),
as well as Fourier-transforming the signals : """

# The shifts (on the x-axis) of integrated versus end charges, saved as a numpy zip file : 
np.savez("figures/" + filepath_skeleton + "/thesis_plots_Anna" + "/data_phaseshifts" + ".npz", 
                phase_shifts_mu=data_phaseshifts[0], 
                phase_shifts_b=data_phaseshifts[1], 
                QR_mu_zeros_and_extrema=data_phaseshifts[2], 
                dEdmu_zeros_and_extrema=data_phaseshifts[3], 
                QR_b_zeros_and_extrema=data_phaseshifts[4], 
                dEdb_zeros_and_extrema=data_phaseshifts[5])

### SHOWING A LINE-CUT IN THE TRIVIAL REGION TO INVESTIGATE FOURIER SPECTRUM OF THE SIGNAL : ###
wf2_E0_mu_sum, wf2_E0_b_sum, spinpol_E0_mu, spinpol_E0_b, \
                                u2_up_mu, u2_down_mu, v2_up_mu, v2_down_mu, u2_up_b, \
                                u2_down_b, v2_up_b, v2_down_b \
                                                = calc_wf_data(wf_mu, wf_b)
# Charges : 
Q_E0_mu = u2_up_mu + u2_down_mu - v2_up_mu - v2_down_mu
Q_E0_b = u2_up_b + u2_down_b - v2_up_b - v2_down_b

# Choose an index to perform the line-cut in the desired signal 
# (wavefunction, spin, or charge) : 
cut_idx = 10    # index of mu for which to cut
cut_idx1 = 10   # index of b for which to cut

# Phase function for estimating 'degree of localization' in thesis : 
Pl_mu_wf2, Pl_b_wf2     = wf2_P_l_A_L_mu_b[0], wf2_P_l_A_L_mu_b[1]
Pl_mu_spin, Pl_b_spin   = spin_P_l_A_L_mu_b[0], spin_P_l_A_L_mu_b[1]
Pl_mu_Q, Pl_b_Q         = Q_P_l_A_L_mu_b[0], Q_P_l_A_L_mu_b[1]

plot_Pl_and_cut_wf_spinpol_Q(
                                wf2_E0_mu_sum, 
                                spinpol_E0_mu, 
                                Q_E0_mu, 
                                Pl_mu_wf2, 
                                Pl_mu_spin, 
                                Pl_mu_Q, 
                                cut_idx = cut_idx)

plot_Pl_and_cut_wf_spinpol_Q(
                                wf2_E0_b_sum, 
                                spinpol_E0_b, 
                                Q_E0_b, 
                                Pl_b_wf2, 
                                Pl_b_spin, 
                                Pl_b_Q, 
                                cut_idx = cut_idx1)


