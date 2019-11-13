"""
Plotting already generated data for self-energy model of three-terminal device,
using functions from the module thesis_plot_library.
"""

from thesis_plot_library.thesis_plot_library import *
import sys
from pathlib import Path
sys.path.append(str(Path('.').absolute().parent))
sys.path.append('../../')

########## Plot Asymmetric differential conductances: ##########
var = "mu"
var_2 = "b"
filepath = "wire_101/sites=800/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"

filename_a = "diff_cond_12_of_mu_in_-600-600_t=47118.61755164704_B=400 (copy)"
filepath_a = (filepath + "/diff_cond_of_" + var)
data_a = np.load(
            "data/" + filepath_a + "/data_" + filename_a +  ".npz"
            )

filename_ = "diff_cond_12_of_mu_in_-150-150_t=47118.61755164704_B=400 (copy)"
filepath_ = (filepath + "/diff_cond_of_" + var)
data_a_ = np.load(
            "data/" + filepath_ + "/data_" + filename_ +  ".npz"
            )

filename_b = "diff_cond_12_of_b_in_0-800_t=47118.61755164704 (copy)"
filepath_b = (filepath + "/diff_cond_of_" + var_2)
data_a_b = np.load(
            "data/" + filepath_b + "/data_" + filename_b +  ".npz"
            )

filename_b_zoom = "diff_cond_12_of_b_in_200-800_t=47118.61755164704 (copy)"
filepath_b_zoom = (filepath + "/diff_cond_of_" + var_2)
data_a_b_zoom = np.load(
            "data/" + filepath_b_zoom + "/data_" + filename_b_zoom +  ".npz"
            )

x_range, e_range, diff_cond_a = \
        data_a["x_range"], data_a["e_range"], data_a["diff_cond"]
x_range_zoom, e_range_zoom, diff_cond_a_zoom = \
        data_a_["x_range"], data_a_["e_range"], data_a_["diff_cond"]
x_range_2, e_range_2, diff_cond_a_2 = \
        data_a_b["x_range"], data_a_b["e_range"], data_a_b["diff_cond"]
x_range_2_zoom, e_range_2_zoom, diff_cond_a_2_zoom = \
        data_a_b_zoom["x_range"], data_a_b_zoom["e_range"], data_a_b_zoom["diff_cond"]

#################### Plot Asymmetric differential conductances:
fig, ax, cbar = plot_pcolormesh_4subplots_wzoom_a(\
                        # 1st variable w zoom:
                        x_range, e_range, diff_cond_a, 
                        x_range_zoom, e_range_zoom, diff_cond_a_zoom,\
                        # 2nd variable w zoom:
                        x_range_2, e_range_2, diff_cond_a_2,\
                        x_range_2_zoom, e_range_2_zoom, diff_cond_a_2_zoom,\
                        # specifications:
                        colormap_log_2_value =0.003,
                        colormap_limits=[False,-0.011,0.011],\
                        colormap_norm=[True, 0], \
                        colormap_name='seismic',#'RdBu_r',
                        colormap_log=[True, 0.0001])

ft = "pdf"
fig.savefig("figures/" + filepath + "/thesis_plots_Anna" + "/fig_" \
                                + filename_a + "." + ft, format=ft, dpi=300)
fig.show()
#################################################################



###########     PLOTTING WAVEFUNCTION/SPIN POLARIZATION/CHARGE 
#               as function of lattice sites in real space      ###########
var = "mu"
var_2 = "b"
filepath = "wire_101/sites=800/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"
filepath = "wire_101/sites=100/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"
filepath = "wire_101/sites=200/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"

filepath_b = (filepath + "/Psi_components_vs_b")
filename_b = "Psi_=b_in_200-800 (another copy)"
filename_b = "Psi_=b_in_0-800"
data_b =np.load(
                "data/" + filepath_b + "/data_" + filename_b + ".npz"
                )

filepath_mu = (filepath + "/Psi_components_vs_mu")
filepath_mu = ("../Jeroen tight binding code/data/")
filepath_mu = "wire_101/sites=100/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180" \
                                        + "/Psi_components_vs_mu/"
filepath_mu = "wire_101/sites=200/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180" \
                                        + "/Psi_components_vs_mu/"


filename_mu = "Psi_=mu_in_-150-150_B=400 (another copy)" # old domain, far inside top region
filename_mu = "dat_190529_114228" # new larger trial domain with Jeroen model
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

# Regions of well-behaved energy 
# (level-crossings, with not-near zero energy modes, have not occured inside said regions) : 
b_domain = 500
b_domain_2 = 975
mu_domain_2 = 792 
mu_domain = 317


""" 
'Cutting and pasting' lowest-energy modes 
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
Calculating the end charges, q, end weights, n, and normalized end charges, Q :
"""
q_E0_R_mu, n_E0_R_mu, q_E0_L_mu, n_E0_L_mu, Q_E0_R_mu, Q_E0_L_mu,\
q_Em0_R_mu, n_Em0_R_mu, q_Em0_L_mu, n_Em0_L_mu, Q_Em0_R_mu, Q_Em0_L_mu = calc_QLR(wf_mu[:])

q_E0_R_b, n_E0_R_b, q_E0_L_b, n_E0_L_b, Q_E0_R_b, Q_E0_L_b,\
q_Em0_R_b, n_Em0_R_b, q_Em0_L_b, n_Em0_L_b, Q_Em0_R_b, Q_Em0_L_b = calc_QLR(wf_b[:])


"""
Calculating the shift between integrated and local charges 
(used to estimate how localized the wavefunction may be, away from the ends of the wire) : 
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
Plotting the probability density, spin polarization, and the charge, as functions of 
        lambda = chemical potential, or magnetic field, on the x-axis, and
        sites (real space) on the y-axis. 
"""
fig, ax, fig_spin, fig_Q, wf2_P_l_A_L_mu_b, spin_P_l_A_L_mu_b, Q_P_l_A_L_mu_b \
                                = calc_plot_wf2_vs_var_and_site(   
                                    mu = 			data_mu["var_range"],
                                    wf_mu = 		wf_mu,
                                    b = 			data_b["var_range"], 
                                    wf_b = 			wf_b,
                                    phase_data=[    False,
                                                    (QR_mu_i+(dEdmu_i))/2., ph_mu_i,
                                                    (QR_b_i+(dEdb_i))/2., ph_b_i
                                                ],
                                    Q_realspace=[	True,
													dEdmu,
													dEdb
												],
                                    Dtop2 = [		True,  
													Dtop_mu_data["var_range"], 
													Dtop_mu_data["D_top"] ,\
                                                    Dtop_b_data["var_range"], 
													Dtop_b_data["D_top"]
											],
                                    level_cross_mu = [True,409,'gray'] )


########        PLOTTING TRACES 
#               of conductances along the two lowest energies   ########
var = "mu"
var_2 = "b"
filepath = "wire_101/sites=800/length=1.5/alpha=0.280204287582202eVA/Delta_ind=180"

# plot gs, ga (traced), E0; QR, dE, gS/gA:
filepath_b = (filepath + "/Psi_components_vs_b")
filename_b = "Psi_=b_in_200-800_gS_and_gA_line_E0_Q__tunnel_t=47118.61755164704 (another copy)"
data_b =np.load(
                "data/" + filepath_b + "/data_" + filename_b + ".npz"
                )

filepath_mu = (filepath + "/Psi_components_vs_mu")
filename_mu = "Psi_=mu_in_-150-150_B=400_gS_and_gA_line_E0_Q__tunnel_t=47118.61755164704 (another copy)"
data_mu = np.load(
                 "data/" + filepath_mu + "/data_" + filename_mu + ".npz"
                 )

fig,ax, data_phaseshifts = plot_E0_gS_gA_and_QR_dE0_gSovergA_4subplots(   
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
QR_mu_zeros_and_extrema=data_phaseshifts[2]
star_point = QR_mu_zeros_and_extrema[2]
triangle_point = QR_mu_zeros_and_extrema[3]

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

filename_a = "E0_gS_gA_and_QR_dE0_gSovergA_4subplots"
print(" - Filepath where conductance trace plot is saved : ", filepath)
ft = "pdf"
fig.savefig("figures/" + filepath + "/thesis_plots_Anna" + "/fig_" \
                                + filename_a + "." + ft, format=ft, dpi=300)
fig.show()























##############################  NOT INCLUDED IN THESIS : #############################################


""" The shifts (on the x-axis) of integrated versus end charges, saved as a numpy zip file : """
np.savez("figures/" + filepath + "/thesis_plots_Anna" + "/data_phaseshifts" + ".npz", 
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

Q_E0_mu = u2_up_mu + u2_down_mu - v2_up_mu - v2_down_mu
Q_E0_b = u2_up_b + u2_down_b - v2_up_b - v2_down_b

cut_idx = 10    # index of mu for which to cut
cut_idx1 = 10   # index of b for which to cut

Pl_mu_wf2, Pl_b_wf2 = wf2_P_l_A_L_mu_b[0], wf2_P_l_A_L_mu_b[1]
Pl_mu_spin, Pl_b_spin = spin_P_l_A_L_mu_b[0], spin_P_l_A_L_mu_b[1]
Pl_mu_Q, Pl_b_Q = Q_P_l_A_L_mu_b[0], Q_P_l_A_L_mu_b[1]

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


