#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''
We use the functions we introduced on the calculations for the hadronic EOS:
    tov_equations(r, y, energy_density_func,d_energy_density_func):
    solve_tov(Pc, energy_density_func,d_energy_density_func):
    compute_k2_y_lambda(y_R, R_star, M_star)
    plot_quantity_vs_quantity(xkey, ykey, xlabel, ylabel, title, fname, xlim=None, ylim=None)
'''

def make_hybrid(b, d, P_tr=30.0, e_gap=200.0): #Maxwell construction for a Hybrid EOS
    '''
    Inputs:

            b and d are CFL “bag constant” and pairing gap parameters.

            P_tr (default 30 MeV·fm⁻³) is the pressure where we switch from hadronic to quark matter.

            e_gap shifts the energy density of the quark phase upward (default 200 MeV·fm⁻³).
            
            for the hadronic eos we use mdi1

    Returns:
        eps_hybrid(P):the energy density at pressure P piecewise-defined (crust → hadronic core → quark core).
        dε/dP, the analytic derivative needed for the tidal calculations.


    '''
    def eps_hybrid(P):
        if P > P_tr:
            return cfl_eos(P, b, d) + e_gap
        elif P >= 0.184:
            return mdi1(P)
        elif P >= 9.34375e-5:
            return h_crust1(P)
        elif P >= 4.1725e-8:
            return h_crust2(P)
        else:
            return h_crust3(P)

    def deps_hybrid(P):
        hbarc = 197.327
        ms    = 95.0
        a     = -ms**2/6 + 2*d**2/3

        if P > P_tr:
            D = np.sqrt((4*np.pi**2/3)*(b + P)*hbarc**3 + (3*a)**2)
            return 3.0 - 6.0*a / D
        elif P >= 0.184:
            return (4.1844 * 0.81449 * P**(0.81449 - 1.0)
                    + 95.00135 * 0.31736 * P**(0.31736 - 1.0))
        elif P >= 9.34375e-5:
            return ((103.17338/0.36227) * np.exp(-P/0.36227)
                    + (7.34979/0.012721) * np.exp(-P/0.012721))
        elif P >= 4.1725e-8:
            return ((0.00203/3448.27) * np.exp(-P/3448.27)
                    + (0.10851/7692.3076) * np.exp(-P/7692.3076))
        else:
            return ((0.000051 * 0.2373e10 * np.exp(-P * 0.2373e10))
                    + (0.00014 * 0.4020e8   * np.exp(-P * 0.4020e8)))

    return eps_hybrid, deps_hybrid


# --- list of the CFL parameter sets ---
parameter_sets = [
    {"b": 60,  "d": 50},  {"b": 60,  "d": 100}, {"b": 60,  "d": 150},
    {"b": 80,  "d": 50},  {"b": 80,  "d": 100}, {"b": 80,  "d": 150},
    {"b":100,  "d": 60},  {"b":100, "d":100},  {"b":100, "d":150},
    {"b":120, "d":100},  {"b":120, "d":150},  {"b":140, "d":120},
    {"b":140, "d":150},  {"b":160, "d":140},
]
labels = [f"CFL {i+1}" for i in range(len(parameter_sets))]

Pc_values  = np.linspace(1.0, 1200.0, 2000)

plt.figure()
#Loop over each hybrid model and solve the TOV
for params, label in zip(parameter_sets, labels):
    eps_hyb, deps_hyb = make_hybrid(params["b"], params["d"])
    R_list, M_list,yr_list = [], [], []
    for Pc in Pc_values:
        R, M, yr = solve_tov(Pc, eps_hyb, deps_hyb)
        R_list.append(R)
        M_list.append(M)
        yr_list.append(yr)
        # convert to arrays
    R_arr  = np.array(R_list)
    M_arr  = np.array(M_list)
    yr_arr = np.array(yr_list)

    save_dict = {}
    # save to .npz for future use
    np.savez_compressed(f"hybrid_{label}.npz",
                        R=R_arr, M=M_arr, yr=yr_arr)
    


# In[ ]:


#Pressure-Energy Density Diagram for Hybrid EOS using Maxwell construction

# --- Given crust & hadron definitions ---
def mdi1(P):
    return 4.1844 * P**0.81449 + 95.00135 * P**0.31736

def cfl_eos(P, b, d):
    hbarc = 197.327
    ms    = 95.0
    a     = -ms**2 / 6 + 2 * d**2 / 3
    mu2   = -3 * a + np.sqrt((4 * np.pi**2 / 3) * (b + P) * hbarc**3 + (3 * a)**2)
    return 3 * P + 4 * b - (9 * a * mu2) / (np.pi**2 * hbarc**3)

# Transition pressure and parameter choices
Ptr = 30.0         # Transition pressure [MeV fm^-3]
b_val = 80         # Bag constant
d_val = 100        # Pairing gap
e_gap_100 = 100    # Energy gap for first case [MeV]
e_gap_200 = 200    # Energy gap for second case [MeV]

# Define pressure ranges
P_nuc   = np.linspace(0.184, Ptr, 300)   # Hadronic region
P_quark = np.linspace(Ptr, 200, 300)     # Quark region

# Compute energy densities
E_nuc = mdi1(P_nuc)  # purely hadronic

E_quark_gap100 = cfl_eos(P_quark, b_val, d_val) + e_gap_100
E_quark_gap200 = cfl_eos(P_quark, b_val, d_val) + e_gap_200

# Maxwell jump lines at Ptr
E_H_Ptr = mdi1(Ptr)
E_Q_Ptr_100 = cfl_eos(Ptr, b_val, d_val) + e_gap_100
E_Q_Ptr_200 = cfl_eos(Ptr, b_val, d_val) + e_gap_200

E_maxwell_100 = np.array([E_H_Ptr, E_Q_Ptr_100])
E_maxwell_200 = np.array([E_H_Ptr, E_Q_Ptr_200])
P_maxwell = np.array([Ptr, Ptr])

# Plot Pressure vs Energy Density (ε horizontal, P vertical)
plt.figure(figsize=(8, 6), dpi=150)

# Pure hadronic curve
plt.plot(E_nuc, P_nuc, label='Hadronic (MDI-1)',
         color='tab:blue', linewidth=2.0)

# Pure quark curves for two gap choices
plt.plot(E_quark_gap100, P_quark, label='CFL4 (gap=100 MeV)',
         color='tab:orange', linewidth=2.0, linestyle='-')
plt.plot(E_quark_gap200, P_quark, label='CFL4 (gap=200 MeV)',
         color='tab:red', linewidth=2.0, linestyle='-.')

# Maxwell construction vertical lines
plt.plot(E_maxwell_100, P_maxwell, label='Maxwell Jump (gap=100)',
         color='black', linewidth=1.5, linestyle='--')
plt.plot(E_maxwell_200, P_maxwell, label='Maxwell Jump (gap=200)',
         color='black', linewidth=1.5, linestyle=':')

# Annotate transition pressure
plt.axhline(Ptr, color='grey', linestyle=':', linewidth=1.0)
plt.text(E_H_Ptr - 50, Ptr + 2, r'$P_{\rm tr}=30\,\mathrm{MeV\,fm^{-3}}$',
         fontsize=10, color='grey')

# Labels and formatting
plt.xlabel(r'Energy Density $\varepsilon$ (MeV fm$^{-3}$)', fontsize=12)
plt.ylabel(r'Pressure $P$ (MeV fm$^{-3}$)', fontsize=12)
plt.title('Pressure vs. Energy Density for Maxwell Hybrid EOS', fontsize=14)

# Axis limits
plt.xlim(0, 650)
plt.ylim(0, 70)

# Tick parameters
plt.tick_params(axis='both', which='major', labelsize=10, length=6)
plt.tick_params(axis='both', which='minor', labelsize=8, length=3)
plt.minorticks_on()

# Grid
plt.grid(which='major', color='grey', linestyle='--', linewidth=0.5, alpha=0.7)
plt.grid(which='minor', color='grey', linestyle=':', linewidth=0.3, alpha=0.5)

# Legend outside plot
plt.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)

plt.tight_layout()

# Save the figure
plt.savefig("HybridEOS_P_vs_Epsilon_Maxwell_Styled.png", dpi=300, bbox_inches="tight")
plt.savefig("HybridEOS_P_vs_Epsilon_Maxwell_Styled.pdf", bbox_inches="tight")

plt.show()


# In[ ]:


#The MR diagram for the hybrid CFL with the MR curve of mdi1

# 1) Load MDI-1 data
mdi_data = np.load('eosresults.npz')
R_mdi = mdi_data['MDI-1_R']
M_mdi = mdi_data['MDI-1_M']

# 2) Load hybrid CFL data
hyb_data = np.load('hybrid_cfl_data.npz')
labels = [f"CFL {i+1}" for i in range(14)]

# 3) Plot everything
plt.figure(figsize=(8, 6), dpi=150)

# Plot all CFL hybrids (dashed lines)
for label in labels:
    R_h = hyb_data[f"{label}_R"]
    M_h = hyb_data[f"{label}_M"]
    plt.plot(R_h, M_h,
             linestyle='--',
             linewidth=1,
             alpha=0.7,
             label=label)

# Plot MDI-1 (solid bold line)
plt.plot(R_mdi, M_mdi,
         color='black',
         linewidth=2.5,
         label='MDI-1')
plt.xlim(right=18)
plt.xlabel(r"$R\ (\mathrm{km})$", fontsize=12)
plt.ylabel(r"$M\ (M_\odot)$",   fontsize=12)
plt.title("MDI-1 & Hybrid CFL M–R Curves", fontsize=14, pad=10)

plt.grid(which='major', linestyle='--', linewidth=0.6, alpha=0.7)
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth=0.4, alpha=0.5)

plt.legend(fontsize=8, ncol=2, frameon=True)
plt.tight_layout()

# 4) Save to disk before showing
plt.savefig("MDI1_HybridCFL_MR.png", dpi=300, bbox_inches="tight")
plt.savefig("MDI1_HybridCFL_MR.pdf", bbox_inches="tight")

# 5) (Optional) Display on screen
plt.show()


# In[ ]:


#The tidal plots for a sample of hybrid eos models

# 2) Load saved results into a dict
labels = ["CFL1","CFL3","CFL5","CFL10"]
models = {}
for lbl in labels:
    data = np.load(f"hybrid_{lbl}.npz")
    R   = data["R"]
    M   = data["M"]
    yR  = data["yr"]
    C   = 1.474 * M / R
    k2, lam = compute_k2_y_lambda(yR, R, M) #the fucntion that calculates k2,λ
    models[lbl] = {"M":M, "R":R, "yR":yR, "C":C, "k2":k2, "lam":lam}


# Now create all  diagrams in a consistent style using our function for ploting: plot_quantity_vs_quantity

plot_quantity_vs_quantity('R',    'k2', r"Radius $R$ (km)",         r"Love number $k_2$",  r"$k_2$ vs. Radius",        "k2_vs_radius.pdf"     )
plot_quantity_vs_quantity('M',    'k2', r"Mass $M\,(M_\odot)$",     r"Love number $k_2$",  r"$k_2$ vs. Mass",          "k2_vs_mass.pdf"       )
plot_quantity_vs_quantity('C',    'k2', r"Compactness $C = M/R$",   r"Love number $k_2$",  r"$k_2$ vs. Compactness",   "k2_vs_compactness.pdf")

plot_quantity_vs_quantity('M',   'yr', r"Mass $M\,(M_\odot)$",      r"Dimensionless $y_R$", r"$y_R$ vs. Mass",         "yR_vs_mass.pdf")
plot_quantity_vs_quantity('R',   'yr', r"Radius $R$ (km)",          r"Dimensionless $y_R$", r"$y_R$ vs. Radius",       "yR_vs_radius.pdf")
plot_quantity_vs_quantity('C',   'yr', r"Compactness $C = M/R$",    r"Dimensionless $y_R$", r"$y_R$ vs. Compactness",  "yR_vs_compactness.pdf")

plot_quantity_vs_quantity('M', 'Lambda', r"Mass $M\,(M_\odot)$",    r"$\lambda$",           r"$\lambda$ vs. Mass",     "Lambda_vs_mass.pdf")
plot_quantity_vs_quantity('R', 'Lambda', r"Radius $R$ (km)",        r"$\lambda$",           r"$\lambda$ vs. Radius",   "Lambda_vs_radius.pdf)
plot_quantity_vs_quantity('C', 'Lambda', r"Compactness $C = M/R$",  r"$\lambda$",           r"$\lambda$ vs. Compactness", "Lambda_vs_compactness.pdf")

