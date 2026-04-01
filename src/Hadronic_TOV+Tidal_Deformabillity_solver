#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Define right -hand system of equations
#Calculates dy/dr which is needed to calculate the tidal deformabilty and love number.
def tov_equations(r, y, energy_density_func,d_energy_density_func):
    P, M, y_r = y
    # get EOS quantities
    eps =  energy_density_func(P)
    deps = d_energy_density_func(P)

    # F, J auxiliary parameters for tidal deformabillity
    F = (1.0
         - 1.474 * 11.2e-6 * r**2 * (eps - P)
        ) / (1.0 - 2.948 * M / r)

    J = (
        1.474 * 11.2e-6 * r**2 * (5*eps + 9*P + (eps+P)*deps)
        / (1.0 - 2.948*M / r)
        - 6.0 / (1.0 - 2.948*M / r)
        - 4.0 * (1.474**2) * (M**2)/(r**2)
          * (1 + 11.2e-6 * r**3 * (P/M))**2
          * (1 - 2.948*M/r)**-2
    )

    # TOV
    dP_dr = (
        -1.474 * eps * M * (1 + P/eps) / r**2
        * (1 + 11.2e-6 * r**3 * P / M)
        / (1 - 2.948*M/r)
    )
    dM_dr = 11.2e-6 * r**2 * eps

    # perturbation
    dy_r_dr = (-y_r**2 - y_r*F - J) / r

    return [dP_dr, dM_dr, dy_r_dr]


# In[ ]:



'''Gives the mass, radius, and parameter yR for a specific central pressure P(0)
 for a specific EOS'''
def solve_tov(Pc, energy_density_func,d_energy_density_func):
    y0    = [Pc, 1e-8, 2.0]   # [P_center, M_small but no equal to zero, y(0)=2]
    r0    = 1e-6
    r_max = 30.0              # km

    # stop when P drops below 1e-4, besically zero
    def event_surface(r, y):
        return y[0] - 1e-4
    event_surface.terminal  = True
    event_surface.direction = -1

    sol = solve_ivp(
        fun=lambda r, y: tov_equations(r, y, energy_density_func,d_energy_density_func),
        t_span=(r0, r_max),
        y0=y0,
        events=event_surface,
        rtol=1e-8,
        atol=1e-10,
        max_step=0.5
    )

    if sol.t_events[0].size > 0:
        R_star    = sol.t_events[0][0]
        P_surf, M_star, y_surf = sol.y_events[0][0]
    else:
        # fallback if the event never fired
        R_star    = sol.t[-1]
        P_surf, M_star, y_surf = sol.y[:, -1]

    return R_star, M_star, y_surf


# In[ ]:


Pc_values = np.linspace(1, 1200, 2000) #the central pressure values(1,1200) we integrate out

EoS_functions = [
    mdi1_crust, mdi2_crust, mdi3_crust, mdi4_crust, nld_crust, 
    hhj1_crust, hhj2_crust, ska_crust, ski4_crust, hlps3_crust, scvbb_crust,
    WFF_1_crust, WFF_2_crust, W_crust, BGP_crust, BL_1_crust, BL_2_crust, DH_crust, APR_1_crust
]# the hadronic EOS we integrate
dEoS_functions = [
    d_mdi1_crust, d_mdi2_crust, d_mdi3_crust, d_mdi4_crust, d_nld_crust, 
    d_hhj1_crust, d_hhj2_crust, d_ska_crust, d_ski4_crust, d_hlps3_crust, d_scvbb_crust,
    d_WFF_1_crust, d_WFF_2_crust, d_W_crust, d_BGP_crust, d_BL_1_crust, d_BL_2_crust, d_DH_crust, d_APR_1_crust
] #the derivatives of the hadronic EOS we integrate
labels = [
    "MDI-1", "MDI-2", "MDI-3", "MDI-4", "NLD", 
    "HHJ-1", "HHJ-2", "SKa", "SKI4", "HLPS3", "SCvBB",
    "WFF-1", "WFF-2",  "W", "BGP", "BL-1", "BL-2", "DH", "APR-1"
]


# --- Dictionary for results of each EOS ---
eos_results = {}
def compute_k2_y_lambda(y_R, R_star, M_star):
    '''Compute the quadrupolar Love number k2 and tidal deformability Lambda
    for a neutron star.
    Parameters
    ----------
    y_R :Value of the metric perturbation at the stellar surface.
    R_km :Stellar radius in kilometers.
    M_solar : Stellar mass in solar masses.
    Returns
    -------
    k2 :Quadrupolar Love number.
    Lambda :Dimensionless tidal deformability
'''
    C =1.474* M_star / R_star  # Compactness
    
    #Build up the four “term” pieces of the analytic k₂ formula
    #(from Hinderer 2008), where y_R is the metric perturbation at the surface.
    term1 = (1 - 2 * C) ** 2 * (2 + 2 * C * (y_R - 1) - y_R)
    term2 = 2 * C * (6 - 3 * y_R + 3 * C * (5 * y_R - 8))
    term3 = 4 * C**3 * (13 - 11 * y_R + C * (3 * y_R - 2) + 2 * C**2 * (1 + y_R))
    term4 = 3 * (1 - 2 * C)**2 * (2 + 2 * C * (y_R - 1) - y_R) * np.log(1 - 2 * C)
    #Gravitational constant in these units (so that everything comes out dimensionless)
    G = 6.674*10**(-33)
    k2 = (8/5) * C**5 * term1 / (term2 + term3 + term4)
    
    #The tidal deformability λ (sometimes called lambda)
    #made dimensionless by the 10⁻³⁶ factor here.
    Lambda = 2/3*k2*((R_star*5)/G)*10**(-36) #λ
    
    return k2, Lambda


def compute_all_params(R_arr, M_arr, yr_arr):
            """
    Compute the quadrupolar Love number k2 and tidal deformability Lambda
    for a set of neutron star models.

        Parameters
    ----------
    R_km_arr :  Array of stellar radii in kilometers.
    M_solar_arr : Array of stellar masses in solar masses.
    y_R_arr :  Array of metric perturbation values at the stellar surface.

    Returns
    -------
    k2_arr :  Quadrupolar Love numbers for each model.
    Lambda_arr : Dimensionless tidal deformabilities for each model.
    C_arr : Compactness values (GM/(Rc^2)) for each model.
    """
    k2_list = []
    lambda_list = []
    compactness_list = []
    for y, R, M in zip(yr_arr, R_arr, M_arr):
        C = 1.474 * M / R
        k2, Lambda = compute_k2_y_lambda(y, R, M)
        k2_list.append(k2)
        lambda_list.append(Lambda)
        compactness_list.append(C)
    return np.array(k2_list), np.array(lambda_list), np.array(compactness_list)

for func, dfunc, label in zip(EoS_functions, dEoS_functions, labels):
        """
    For each equation of state:
      - Solve the TOV equations over a range of central pressures (Pc_values)
      - Collect valid (R_star, M_star, y_R) outputs
      - Compute k2, Lambda, and C arrays
      - Store everything under the eos_label key
    """
    R_values, M_values, yr_values = [], [], []
    for Pc in Pc_values:
        ## Solve the TOV equations
        R_star, M_star, y_R = solve_tov(Pc, func, dfunc)
        if R_star is not None and M_star is not None:
            R_values.append(R_star)
            M_values.append(M_star)
            yr_values.append(y_R)

    R_arr = np.array(R_values)
    M_arr = np.array(M_values)
    yr_arr = np.array(yr_values)
    # Compute all tidal parameters
    k2_arr, lambda_arr, C_arr = compute_all_params(R_arr, M_arr, yr_arr)
  # Store in results dict
    eos_results[label] = {
        'R': R_arr,
        'M': M_arr,
        'yr': yr_arr,
        'k2': k2_arr,
        'Lambda': lambda_arr,
        'C': C_arr
    }

# Save ALL results at once after the loop:
np.savez('eos_results.npz', **eos_results)

print("All results saved to eos_results.npz")


# In[ ]:


import matplotlib.pyplot as plt
#The MR diagram for our range of hadronic EOS
colors = plt.cm.tab10.colors
line_styles = ['-', '--', '-.', ':']

plt.figure(figsize=(8,6), dpi=300)
ax = plt.gca()
#Loop over each EoS model and plot Radius vs. Mass
for i, (label, data) in enumerate(eos_results.items()):
    ax.plot(
        data['R'], data['M'],
        color=colors[i % 10],
        linestyle=line_styles[i % len(line_styles)],
        linewidth=2,
        label=label,
    )

# Format plot
ax.set_xlabel(r"$R\ (\mathrm{km})$", fontsize=14)
ax.set_ylabel(r"$M\ (M_\odot)$", fontsize=14)
ax.set_title(r"Mass–Radius Diagram for Hadronic EoS Models", fontsize=16, pad=15)
ax.set_xlim(8, 20)
ax.set_ylim(bottom=0)
ax.grid(which='major', linestyle='--', linewidth=0.6)
ax.grid(which='minor', linestyle=':',  linewidth=0.4)
ax.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=12, length=6)
ax.tick_params(axis='both', which='minor', labelsize=10, length=3)
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10,
          frameon=True, title="EoS Models", title_fontsize=12)
plt.tight_layout()
plt.savefig('MR_Diagram.pdf', bbox_inches='tight')
plt.show()


# In[ ]:


#--------THE TIDAL PARAMS. DIAGRAMS-------------

colors = plt.cm.tab10.colors # Global style settings
line_styles = ['-', '--', '-.', ':']

def plot_quantity_vs_quantity(xkey, ykey, xlabel, ylabel, title, fname, xlim=None, ylim=None):
    '''Parameters:

        -xkey, ykey: strings that index into each data dict (e.g. 'R', 'k2').

        -xlabel, ylabel, title: the text for axes and plot title .
    
        -filename: output filename (e.g. "k2_vs_radius.pdf").
    
        -xlim, ylim: optional 2‐tuples like (xmin, xmax)
    We loop over the eos_results dict, plotting the chosen quantities with 
    a unique color/style and adding a legend entry, while labeling everything
        '''
    plt.figure(figsize=(8,6), dpi=300)
    ax = plt.gca()
    for i, (label, data) in enumerate(eos_results.items()):
        ax.plot(
            data[xkey], data[ykey],
            color=colors[i % 10],
            linestyle=line_styles[i % len(line_styles)],
            linewidth=2,
            label=label,
        )
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16, pad=15)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(which='major', linestyle='--', linewidth=0.6)
    ax.grid(which='minor', linestyle=':', linewidth=0.4)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=12, length=6)
    ax.tick_params(axis='both', which='minor', labelsize=10, length=3)
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10,
              frameon=True, title="EoS Models", title_fontsize=12)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    plt.show()

# Now create all  diagrams in this consistent style:
plot_quantity_vs_quantity('R',    'k2', r"Radius $R$ (km)",         r"Love number $k_2$",  r"$k_2$ vs. Radius",        "k2_vs_radius.pdf",      xlim=(8, 20))
plot_quantity_vs_quantity('M',    'k2', r"Mass $M\,(M_\odot)$",     r"Love number $k_2$",  r"$k_2$ vs. Mass",          "k2_vs_mass.pdf",        xlim=(0, 3))
plot_quantity_vs_quantity('C',    'k2', r"Compactness $C = M/R$",   r"Love number $k_2$",  r"$k_2$ vs. Compactness",   "k2_vs_compactness.pdf", xlim=(0, 0.4))

plot_quantity_vs_quantity('M',   'yr', r"Mass $M\,(M_\odot)$",      r"Dimensionless $y_R$", r"$y_R$ vs. Mass",         "yR_vs_mass.pdf",        xlim=(0, 3))
plot_quantity_vs_quantity('R',   'yr', r"Radius $R$ (km)",          r"Dimensionless $y_R$", r"$y_R$ vs. Radius",       "yR_vs_radius.pdf",      xlim=(8, 20))
plot_quantity_vs_quantity('C',   'yr', r"Compactness $C = M/R$",    r"Dimensionless $y_R$", r"$y_R$ vs. Compactness",  "yR_vs_compactness.pdf", xlim=(0, 0.4))

plot_quantity_vs_quantity('M', 'Lambda', r"Mass $M\,(M_\odot)$",    r"$\lambda$",           r"$\lambda$ vs. Mass",     "Lambda_vs_mass.pdf",    xlim=(0, 3))
plot_quantity_vs_quantity('R', 'Lambda', r"Radius $R$ (km)",        r"$\lambda$",           r"$\lambda$ vs. Radius",   "Lambda_vs_radius.pdf",  xlim=(8, 20))
plot_quantity_vs_quantity('C', 'Lambda', r"Compactness $C = M/R$",  r"$\lambda$",           r"$\lambda$ vs. Compactness", "Lambda_vs_compactness.pdf", xlim=(0, 0.4))

