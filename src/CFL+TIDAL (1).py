#!/usr/bin/env python
# coding: utf-8

# In[1]:



import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from scipy.misc import derivative
from scipy.integrate import solve_ivp


# In[2]:


#TOV_solver_from prev
def tov_equations(r, y,  energy_density_func) :
    P, M, yr = y
    epsilon_r = energy_density_func(P)# Energy density from the equation of state
    dedp = derivative(energy_density_func, P, dx=1e-6) 
    dedp = derivative(energy_density_func, P, dx=1e-6) 
    # Modified TOV equations
    F = (1-1.474 * 11.2 * (10 **(-6))*(r**2) *(epsilon_r-P))*(1/(1-2.948 *M/r))
       
    J = (
    1.474 * 11.2e-6 * (r ** 2) * (5 * epsilon_r + 9 * P + (epsilon_r + P) * dedp) * (1 / (1 - 2.948 * M / r))
    - 6 * (1 / (1 - 2.948 * M / r))
    - 4 * ((1.474 ** 2) * (M ** 2) / (r ** 2))
      * ((1 + 11.2e-6 * (r ** 3) * (P / M)) ** 2)
      * ((1 - 2.948 * M / r) ** -2)
    )  ##J=r^2Q


    # Modified TOV equations
    dP_dr = -1.474 * (epsilon_r * M * (1 + P / epsilon_r)) / r**2 *             (1 + 11.2e-6 * r**3 * P / M) * (1 - 2.948 * M / r)**-1
    dM_dr = 11.2e-6*r**2 * epsilon_r

    dyrdt= (-yr**2-yr*F-J)/r

    return [dP_dr, dM_dr,dyrdt]


# In[3]:



def tov_initial_conditions(Pc):
    M_c = 1e-3  # Small initial mass at the center
    yc = 2
    return [Pc, M_c, yc]


# In[4]:


# Define the CFL EoS function
def cfl_eos(P, b, d):#6
    """
    Compute the equation of state for CFL quark stars with given parameters.
    
    Parameters:
    P : 
        Pressure values.
    b :
        Bag constant.
    d :
        Gap parameter.
    ms : 
        Strange quark mass = 95MeV

    Returns:
    
        Energy density ε(P).
    """
    hbarc = 197.327
    ms = 95
    a = -ms**2 / 6 +  2 * d**2 / 3
    mm =  -3 * a + np.sqrt(4/3 * np.pi**2 * (b+P)*hbarc**3  + 9 * a**2) #7μ^2
    return 3 * P + 4 * b - (9 * a * mm) / ( np.pi**2*hbarc**3) #we multiply by hbarc^3 so the values are in MeVfm_3

# Define parameter sets for different EoS models table 1
parameter_sets = [
    {"b": 60, "d": 50},
    {"b": 60, "d": 100},
    {"b": 60, "d": 150},
    {"b": 80, "d": 50,},
    {"b": 80, "d": 100,},
    {"b": 80, "d": 150,},
    {"b": 100, "d": 60,},
    {"b": 100, "d": 100},
    {"b": 100, "d": 150},
    {"b": 120, "d": 100},
    {"b": 120, "d": 150},
    {"b": 140, "d": 120},
    {"b": 140, "d": 150},
    {"b": 160, "d": 140},
    ]

# Generate labels for each parameter set
labels = [f"cfl{i+1}" for i in range(len(parameter_sets))]

# Generate line styles for plotting
line_styles = ['-', '--']

# Create a function for each EoS using lambda
EoS_functions = [
    lambda P, b=b, d=d: cfl_eos(P, b, d)
    for params in parameter_sets
    for b, d in [(params['b'], params['d'])]
]


# In[5]:



# Set up central pressures and prepare for solving TOV equations
Pc_values = np.linspace(1, 1200, 2000)  # Central pressures in MeV/fm^3


# In[6]:


# --- Constants  ---
hbarc = 197.327  # MeV·fm
ms    = 95.0     # MeV
G     = 1.474    # km * c^2 / M_sun

# Analytic sound-speed derivative (used in dydr)
def dpdE(E, B, Delta):
    a = -ms**2/6 + 2*Delta**2/3
    return 1/3 + (2/3)*a/np.sqrt(a*a + (4/9)*np.pi**2*(E - B)*hbarc**3)

#We use the tov_equations functions from previus

# Event to stop integration when pressure crosses zero (surface)
def pressure_zero(r, y):
    return y[0]  # y[0] is pressure P(r)
pressure_zero.terminal = True
pressure_zero.direction = -1
"""The solve_tov function is like the previus one but now due to the surface jump for Self-bound configurations 
    we have to add a correction term on yr
"""


def solve_tov(Pc, energy_density_func):
    # Initial conditions at small r0
    r0 = 1e-6
    epsilon_c = energy_density_func(Pc)
    M0 = 11.2e-6 * epsilon_c * r0**3
    y0 = 2.0
    y_init = [Pc, M0, y0]

    # Integrate until P(r)=0
    sol = solve_ivp(
        fun=lambda r, y: tov_equations(r, y, energy_density_func),
        t_span=(r0, 30.0),
        y0=y_init,
        events=pressure_zero,
        max_step=0.05
    )

    # Extract surface values from event if triggered
    if sol.t_events[0].size > 0:
        R_star = sol.t_events[0][0]
        P_surf, M_star, y_R = sol.y_events[0][0]
    else:
        # fallback to last integration point
        R_star = sol.t[-1]
        P_surf, M_star, y_R = sol.y[:, -1]

    # Compute corrected y_R and tidal Love number as before...
    epsilon_surf = energy_density_func(0.0)
    yR_corr = y_R - 11.2e-6 * epsilon_surf * R_star**3 / M_star
    C = G * M_star / R_star
    term1 = (1 - 2*C)**2 * (2 + 2*C*(yR_corr - 1) - yR_corr)
    term2 = 2*C*(6 - 3*yR_corr + 3*C*(5*yR_corr - 8))
    term3 = 4*C**3*(13 - 11*yR_corr + C*(3*yR_corr - 2) + 2*C**2*(1 + yR_corr))
    term4 = 3*(1 - 2*C)**2*(2 - yR_corr + 2*C*(yR_corr - 1)) * np.log(1 - 2*C)
    k2 = (8/5)*C**5 * term1 / (term2 + term3 + term4)
    Lambda = (2/3) * k2 * R_star**5*10**(-4) #λ multiplied byt 10^(-4) 

    return R_star, M_star, yR_corr, k2, C, Lambda


# In[7]:


# --- Store M-R data for each CFL EoS in a dictionary ---
cfl_results = {}
for func, label in zip(EoS_functions, labels):
    R_values, M_values,yr_values, k2_values, beta_values, Lambda_values = [], [],[],[],[],[]
    for Pc in Pc_values:
        R_star, M_star, y_Rf, k2, beta, Lambda = solve_tov(Pc, func)
        if R_star is not None and M_star is not None:
            R_values.append(R_star)
            M_values.append(M_star)
            yr_values.append(y_Rf)
            k2_values.append(k2)
            beta_values.append(beta)
            Lambda_values.append(Lambda)
    # Store the results for plotting later.
    cfl_results[label] = {
        'R': np.array(R_values),
        'M': np.array(M_values),
        'yr': np.array(yr_values),
        'k2': np.array(k2_values),
        'beta': np.array(beta_values),
        'Lambda': np.array(Lambda_values)
    }
    


# In[ ]:


"""We now plot the MR diagram for the CFL stars alongside
    with the forbiden regions and observational constraints
"""
# Observational bands: (label, color, r_min, r_max, m_min, m_max)
obs_bands = [
    ('PSR J0740+6620', 'green', 12.0, 14.0, 2.00, 2.14),
    ('PSR J0952−0607', 'yellow', 6.0, 16.0, 2.18, 2.52),
    ('GW170817',       'purple', 10.5, 13.4, 1.20, 1.60),
    ('HESS J1731−347', 'cyan',   9.5, 11.3, 0.60, 1.00),
]

# Create plot
fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

# Plot CFL model curves (thick lines)
for label, data in cfl_results.items():
    ax.plot(data['R'], data['M'], lw=2.0, label=label)

# Overlay observational bands with light shading (no individual legend entry)
for name, color, rmin, rmax, mmin, mmax in obs_bands:
    ax.fill_betweenx([mmin, mmax], rmin, rmax, color=color, alpha=0.2)

# Forbidden regions: Black hole, Buchdahl, Causality
R_vals = np.linspace(0, 16, 500)
conv = 1.47664  # Solar mass in km (G=c=1)
M_bh   = R_vals / (2 * conv)           # R = 2M
M_buch = (4 * R_vals) / (9 * conv)     # R = 9/4 M
M_caus = (0.354 * R_vals) / conv       # M/R = 0.354

ymax = max(data['M'].max() for data in cfl_results.values()) * 1.05
ax.fill_between(R_vals, M_bh, ymax, color='black',    alpha=0.15)
ax.fill_between(R_vals, M_buch, M_bh, color='dimgray', alpha=0.15)
ax.fill_between(R_vals, M_caus, M_buch, color='lightgray', alpha=0.15)

# Formatting
ax.set_xlim(0, 16)
ax.set_ylim(0, ymax)
ax.set_xlabel('Radius (km)', fontsize=14)
ax.set_ylabel(r'Mass ($M_\odot$)', fontsize=14)
ax.set_title('CFL Quark Star Mass–Radius Diagram with Observational Bands', fontsize=16)

ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth=0.8, alpha=0.6)
ax.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.4)

# Clean spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend: curves + patches for bands and limits
handles, labels = ax.get_legend_handles_labels()
legend_patches = [
    Patch(facecolor='green', alpha=0.2, label='PSR J0740+6620'),
    Patch(facecolor='yellow', alpha=0.2, label='PSR J0952−0607'),
    Patch(facecolor='purple', alpha=0.2, label='GW170817'),
    Patch(facecolor='cyan', alpha=0.2, label='HESS J1731−347'),
    Patch(facecolor='black', alpha=0.15, label='Black hole limit'),
    Patch(facecolor='dimgray', alpha=0.15, label='Buchdahl limit'),
    Patch(facecolor='lightgray', alpha=0.15, label='Causality limit'),
]
all_handles = handles + legend_patches
ax.legend(handles=all_handles, loc='upper left', fontsize=10, ncol=2, frameon=True, framealpha=0.8)

plt.tight_layout()
fig.savefig('mr_diagram_styled.png', dpi=300, bbox_inches='tight')
plt.show()


# In[23]:


#The plots for the tidal deformability parameters
colors = plt.cm.tab10.colors

# Styling 
def style_axes(ax):
    ax.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth=0.8, alpha=0.6)
    ax.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# 1) Love number vs Radius
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['R'], data['k2'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Radius $R$ (km)', fontsize=14)
plt.ylabel('Love number $k_2$', fontsize=14)
plt.title('Love Number vs. Radius', fontsize=16)
plt.xlim(0,16)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('love_vs_radius.png', bbox_inches='tight')
plt.show()

# 2) Love number vs Mass
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['M'], data['k2'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Mass $M\\,(M_\\odot)$', fontsize=14)
plt.ylabel('Love number $k_2$', fontsize=14)
plt.title('Love Number vs. Mass', fontsize=16)
plt.xlim(0,3)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('love_vs_mass.png', bbox_inches='tight')
plt.show()

# 3) Love number vs Compactness
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['beta'], data['k2'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Compactness $\\beta = M/R$', fontsize=14)
plt.ylabel('Love number $k_2$', fontsize=14)
plt.title('Love Number vs. Compactness', fontsize=16)
plt.xlim(0,0.4)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('love_vs_compactness.png', bbox_inches='tight')
plt.show()

# 4) y_R vs Mass
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['M'], data['yr'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Mass $M\\,(M_\\odot)$', fontsize=14)
plt.ylabel('Dimensionless $y_R$', fontsize=14)
plt.title('$y_R$ vs. Mass', fontsize=16)
plt.xlim(0,3)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('yr_vs_mass.png', bbox_inches='tight')
plt.show()

# 5) y_R vs Radius
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['R'], data['yr'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Radius $R$ (km)', fontsize=14)
plt.ylabel('Dimensionless $y_R$', fontsize=14)
plt.title('$y_R$ vs. Radius', fontsize=16)
plt.xlim(0,16)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('yr_vs_radius.png', bbox_inches='tight')
plt.show()

# 6) y_R vs Compactness
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['beta'], data['yr'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Compactness $\\beta = M/R$', fontsize=14)
plt.ylabel('Dimensionless $y_R$', fontsize=14)
plt.title('$y_R$ vs. Compactness', fontsize=16)
plt.xlim(0,0.4)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('yr_vs_compactness.png', bbox_inches='tight')
plt.show()

# 7) Lambda vs Mass
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['M'], data['Lambda'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Mass $M\\,(M_\\odot)$', fontsize=14)
plt.ylabel('$\\lambda\\times10^{-4}\\,[\\mathrm{km}^5]$', fontsize=14)
plt.title('$\\lambda\\times10^{-4}$ vs. Mass', fontsize=16)
plt.xlim(0,3)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('Lambda_vs_mass.png', bbox_inches='tight')
plt.show()

# 8) Lambda vs Radius
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['R'], data['Lambda'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Radius $R$ (km)', fontsize=14)
plt.ylabel('$\\lambda\\times10^{-4}\\,[\\mathrm{km}^5]$', fontsize=14)
plt.title('$\\lambda\\times10^{-4}$ vs. Radius', fontsize=16)
plt.xlim(0,16)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('Lambda_vs_radius.png', bbox_inches='tight')
plt.show()

# 9) Lambda vs Compactness
plt.figure(figsize=(8,6), dpi=300)
for i, (label, data) in enumerate(cfl_results.items()):
    plt.plot(data['beta'], data['Lambda'], color=colors[i%10], lw=2, label=label)
plt.xlabel('Compactness $\\beta = M/R$', fontsize=14)
plt.ylabel('$\\lambda\\times10^{-4}\\,[\\mathrm{km}^5]$', fontsize=14)
plt.title('$\\lambda\\times10^{-4}$ vs. Compactness', fontsize=16)
plt.xlim(0,0.4)
style_axes(plt.gca())
plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('Lambda_vs_compactness.png', bbox_inches='tight')
plt.show()


# In[ ]:




