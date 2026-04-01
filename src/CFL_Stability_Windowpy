#!/usr/bin/env python
# coding: utf-8
#CFL Stability Window
# In[7]:


import numpy as np
import matplotlib.pyplot as plt

# Constants
m_s    = 95.0        # MeV,Strange quark mass
m_n    = 939.565     # MeV, Nucleon mass
B_min  = 57.0        # MeV fm^-3,  Minimum bag constant for stability (MeV fm⁻³)
hbarc  = 197.326     # MeV·fm, to convert natural units

# Δ range
Delta = np.linspace(40, 160, 400)

# Natural-unit B(Δ) from Eq. (11) [MeV^4]
B_nat4 = (
    - (m_s**2 * m_n**2) / (12 * np.pi**2)
    + (Delta**2 * m_n**2) / (3 * np.pi**2)
    + (m_n**4)        / (108 * np.pi**2)
)

# Convert to MeV fm^-3 from MeV⁴ by dividing by (ℏc)³
B_curve = B_nat4 / (hbarc**3)

# CFL models from Table I
models = {
    'CFL-1':   (50,  60),
    'CFL-2':   (100, 60),
    'CFL-3':   (150, 60),
    'CFL-4':   (50,  80),
    'CFL-5':   (100, 80),
    'CFL-6':   (150, 80),
    'CFL-7':   (60,  100),
    'CFL-8':   (100, 100),
    'CFL-9':   (150, 100),
    'CFL-10':  (100, 120),
    'CFL-11':  (150, 120),
    'CFL-12':  (120, 140),
    'CFL-13':  (150, 140),
    'CFL-14':  (140, 160),
}

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 6), dpi=100)

# Plot stability boundary
ax.plot(Delta, B_curve, lw=2, color='teal', label=r'$m_s = 95\ \mathrm{MeV}$')

# Plot bag constant line
ax.axhline(B_min, ls='--', lw=1.5, color='red', label=r'$B = 57\ \mathrm{MeV\,fm}^{-3}$')

# Shade stable region
ax.fill_between(Delta, B_min, B_curve, where=(B_curve >= B_min),
                color='lavender', alpha=0.5, label='Stable CFL quark matter')

# Scatter & annotate models
for name, (d, b) in models.items():
    ax.plot(d, b, 'o', markersize=6)
    ax.text(d + 2, b + 2, name, fontsize=8)

# Improve appearance
ax.set_xlabel(r'$\Delta\ (\mathrm{MeV})$', fontsize=12)
ax.set_ylabel(r'$B\ (\mathrm{MeV\,fm}^{-3})$', fontsize=12)
ax.set_title('CFL Stability Window', fontsize=14)
ax.grid(which='major', linestyle='-', linewidth=0.5)
ax.grid(which='minor', linestyle=':', linewidth=0.5)
ax.minorticks_on()

# Show full frame
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

# Legend with frame
ax.legend(frameon=True, framealpha=0.9, edgecolor='gray', fontsize=10, loc='upper left')

# Save to file
fig.savefig('cfl_stability_window.png', dpi=300, bbox_inches='tight')

plt.tight_layout()
plt.show()


# In[ ]:




