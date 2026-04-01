#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp


# In[2]:




# MDI Equations
def mdi1(P):
    if P>=0.184:
        return 4.1844 * P**0.81449 + 95.00135 * P**0.31736
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

def mdi2(P):
    if P>=0.184:
        return 5.97365 * P**0.77374 + 89.24 * P**0.30993
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

def mdi3(P):
    if P>=0.184:
        return 15.55 * P**0.666 + 76.71 * P**0.247
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

def mdi4(P):
    if P>=0.184:
        return 25.99587 * P**0.61209 + 65.62193 * P**0.15512
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# NLD Equation
def nld(P):
    if P>=0.184:
        return (119.05736  + 304.80445 * (1 - np.exp(-P / 48.61465))+ 33722.34448 * (1 - np.exp(-P / 17499.47411)))
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# HHJ Equations
def hhj1(P):
    if P>=0.184:
        return 1.78429 * P**0.93761 + 106.93652 * P**0.31715
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

def hhj2(P):
    if P>=0.184:
        return 1.18961 * P**0.96539 + 108.40302 * P**0.31264
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# Ska Equation
def ska(P):
    if P>=0.184:
        return 0.53928 * P**1.01394 + 94.31452 * P**0.35135
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# SkI4 Equation
def ski4(P):
    if P>=0.184:
        return 4.75668 * P**0.76537 + 105.7220 * P**0.2745
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# HLPS Equations


def hlps3(P):
    if P>=0.184:
        return (131.811 * (1 - np.exp(-P / 4.41577))+ 924.143 * (1 - np.exp(-P / 523.736))+ 81.5682)
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# SCVBB Equation
def scvbb(P):
    if P>=0.184:
        return 0.3714141 * P**1.08004 + 109.258 * P**0.351019
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))


# In[3]:


def energy_density_HLPS_3(P):
    if P>=0.184:
         return 131.811*(1-np.exp(-P/4.41577))+924.143*(1-np.exp(-P/523.736))+81.5682 #HLPS-3
       
    
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))


# In[4]:


import numpy as np

# WFF-1
def WFF_1(P):
    if P>=0.184:
        return 0.00127717 * P**1.69617 + 135.233 * P**0.331471
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# WFF-2
def WFF_2(P):
    if P>=0.184:
        return 0.00244523 * P**1.62962 + 122.076 * P**0.304401
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# PS
def PS(P):
    if P>=0.184:
        return 9805.95 * (1 - np.exp(-0.000193624 * P)) + 212.072 * (1 - np.exp(-0.401508 * P)) + 1.69483
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# W
def W(P):
    if P>=0.184:
        return 0.261822 * P**1.16851 + 92.4893 * P**0.307728
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# BGP
def BGP(P):
    if P>=0.184:
        return 0.0112475 * P**1.59689 + 102.302 * P**0.335526
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# BL-1
def BL_1(P):
    if P>=0.184:
        return 0.488686 * P**1.01457 + 102.26 * P**0.355095
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# BL-2
def BL_2(P):
    if P>=0.184:
        return 1.34241 * P**0.910079 + 100.756 * P**0.354129
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))
# DH
def DH(P):
    if P>=0.184:
        return 39.5021 * P**0.541485 + 96.0528 * P**0.00401285
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))

# APR-1
def APR_1(P):
    if P>=0.184:
        return 0.000719964 * P**1.85898 + 108.975 * P**0.340074
    elif P>= 9.34375e-5:
        # Equation (1.41)
        return 0.00873 + 103.17338 * (1 - np.exp(-P / 0.36227)) +7.34979 * (1 - np.exp(-P / 0.012721))
    
    elif P>=4.1725e-8:
        # Equation (1.42)
        return (0.00015 + 0.00203 * (1 - np.exp(-P / 3448.27)) + 0.10851 * (1 - np.exp(-P / 7692.3076)))
    
    elif P>=1.44875e-11:
        # Equation (1.43)
        return (0.000051 * (1 - np.exp(-P / (0.2373 * 1e10))) +
                0.00014 * (1 - np.exp(-P / (0.4020 * 1e8))))
    else:
       
       
        return 10 **(31.93753+10.82611 *np.log10(P)+1.29312 *(np.log10(P)**2)+0.08014 *(np.log10(P) **3)+0.00242 *(np.log10(P) **4)+0.000028 *(np.log10(P) **5))


# In[48]:


from scipy.misc import derivative

def der_energy_density(P, energy_density_func):
    """Compute dε/dP using SciPy's derivative function."""
    
    return derivative(energy_density_func, P, dx=1e-16)


# In[55]:


def tov_equations(r, y,  energy_density_func) :
    P, M, yr = y
    epsilon_r = energy_density_func(P)# Energy density from the equation of state
    dedp = derivative(energy_density_func, P, dx=1e-16) 
    dedp = derivative(energy_density_func, P, dx=1e-16) 
    # Modified TOV equations
    F = (1-1.474 * 11.2 * (10 **(-6))*(r**2) *(epsilon_r-P))*(1/(1-2.948 *M/r))
       
    J = (
    1.474 * 11.2e-6 * (r ** 2) * (5 * epsilon_r + 9 * P + (epsilon_r + P) * dedp) * (1 / (1 - 2.948 * M / r))
    - 6 * (1 / (1 - 2.948 * M / r))
    - 4 * ((1.474 ** 2) * (M ** 2) / (r ** 2))
      * ((1 + 11.2e-6 * (r ** 3) * (P / M)) ** 2)
      * ((1 - 2.948 * M / r) ** -2)
    )  ##J=r^2Q


    if P <= 1e-8:
        return [0, 0, 0]#op when pressure becomes zero

    # Modified TOV equations
    dP_dr = -1.474 * (epsilon_r * M * (1 + P / epsilon_r)) / r**2 *             (1 + 11.2e-6 * r**3 * P / M) * (1 - 2.948 * M / r)**-1
    dM_dr = 11.2e-6*r**2 * epsilon_r

    dyrdt= (-yr**2-yr*F-J)/r

    return [dP_dr, dM_dr,dyrdt]


# In[56]:


def tov_initial_conditions(Pc):
    M_c = 1e-3  # Small initial mass at the center
    yc = 2
    return [Pc, M_c, yc]


# In[ ]:





# In[64]:


def solve_tov(Pc,energy_density_func):
    # Initial conditions: central pressure and mass

    r0 = 1e-6
    epsilon_c = energy_density_func(Pc)
    M0 = 11.2e-6 * epsilon_c * r0**3
    y0 = 2.0
    y_init = [Pc, M0, y0]


    # Define the radial range for integration (we'll stop when pressure drops to ~0)
    r_max = 30  # Maximum radius to integrate out to [km]
    r_eval = np.linspace(0.1, r_max, 1000)

    # Solve TOV equations using solve_ivp ()
    solution = solve_ivp(lambda r, y: tov_equations(r, y, energy_density_func), 
                         t_span=(r0, 30.0),
        y0=y_init,
        max_step=0.05
    )
    
    radius = solution.t
    pressure = solution.y[0]
    mass = solution.y[1]
    y_values = solution.y[2] 
     # … after extracting pressure, radius, mass, y_values …

    # threshold in dyn/cm²
    threshold = 1e-4

    # find all radii where pressure < threshold
    surface_idxs = np.where(pressure < threshold)[0]

    if surface_idxs.size > 0:
        surface_index = surface_idxs[0]
    else:
        # fallback to last point if we never dropped below threshold
        surface_index = len(pressure) - 1

    R_star = radius[surface_index]
    M_star = mass[surface_index]
    y_R     = y_values[surface_index]
    beta =  (1.474*M_star) / R_star  # Compactness
    # Compute k2
    term1 = (1 - 2 * beta) ** 2 * (2 + 2 * beta * (y_R - 1) - y_R)
    term2 = 2 * beta * (6 - 3 * y_R + 3 * beta * (5 * y_R - 8))
    term3 = 4 * beta**3 * (13 - 11 * y_R + beta * (3 * y_R - 2) + 2 * beta**2 * (1 + y_R))
    term4 = 3 * (1 - 2 * beta)**2 * (2 + 2 * beta * (y_R - 1) - y_R) * np.log(1 - 2 * beta)
    k2 = (8/5) * beta**5 * term1 / (term2 + term3 + term4)
    Lambda = (2/3) * k2 * R_star**5*10**(-4) #converted lmabda
   

    return R_star, M_star, y_R, k2,beta, Lambda


# In[65]:


solve_tov(600,mdi1)


# In[ ]:



Pc_values = np.linspace(1, 1200, 2200)

# List of EoS functions and labels
EoS_functions = [mdi1]
labels = ["MDI-1"]

# Create a dictionary to store the results for each EoS
eos_results1 = {}

# Loop over each EoS function
for func, label in zip(EoS_functions, labels):
   # Initialize lists for the results
   R_values, M_values, yr_values, k2_values, beta_values, Lambda_values = [], [], [], [], [],[]
   
   # Loop over central pressure values
   for Pc in Pc_values:
       # solve_tov returns a tuple: (R_star, M_star, y_R, k2, beta, Lambda)
       R_star, M_star, y_R, k2, beta, Lambda = solve_tov(Pc, func)
       if R_star is not None and M_star is not None:
           R_values.append(R_star)
           M_values.append(M_star)
           yr_values.append(y_R)
           k2_values.append(k2)
           beta_values.append(beta)
           Lambda_values.append(Lambda)
   
   # Store the results for this EoS in the dictionary (converted to numpy arrays)
   eos_results1[label] = {
       'R': np.array(R_values),
       'M': np.array(M_values),
       'yr': np.array(yr_values),
       'k2': np.array(k2_values),
       'beta': np.array(beta_values),
       'Lambda': np.array(Lambda_values)
   }


# In[ ]:


# --- Plot 1: Mass vs. Radius (M-R) ---
plt.figure(figsize=(10, 6))
# Plot curves from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['R'], data['M'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Radius (km)")
plt.ylabel("Mass (M☉)")
plt.title("Mass vs. Radius")
plt.xlim(8, 20)  # Adjust as needed
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 2: k₂ vs. Radius (k2-R) ---
plt.figure(figsize=(10, 6))
# Plot eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['R'], data['k2'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Radius (km)")
plt.ylabel("k₂ (Love Number)")
plt.title("k₂ vs. Radius")
plt.xlim(8, 20)
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 3: k₂ vs. Mass (k2-M) ---
plt.figure(figsize=(10, 6))
# Plot eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['M'], data['k2'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Mass (M☉)")
plt.ylabel("k₂ (Love Number)")
plt.title("k₂ vs. Mass")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 4: k₂ vs. Compactness (k2-β) ---
plt.figure(figsize=(10, 6))
# Plot eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['beta'], data['k2'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Compactness (β = M/R)")
plt.ylabel("k₂ (Love Number)")
plt.title("k₂ vs. Compactness")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()



# --- Plot A: y_R vs. Mass ---
plt.figure(figsize=(10, 6))
# Plot from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['M'], data['yr'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Mass (M☉)")
plt.ylabel("y_R")
plt.title("y_R vs. Mass")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot B: y_R vs. Radius ---
plt.figure(figsize=(10, 6))
# Plot from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['R'], data['yr'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Radius (km)")
plt.ylabel("y_R")
plt.title("y_R vs. Radius")
plt.xlim(8, 20)  # Optional: limit Radius axis if desired
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot C: y_R vs. Compactness (β) ---
plt.figure(figsize=(10, 6))
# Plot from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['beta'], data['yr'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Compactness (β = M/R)")
plt.ylabel("y_R")
plt.title("y_R vs. Compactness")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 1: Lambda vs. Mass (Λ-M) ---
plt.figure(figsize=(10, 6))
# Plot curves from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['M'], data['Lambda'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Mass (M☉)")
plt.ylabel("λ")
plt.title("λ vs. M")
plt.xlim(0, 3)  # Adjust as needed
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 2: Lambda vs. Radius (Λ-M) ---
plt.figure(figsize=(10, 6))
# Plot curves from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['R'], data['Lambda'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("Radius (km)")
plt.ylabel("λ")
plt.title("λ vs. R")
plt.xlim(8, 20)  # Adjust as needed
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 2: Lambda vs. β (Λ-M) ---
plt.figure(figsize=(10, 6))
# Plot curves from eos_results1 (set 1)
for label, data in eos_results1.items():
    plt.plot(data['beta'], data['Lambda'], label=f"{label} (set 1)", linestyle='-')
plt.xlabel("beta")
plt.ylabel("λ")
plt.title("λ vs. β")
plt.xlim(0, 0.4)  # Adjust as needed
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.show()


# In[ ]:




