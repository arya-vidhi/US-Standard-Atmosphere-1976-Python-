"""
Created on Wed Jun  5 17:25:17 2024 2024

@author: Vidhi Arya
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
# Defining constants
g0 = 9.80665
r0 = 6.3567 * (10**11)  # Earth's radius in meters
M0 = 28.9644*(10**-3)  # Molar mass of Earth's air
R0 = 8.3143 * (10**-3)  # Universal gas constant in kJ/(kmol*K)
Tb = np.array([288.15, 216.6500012372771, 216.65000123727717, 228.6499996263785,270.6499898961733, 270.6499898961733, 214.6500121007738, 184.65003537070632, 186.8673, 186.8691230926484, 306.86912309264835])# Molecular scale temperature gradient
Pb = [101325, 22631.983047306938, 5474.850990135241, 868.0090831541909, 110.90457815919828, 66.93778253925619, 3.9563335745933146,0.30232174031073117, 0.11987083895166374, 0.0037166463303773345, 0.0009054783416967828]

#Z = int(input("enter the geometric height : ")) # Altitude array from 0 to 20 km, in meters

# Defining the reference table for gradients of the linearly segmented temperature-height profile
Hb = np.array([0, 11, 20, 32, 47, 51, 71,86,91,110,120,1000])*1000  # Base heights for atmospheric layers in meters
Zb = np.array([0, 11, 20, 32, 47, 51, 71,86,91,110,120,1000])*1000  # geometric Base heights for atmospheric layers in meters
Lmb = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0,0.0,0.0,12.0,0.0])  # Lapse rates in degrees per kilometer


def calculate_density(Z):
    H = (r0 * Z) / (r0 + Z)
    
    if Z <= Zb[1]:   
        T =  Tb[0] + (Lmb[0] * (H - Hb[0])/1000)
        P =  Pb[0] * ((Tb[0] / (Tb[0] + (Lmb[0] * (H - Hb[0])/1000)))**((g0 * M0) / (R0 * Lmb[0])))
        D = (P * M0) / (R0*(10**3) * T) 
 
    elif (Z > Zb[1]) and (Z <= Zb[2]):   
        T =  Tb[1] + (Lmb[1] * (H - Hb[1])/1000)
        P =  Pb[1] * (np.exp((-g0 * M0 * (H - Hb[1])) / (R0*(10**3) * Tb[1])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[2]) and (Z <= Zb[3]):   
        T =  Tb[2] + (Lmb[2] * (H - Hb[2])/1000)
        P = Pb[2] * ((Tb[2] / (Tb[2]  + (Lmb[2] * (H - Hb[2])/1000)))**((g0 * M0) / (R0 * Lmb[2])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[3]) and (Z <= Zb[4]):   
        T =  Tb[3] + (Lmb[3] * (H - Hb[3])/1000)
        P =  Pb[3] * ((Tb[3] / (Tb[3]  + (Lmb[3] * (H - Hb[3])/1000)))**((g0 * M0) / (R0 * Lmb[3])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[4]) and (Z <= Zb[5]): 
        T =  Tb[4] + (Lmb[4] * (H - Hb[4])/1000)
        P =  Pb[4] * (np.exp((-g0 * M0 * (H - Hb[4])) / (R0*(10**3) * Tb[4])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[5]) and (Z <= Zb[6]):  
        T =  Tb[5] + (Lmb[5] * (H - Hb[5])/1000)       
        P =  Pb[5]* ((Tb[5] / (Tb[5]  + (Lmb[5] * (H - Hb[5])/1000)))**((g0 * M0) / (R0 * Lmb[5])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[6]) and (Z <= Zb[7]):   
        T =  Tb[6] + (Lmb[6] * (H - Hb[6])/1000)
        P =  Pb[6] * ((Tb[6] / (Tb[6]  + (Lmb[6] * (H - Hb[6])/1000)))**((g0 * M0) / (R0 * Lmb[6])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[7]) and (Z <= Zb[8]):
        T = Tb[7] + (Lmb[7] * (Z - Zb[7])/1000)
        P = Pb[7] * (np.exp((-g0 * M0 * (H - Hb[7])) / (R0*(10**3) * Tb[7])))
        D = (P * M0) / (R0*(10**3) * T)
    
    elif (Z > Zb[8]) and (Z <= Zb[9]):
        T = 263.1905 + -76.3232 * (1-((Z - Zb[8])/(-19.9429*1000)**2))**0.5
        P = Pb[8] * (np.exp((-g0 * M0 * (H - Hb[8])) / (R0*(10**3) * Tb[8])))
        D = (P * M0) / (R0*(10**3) * T)
        
    elif (Z > Zb[9]) and (Z <= Zb[10]):
        T = Tb[9] + (Lmb[9] * (Z - Zb[9])/1000)
        P = Pb[9] * ((Tb[9] / (Tb[9]  + (Lmb[9] * (H - Hb[9])/1000)))**((g0 * M0) / (R0 * Lmb[9])))
        D = (P * M0) / (R0*(10**3) * T)

    else:
        zai = (Z-Zb[10])*((r0+Zb[10])/(r0 +Z))
        T = 1000 - (1000 - Tb[10])*(np.exp(-0.01875*zai)) 
        P = Pb[10] * (np.exp((-g0 * M0 * (H - Hb[10])) / (R0*(10**3) * Tb[10])))
        D = (P * M0) / (R0*(10**3) * T)
        
    return D

#print(calculate_density(5000))
Z = np.arange(1, 1000, 1)*1000
n = np.size(Z)
rhodist = np.zeros(n)

den_data = np.array([1.1117, 1.0066, 0.90912, 0.81913, 0.73612, 0.6597, 0.5895, 0.52517, 0.46635, 0.41271, 0.36392])
alt = np.arange(1,12,1)*1000
for j in range(n):
    rhodist[j] = calculate_density(Z[j])
    
df = pd.DataFrame({'alt(km)': Z, 'Density': rhodist})
plt.plot(rhodist, Z)
plt.plot(den_data, alt)