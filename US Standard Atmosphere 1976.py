"""
Created on Wed Jun  5 11:30:58 2024

@author: Vidhi Arya
"""
import numpy as np
import matplotlib.pyplot as plt

# Defining constants
g0 = 9.80665
r0 = 6.3567 * (10**11)  # Earth's radius in meters
M0 = 28.9644*(10**-3)  # Molar mass of Earth's air
R0 = 8.3143 * (10**-3)  # Universal gas constant in kJ/(kmol*K)
P0 = 101325.0  # Standard atmospheric pressure in Pascals
Tmb = 288.15  # Standard temperature at sea level in Kelvin

Z = np.arange(0, 1001, 1)*1000  # Altitude array from 0 to 20 km, in meters

# Defining the reference table for gradients of the linearly segmented temperature-height profile
Hb = np.array([0, 11, 20, 32, 47, 51, 71,86,91,110,120,1000])*1000  # Base heights for atmospheric layers in meters
Lmb = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0,0.0,0.0,12.0,0.0])  # Lapse rates in degrees per kilometer

# Initializing arrays
height = np.zeros(len(Z))
Temp = np.zeros(len(Z))
press = np.zeros(len(Z))
dens = np.zeros(len(Z))

# Defining gravity as a function of geometric altitude
g = g0 * (r0 / (r0 + Z))**2

# Defining geopotential altitude
H = (r0 * Z) / (r0 + Z)
height = H  # Assigning calculated geopotential altitude to height array

for i in range(12):
    # Defining temp for 0-11km
    T = Tmb + (Lmb[0] * (H[i] - Hb[0])/1000)
    Temp[i] = T  # Assigning value to Temp array
    
    # Defining pressure for 0-11km
    P = P0 * ((Tmb / (Tmb + (Lmb[0] * (H[i] - Hb[0])/1000)))**((g0 * M0) / (R0 * Lmb[0])))
    press[i] = P  # Assigning value to press array


    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D
    
for i in range(12, 21): 
    # Defining temp for altitudes above 11km
    T = Temp[11] + (Lmb[1] * (H[i] - Hb[1])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[11] * (np.exp((-g0 * M0 * (H[i] - Hb[1])) / (R0*(10**3) * Temp[11])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D 

for i in range(21, 33): 
    # Defining temp for altitudes above 11km
    T = Temp[20] + (Lmb[2] * (H[i] - Hb[2])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[20] * ((Temp[20] / (Temp[20]  + (Lmb[2] * (H[i] - Hb[2])/1000)))**((g0 * M0) / (R0 * Lmb[2])))
    press[i] = P

    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D  

for i in range(33, 48): 
    # Defining temp for altitudes above 11km
    T = Temp[32] + (Lmb[3] * (H[i] - Hb[3])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[32] * ((Temp[32] / (Temp[32]  + (Lmb[3] * (H[i] - Hb[3])/1000)))**((g0 * M0) / (R0 * Lmb[3])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D

for i in range(48, 52): 
    # Defining temp for altitudes above 11km
    T = Temp[47] + (Lmb[4] * (H[i] - Hb[4])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[47] * (np.exp((-g0 * M0 * (H[i] - Hb[4])) / (R0*(10**3) * Temp[47])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D 

for i in range(52, 72): 
    # Defining temp for altitudes above 11km
    T = Temp[51] + (Lmb[5] * (H[i] - Hb[5])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[51]* ((Temp[51] / (Temp[51]  + (Lmb[5] * (H[i] - Hb[5])/1000)))**((g0 * M0) / (R0 * Lmb[5])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D

for i in range(72, 87): 
    # Defining temp for altitudes above 11km
    T = Temp[71] + (Lmb[6] * (H[i] - Hb[6])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[71] * ((Temp[71] / (Temp[71]  + (Lmb[6] * (H[i] - Hb[6])/1000)))**((g0 * M0) / (R0 * Lmb[6])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D

for i in range(87, 92): 
    # Defining temp for altitudes above 11km
    T = Temp[86] + (Lmb[7] * (Z[i] - Z[86])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[86] * (np.exp((-g0 * M0 * (H[i] - Hb[7])) / (R0*(10**3) * Temp[86])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D

for i in range(91, 111): 
    # Defining temp for altitudes above 11km
    T = 263.1905 + -76.3232 * (1-((Z[i] - Z[91])/(-19.9429*1000)**2))**0.5
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[91] * (np.exp((-g0 * M0 * (H[i] - Hb[8])) / (R0*(10**3) * Temp[91])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D

for i in range(110, 121): 
    # Defining temp for altitudes above 11km
    T = Temp[110] + (Lmb[9] * (Z[i] - Z[110])/1000)
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[110] * ((Temp[110] / (Temp[110]  + (Lmb[9] * (H[i] - Hb[9])/1000)))**((g0 * M0) / (R0 * Lmb[9])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D
    
for i in range(120, 1001): 
    # Defining temp for altitudes above 11km
    zai = (Z[i]-Z[120])*((r0+Z[120])/(r0 +Z[i]))
    T = 1000 - (1000 - Temp[120])*(np.exp(-0.01875*zai))
    Temp[i] = T
    
    # Defining pressure for altitudes above 11km
    P = press[120] * (np.exp((-g0 * M0 * (H[i] - Hb[10])) / (R0*(10**3) * Temp[120])))
    press[i] = P
    
    # Calculating density
    D = (P * M0) / (R0*(10**3) * T) 
    dens[i] = D


figure, axis = plt.subplots(2, 2)
plt.suptitle('US Standard altmosphere',fontsize=16, weight = 'bold')

axis[0, 0].plot(dens,Z, color = 'k',marker = 'o',markersize=4)
axis[0, 0].set_title("Altitude vs density", weight = 'bold') 
axis[0, 0].set_xlabel('Density')
axis[0, 0].set_ylabel('Altitude')

axis[0, 1].plot(Temp,Z, color = 'k',marker = 'o',markersize=4)
axis[0, 1].set_title("Altitude vs Temprature", weight = 'bold') 
axis[0, 1].set_xlabel('Temprature')
axis[0, 1].set_ylabel('Altitude')

axis[1, 0].plot(press,Z, color = 'k',marker = 'o',markersize=4)
axis[1, 0].set_title("Altitude vs Pressure", weight = 'bold') 
axis[1, 0].set_xlabel('Pressure')
axis[1, 0].set_ylabel('Altitude')

axis[1, 1].plot(Temp,Z, color = 'k',marker = 'o',markersize=4)
axis[1, 1].plot(dens,Z, color = 'r',marker = 'o',markersize=4)
axis[1, 1].plot(press,Z, color = 'b',marker = 'o',markersize=4)
axis[1, 1].set_title("Altitude vs parameters", weight = 'bold') 
axis[1, 1].set_xlabel('parameters')
axis[1, 1].set_ylabel('Altitude')

plt.show() 
