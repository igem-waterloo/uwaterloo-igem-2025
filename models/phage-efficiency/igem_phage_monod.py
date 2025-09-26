"""
Created on Sun Aug 17 23:21:02 2025

@author: angelinachen. Revised by Emily Wang, August 18, 2025
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.DataFrame({
    
    "Strains"
: ["4°C", "18°C", "20°C", "26°C", "30°C"],
    "α": [4.2020e-01, 4.3590e-01, 1.8510e-01, 3.7200e-01, 3.7800e-02],
    "μ*": [1.5803e+01, 2.8733e+01, 1.9318e+01, 2.0395e+01, 2.4604e+02],
    "K1": [4.9630e+02, 5.4636e+01, 1.2850e+02, 3.8677e+02, 4.9429e+01],
    "K2": [8.2411e+00, 6.6885e+00, 3.7309e+01, 2.7090e+01, 4.0856e+01],
    "s0": [5.6840e-02, 4.4870e-01, 1.9650e-01, 7.9280e-01, 7.4600e-02],
    "x0": [3.5400e-02, 1.2200e-02, 1.1000e-02, 1.9500e-02, 1.8500e-02]
})
temperatures = ["4°C", "18°C", "20°C", "26°C", "30°C"]

def monod(s, x, a, mu_max, mu_temp, K_m, V_max):
    # model function
    K_m = s*(mu_max-mu_temp)/mu_temp
    return K_m

def main():
    K_m = []
    # model parameters
    for i, row in data.iterrows(): 
       mu_max = row["μ*"]
       #assume mu is approximately a
       mu_temp = row["α"]
       # K1 = row["K1"]
       # K2 = row["K2"] 
       s0 = row["s0"] 
       # x0  = row["x0"] 
       new_KM = s0*(mu_max-mu_temp)/mu_temp
       K_m.append(new_KM) 
    return K_m 

results = main() 


# convert temperatures to numbers for plotting
temperatures = [int(t.replace("°C","")) for t in data["Strains"]]

#plot K_m vs temperature
x = temperatures
y = results
plt.plot(x,y)
plt.xlabel("Temperature(°C)")
plt.ylabel("K_m")
plt.show()
