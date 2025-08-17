# -*- coding: utf-8 -*-
#Created on Tue Aug 12 22:21:34 2025

#@author: emily
#credits to 
# https://www.researchgate.net/publication/349961642_Mathematical_Models_Describing_Biological_Systems_Under_Inhibitory_Conditions

import numpy as np
from scipy.optimize import curve_fit

temperatures = []
K_m = []
def monod(s, x, a, mu_max, mu_temp, K_m, V_max):
    # model function
    return (s*(mu_max-mu_temp))/mu_temp

def main():
    # model parameters 
    V_max= 2.3982e+001
    K_m= 1.2706e+002

#dataset



