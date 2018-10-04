
"""
This script loops through the python dictionary and 
calculates zb and dune erosion
volumes for all profiles. 

The output is a time series of profile evolution and erosion volumes.
This is compared to measured zb and dV

"It is easier to write a new code than to understand an old one"
-John von Neumann to Marston Morse, 1952

EBG Aug. 2018 (Matlab); Oct 2018 (Python)

    
"""

import pickle

with open('DIM_data.pkl', 'rb') as f:
    data = pickle.load(f)
    
#loop through those indicies
for k in dict:
    profile=data[k]
    
    #pull all relevant data out of the dictionary
    dv=profile['dv']
    zb=profile['zb']
    T=profile['Tp']
    R_st=profile['R_st']
    R_gp=profile['R_gp']
    R_gp_draws=profile['R_gp_draws']
