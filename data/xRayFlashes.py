# Open csv file with flash data
# For each flash
    # Save the coordinates
    # Look at direction and open corresponding galaxy-csv file
    # Find the galaxy with coordinates closest to the flash's coords

import os 
import numpy as np
import pandas as pd

# Open csv file with flash data
datapath = os.getcwd() + "/data"
flashdata = pd.read_csv(datapath + f'/Flash_Data.csv', delimiter=',') 
flash_velocity_data = []
# For each flash
for flash in flashdata.values:
    # Save the coordinates
    x = flash[2]
    y = flash[3]
    # Look at direction and open corresponding galaxy-csv file
    if flash[1] == 'Top':
        galaxydata = pd.read_csv(datapath + f'/Top/Distant_Galaxy_Data.csv', delimiter=',') 
    elif flash[1] == 'Bottom':
        galaxydata = pd.read_csv(datapath + f'/Bottom/Distant_Galaxy_Data.csv', delimiter=',') 
    elif flash[1] == 'Right':
        galaxydata = pd.read_csv(datapath + f'/Right/Distant_Galaxy_Data.csv', delimiter=',') 
    elif flash[1] == 'Left':
        galaxydata = pd.read_csv(datapath + f'/Left/Distant_Galaxy_Data.csv', delimiter=',') 
    elif flash[1] == 'Back':
        galaxydata = pd.read_csv(datapath + f'/Back/Distant_Galaxy_Data.csv', delimiter=',') 
    elif flash[1] == 'Front':
        galaxydata = pd.read_csv(datapath + f'/Front/Distant_Galaxy_Data.csv', delimiter=',') 
    else:
        print(f"No match. flash direction = {flash[1]}")
    # Find the galaxy with coordinates closest to the flash's coords
    closest_galaxy_velocity = galaxydata['RadialVelocity'][0]
    min_dist = np.sqrt(np.power((galaxydata['X'][0]-x), 2)+np.power((galaxydata['Y'][0]-y), 2))
    for galaxy in galaxydata.values:
        if np.sqrt(np.power((galaxy[1]-x), 2)+np.power((galaxy[2]-y), 2)) < min_dist:
            min_dist = np.sqrt(np.power((galaxy[1]-x), 2)+np.power((galaxy[2]-y), 2))
            closest_galaxy_velocity = galaxy[7]
            #print(f"Its x: {galaxy[1]}")
            #print(f"Its y: {galaxy[2]}")
            #print(f"flash x: {x}")
            #print(f"flash y: {y}")
    flash_velocity_data.append([flash[0], closest_galaxy_velocity])
    print(flash_velocity_data)
           
     
    

