# Open csv file with flash data
# For each flash
    # Save the coordinates
    # Look at direction and open corresponding galaxy-csv file
    # Find the galaxy with coordinates closest to the flash's coords

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_Hubble(velocity, distance):
    #plt.xlim([0,0.0003])
    plt.plot(distance,velocity,'.',color='C2', label='Data')
    A = np.vander(distance,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
    b, residuals, rank, s = np.linalg.lstsq(A,velocity, rcond=None)
    reconstructed = A @ b # @ is shorthand for matrix multiplication in python
    plt.plot(distance,reconstructed,'-r',label='Slope')
    plt.legend()
    plt.xlabel('Distance')
    plt.ylabel('Velocity')
    plt.grid()
    plt.show()

# Open csv file with flash data
datapath = os.getcwd() + "\\data"
flashdata = pd.read_csv(datapath + f'\\Flash_Data.csv', delimiter=',') 
flash_velocity_data = []
# For each flash
for flash in flashdata.values:
    photon_count = flash[4]
    # Save the coordinates
    x = flash[2]
    y = flash[3]
    #print(f"1: {np.power(d,2)}")
    #print(f"2: {flash[4]}")
    #print(f"3: {np.sqrt(1/max_xRay_photoncount)}")
    dist_list.append(np.sqrt(np.power(d,2)*flash[4])/max_xRay_photoncount)

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
    closest_galaxy_x = galaxydata['X'][0]
    closest_galaxy_y = galaxydata['Y'][0]
    closest_galaxy_name = galaxydata['Name'][0]
    
    closest_galaxy_velocity = galaxydata['RadialVelocity'][0]
    min_dist = np.sqrt(np.power((galaxydata['X'][0]-x), 2)+np.power((galaxydata['Y'][0]-y), 2))
    for galaxy in galaxydata.values:
        if np.sqrt(np.power((galaxy[1]-x), 2)+np.power((galaxy[2]-y), 2)) < min_dist:
            min_dist = np.sqrt(np.power((galaxy[1]-x), 2)+np.power((galaxy[2]-y), 2))
            closest_galaxy_velocity = galaxy[7]
            closest_galaxy_name = galaxy[0]
            closest_galaxy_x = galaxy[1]
            closest_galaxy_y = galaxy[2]
            #print(f"Its x: {galaxy[1]}")
            #print(f"Its y: {galaxy[2]}")
            #print(f"flash x: {x}")
            #print(f"flash y: {y}")
    flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1]}
    #flash_velocity_df.append([flash[0], closest_galaxy_velocity, closest_galaxy_x, closest_galaxy_y, photon_count])
    flash_velocity_df.sort_values(by=['Photon_count'], inplace=True)
    flash_velocity_df['Distance'] = dist_list

velocity = flash_velocity_df['Radial_velocity'].tolist()
distance = flash_velocity_df['Distance'].tolist()
#print(velocity)
print(distance)
plot_Hubble(velocity, distance)
#print(flash_velocity_df)
#print(len(flash_velocity_df))
           
     
    

