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
    #plt.xlim([0.8,1.4])
    #velocity = [abs(v) for v in velocity]
    #distance = [abs(d) for d in distance]

    plt.plot(distance,velocity,'.',color='C2', label='Data')
    A = np.vander(distance,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
    b, residuals, rank, s = np.linalg.lstsq(A,velocity, rcond=None)
    reconstructed = A @ b # @ is shorthand for matrix multiplication in python
    print(f"b = {b}")
    plt.plot(distance,reconstructed,'-r',label='Slope')
    #plt.plot(np.log10(distance), 1421*np.log10(distance)-2125, '-r',label='test')
    #plt.plot(np.log10(distance), -2125*np.log10(distance)+1421, '-r',label='test')
    plt.legend()
    plt.xlabel('Distance [Mpc]')
    plt.ylabel('Velocity [km/s]')
    plt.savefig('hubble.png') 
    plt.grid()
    #plt.savefig('hubble.png')  
    plt.show()

# Open csv file with flash data
datapath = os.getcwd() + "/data"
flashdata = pd.read_csv(datapath + f'/Flash_Data.csv', delimiter=',') 
flash_velocity_df = pd.DataFrame(columns=['Flash_name','Radial_velocity','Galaxy_name','Galaxy_X','Galaxy_Y','Photon_count', 'Direction', 'Distance'])
#dist_list = []
#d= 1778#0.00018073525133440103# From variable_stars.py (= dist to galaxy with largest x-ray)
#max_xRay_photoncount=20089120
d=0.374630
max_xRay_photoncount=331
# For each flash
for flash in flashdata.values:
    photon_count = flash[4]
    # Save the coordinates
    x = flash[2]
    y = flash[3]
    #dist_list.append(np.sqrt((float(np.power(d,2))*flash[4])/max_xRay_photoncount))

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
    print(f"A: {float(np.power(d,2))}")
    print(f"B: {float(flash[4])}")
    print(f"C: {float(max_xRay_photoncount)}")
    print(f"D: {(np.sqrt(((np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount)))}")
    np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
    flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':np.sqrt(2.165*(10**(33))/((4*np.pi*(1.9865*(10**(-16))*flash[4]))/0.051))*3.24077929*10**(-23)} #np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
    #flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))} #np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
    
    
    #flash_velocity_df.append([flash[0], closest_galaxy_velocity, closest_galaxy_x, closest_galaxy_y, photon_count])
    flash_velocity_df.sort_values(by=['Photon_count'], inplace=True)



# Dropping last n rows using drop
flash_velocity_df.drop(flash_velocity_df.tail(9).index,inplace = True)
print(flash_velocity_df)
velocity = flash_velocity_df['Radial_velocity'].tolist()
distance = flash_velocity_df['Distance'].tolist()
#print(velocity)
#print(distance)
plot_Hubble(velocity, distance)
#print(flash_velocity_df)
#print(len(flash_velocity_df))
           
     
    

