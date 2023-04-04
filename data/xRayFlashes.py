# Open csv file with flash data
# For each flash
    # Save the coordinates
    # Look at direction and open corresponding galaxy-csv file
    # Find the galaxy with coordinates closest to the flash's coords

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from numpy import *
from astropy.table import Table
import os 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.metrics import r2_score

def plotHubble(velocity, distance, uncertainty):
    #x = np.linspace(0,5.5,10)
    #y = 10*np.exp(-x)
    #xerr = np.random.random_sample(10)
    #yerr = np.random.random_sample(10)
    #fig, ax = plt.subplots()            #initialize axes
    #ax.set_ylabel('Radial Velocity (km/s)')
    #ax.set_xlabel('Distance from X-Ray Source (pc)')
    #ax.invert_xaxis()
    #plt.scatter(distance, velocity, s=3)
    """
    z,cov = np.polyfit(distance, velocity, 1, cov=True)     #this finds the linear fit for the data
    p = np.poly1d(z)
    gradUnc, intUnc = np.sqrt(np.diag(cov))       #this is the uncertainty in the fit
    #the following plots an uncertainty "area" for the trendline
    x = np.arange(min(distance), max(distance), 1000)
    upper = (z[0] + gradUnc) * x + (z[1] - intUnc)
    lower = (z[0] - gradUnc) * x + (z[1] + intUnc)
    plt.fill_between(x, lower, upper, color='r', alpha=0.2)

    plt.plot(x,p(x),"r--", linewidth=0.5)     #plot the trendline on top of the data
    plt.errorbar(distance, velocity, xerr=distUnc, yerr=0.5, fmt=',', linewidth=0.5)"""


    fig, ax = plt.subplots()
    plt.xlim([0,0.55])
    plt.ylim([-1100, 0])
    A = np.vander(distance,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
    b, residuals, rank, s = np.linalg.lstsq(A,velocity, rcond=None)
    reconstructed = A @ b # @ is shorthand for matrix multiplication in python
    ax.plot(distance, reconstructed, '-r')
    ax.errorbar(distance, velocity,
                xerr=[0]*len(distance),
                yerr=[0.3]*len(velocity),
                fmt='.', elinewidth=1)
    ax.set_xlabel('Distance [Mpc]')
    ax.set_ylabel('Radial Velocity [km/s]')
    ax.set_title('Hubble constant')
    plt.gca().legend(('LS Estimation','Galaxy Data'))
    plt.grid()
    plt.savefig('hubble.png') 
    plt.show()

    print(f"Hubble constant, m = {b}")

def plot_Hubble_old(velocity, distance, uncertainty):
    #plt.xlim([0.8,1.4])
    #velocity = [abs(v) for v in velocity]
    #distance = [abs(d) for d in distance]


    #plt.plot(distance,velocity,'.',color='C2', label='Data')
    A = np.vander(distance,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
    b, residuals, rank, s = np.linalg.lstsq(A,velocity, rcond=None)
    reconstructed = A @ b # @ is shorthand for matrix multiplication in python
    print(f"b = {b}")
    plt.plot(distance,reconstructed,'-r',label='Slope')
    #plt.plot(np.log10(distance), 1421*np.log10(distance)-2125, '-r',label='test')
    #plt.plot(np.log10(distance), -2125*np.log10(distance)+1421, '-r',label='test')
    plt.legend()
    #plt.errorbar(distance, reconstructed, uncertainty, [0]*len(velocity))
    plt.xlabel('Distance [Mpc]')
    plt.ylabel('Velocity [km/s]')
    plt.savefig('hubble?old.png') 
    plt.grid()
    #plt.savefig('hubble.png')  
    plt.show()

# Open csv file with flash data
datapath = os.getcwd() + "/data"
flashdata = pd.read_csv(datapath + f'/Flash_Data.csv', delimiter=',') 
flash_velocity_df = pd.DataFrame(columns=['Flash_name','Radial_velocity','Galaxy_name','Galaxy_X','Galaxy_Y','Photon_count', 'Direction', 'Distance'])
#dist_list = []
d= 1778.27941*0.000001#1778#0.00018073525133440103# From variable_stars.py (= dist to galaxy with largest x-ray)
max_xRay_photoncount=20089120
uncertaintyList = []
#d=0.374630
#max_xRay_photoncount=331
# For each flash
for i, flash in enumerate(flashdata.values):
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
    #np.sqrt((float(np.power(d,2))*float(max_xRay_photoncount))/float(flash[4]))
  
    dist = np.sqrt((float(np.power(d,2))*float(max_xRay_photoncount))/float(flash[4]))
    flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':dist} 
    
    #uncertaintyList.append(dist/2*np.sqrt(sigma_m*))
    #####flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':np.sqrt((float(np.power(d,2))*float(max_xRay_photoncount))/float(flash[4])), 'Uncertainty':float(0.5)*(np.sqrt(float((max_xRay_photoncount+174)*(d**2))/float(flash[4]-174))-np.sqrt(((float(max_xRay_photoncount-174)*(d**2)))/float(flash[4]+174)))} #np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
    #print(f"len = {len(unc)}")
    #print(f"i = {i}")
    #if i < len(unc):
    ###flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':np.sqrt((float(np.power(d,2))*float(max_xRay_photoncount))/float(flash[4]))} #np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
   
    #print(f"Uncertainty: {np.sqrt(float((max_xRay_photoncount+174)*(d**2))/float(flash[4]-174))-np.sqrt(((float(max_xRay_photoncount-174)*(d**2)))/float(flash[4]+174))}")
    #flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':np.sqrt(2.165*(10**(33))/((4*np.pi*(1.9865*(10**(-16))*flash[4]))/0.051))*3.24077929*10**(-23)} #np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
    #flash_velocity_df.loc[len(flash_velocity_df)] = {'Flash_name':flash[0],'Radial_velocity':closest_galaxy_velocity,'Galaxy_name':closest_galaxy_name,'Galaxy_X':closest_galaxy_x,'Galaxy_Y':closest_galaxy_y,'Photon_count':photon_count, 'Direction':flash[1], 'Distance':np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))} #np.sqrt((float(np.power(d,2))*float(flash[4]))/float(max_xRay_photoncount))
    
    
    #flash_velocity_df.append([flash[0], closest_galaxy_velocity, closest_galaxy_x, closest_galaxy_y, photon_count])
    flash_velocity_df.sort_values(by=['Photon_count'], inplace=True)
    #flash_velocity_df['Uncertainty'] = [5]*len()



# Dropping last n rows using drop
flash_velocity_df.drop(flash_velocity_df.tail(9).index,inplace = True)
print(flash_velocity_df)
velocity = flash_velocity_df['Radial_velocity'].tolist()
distance = flash_velocity_df['Distance'].tolist()

unc = [0.149576812,
0.134774295,
0.115130257,
0.105883771,
0.104982635,
0.10409421,
0.10409421,
0.100662708,
0.099424215,
0.09663113,
0.079976573,
0.079693136,
0.063337525,
0.062198605,
0.055404918,
0.053026381,
0.051488094,
0.05121639,
0.050947066,
0.050022771,
0.048378617,
0.044265968,
0.041933819,
0.040704693,
0.040429558,
0.039272353,
0.036570177,
0.036263724,
0.033553591,
0.030186616,
0.027780794,
0.025721242,
0.025297968,
0.025131899,
0.024764815,
0.022325236,
0.022056798,
0.017639843,
0.012767555,
0.012313989,
0.012152319,
0.011792443,
0.010564529,
0.007875221,
0.007072722,
0.006176421]
unc2 = [
    0.011339389,
0.010217213,
0.008728002,
0.008027028,
0.007958713,
0.007891362,
0.007891362,
0.00763122,
0.00753733,
0.007325587,
0.006063008,
0.006041521,
0.004801605,
0.004715264,
0.004200236,
0.00401992,
0.003903302,
0.003882705,
0.003862287,
0.003792216,
0.003667573,
0.003355794,
0.003178995,
0.003085815,
0.003064957,
0.002977229,
0.002772378,
0.002749146,
0.002543691,
0.002288441,
0.002106057,
0.001949922,
0.001917834,
0.001905244,
0.001877416,
0.001692472,
0.001672122,
0.001337273,
0.000967906,
0.000933521,
0.000921265,
0.000893983,
0.000800895,
0.000597019,
0.000536182,
0.000468233,
]

unc3 = []
###uncertainty = flash_velocity_df['Uncertainty'].tolist() #unc
#print(velocity)
#print(distance)
###myTest(velocity, distance, uncertainty)
plotHubble(velocity, distance, unc3)
#plot_Hubble_old(velocity, distance, uncertainty)
#print(flash_velocity_df)
#print(len(flash_velocity_df))
           
     
    

