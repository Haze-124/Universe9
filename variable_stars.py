import numpy as np # for maths 
import matplotlib # for plotting 
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob # this package lets you search for filenames
from pathlib import Path

from tqdm import tqdm # tqdm is a package that lets you make progress bars to see how a loop is going

import os 
from sklearn.metrics import pairwise_distances

import pandas as pd # pandas is a popular library in industry for manipulating large data tables
from astropy.timeseries import LombScargle

# configure notebook for plotting
#matplotlib inline

mpl.style.use('seaborn-colorblind') # colourblind-friendly colour scheme

# subsequent lines default plot settings
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams['figure.figsize']=(8.0,6.0)   
matplotlib.rcParams['font.size']=16              
matplotlib.rcParams['savefig.dpi']= 300             

import warnings
warnings.filterwarnings('ignore')

def load_and_get_nyquist(fname):
    current_dir = os.getcwd()
    ddir = str(current_dir)+'\\data\\Variable_Star_Data\\'
    #fname = 'BackS023442.csv' # put your filename here
    data = pd.read_csv(ddir+fname) # load in CSV data as a Pandas object
    print(data.keys()) # see what's in it
    time, flux = data.Time, data.NormalisedFlux # just extract the columns as variables
    dt = np.median(np.diff(time))
    print('Nyquist Limit',0.5/dt,'cycles per hour') # can't get frequencies higher than the Nyquist limit
    return time, flux

def plot(time, flux):
    plt.plot(time,flux,'.',markersize=16)
    plt.xlabel('Time (h)')
    plt.ylabel('Relative Flux')
    plt.show()

def get_period_and_freqPlot(time, flux):
    LS = LombScargle(time,flux) # initialize a Lomb-Scargle algorithm from Astropy
    freqs = np.linspace(1/100,0.45,10000) # frequency grid shouldn't go higher than Nyquist limit
    power = LS.power(freqs) # calculate LS power
    print('Best period: %.2f h' % (1/freqs[np.argmax(power)]))
    plt.plot(freqs,power)
    plt.xlabel('Frequency (c/h)')
    plt.ylabel('LS Power') # LS stands for Lomb-Scargle
    plt.show()

def loop_through_all_stars():
    current_dir = os.getcwd()
    ddir = str(current_dir)+'\\data\\Variable_Star_Data\\'
    fnames = glob.glob(ddir+'*.csv')
    #print(fnames[:10])
    freqs = np.linspace(1/100,0.45,10000) # frequency grid shouldn't go higher than Nyquist limit
    
    periods = [] # start an empty list to hold the period 
    names = []
    #fig, axes = plt.subplots(3,4,figsize=(18,12))
    #for fname, ax in zip(fnames[:12], axes.ravel()): # you can loop over two things
    for fname in tqdm(fnames): # tqdm is a package that gives you a progress bar - neat! 
        data = pd.read_csv(fname) # load in CSV data as a Pandas object
        time, flux = data.Time, data.NormalisedFlux # just extract the columns as variables
        LS = LombScargle(time,flux) # initialize a Lomb-Scargle
        power = LS.power(freqs) # calculate LS power 
        bestfreq = freqs[np.argmax(power)] # which frequency has the highest Lomb-Scargle power?
        #pred = LS.model(time,bestfreq) # make a sine wave prediction at the best frequency
        #ax.plot(time,flux,'.')
        #ax.plot(time,pred) # plot the model over the data
        periods.append(1/bestfreq) # add each period to the list
        names.append(Path(fname).stem)
        
    periods = np.array(periods) # turn it from a list to an array
    #plt.show()
    return names, periods

def PeriodLuminosity(varNames, varPeriods):
    variables = pd.DataFrame({'Name':varNames,
              'Period':varPeriods}) # you can turn a dictionary into a dataframe like this
    variables.Name = variables.Name.astype('|S') # have to do this so that it knows the names are strings

    current_dir = os.getcwd()
    #all_star_files = glob.glob(current_dir+'/data/Converted_Star_Data.csv')
    all_star_files = glob.glob(current_dir+'/data/*/Star_Data.csv')

    all_stars = pd.concat([pd.read_csv(table) for table in all_star_files]) # we are concatenating a list of dataframes; 
    #we generate this list with a "list comprehension", a loop you write inside a list bracket 

    all_stars.Name = all_stars.Name.astype('|S') # have to do this so that it knows the names are strings
    all_stars = all_stars[all_stars.Parallax > 0.01] # 10 mas parallax cut
    print(len(all_stars),'stars above 10 mas parallax') # check how many stars there are total with good parallax

    variables = pd.merge(all_stars,variables,on='Name') # merge these two arrays according to the keyword 'name'
    print('Of which',len(variables),'variables') # cut down to a small list
    print(variables.head())
    return all_stars, variables

def plot_HR(all_stars, variables):
    m0, m1, m2 = np.log10(all_stars['BlueF']), np.log10(all_stars['GreenF']), np.log10(all_stars['RedF']) 
    colour = m2-m0
    abs_mag = m1 + 2*np.log10(1./all_stars.Parallax) 

    v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
    variable_colour = v2-v0
    abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)

    s = plt.plot(colour,abs_mag,'.C0')
    h = plt.plot(variable_colour,abs_mag_v,'.C2',marker='*',markersize=10)  
    plt.legend([s, h],['Steady','Variable'])
    plt.ylabel('Log Flux 1')
    plt.xlabel('Log Flux 2 - Log Flux 0')
    plt.grid()
    plt.show()

    plt.plot(variables.Period,abs_mag_v,'.',color='C2')
    plt.xlabel('Period (h)')
    plt.ylabel('Log Flux')
    plt.grid()
    plt.show()

def plot_zoom_in(starClass, variables):
    #print(variables.head())

    if starClass == 1:
        variables = variables.loc[(variables['Period']!=0) & (variables['Period']>15) & (variables['Period']<25)]
        v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
        abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)
        plt.xlim([18,24])
        plt.plot(variables.Period,abs_mag_v,'.',color='C2',label='Data')
        A = np.vander(variables.Period,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
        b, residuals, rank, s = np.linalg.lstsq(A,abs_mag_v)
        reconstructed = A @ b # @ is shorthand for matrix multiplication in python
        plt.plot(variables.Period,reconstructed,'-r',label='LS Estimation')
        plt.legend()
        plt.xlabel('Period (h)')
        plt.ylabel('Log Flux')
        plt.grid()
        plt.show()

    elif starClass == 2:
        variables = variables.loc[(variables['Period']!=0) & (variables['Period']>40) & (variables['Period']<50)]
        v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
        abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)
        plt.xlim([43,48])
        plt.plot(variables.Period,abs_mag_v,'.',color='C2', label='Data')
        A = np.vander(variables.Period,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
        b, residuals, rank, s = np.linalg.lstsq(A,abs_mag_v)
        reconstructed = A @ b # @ is shorthand for matrix multiplication in python
        plt.plot(variables.Period,reconstructed,'-r',label='LS estimation')
        plt.legend()
        plt.xlabel('Period (h)')
        plt.ylabel('Log Flux')
        plt.grid()
        plt.show()
    else:
        print("You didn't specify a star class. Input 1 or 2.")

def find_dist_to_galaxy(variables, galaxyX, galaxyY, direction):
    ## Set default values
    min_dist = np.sqrt(np.power((variables['X'][0]-galaxyX), 2)+np.power((variables['Y'][0]-galaxyY), 2))
    closest_variable_star = variables['Name'][0]
    closest_variable_starX = variables['X'][0]
    closest_variable_starY = variables['Y'][0]
    dist=0
    ## Go through all variable stars
    for i, var in enumerate(variables.values):
        ## If the distance between the variable star and the galaxy (and looking in correct direction) is the smallest one yet
        if (np.sqrt(np.power((var[1]-galaxyX), 2)+np.power((var[2]-galaxyY), 2)) < min_dist) & (direction in bytes.decode(var[0], 'utf-8')):
            ## Set this as new minimum distance
            min_dist = np.sqrt(np.power((var[1]-galaxyX), 2)+np.power((var[2]-galaxyY), 2))
            closest_variable_star = var[0]
            closest_variable_starX = var[1]
            closest_variable_starY = var[2]

            v0, v1, v2 = np.log10(variables['BlueF'][i]), np.log10(variables['GreenF'][i]), np.log10(variables['RedF'][i]) 
            luminosity = v1 + 2*np.log10(1./variables['Parallax'][i])
            brightness = 20089120 # photon count/flux
            print(f"L = {luminosity}")
            print(f"b = {brightness}")
            dist = np.sqrt(abs(luminosity)/(4*np.pi*brightness)) # [distance units]

    print(f"Closest variable star: {closest_variable_star}")
    print(f"star X: {closest_variable_starX}")
    print(f"star Y: {closest_variable_starY}")

    variables = variables.loc[(variables['Period']!=0) & (variables['Period']>15) & (variables['Period']<25)]
    print('DISTANCE:')
    print(dist)
    return dist
        




#time, flux = load_and_get_nyquist('BackS023442.csv')
#plot(time, flux)
#get_period_and_freqPlot(time, flux)
names, periods = loop_through_all_stars()
#print(f"NAMES = {names}")
#print(f"PERIODS = {periods}")
all_stars, variables = PeriodLuminosity(names, periods)
#plot_HR(all_stars, variables)
#plot_zoom_in(2, variables)
#find_dist_to_galaxy(variables, 22.3940, 13.1841, 'Top')
dist = find_dist_to_galaxy(variables, -4.3630,9.2000, 'Right')
# GALAXY NAME       EQUAT       POLAR           X           Y
# RightDG0412     264.4622      78.4979      -4.3630      9.2000
# TopDG1186       210.4868	    30.0057       22.3940    13.1841
print(f"Distance to galaxy: {str(dist)} DU (dist units)")

