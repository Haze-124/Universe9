import numpy as np # for maths 
import matplotlib # for plotting 
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob # this package lets you search for filenames
from pathlib import Path
import random
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
    data = pd.read_csv(ddir+fname) # load in CSV data as a Pandas object
    print(data.keys()) # see what's in it
    time, flux = data.Time, data.NormalisedFlux # just extract the columns as variables
    dt = np.median(np.diff(time))
    print('Nyquist Limit',0.5/dt,'cycles per hour') # can't get frequencies higher than the Nyquist limit
    return time, flux

def plot(time, flux):
    plt.plot(time,flux,'.',markersize=10)
    plt.xlabel('Time (h)')
    plt.ylabel('Normalized Flux')
    plt.savefig('normalizedFlux.png') 
    plt.show()

def get_period_for_MC(time, flux):
    LS = LombScargle(time,flux) # initialize a Lomb-Scargle algorithm from Astropy
    freqs = np.linspace(1/100,0.45,10000) # frequency grid shouldn't go higher than Nyquist limit
    power = LS.power(freqs) # calculate LS power
    return (1/freqs[np.argmax(power)])

def get_period_uncertainty(flux, time, stdTime, stdFlux):
    periodList = []
    for _ in range(10):
        flux2 = []
        time2 = []
        if random.randint(0, 9) <=4:
            flux2 = [f*(1+random.random()*stdFlux) for f in flux]#flux + random.random()*stdFlux
        else:
            flux2 = [f*(1-random.random()*stdFlux) for f in flux]#flux + random.random()*stdFlux
        time2 = [t + random.random()*stdTime for t in time]
        periodList.append(get_period_for_MC(time2.copy(), flux2.copy()))
    return np.std(periodList)

def get_uncertainty_for_luminositySlope(variables, starclass):
    stdPeriod=np.std(variables.Period)
    if starclass == 1:
        m_list = []
        c_list = []
        variables = variables.loc[(variables['Period']!=0) & (variables['Period']>15) & (variables['Period']<25)]
        v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
        abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)
        dp=0.001 # Uncertainty for parallax [arcsec]
        abs_mag_v_uncertainty = np.log10(1/(variables.Parallax-dp))-np.log10(1/(variables.Parallax+dp))#2*np.log10(1./(variables.Parallax*0.001))
        for _ in range(20):
            A = np.vander(variables.Period+random.uniform(-stdPeriod, +stdPeriod),2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
            b, residuals, rank, s = np.linalg.lstsq(A,abs_mag_v+random.uniform(-abs_mag_v_uncertainty, +abs_mag_v_uncertainty))
            m_list.append(b[0])
            c_list.append(b[1])
        return np.std(m_list), np.std(c_list)
    elif starclass == 2:
        m_list = []
        c_list = []
        variables = variables.loc[(variables['Period']!=0) & (variables['Period']>40) & (variables['Period']<50)]
        v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
        abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)
        dp=0.001 # Uncertainty for parallax [arcsec]
        abs_mag_v_uncertainty = np.log10(1/(variables.Parallax-dp))-np.log10(1/(variables.Parallax+dp))#2*np.log10(1./(variables.Parallax*0.001))
        for _ in range(20):
            A = np.vander(variables.Period+random.uniform(-stdPeriod, +stdPeriod),2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
            b, residuals, rank, s = np.linalg.lstsq(A,abs_mag_v+random.uniform(-abs_mag_v_uncertainty, +abs_mag_v_uncertainty))
            m_list.append(b[0])
            c_list.append(b[1])
        return np.std(m_list), np.std(c_list)

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
    freqs = np.linspace(1/100,0.45,10000) # frequency grid shouldn't go higher than Nyquist limit
    
    periods = [] # start an empty list to hold the period 
    names = []
    periodUnc = []
    for fname in tqdm(fnames): # tqdm is a package that gives you a progress bar - neat! 
        data = pd.read_csv(fname) # load in CSV data as a Pandas object
        time, flux = data.Time, data.NormalisedFlux # just extract the columns as variables
        LS = LombScargle(time,flux) # initialize a Lomb-Scargle
        power = LS.power(freqs) # calculate LS power 
        bestfreq = freqs[np.argmax(power)] # which frequency has the highest Lomb-Scargle power?
        periods.append(1/bestfreq) # add each period to the list
        names.append(Path(fname).stem)
        periodUnc.append(get_period_uncertainty(flux, time, 0.3, 0.015))
        
    periodUnc = np.array(periodUnc)
    periods = np.array(periods) # turn it from a list to an array
    #plt.show()
    return names, periods, periodUnc

def PeriodLuminosity(varNames, varPeriods, varPeriodUnc):
    variables = pd.DataFrame({'Name':varNames,
              'Period':varPeriods, 'periodUnc':varPeriodUnc}) # you can turn a dictionary into a dataframe like this
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

    plt.plot(colour,abs_mag,'.C0', label='Benchmark')
    plt.plot(variable_colour,abs_mag_v,'.C2',marker='*',markersize=8, label='Variable Stars')  
    plt.legend()
    plt.ylabel('Log($M_G$)')
    plt.xlabel('Log($M_R$)- Log($M_B$)')
    plt.grid()
    plt.savefig('variables_HR.png')
    
    plt.show()

    plt.plot(variables.Period,abs_mag_v,'.',color='C2')
    plt.xlabel('Variable Period (h)')
    plt.ylabel('Log($M_{G}$)')
    plt.grid()
    plt.savefig('bothClasses.png') 
    plt.show()

def plot_zoom_in(starClass, variables):

    if starClass == 1:
        variables = variables.loc[(variables['Period']!=0) & (variables['Period']>15) & (variables['Period']<25)]
        v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
        abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)
        dp=0.001 # Uncertainty for parallax [arcsec]
        abs_mag_v_uncertainty = np.log10(1/(variables.Parallax-dp))-np.log10(1/(variables.Parallax+dp))#2*np.log10(1./(variables.Parallax*0.001))
        #abs_mag_v_uncertainty = 2*np.log10(1./(variables.Parallax*0.001))
        A = np.vander(variables.Period,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
        b, residuals, rank, s = np.linalg.lstsq(A,abs_mag_v)
        reconstructed = A @ b # @ is shorthand for matrix multiplication in python
        
        fig, ax = plt.subplots()
        plt.xlim([18,24])
        ax.plot(variables.Period,reconstructed,'-r')#,label='LS Estimation')
        ax.errorbar(variables.Period, abs_mag_v,
                    xerr=variables.periodUnc,
                    yerr=abs_mag_v_uncertainty,
                    fmt='.', elinewidth=1) #abs_mag_v_uncertainty,
        ax.set_xlabel('Variable Period (h)')
        ax.set_ylabel('Log($M_G$)')
        plt.gca().legend(('LS Estimation','Class A Data'))
        plt.grid()
        plt.savefig('zoomin1.png') 
        plt.show()


    elif starClass == 2:
        variables = variables.loc[(variables['Period']!=0) & (variables['Period']>40) & (variables['Period']<50)]
        v0, v1, v2 = np.log10(variables['BlueF']), np.log10(variables['GreenF']), np.log10(variables['RedF']) 
        abs_mag_v = v1 + 2*np.log10(1./variables.Parallax)
        dp = 0.001# Uncertainty for parallax [arcsec]
        abs_mag_v_uncertainty = np.log10(1/(variables.Parallax-dp))-np.log10(1/(variables.Parallax+dp))#2*np.log10(1./(variables.Parallax*0.001))
        A = np.vander(variables.Period,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
        b, residuals, rank, s = np.linalg.lstsq(A,abs_mag_v)
        reconstructed = A @ b # @ is shorthand for matrix multiplication in python
        
        fig, ax = plt.subplots()
        plt.xlim([43,48])
        ax.plot(variables.Period,reconstructed,'-r')#,label='LS Estimation')
        ax.errorbar(variables.Period, abs_mag_v,
                    xerr=variables.periodUnc,
                    yerr=abs_mag_v_uncertainty,
                    fmt='.', elinewidth=1) #abs_mag_v_uncertainty,
        ax.set_xlabel('Variable Period (h)')
        ax.set_ylabel('Log($M_G$)')
        plt.gca().legend(('LS Estimation','Class B Data'))
        plt.grid()
        plt.savefig('zoomin2.png') 
    
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

            ## Open corresponding variable-star file
            current_dir = os.getcwd()
            ddir = str(current_dir)+'\\data\\Variable_Star_Data\\'
            flux_data = pd.read_csv(ddir+bytes.decode(var[0], 'utf-8')+'.csv') # load in CSV data as a Pandas object
            var = np.append(var, flux_data['NormalisedFlux'])
            
            if var[9] < 30: # Class A
                m=-0.09357367
                c=-6.10646631
                dist = 10**((10**(m*var[9]+c))/2)
            else: # Class B
                print("class B")
                m=-0.06474778
                c=-1.51789782
                dist=np.sqrt((10**((m*var[9]+c)/2))/(var[10]*4*np.pi))#1./variables['Parallax'][i]
            print(f'DIST = {dist}')

            

    print(f"Closest variable star: {closest_variable_star}")
    print(f"star X: {closest_variable_starX}")
    print(f"star Y: {closest_variable_starY}")

    variables = variables.loc[(variables['Period']!=0) & (variables['Period']>15) & (variables['Period']<25)]
    return dist
        



## RUN FUNCTIONS HERE ##

time, flux = load_and_get_nyquist('BackS023442.csv')
#plot(time, flux)
#get_period_and_freqPlot(time, flux)
names, periods, periodUnc = loop_through_all_stars()
#get_period_uncertainty(time, flux, 0.3, 0.015)
all_stars, variables = PeriodLuminosity(names, periods, periodUnc)
plot_HR(all_stars, variables)
plot_zoom_in(1, variables)
plot_zoom_in(2, variables)

## Fining uncertainty for m and c (y = mx + c) in periodLuminosity graph
m1_std, c1_std = get_uncertainty_for_luminositySlope(variables, 1)
print(f"class A: std(m) = {m1_std}, std(c) = {c1_std}")
m2_std, c2_std = get_uncertainty_for_luminositySlope(variables, 2)
print(f"class B: std(m) = {m2_std}, std(c) = {c2_std}")

#find_dist_to_galaxy(variables, 22.3940, 13.1841, 'Top')
#dist = find_dist_to_galaxy(variables, -4.3630,9.2000, 'Right')
# GALAXY NAME       EQUAT       POLAR           X           Y
# RightDG0412     264.4622      78.4979      -4.3630      9.2000
# TopDG1186       210.4868	    30.0057       22.3940    13.1841
#print(f"Distance to galaxy: {str(dist)} pc (dist units)")

