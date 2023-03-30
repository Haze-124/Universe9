import numpy as np # for maths 
import matplotlib # for plotting 
import matplotlib as mpl
import matplotlib.pyplot as plt

import os 
import glob # this package lets you search for filenames

import pandas as pd # pandas is a popular library in industry for manipulating large data tables

# configure notebook for plotting
#matplotlib inline 
mpl.style.use('seaborn-colorblind') # colourblind-friendly colour scheme

# define default plot settings
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams['figure.figsize']=(8.0,6.0)    #(6.0,4.0)
matplotlib.rcParams['font.size']=16              #10 
matplotlib.rcParams['savefig.dpi']= 300             #72 

import warnings
warnings.filterwarnings('ignore')

def load_data(direction):
    current_dir = os.getcwd()
    #stars = pd.read_csv(str(current_dir)+'/data/Right/Star_Data.csv') 
    stars = pd.read_csv(str(current_dir)+'/data/' + direction + '/Star_Data.csv') 
    print(stars.keys()) # this tells us what column names we have
    return stars

def plot_all_star_positions(stars):
    plt.scatter(stars.X,stars.Y)
    plt.xlabel('x (pix)')
    plt.ylabel('y (pix)')
    plt.show()

def plot_zoom_in(centre, max_dist, stars):
    d = np.sqrt((stars.X-centre[0])** 2 + (stars.Y - centre[1])**2)
    galaxy = stars[d<max_dist] # filter to only close ones
    plt.scatter(galaxy.X,galaxy.Y,c=galaxy.RadialVelocity,cmap=mpl.cm.seismic) # let's overplot the radial velocities
    plt.colorbar()
    plt.scatter(*centre,color='C2',marker='X') # * expands the elements of a list 
    plt.xlabel('x (pix)')
    plt.ylabel('y (pix)')
    plt.show()
    return galaxy

def plot_HR_diagram(galaxy):
    m0, m1, m2 = (np.log10(galaxy['BlueF']), 
                np.log10(galaxy['GreenF']), 
                np.log10(galaxy['RedF'])) 
    colour = m2-m0
    s = plt.scatter(colour,m1)
    plt.ylabel('Log Flux 1')
    plt.xlabel('Log Flux 2 - Log Flux 0')
    plt.show()
    print('Parallaxes: mean %.3f, sd %.3f' % (np.mean(galaxy['Parallax']),np.std(galaxy['Parallax'])))
    return max(m1)

def plot_HR_diagram_for_nearby_stars():
    current_dir = os.getcwd()
    all_stars = glob.glob(current_dir+'/data/*/Star_Data.csv')

    fig, ax1 = plt.subplots(1,1)
    for j, catalog in enumerate(all_stars):
        try:
            this = pd.read_csv(catalog)
            
            thispar = this.Parallax
            thism0, thism1, thism2 = (np.log10(this.BlueF), 
                                    np.log10(this.GreenF), 
                                    np.log10(this.RedF))
            thiscolour = thism2-thism0
            dist = 1/thispar
            abs_mag = thism1 + 2*np.log10(dist) 
            mm = thispar>0.010 # only pick the ones with good signal-to-noise - 10 mas is ok 
            
            ax1.scatter(thiscolour[mm],abs_mag[mm],color='C1')
        except:
            pass

    #plt.ylabel('Log Flux 1')
    #plt.xlabel('Log Flux 2 - Log Flux 0')
    plt.ylabel('Log($M_1$)')
    plt.xlabel('Log($M_2$) - Log($M_0$)')
    plt.grid()
    plt.show()
    return all_stars

def plot_Benchmark_and_Cluster(max_cluster):
    current_dir = os.getcwd()
    all_stars = glob.glob(current_dir+'/data/*/Star_Data.csv')
    fig = plt.figure()
    max_benchmark = -2000
    for j, catalog in enumerate(all_stars):
        try:
            this = pd.read_csv(catalog)
            
            thispar = this.Parallax
            thism0, thism1, thism2 = (np.log10(this.BlueF), 
                                    np.log10(this.GreenF), 
                                    np.log10(this.RedF))
            thiscolour = thism2-thism0
            dist = 1/thispar
            abs_mag = thism1 + 2*np.log10(dist) 
            mm = thispar>0.010 # only pick the ones with good signal-to-noise - 10 mas is ok 
            h = plt.scatter(thiscolour[mm],abs_mag[mm], 5,color='C1')
            if max(abs_mag[mm]) > max_benchmark:
                max_benchmark = max(abs_mag[mm])
        except:
            pass

    m0, m1, m2 = (np.log10(galaxy['BlueF']), 
                np.log10(galaxy['GreenF']), 
                np.log10(galaxy['RedF'])) 
    colour = m2-m0
    max_cluster = max(m1)
    print(f"bench max = {max_benchmark}")
    print(f"cluster max = {max_cluster}")
    dist=max_cluster-max_benchmark
    print(f"dist = {dist}")
    #s = plt.scatter(colour,m1+np.log10(dist),color='C0')
    #s = plt.scatter(colour,m1-dist,color='C0')
    s = plt.scatter(colour,m1+6.5, 5,color='C0')
    plt.ylabel('Log($M_1$)')
    plt.xlabel('Log($M_2$) - Log($M_0$)')
    plt.legend([h,s],['Benchmark','RightDG0412']) # 'Cluster'
    plt.grid()
    plt.savefig('benchmark.png') 
    plt.show()

    return 10**abs(dist/2)


# Identify each galaxy (clusters of stars) and return its center point
# Then for each galaxy  
    # Zoom into it
    # Find its maximum value in log Flux 1 scale

#plt.close('all')
stars = load_data('Right')
#plot_all_star_positions(stars)
#centre = ( 0.6622, 0.2589) # Choose point to zoom in around
centre = ( -4.3630, 9.2000) # Choose point to zoom in around
zoom_in_distance=10
galaxy = plot_zoom_in(centre, zoom_in_distance, stars)
max_cluster = plot_HR_diagram(galaxy)
#plot_HR_diagram_for_nearby_stars() # includes all stars in the galaxy
dist_to_galaxy = plot_Benchmark_and_Cluster(max_cluster)
#print(f"Distance to galaxy: {(int) (dist_to_galaxy)} pc")






