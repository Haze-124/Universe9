import os 
import numpy as np
import pandas as pd

def load_data():
    """ Loads data """

def cube_to_equirect(direction, u, v):
    # convert range -45 to 45 to -1 to 1
    uc = u / 45
    vc = v / 45
    if direction == "Front": # POSITIVE X
        x = 1
        y = vc
        z = -uc 
    elif direction == "Back":  # NEGATIVE X
        x = -1
        y = vc
        z = uc
    elif direction == "Top": # POSITIVE Y
        x = uc
        y = 1
        z = -vc
    elif direction == "Bottom": # NEGATIVE Y
        x = uc
        y = -1
        z = vc
    elif direction == "Left": # POSITIVE Z
        x = uc
        y = vc
        z = 1
    else: # direction == "Right": # NEGATIVE Z
        x = -uc
        y = vc
        z = -1 
    # now to convert the XYZ to spherical coordinates
    # this is using the physics convention of spherical coords!
    r = np.sqrt(x**2 + y**2 + z**2)
    azimuth = np.arctan2(z, x)
    theta = np.arccos(y / r)

    theta = theta * 180 / np.pi
    azimuth = (- azimuth + np.pi) * 360 / (2 * np.pi)
    
    return azimuth, theta

def un_cube_plot():
    """ Code taken from the instruction page: https://astrouq.github.io/ladder/tutorials/uncubemapping/
    Shows an image of the un_cubed sides together """
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    
    datapath = os.getcwd() + "/data"

    fig, axes = plt.subplots(3, 4, figsize=(12, 9)) # generate a figure to fit 3 high by 4 wide square images
    fig.subplots_adjust(wspace=0, hspace=0) # we want these squares to be adjacent to each other, with no gap
    # now to iterate over each subplot and remove the axis bars and ticks/labels
    for row in axes:
        for ax in row:
            for side in ['top','right','bottom','left']:
                ax.spines[side].set_visible(False)
            ax.tick_params(axis='both', which='both', labelbottom=False, bottom=False, left=False, labelleft=False)

    # now we load in the images and put them at their correct location
    for i, direct in enumerate(["Back", "Left", "Front", "Right", "Top", "Bottom"]): # one loop for each direction
        img = mpimg.imread(datapath + f'/{direct}/{direct}.png') # this loads in the image from the corresponding folder
        img_cropped = img[:-700, 900:] # crop the image to remove the axis labels/ticks
        if i == 4: # if the Top image
            imgplot = axes[0][1].imshow(img_cropped) # image needs to go at the top
        elif i == 5: # if the Bottom image
            imgplot = axes[2][1].imshow(img_cropped) # image needs to go at the bottom
        else:
            imgplot = axes[1][i].imshow(img_cropped) # put the image in the middle row at the correct column
    plt.show()
    return plt

def un_cube_stars():
    plt = un_cube_plot()
    datapath = os.getcwd() + "\\data"
    for i, direct in enumerate(["Front", "Back", "Left", "Right", "Top", "Bottom"]):
        # read the data from the .txt file into a dataframe
        stardata = pd.read_csv(datapath + f'\\{direct}\\Star_Data.csv', delimiter=',') 
        u = stardata['X'].to_numpy()
        v = stardata['Y'].to_numpy() # convert X and Y data to "U" and "V" data
        azimuth, theta = cube_to_equirect(direct, u, v) # perform the coordinate transform
        azimuth = np.around(azimuth, decimals=4); theta = np.around(theta, decimals=4) # round to appropriate decimals
        
        df = pd.DataFrame({"Equat": azimuth, "Polar": theta}) # make a temporary DataFrame object with new coordinates
        # now overwrite the old coordinates with the new ones
        stardata['X'] = df['Equat']
        stardata['Y'] = df['Polar']
        stardata = stardata.rename(columns={"X": "Equat", "Y": "Polar"}) # and finally change the name of the columns 
        if i == 0:
            # if this is the first iteration, write to a new DataFrame that will store all of the star data
            all_stardata = stardata
        else:
            all_stardata = pd.concat([all_stardata, stardata]) # add this face stardata to the rest of the data

    # now let's plot the data to see if it's worked!
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(all_stardata["Equat"].to_numpy(), all_stardata["Polar"].to_numpy(), s=0.1, c='k', lw=0)
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.invert_yaxis()
    #plt.show()
    all_stardata.to_csv(datapath + "/Converted_Star_Data.csv", index=False, sep=',') # Save data to csv file
    # dont want to save the 'indices' of the data, and I want a space character to separate the data

def un_cube_galaxies():
    plt = un_cube_plot()
    datapath = os.getcwd() + "\\data"
    for i, direct in enumerate(["Front", "Back", "Left", "Right", "Top", "Bottom"]):
        # read the data from the .txt file into a dataframe
        galaxdata = pd.read_csv(datapath + f'\\{direct}\\Distant_Galaxy_Data.csv', delimiter=',')  
        u = galaxdata["X"].to_numpy(); v = galaxdata["Y"].to_numpy() # convert X and Y data to "U" and "V" data
        azimuth, theta = cube_to_equirect(direct, u, v) # perform the coordinate transform
        azimuth = np.around(azimuth, decimals=4); theta = np.around(theta, decimals=4) # round to appropriate decimals
        
        df = pd.DataFrame({"Equat": azimuth, "Polar": theta}) # make a temporary DataFrame object with new coordinates
        # now overwrite the old coordinates with the new ones
        galaxdata['X'] = df['Equat']
        galaxdata["Y"] = df["Polar"]
        galaxdata = galaxdata.rename(columns={"X": "Equat", "Y": "Polar"}) # and finally change the name of the columns 
        if i == 0:
            # if this is the first iteration, write to a new DataFrame that will store all of the star data
            all_galaxdata = galaxdata
        else:
            all_galaxdata = pd.concat([all_galaxdata, galaxdata]) # add this face stardata to the rest of the data

    # now let's plot the data to see if it's worked!
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(all_galaxdata["Equat"].to_numpy(), all_galaxdata["Polar"].to_numpy(), s=0.5, c='k', lw=0);
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.invert_yaxis()
    plt.show()
    all_galaxdata.to_csv(datapath + "/Converted_Galaxy_Data.csv", index=False, sep=',') # Save data to csv file
    # dont want to save the 'indices' of the data, and I want a space character to separate the data

def identify_galaies(data):
    """ Returns all the galaxies in the input data """

def find_distances(galaxies):
    """ Returns distances to all input galaxies """

#un_cube_stars()
un_cube_galaxies()

