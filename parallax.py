import os
import pandas as pd
import numpy as np

def load_data(direction):
    current_dir = os.getcwd()
    #stars = pd.read_csv(str(current_dir)+'/data/Converted_Star_Data.csv') 
    stars = pd.read_csv(str(current_dir)+'/data/' + direction + '/Star_Data.csv') 
    print(stars.keys()) # this tells us what column names we have
    return stars

def get_dist_to_stars(stars):
    stars['Distance [pc]'] = np.where(abs(stars['Parallax']) > 0.01, 1/stars['Parallax'], 10000)
    stars.drop(columns=['BlueF', 'GreenF', 'RedF', 'RadialVelocity',"Variable?"], inplace=True)
    stars.drop(stars[stars['Distance [pc]'] == 10000].index, inplace = True)
    return stars

def get_dist_to_galaxies(star_dist_df):
    for star in star_dist_df.values:
        print(star)
    


stars = load_data('Right')
star_dist_df = get_dist_to_stars(stars)

star_dist_df.sort_values(by=['Distance [pc]'], inplace=True)
print(len(star_dist_df))
#print(max(star_dist_df['Distance [pc]']))
get_dist_to_galaxies(star_dist_df)
