# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 09:15:39 2023

@author: lbr20
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# read CSV file into a pandas DataFrame
df = pd.read_csv('lbr0024_A549-EGFRmEos_NoVirus_005.csv', skiprows = 3)

# define the frame ranges for plotting
frame_ranges = [(0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 2500), (2500, 3000), (3000, 3500), (3500, 4000), (4000, 4500), (4500, 5000), 
                (5000, 5500), (5500, 6000), (6000, 6500), (6500, 7000), (7000, 7500), (7500, 8000), (8000, 8500), (8500, 9000), (9000, 9500), (9500, 9999)]

# set the style to use white circles on black background
plt.style.use('dark_background')

# loop through the frame ranges and plot the corresponding data
for frame_range in frame_ranges:
    # filter the DataFrame to only include rows within the current frame range
    frame_min, frame_max = frame_range
    df_filtered = df[(df['frame_ix'] >= frame_min) & (df['frame_ix'] < frame_max)]
    
    # create a 256x256 pixel figure and axis object
    fig, ax = plt.subplots(figsize=(256/100, 256/100))
    x=np.linspace(0,256)
    y=np.linspace(0,256)
    plt.plot(x, y, linewidth=0)
    
    # create scatter plot using the 'x' and 'y' columns from the filtered DataFrame
    ax.scatter(df_filtered['x'], df_filtered['y'], marker='o', s=0.1, facecolor='white', edgecolor='white', linewidths=0.1)

    # set title and axis labels for the current plot
    ax.set_title(f"Frame range: {frame_min} - {frame_max}")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    # save the scatter plot as a PNG file with a name based on the current frame range
    filename = f"scatterplot_{frame_min}_{frame_max}.png"
    fig.savefig(filename, dpi=1200)
    
    # clear the plot for the next iteration
    plt.clf()
