# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 13:45:13 2023

@author: lukas Broich
"""


# %%

## Import all necessary packages
## For compatibility with Python 2 and 3
from __future__ import division, unicode_literals, print_function  

## Prevent some error messsages
import warnings
warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", RuntimeWarning)
warnings.simplefilter("ignore", FutureWarning)

## Import packages
import numpy as np
import pandas as pd
#from pandas import DataFrame, Series
import pims
import trackpy as tp
import os
#import scipy.misc
from scipy import ndimage
#import functools
#import operator

## Import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

## Optionally, tweak matplot styles.
plt.style.use('default')
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

# %%

## Add necessary data
folderpath = r'C:/Users/lbr20/Desktop/WD/lbr0026/013/' ## Parentfolder of your data (Where Image and localizationfile are stored)
imagename = r'lbr0026_A549_EGFR-mEos_PR8_Sialidase_Virus_013.nd2 - C=0_2.tif' ## Virus Image
locfilename = r'lbr0026_A549_EGFR-mEos_PR8_Sialidase_imageStack_013.csv' ## Decode localizations

videofps = 31.5 ## Frames per second of your rawdata
plot_single = False ## Do you want to have plots for individual Tracks?
conf_tresh_D = 0.01 ## Set diffusivity threshold to consider confined particles (µm²/s)
conf_tresh_a = 0.6 ## Set power-law threshold to consider confined particles


## Calculate first and third tau for Diffusion coefficient calculation
tau1 = 1/videofps
tau3 = (1/videofps)*3

## Modify image to fit to decode localizations
frame = pims.open(folderpath + imagename) ## Load in virus image
frame = ndimage.rotate(frame[0], 90) ## Rotate virus image (Needs to be done because of decode)
frame = np.flipud(frame) ## Mirror virus image (Needs to be done because of decode)

##Load in decode localizations and make them trackpy suitable
data_raw = pd.read_csv(folderpath + locfilename, skiprows = 3) ## Load in Localizations
data_1 = data_raw.rename(columns={'frame_ix': 'frame'}) ## Rename 'frame' column so Trackpy recognizes it
data_2 = data_1[['x','y','phot','frame','prob','bg','phot_sig','x_sig','y_sig']] ## Only keep the columns with relevant information to reduce dataframe

## Create an outputfolder to save the data
outfolder_path = folderpath + 'Analysis/' ## Add '/' to the end!
outfolder_path = outfolder_path.replace("\\", "/")
if not os.path.exists(outfolder_path):
    os.makedirs(outfolder_path)
# %%
    
## Filter the Decode localizations
## Cut out first 10% of frames
framecutoff = int(round((max(data_2['frame']) * 0.1)-1, 0)) ## Change 0.1 to desired % of frames removed from the beginning
data_filt_frame = data_2[data_2['frame'] >= framecutoff] ## This is done to get rid of eventual bleedthrough effects

## Filter the dimmest and brightest 2.5% of photons
photdnlim = np.percentile(data_2['phot'], 2.5) ## Change 2.5 to desired lower cut off
photuplim = np.percentile(data_2['phot'], 97.5) ## Change 97.5 to desired upper cut off
data_filt_phot = data_filt_frame[(data_filt_frame['phot'] >= photdnlim) & (data_filt_frame['phot'] <= photuplim)]

## Filter the 2.5 most inaccurate localizations made by decode
xsiguplim = np.percentile(data_2['x_sig'], 97.5)
ysiguplim = np.percentile(data_2['y_sig'], 97.5)
data_filt_sig = data_filt_phot[(data_filt_phot['x_sig'] <= xsiguplim) & (data_filt_phot['y_sig'] <= ysiguplim)]

## Save the filtered localizations
data_filt = data_filt_sig ## Rename the data frame
filt_loc_name = outfolder_path + locfilename ## Adds directory and base name of localizations
filt_loc_name = filt_loc_name.replace('.csv', '_filtered_localizations.csv') ## Adds 'filtered_localizations.csv' to the basefilename
data_filt.to_csv(filt_loc_name, index=False) ## Saves


#%%

# filter dataframe
#data_filt = data_filt.drop(data_filt[((data_filt['x'] >= 174) & (data_filt['x'] <= 181)) & ((data_filt['y'] >= 18) & (data_filt['y'] <= 26))].index)
#data_filt.to_csv(filt_loc_name, index=False) ## Saves

# %%

## Link localizations between frames to get trajectories
all_traj = tp.link(data_filt, search_range = 0.8, memory= 0, adaptive_stop = 0.1, adaptive_step = 0.9)

## Save the trajectories in a .csv
all_traj_name = outfolder_path + locfilename ## Adds directory and base name of localizations
all_traj_name = all_traj_name.replace('.csv', '_all_tracks.csv') ## Adds 'all_tracks.csv' to the basefilename
all_traj.to_csv(all_traj_name, index=False) ## Saves


#%%

## Filter Trajectories
non_stub_traj = tp.filter_stubs(all_traj, 10) ## Filter tracks short than 0.3 secs (10 frames)
filt_traj = non_stub_traj.groupby('particle').filter(lambda x : len(x)<=133) ## Filter tracks longer than 4secs (133 frames)
work_traj = filt_traj

## Print the before and after values
print('Trajectories before filtering:', all_traj['particle'].nunique()) ## Prints before
print('Trajectories after filtering:', filt_traj['particle'].nunique()) ## Prints after

## Save the filtered trajectories
filt_traj_name = outfolder_path + locfilename ## Adds directory and base name of localizations
filt_traj_name = filt_traj_name.replace('.csv', '_filtered_tracks.csv') ## Adds 'filtered_tracks.csv' to the basefilename
filt_traj.to_csv(filt_traj_name, index=False) ## Saves

#%%

##############################
## Optionally drift correct ##
##############################

drift = tp.compute_drift(work_traj)
drift_corr_traj = tp.subtract_drift(work_traj.copy(), drift)
work_traj = drift_corr_traj

## Save the drift corrected trajectories
drift_traj_name = outfolder_path + locfilename ## Adds directory and base name of localizations
drift_traj_name = drift_traj_name.replace('.csv', '_drift corrected tracks.csv') ## Adds 'filtered_tracks.csv' to the basefilename
drift_raw_name = drift_traj_name.replace('_drift corrected tracks.csv', '_drift values.csv')
#drift_corr_traj.to_csv(drift_traj_name, index=False) ## Saves
drift.to_csv(drift_raw_name)

## Plot the drift
fig, ax = plt.subplots()
ax.plot(drift[drift.columns[0]], label='y')
ax.plot(drift[drift.columns[1]], label='x')
ax.set(ylabel='Drift [px]',
       xlabel='Frame')
plt.legend();

##Save the plot
plt.savefig(drift_traj_name.replace('_drift corrected tracks.csv', '_drift.png'))

#%%

## Calculate individual MSDs
iMSD = tp.imsd(work_traj, mpp = 0.160, fps = videofps)

## For some particles the drift correction adds an additional value at the end or beginning for unknown reasons. 
## This is to filter them out
for col in iMSD.columns:
    # check if first MSD value is greater than 0.1
    if iMSD.iloc[0][col] > 0.1:
        del iMSD[col] ## removes particle from iMSD
        work_traj = work_traj[work_traj.particle != col] ## removes particle from trajectories

iMSD_name = outfolder_path + locfilename ## Adds directory and base name of localizations
iMSD_name = iMSD_name.replace('.csv', '_iMSDs of all particles.csv') ## Adds 'filtered_tracks.csv' to the basefilename
iMSD.to_csv(iMSD_name)
work_traj.to_csv(drift_traj_name, index=False) ## Saves

## Plot individual MSDs
fig, ax = plt.subplots()
ax.plot(iMSD.index, iMSD, 'k-', alpha = 0.1)  # black lines, semitransparent
ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
       xlabel='lag time $t$ [s]')
ax.set_xscale('log')
ax.set_yscale('log')

## Save the Plot
plt.savefig(iMSD_name.replace('_iMSDs of all particles.csv', '_IMSDs plot.png'))

#%%

## Calculate all diffusion and power-law exponents
work_traj['Init D'] = 0
work_traj['Confined D'] = 0
count_all = 0
count_conf = 0

work_traj['Power-law'] = 0
work_traj['Confined a'] = 0
count_all_a = 0
count_conf_a = 0
## Calculate all

for i in iMSD:
    iCalc = iMSD.loc[tau1:tau3,i] ## Calculate diffusion coefficient from the first 3 taus
    iFIT = tp.utils.fit_powerlaw(iCalc, plot = False)
    iDiffu = iFIT.iloc[0]['A'] / 4
    work_traj.loc[work_traj['particle'] == i,'Init D'] = iDiffu
    count_all = count_all + 1
    
    if iDiffu <= conf_tresh_D:
        work_traj.loc[work_traj['particle'] == i,'Confined D'] = 1
        count_conf = count_conf + 1
        
    ### Power-law calculation
    length = round(iMSD[i].notna().sum() / 2) ## Bestimmt wie viele Datenpunkte es für Partikel i gibt. Ersetzt man die 1 durch 2 wird später nur aus der ersten Hälfte der Daten alpha berechnet.
    tauend = (1/videofps*length) ## Setzt den Endwert für die lagtime beim berechnen.
    iCalc = iMSD.loc[tau1:tauend,i] ## Macht dataframe aus allen vorher festgelegten Datenpunkten
    iFIT = tp.utils.fit_powerlaw(iCalc, plot = False)
    iAlpha = iFIT.iloc[0]['n']
    work_traj.loc[work_traj['particle'] == i,'Power-law'] = iAlpha
    count_all_a = count_all_a + 1
    
    if iAlpha <= conf_tresh_a:
        work_traj.loc[work_traj['particle'] == i,'Confined a'] = 1
        count_conf_a = count_conf_a + 1    

## Save file
traj_and_diffu_name = outfolder_path + locfilename ## Adds directory and base name of localizations
traj_and_diffu_name = traj_and_diffu_name.replace('.csv', '_Tracks with diffusion and power-law.csv') ## Adds 'filtered_tracks.csv' to the basefilename
work_traj.to_csv(traj_and_diffu_name, index=False) ## Saves

## Calculate percent of confined particles
p_d = round((count_conf/count_all)*100 , 2)
p_d_str = str(p_d)
t_d_str = str(conf_tresh_D)

out_str_d = str(p_d_str + '% of your particles are considered confined.\nTreshold for confinement was: ' + t_d_str +' µm²/s')
conf_perc_name = outfolder_path + locfilename ## Adds directory and base name of localizations
conf_perc_name = conf_perc_name.replace('.csv', '_percent of confined particles by diffusion.txt')
text_file = open(conf_perc_name, "w")
text_file.write(out_str_d)
text_file.close()

## Calculate percent of confined particles
p_a = round((count_conf_a/count_all_a)*100 , 2)
p_a_str = str(p_a)
t_a_str = str(conf_tresh_a)

out_str_a = str(p_a_str + '% of your particles are considered confined.\nTreshold for confinement was: ' + t_a_str)
conf_perc_name = outfolder_path + locfilename ## Adds directory and base name of localizations
conf_perc_name = conf_perc_name.replace('.csv', '_percent of confined particles by power-law.txt')
text_file = open(conf_perc_name, "w")
text_file.write(out_str_a)
text_file.close()
print(out_str_d)
print(out_str_a)

# %%

group_confined = []
group_free = []

for col in iMSD.columns:
    confined = work_traj.loc[work_traj['particle'] == col, 'Confined D'].iloc[0]
    if confined == 1:
        group_confined.append(col)
    else:
        group_free.append(col)

## Create new dataframes for each group
iMSD_group_confined = iMSD[group_confined]
iMSD_group_free = iMSD[group_free]

## Save names

outname = outfolder_path + locfilename ## Adds directory and base name of localizations
outname_conf = outname.replace('.csv', '_iMSDs of confined particles.csv')
outname_free = outname.replace('.csv', '_iMSDs of free particles.csv')

iMSD_group_confined.to_csv(outname_conf)
iMSD_group_free.to_csv(outname_free)

## Print the results
print(iMSD_group_confined)
print(iMSD_group_free)
#%%

## Plot Trajectories OVER WHOLE IMAGE

######################################
## If you dont have a Virus picture ##
######################################
#x=np.linspace(0,256)
#y=np.linspace(0,256)
#plt.plot(x, y, linewidth=0)
#plt.axis('square')
#tp.plot_traj(work_traj);

##################################
## If you  have a Virus picture ##
##################################
tp.plot_traj(work_traj, superimpose = frame); ## With drift correction

###################################
## Always execute to empty cache ##
###################################
overall_alpha = []
overall_diffu = []

squaresize = 5 ## Side length of the square you want to analyze (in px)

#%%
virusx = 150
virusy = 97

diffudata = []
alphadata = []
group_confined_spot = pd.DataFrame()
group_free_spot = pd.DataFrame()

## Create Outputfolder
p1 = 'x' + str(virusx) ## Converts X-coordinate to string
p2 = 'y' + str(virusy) ## Converts Y-coordinate to string
outfolder_path_calc = outfolder_path + '/Squares/' + p1 + p2 ## Creates the folder name
if not os.path.exists(outfolder_path_calc):
    os.makedirs(outfolder_path_calc) ## Creates the actual folder

## Set limits    
xlimup = virusx + squaresize/2 ## Upper X
xlimdn = virusx - squaresize/2 ## Lower X
ylimup = virusy + squaresize/2 ## Upper Y
ylimdn = virusy - squaresize/2 ## Lower Y

## Make new dataframe containing only trajectories in the square
point_traj = work_traj[(work_traj['x'] > xlimdn) & (work_traj['x'] < xlimup) & (work_traj['y'] > ylimdn) & (work_traj['y'] < ylimup)] ## Look which particles are in the square
particle_unique = point_traj['particle'].unique() ## Get every particles ID once
particle_list = particle_unique.tolist() ## Put particles ID in a list
traj_calc = work_traj[work_traj['particle'].isin(particle_list)] ## Get ALL values for the particles in the list from the original data

## Calculate initial diffusivities
for i in (point_traj['particle'].unique()): ## Go through all particles in frame
    iCalc = iMSD.loc[tau1:tau3,i] ## Diffu aus den ersten 3 Spots
    iFIT = tp.utils.fit_powerlaw(iCalc, plot = False)
    iDiffu = iFIT.iloc[0]['A'] / 4
    diffudata.append(iDiffu)
    overall_diffu.append(iDiffu)
    print('The initial diffusivity of particle number',i,'is: ',iDiffu,'µm^2/s.')
 ###
    if iDiffu < conf_tresh_D:
        # move column to new dataframe
        group_confined_spot[i] = iMSD[i]
    else:
        group_free_spot[i] = iMSD[i]
##
    ## Save iMSDs

    group_confined_spot.to_csv(outfolder_path_calc + '/' + p1 + p2 + '_confined iMSDs on spot.csv', sep = ',')
    group_free_spot.to_csv(outfolder_path_calc + '/' + p1 + p2 + '_Free iMSDs on spot.csv', sep = ',')
    
    ## Save the diffusivites
    iDiffudata = pd.DataFrame(diffudata, columns = ['Initial diffusion coefficient']) ## Save Diffusioncoefficients
    iDiffudata.to_csv(outfolder_path_calc + '/' + p1 + p2 + '_Initial diffusion coefficients.csv', sep = ',', index = False)
    
## Calculate power-law exponents
for i in (point_traj['particle'].unique()): ## Go through all particles in frame
    length = round(iMSD[i].notna().sum() / 2) ## Bestimmt wie viele Datenpunkte es für Partikel i gibt. Ersetzt man die 1 durch 2 wird später nur aus der ersten Hälfte der Daten alpha berechnet.
    calc = (1/videofps*length) ## Setzt den Endwert für die lagtime beim berechnen.
    
    iCalc = iMSD.loc[tau1:calc,i] ## Macht dataframe aus allen vorher festgelegten Datenpunkten
    iFIT = tp.utils.fit_powerlaw(iCalc, plot = False)
    iAlpha = iFIT.iloc[0]['n']
    alphadata.append(iAlpha)
    overall_alpha.append(iAlpha)
    print('The power-law exponent of particle number',i,'is: ',iAlpha)
    
    ## Save the power-laws
    iAlphadata = pd.DataFrame(alphadata, columns = ['Power-law exponents']) ## Save power-law exponents
    iAlphadata.to_csv(outfolder_path_calc + '/' + p1 + p2 + '_Power-law exponents.csv', sep = ',', index = False)
    
    ## Plot the MSD plots
    if plot_single == True:
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5)) ## Create figure 'canvas'
        tp.utils.fit_powerlaw(iCalc, ax=ax1, plot = True) ## Plot the Powerlawfit
        tp.plot_traj(traj = point_traj.loc[point_traj['particle'] == i], ax=ax2, label = True, plot = True) ## Plot the single trajectory

        figure_path = outfolder_path_calc + '/Graphs/' ## Output and save
        if not os.path.exists(figure_path):
            os.makedirs(figure_path) ## Creates the actual folder
        plt.savefig(figure_path + 'Particle_ID_' + str(i) + '_Figure.png')
    
## Calculate all MSDs over the virus spot
iMSDspot = tp.imsd(traj_calc, mpp = 0.160, fps = videofps)
iMSDspot.to_csv(outfolder_path_calc + '/' + p1 + p2 + '_MSDs on spot.csv', sep = ',')

## Plot the MSDs over the virus spot
fig, ax = plt.subplots()
ax.plot(iMSDspot.index, iMSDspot, 'k-', alpha = 0.1)  # black lines, semitransparent
ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
xlabel='lag time $t$ [s]')
ax.set_xscale('log')
ax.set_yscale('log')

## Save the Plot
plt.savefig(outfolder_path_calc + '/' + p1 + p2 + '_MSDs on spot.png')

# %%

### Save the overall files for convenience

overall_diffu_name = outfolder_path + locfilename
overall_diffu_name = overall_diffu_name.replace('.csv', '_All diffusioncoefficients on spots.csv')

overall_diffu = pd.DataFrame(overall_diffu, columns = ['Initial diffusion coefficients']) ## Save Diffusioncoefficients
overall_diffu.to_csv(overall_diffu_name, sep = ',', index = False)

overall_alpha_name = overall_diffu_name.replace('_All diffusioncoefficients on spots.csv','_All power-law exponents on spots.csv')
overall_alpha = pd.DataFrame(overall_alpha, columns = ['Power-law exponents']) ## Save Diffusioncoefficients
overall_alpha.to_csv(overall_alpha_name, sep = ',', index = False)

#%%
## Only if you want to plot tracks again

tracks_only_plot = pd.read_csv('C:/Users/lbr20/Desktop/WD/lbr0026/013/Analysis/lbr0026_A549_EGFR-mEos_PR8_Sialidase_imageStack_013_drift corrected tracks.csv')
tracks_only_image = r'C:\Users\lbr20\Desktop\WD\lbr0026\013\lbr0026_A549_EGFR-mEos_PR8_Sialidase_Virus_013.nd2 - C=0_2.tif'
only_image = pims.open(tracks_only_image) ## Load in virus image
only_image = ndimage.rotate(only_image[0], 90) ## Rotate virus image (Needs to be done because of decode)
only_image = np.flipud(only_image) ## Mirror virus image (Needs to be done because of decode)
tp.plot_traj(tracks_only_plot, superimpose = only_image); ## With drift correction