# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 14:43:20 2023

@author: lukas
"""
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('default')
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

#%%

master_folder = r'C:/Users/lbr20/Desktop/WD/lbr0023 006/'

data_confined = pd.read_csv(master_folder + 'x194y149_confined iMSDs on spot.csv', index_col = 0)
n1 = int(len(data_confined)/2) - 0
data_confined_half = data_confined.iloc[:n1,:]

data_free = pd.read_csv(master_folder + 'x194y149_Free iMSDs on spot.csv', index_col = 0)
n2 = int(len(data_free)/2) - 0
data_free_half = data_free.iloc[:n2,:]

mean_confined = data_confined_half.mean(axis=1)
std_confined = data_confined_half.std(axis=1)

mean_free = data_free_half.mean(axis=1)
std_free = data_free_half.std(axis=1)

#%%
data_confined_2 = pd.read_csv(master_folder + 'x194y168_confined iMSDs on spot.csv', index_col = 0)
n1 = int(len(data_confined)/2) - 0
data_confined_half_2 = data_confined_2.iloc[:n1,:]

data_free_2 = pd.read_csv(master_folder + 'x194y168_Free iMSDs on spot.csv', index_col = 0)
n2 = int(len(data_free)/2) - 0
data_free_half_2 = data_free_2.iloc[:n2,:]

mean_confined_2 = data_confined_half_2.mean(axis=1)
std_confined_2 = data_confined_half_2.std(axis=1)

mean_free_2 = data_free_half_2.mean(axis=1)
std_free_2 = data_free_half_2.std(axis=1)


# %%
fig, ax = plt.subplots()
ax.plot(mean_confined.index, mean_confined, 'k-', alpha = 0.9, label = 'Immobile over Virus')
#ax.fill_between(mean_confined.index, mean_confined - std_confined, mean_confined + std_confined,
               #alpha = 0.25, facecolor='#000000')

ax.plot(mean_free.index, mean_free, 'r-', alpha = 0.9, label = 'Mobile over Virus')
#ax.fill_between(mean_free.index, mean_free - std_free, mean_free + std_free,
               #alpha = 0.25, facecolor='#FF0000')

ax.plot(mean_confined_2.index, mean_confined_2, 'b-', alpha = 0.9, label = 'Immobile Background')
#ax.fill_between(mean_confined_2.index, mean_confined_2 - std_confined_2, mean_confined_2 + std_confined_2,
               #alpha = 0.25, facecolor='#F0FFFF')

ax.plot(mean_free_2.index, mean_free_2, 'c-', alpha = 0.9, label = 'Mobile Background')
#ax.fill_between(mean_free_2.index, mean_free_2 - std_free_2, mean_free_2 + std_free_2,
               #alpha = 0.25, facecolor='#00FFFF')

ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
       xlabel='lag time $t$ [s]')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlim(left = 0.03, right = 0.2)
ax.set_ylim(bottom = 0, top = 0.03)
ax.legend()
plt.savefig(master_folder + 'Overall_iMSD plot_Immobile vs Mobile.png', format ='png', dpi=1200);