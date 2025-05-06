'''
script purpose:             visualization of NED boxplots splitted into 4 time periods
outputs:                    time_period_NED_separated.png
modifications required:     commented with "MODIFY"

codeauthor: Fabian Weibel (2022)
adapted for publication by: Fiona Federer (2025)
'''

'''Import libraries'''
import os
import os.path
import random
from operator import add
from datetime import datetime, date, timedelta
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import seaborn as sns
from decimal import Decimal, ROUND_DOWN, ROUND_UP
from itertools import groupby
import ema_workbench
import csv as csv
import seaborn as sns

'''Load the data for CH2018: uncomment either the Basel or Zurich code block, change paths (MODIFY)'''
# Define the paths for CH2018 Basel
# city = 'Basel' # name used for directories
# station = 'Basel-Binningen' # name used for plots
# rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/26_30_bernard' 

# Define the paths for CH2018 Basel
city = 'Zurich' # name used for directories
station = 'Zurich-Fluntern' # name used for plots
rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/26_30_bernard' 

# Open time periods files that have been created in "time_periods_data_preparation.py"
df_26=pd.read_csv(os.path.join(rootOut, f'time_periods_NED_RCP26_{city}_2.csv'), sep=',') 
df_45=pd.read_csv(os.path.join(rootOut, f'time_periods_NED_RCP45_{city}_2.csv'), sep=',')
df_85=pd.read_csv(os.path.join(rootOut, f'time_periods_NED_RCP85_{city}_2.csv'), sep=',')


'''Boxplots number of days, 1981-2020, 2021-2050, 2051-2080, 2081-2099, separated'''
fig1, axs = plt.subplots(3,1, figsize=(6, 10), dpi=1000)
sns.set_style("darkgrid")

data26 = [df_26.iloc[:, 2:40].mean(), df_26.iloc[:, 40:70].mean(), df_26.iloc[:, 70:100].mean(), df_26.iloc[:, 100:].mean()]
data45 = [df_45.iloc[:, 2:40].mean(), df_45.iloc[:, 40:70].mean(), df_45.iloc[:, 70:100].mean(), df_45.iloc[:, 100:].mean()]
data85 = [df_85.iloc[:, 2:40].mean(), df_85.iloc[:, 40:70].mean(), df_85.iloc[:, 70:100].mean(), df_85.iloc[:, 100:].mean()]
labels = ['1981-2020', '2021-2050', '2051-2080', '2081-2099']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box26 = axs[0].boxplot(data26, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
box45 = axs[1].boxplot(data45, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
box85 = axs[2].boxplot(data85, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
axs[0].set_title(f'Average number of extreme days per year in {station} \n for four different time periods', size=12, pad=10)
axs[0].set_ylim(bottom=0, top =41) # Zurich: 41 / Basel: 51 (MODIFY)
axs[1].set_ylim(bottom=0, top =41) # Zurich: 41 / Basel: 51 (MODIFY)
axs[2].set_ylim(bottom=0, top =41) # Zurich: 41 / Basel: 51 (MODIFY)
axs[0].yaxis.set_ticks(np.arange(0, 41, 10)) # Zurich: 41 / Basel: 51 (MODIFY)
axs[1].yaxis.set_ticks(np.arange(0, 41, 10)) # Zurich: 41 / Basel: 51 (MODIFY)
axs[2].yaxis.set_ticks(np.arange(0, 41, 10)) # Zurich: 41 / Basel: 51 (MODIFY)
axs[2].set_xlabel('time period', size=10, labelpad=10)
axs[1].set_ylabel('# of days', size=10, labelpad=10)
axs[0].tick_params(axis='x', labelsize=8)
axs[1].tick_params(axis='x', labelsize=8)
axs[2].tick_params(axis='x', labelsize=8)
axs[0].tick_params(axis='y', labelsize=8)
axs[1].tick_params(axis='y', labelsize=8)
axs[2].tick_params(axis='y', labelsize=8)
axs[0].set_facecolor('#eaeaf2')
axs[1].set_facecolor('#eaeaf2')
axs[2].set_facecolor('#eaeaf2')
for pos in ['right', 'top', 'bottom', 'left']:
   axs[0].spines[pos].set_visible(False)
for pos in ['right', 'top', 'bottom', 'left']:
   axs[1].spines[pos].set_visible(False)
for pos in ['right', 'top', 'bottom', 'left']:
   axs[2].spines[pos].set_visible(False)
axs[0].annotate('RCP 2.6', fontsize=10, xy=(1, 45), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
axs[1].annotate('RCP 4.5', fontsize=10, xy=(1, 45), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
axs[2].annotate('RCP 8.5', fontsize=10, xy=(1, 45), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# Loop through all subplots (axs) and adjust gridlines
for ax in axs:
    ax.yaxis.grid(True)  # Enable horizontal gridlines
    ax.xaxis.grid(False)  # Disable vertical gridlines

# fill with colors
colors26 = ['pink', 'pink', 'pink', 'pink']
for patch, color in zip(box26['boxes'], colors26):
    patch.set_facecolor(color)
for element in ['boxes', 'whiskers', 'caps', 'medians']: 
    for line in box26[element]:
        line.set(linewidth=0.5)
for median in box26['medians']:
    median.set(color ='black',
               linewidth = 0.5)
colors45 = ['orange', 'orange', 'orange', 'orange']
for patch, color in zip(box45['boxes'], colors45):
    patch.set_facecolor(color)
for element in ['boxes', 'whiskers', 'caps', 'medians']: 
    for line in box45[element]:
        line.set(linewidth=0.5)
for median in box45['medians']:
    median.set(color ='black',
               linewidth = 0.5)
colors85 = ['red', 'red', 'red', 'red']
for patch, color in zip(box85['boxes'], colors85):
    patch.set_facecolor(color)
for element in ['boxes', 'whiskers', 'caps', 'medians']: 
    for line in box85[element]:
        line.set(linewidth=0.5)
for median in box85['medians']:
    median.set(color ='black',
               linewidth = 0.5)


fig1.savefig(os.path.join(rootOut, 'z_figures/time_period_NED_separated.png'), format='png', dpi=1000)

plt.show()
