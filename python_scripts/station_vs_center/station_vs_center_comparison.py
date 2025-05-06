'''
script purpose:             modification by station-center correction factor
                            visualization of the impacts on the heat indicators NED, NEE, LEE
outputs:                    'Basel_vs_Center_LEE_boxplot.png'
                            'Basel_vs_Center_NEE_boxplot.png'
                            'Basel_vs_Center_NED_boxplot.png'

modifications required:     commented with "MODIFY"

codeauthor: Fabian Weibel (2022)
adapted for publication by: Fiona Federer (2025)
'''

'''Import libraries'''
from itertools import groupby
import random
from operator import add
from datetime import datetime, date, timedelta
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import time
from numpy.polynomial.polynomial import polyfit
import os
import os.path
from decimal import Decimal, ROUND_DOWN, ROUND_UP
from itertools import groupby
from ema_workbench import load_results, ema_logging

'''Import the ema workbench outputs'''
from ema_workbench.analysis.pairs_plotting import (pairs_lines, pairs_scatter,
                                                   pairs_density)
ema_logging.log_to_stderr(level=ema_logging.DEFAULT_LEVEL)

###### Data preparation ######

'''Load the data for CH2018: uncomment either the Basel or station code block, change paths (MODIFY)'''
# Paths Basel
# station = 'Basel-Binningen'
# rootOut_station = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/26_30_bernard'
# rootOut_Center= r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/Corrfactor_Tmp_1_3_5'

# rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/sensitivity_analysis_figures'

# Paths Zurich
station = 'Zurich-Fluntern'
rootOut_station = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/26_30_bernard'
rootOut_Center = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Corrfactor_Tmp_1_3_5'

rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/sensitivity_analysis_figures'

out_zip_file_name_station = '20250301_runs.tar.gz'
out_zip_file_name_Center = '20250325_runs.tar.gz'

fh_station = os.path.join(rootOut_station, out_zip_file_name_station)
experiments_station, outcomes_station = load_results(fh_station)
fh_Center = os.path.join(rootOut_Center, out_zip_file_name_Center)
experiments_Center, outcomes_Center = load_results(fh_Center)

'''Create data frames for input data station'''
import tarfile
with tarfile.open(os.path.join(rootOut_station, out_zip_file_name_station),"r") as zip_ref_station:
    zip_ref_station.extractall(os.path.join(rootOut_station, out_zip_file_name_station[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_station, out_zip_file_name_station[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_station, out_zip_file_name_station[:-7]), file),
                os.path.join(os.path.join(rootOut_station, out_zip_file_name_station[:-7]), file[:-4] + ".csv"))

# We have 6 types of outputs
outDaily_station = os.path.join(rootOut_station, 'Outputs_py')
#M1
outSeasonTippingPoint_station = os.path.join(rootOut_station, 'outSeason')
outSeason_Length_station = os.path.join(rootOut_station, 'outSeason_length')
outSeason_StartEvents_station = os.path.join(rootOut_station, 'outSeason_StartEvents')
outSeason_EndEvents_station = os.path.join(rootOut_station, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_stationh = os.path.join(rootOut_station, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_station = os.path.join(rootOut_station, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_station = os.path.join(rootOut_station, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_station = os.path.join(rootOut_station, 'outSeason_Ave_StartEvents')
out_ema_station = os.path.join(rootOut_station, out_zip_file_name_station[:-7])

df4_ema_experiment_station = pd.read_csv(os.path.join(out_ema_station, 'experiments.csv'))
#M1 
df4_ema_y0_station = pd.read_csv(os.path.join(out_ema_station, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_station.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_station = pd.read_csv(os.path.join(out_ema_station, 'GCM_RCM.csv'), header=None)
df4_ema_y1_station.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_station = pd.read_csv(os.path.join(out_ema_station, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_station.columns = ["Yout2_Ave_Length"]
df4_ema_y3_station = pd.read_csv(os.path.join(out_ema_station, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_station.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_station = pd.read_csv(os.path.join(out_ema_station, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_station.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_station = pd.read_csv(os.path.join(out_ema_station, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_station.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_station = pd.concat((df4_ema_experiment_station, df4_ema_y0_station, df4_ema_y1_station, df4_ema_y2_station, df4_ema_y3_station, df4_ema_y4_station, df4_ema_y5_station), axis = 1)
# Taking the values and storing in separate variables (CH2018)
xRandomness_station = df_final_ema_station['Xfactor1'].values
xClimateModel_station = df_final_ema_station['xClimateModel'].values
xRCP_station = df_final_ema_station['xRCP'].values
x7Wbgtthreshold_station = df_final_ema_station['x7Wbgtthreshold'].values
#M1
y0_station = df_final_ema_station['Yout0_S_Ave_HotDay'].values
y1_station = df_final_ema_station['Yout1_GCM_RCM'].values
#M3
y2_station = df_final_ema_station['Yout2_Ave_Length'].values
y3_station = df_final_ema_station['Yout3_Ave_StartEvent'].values
y4_station = df_final_ema_station['Yout4_Ave_EndEvent'].values
#M2
y5_station = df_final_ema_station['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_station['xRCP']) == 1.0)
filt2 = (round(df_final_ema_station['xRCP']) == 2.0)
filt3 = (round(df_final_ema_station['xRCP']) == 3.0)
df_final_ema_f1_day_station = df_final_ema_station.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_station = df_final_ema_station.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_station = df_final_ema_station.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_station = df_final_ema_station.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_station = df_final_ema_station.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_station = df_final_ema_station.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_station = df_final_ema_station.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_station = df_final_ema_station.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_station = df_final_ema_station.loc[filt3, 'Yout2_Ave_Length']

'''Create data frames for input data Center'''
import tarfile
with tarfile.open(os.path.join(rootOut_Center, out_zip_file_name_Center),"r") as zip_ref_Center:
    zip_ref_Center.extractall(os.path.join(rootOut_Center, out_zip_file_name_Center[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_Center, out_zip_file_name_Center[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_Center, out_zip_file_name_Center[:-7]), file),
                os.path.join(os.path.join(rootOut_Center, out_zip_file_name_Center[:-7]), file[:-4] + ".csv"))

# We have 6 types of outputs
outDaily_Center = os.path.join(rootOut_Center, 'Outputs_py')
#M1
outSeasonTippingPoint_Center = os.path.join(rootOut_Center, 'outSeason')
outSeason_Length_Center = os.path.join(rootOut_Center, 'outSeason_length')
outSeason_StartEvents_Center = os.path.join(rootOut_Center, 'outSeason_StartEvents')
outSeason_EndEvents_Center = os.path.join(rootOut_Center, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Center = os.path.join(rootOut_Center, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Center = os.path.join(rootOut_Center, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Center = os.path.join(rootOut_Center, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Center = os.path.join(rootOut_Center, 'outSeason_Ave_StartEvents')
out_ema_Center = os.path.join(rootOut_Center, out_zip_file_name_Center[:-7])

df4_ema_experiment_Center = pd.read_csv(os.path.join(out_ema_Center, 'experiments.csv'))
#M1 
df4_ema_y0_Center = pd.read_csv(os.path.join(out_ema_Center, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Center.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Center = pd.read_csv(os.path.join(out_ema_Center, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Center.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Center = pd.read_csv(os.path.join(out_ema_Center, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Center.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Center = pd.read_csv(os.path.join(out_ema_Center, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Center.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Center = pd.read_csv(os.path.join(out_ema_Center, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Center.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Center = pd.read_csv(os.path.join(out_ema_Center, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Center.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Center = pd.concat((df4_ema_experiment_Center, df4_ema_y0_Center, df4_ema_y1_Center, df4_ema_y2_Center, df4_ema_y3_Center, df4_ema_y4_Center, df4_ema_y5_Center), axis = 1)
# Taking the values and storing in separate variables (CH2018)
xRandomness_Center = df_final_ema_Center['Xfactor1'].values
xClimateModel_Center = df_final_ema_Center['xClimateModel'].values
xRCP_Center = df_final_ema_Center['xRCP'].values
x7Wbgtthreshold_Center = df_final_ema_Center['x7Wbgtthreshold'].values
#M1
y0_Center = df_final_ema_Center['Yout0_S_Ave_HotDay'].values
y1_Center = df_final_ema_Center['Yout1_GCM_RCM'].values
#M3
y2_Center = df_final_ema_Center['Yout2_Ave_Length'].values
y3_Center = df_final_ema_Center['Yout3_Ave_StartEvent'].values
y4_Center = df_final_ema_Center['Yout4_Ave_EndEvent'].values
#M2
y5_Center = df_final_ema_Center['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Center['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Center['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Center['xRCP']) == 3.0)
df_final_ema_f1_day_Center = df_final_ema_Center.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Center = df_final_ema_Center.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Center = df_final_ema_Center.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_Center = df_final_ema_Center.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_Center = df_final_ema_Center.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_Center = df_final_ema_Center.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_Center = df_final_ema_Center.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_Center = df_final_ema_Center.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_Center = df_final_ema_Center.loc[filt3, 'Yout2_Ave_Length']

###### Visualization ######

'''Create Boxplots number extreme days'''
fig1, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data1 = [df_final_ema_f1_day_station, df_final_ema_f1_day_Center, df_final_ema_f2_day_station, df_final_ema_f2_day_Center, df_final_ema_f3_day_station, df_final_ema_f3_day_Center]
labels = ['RCP2.6', 'RCP2.6', 'RCP4.5', 'RCP4.5', 'RCP8.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['blue', 'red', 'blue', 'red', 'blue', 'red']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for element in ['boxes', 'whiskers', 'caps', 'medians']: 
    for line in box[element]:
        line.set(linewidth=0.5)
for spine in axs.spines.values(): 
    spine.set_linewidth(0.5)
for median in box['medians']:
    median.set(color ='black',
               linewidth = 0.5)
axs.set_title(f'Average number of extreme days per year. Station {station} vs. city center \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =60)
axs.set_ylabel('# of days', size=10, labelpad=10)
legend_elements = [Patch(facecolor='blue', label=station), Patch(facecolor='red', label='city center')]
axs.legend(handles=legend_elements, prop={'size': 8})
os.makedirs(rootOut, exist_ok=True)
fig1.savefig(os.path.join(rootOut, 'Station_vs_Center_NED_boxplot.png'), format='png', dpi=1000)


'''Create Boxplots number extreme events'''
fig2, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data1 = [df_final_ema_f1_station, df_final_ema_f1_Center, df_final_ema_f2_station, df_final_ema_f2_Center, df_final_ema_f3_station, df_final_ema_f3_Center]
labels = ['RCP2.6', 'RCP2.6', 'RCP4.5', 'RCP4.5', 'RCP8.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['blue', 'red', 'blue', 'red', 'blue', 'red']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for element in ['boxes', 'whiskers', 'caps', 'medians']: 
    for line in box[element]:
        line.set(linewidth=0.5)
for spine in axs.spines.values(): 
    spine.set_linewidth(0.5)
for median in box['medians']:
    median.set(color ='black',
               linewidth = 0.5)
axs.set_title(f'Average number of extreme events per year. Station {station} vs. city center \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =17)
axs.set_ylabel('# of events', size=10, labelpad=10)
legend_elements = [Patch(facecolor='blue', label=station), Patch(facecolor='red', label='City center')]
axs.legend(handles=legend_elements, prop={'size': 8})
fig2.savefig(os.path.join(rootOut, 'Station_vs_Center_NEE_boxplot.png'), format='png', dpi=1000)

'''Create Boxplots length extreme events'''
fig3, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data1 = [df_final_ema_f1_L_station, df_final_ema_f1_L_Center, df_final_ema_f2_L_station, df_final_ema_f2_L_Center, df_final_ema_f3_L_station, df_final_ema_f3_L_Center]
labels = ['RCP2.6', 'RCP2.6', 'RCP4.5', 'RCP4.5', 'RCP8.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['blue', 'red', 'blue', 'red', 'blue', 'red']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for element in ['boxes', 'whiskers', 'caps', 'medians']: 
    for line in box[element]:
        line.set(linewidth=0.5)
for spine in axs.spines.values(): 
    spine.set_linewidth(0.5)
for median in box['medians']:
    median.set(color ='black',
               linewidth = 0.5)
axs.set_title(f'Average length of extreme events per year. Station {station} vs. city center \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =6)
axs.set_ylabel('# of days', size=10, labelpad=10)
legend_elements = [Patch(facecolor='blue', label=station), Patch(facecolor='red', label='City center')]
axs.legend(handles=legend_elements, prop={'size': 8})
fig3.savefig(os.path.join(rootOut, 'Station_vs_Center_LEE_boxplot.png'), format='png', dpi=1000)

'''Show plots'''
plt.show()

print(df_final_ema_f3_L_station.mean())
print(df_final_ema_f3_L_Center.mean())
