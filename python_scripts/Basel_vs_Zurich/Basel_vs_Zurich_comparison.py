'''
script purpose:             visualization of three different between the heat indicators NED, NEE, LEE in Basel and Zurich
outputs:                    'Basel_vs_Zurich_NED_boxplot.png'
                            'Basel_vs_Zurich_NEE_boxplot.png'
                            'Basel_vs_Zurich_LEE_boxplot.png'
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
import tarfile

'''Import the ema workbench outputs'''
from ema_workbench.analysis.pairs_plotting import (pairs_lines, pairs_scatter,pairs_density)
ema_logging.log_to_stderr(level=ema_logging.DEFAULT_LEVEL)

###### Data preparation ######
'''Load the data for CH2018'''
rootOut_Basel = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/26_30_bernard'
rootOut_Zurich = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/26_30_bernard'
rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/Zurich_vs_Basel_figures'

out_zip_file_name_Basel = '20250301_runs.tar.gz'
out_zip_file_name_Zurich = '20250301_runs.tar.gz'

fh_Basel = os.path.join(rootOut_Basel, out_zip_file_name_Basel)
experiments_Basel, outcomes_Basel = load_results(fh_Basel)
fh_Zurich = os.path.join(rootOut_Zurich, out_zip_file_name_Zurich)
experiments_Zurich, outcomes_Zurich = load_results(fh_Zurich)

'''Create data frames for input data Basel'''
with tarfile.open(os.path.join(rootOut_Basel, out_zip_file_name_Basel),"r") as zip_ref_Basel:
    zip_ref_Basel.extractall(os.path.join(rootOut_Basel, out_zip_file_name_Basel[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for filename in os.listdir(os.path.join(rootOut_Basel, out_zip_file_name_Basel[:-7])):
        if filename.endswith('.cls'):
            os.rename(os.path.join(rootOut_Basel, out_zip_file_name_Basel[:-7], filename), os.path.join(rootOut_Basel, out_zip_file_name_Basel[:-7], filename[:-4]+'.csv'))
# We have 6 types of outputs
outDaily_Basel = os.path.join(rootOut_Basel, 'Outputs_py')
#M1
outSeasonTippingPoint_Basel = os.path.join(rootOut_Basel, 'outSeason')
outSeason_Length_Basel = os.path.join(rootOut_Basel, 'outSeason_length')
outSeason_StartEvents_Basel = os.path.join(rootOut_Basel, 'outSeason_StartEvents')
outSeason_EndEvents_Basel = os.path.join(rootOut_Basel, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Basel = os.path.join(rootOut_Basel, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Basel = os.path.join(rootOut_Basel, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Basel = os.path.join(rootOut_Basel, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Basel = os.path.join(rootOut_Basel, 'outSeason_Ave_StartEvents')
out_ema_Basel = os.path.join(rootOut_Basel, out_zip_file_name_Basel[:-7])

df4_ema_experiment_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'experiments.csv'))
#M1 
df4_ema_y0_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Basel.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Basel.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Basel.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Basel.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Basel.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Basel = pd.read_csv(os.path.join(out_ema_Basel, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Basel.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Basel = pd.concat((df4_ema_experiment_Basel, df4_ema_y0_Basel, df4_ema_y1_Basel, df4_ema_y2_Basel, df4_ema_y3_Basel, df4_ema_y4_Basel, df4_ema_y5_Basel), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Basel = df_final_ema_Basel['Xfactor1'].values
xClimateModel_Basel = df_final_ema_Basel['xClimateModel'].values
xRCP_Basel = df_final_ema_Basel['xRCP'].values
x7Wbgtthreshold_Basel = df_final_ema_Basel['x7Wbgtthreshold'].values
#M1
y0_Basel = df_final_ema_Basel['Yout0_S_Ave_HotDay'].values
y1_Basel = df_final_ema_Basel['Yout1_GCM_RCM'].values
#M3
y2_Basel = df_final_ema_Basel['Yout2_Ave_Length'].values
y3_Basel = df_final_ema_Basel['Yout3_Ave_StartEvent'].values
y4_Basel = df_final_ema_Basel['Yout4_Ave_EndEvent'].values
#M2
y5_Basel = df_final_ema_Basel['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Basel['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Basel['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Basel['xRCP']) == 3.0)
df_final_ema_f1_day_Basel = df_final_ema_Basel.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Basel = df_final_ema_Basel.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Basel = df_final_ema_Basel.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_Basel = df_final_ema_Basel.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_Basel = df_final_ema_Basel.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_Basel = df_final_ema_Basel.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_Basel = df_final_ema_Basel.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_Basel = df_final_ema_Basel.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_Basel = df_final_ema_Basel.loc[filt3, 'Yout2_Ave_Length']

'''Create data frames for input data Zurich'''
with tarfile.open(os.path.join(rootOut_Zurich, out_zip_file_name_Zurich),"r") as zip_ref_Zurich:
    zip_ref_Zurich.extractall(os.path.join(rootOut_Zurich, out_zip_file_name_Zurich[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for filename in os.listdir(os.path.join(rootOut_Zurich, out_zip_file_name_Zurich[:-7])):
        if filename.endswith('.cls'):
            os.rename(os.path.join(rootOut_Zurich, out_zip_file_name_Zurich[:-7], filename), os.path.join(rootOut_Zurich, out_zip_file_name_Zurich[:-7], filename[:-4]+'.csv'))
# We have 6 types of outputs
outDaily_Zurich = os.path.join(rootOut_Zurich, 'Outputs_py')
#M1
outSeasonTippingPoint_Zurich = os.path.join(rootOut_Zurich, 'outSeason')
outSeason_Length_Zurich = os.path.join(rootOut_Zurich, 'outSeason_length')
outSeason_StartEvents_Zurich = os.path.join(rootOut_Zurich, 'outSeason_StartEvents')
outSeason_EndEvents_Zurich = os.path.join(rootOut_Zurich, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Zurich = os.path.join(rootOut_Zurich, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Zurich = os.path.join(rootOut_Zurich, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Zurich = os.path.join(rootOut_Zurich, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Zurich = os.path.join(rootOut_Zurich, 'outSeason_Ave_StartEvents')
out_ema_Zurich = os.path.join(rootOut_Zurich, out_zip_file_name_Zurich[:-7])

df4_ema_experiment_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'experiments.csv'))
#M1 
df4_ema_y0_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Zurich.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Zurich.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Zurich.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Zurich.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Zurich.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Zurich = pd.read_csv(os.path.join(out_ema_Zurich, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Zurich.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Zurich = pd.concat((df4_ema_experiment_Zurich, df4_ema_y0_Zurich, df4_ema_y1_Zurich, df4_ema_y2_Zurich, df4_ema_y3_Zurich, df4_ema_y4_Zurich, df4_ema_y5_Zurich), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Zurich = df_final_ema_Zurich['Xfactor1'].values
xClimateModel_Zurich = df_final_ema_Zurich['xClimateModel'].values
xRCP_Zurich = df_final_ema_Zurich['xRCP'].values
x7Wbgtthreshold_Zurich = df_final_ema_Zurich['x7Wbgtthreshold'].values
#M1
y0_Zurich = df_final_ema_Zurich['Yout0_S_Ave_HotDay'].values
y1_Zurich = df_final_ema_Zurich['Yout1_GCM_RCM'].values
#M3
y2_Zurich = df_final_ema_Zurich['Yout2_Ave_Length'].values
y3_Zurich = df_final_ema_Zurich['Yout3_Ave_StartEvent'].values
y4_Zurich = df_final_ema_Zurich['Yout4_Ave_EndEvent'].values
#M2
y5_Zurich = df_final_ema_Zurich['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Zurich['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Zurich['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Zurich['xRCP']) == 3.0)
df_final_ema_f1_day_Zurich = df_final_ema_Zurich.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Zurich = df_final_ema_Zurich.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Zurich = df_final_ema_Zurich.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_Zurich = df_final_ema_Zurich.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_Zurich = df_final_ema_Zurich.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_Zurich = df_final_ema_Zurich.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_Zurich = df_final_ema_Zurich.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_Zurich = df_final_ema_Zurich.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_Zurich = df_final_ema_Zurich.loc[filt3, 'Yout2_Ave_Length']


###### Visualization ######

'''Create Boxplots number extreme days'''
fig1, axs = plt.subplots(figsize=(8,5), dpi=1000)
data1 = [df_final_ema_f1_day_Basel, df_final_ema_f1_day_Zurich, df_final_ema_f2_day_Basel, df_final_ema_f2_day_Zurich, df_final_ema_f3_day_Basel, df_final_ema_f3_day_Zurich]
labels = ['RCP2.6', 'RCP2.6', 'RCP4.5', 'RCP4.5', 'RCP8.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['blue', 'green', 'blue', 'green', 'blue', 'green']
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
axs.set_title(f'Average number of extreme days per year. Basel-Binningen vs. Zurich-Fluntern \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =40)
axs.set_ylabel('# of days', size=10, labelpad=10)
legend_elements = [Patch(facecolor='blue', label='Basel-Binningen'), Patch(facecolor='green', label='Zurich-Fluntern')]
axs.legend(handles=legend_elements, prop={'size': 10})
fig1.savefig(os.path.join(rootOut, 'Basel_vs_Zurich_NED_boxplot.png'), format='png', dpi=1000)

'''Create Boxplots number extreme events'''
fig2, axs = plt.subplots(figsize=(8,5), dpi=1000)
data1 = [df_final_ema_f1_Basel, df_final_ema_f1_Zurich, df_final_ema_f2_Basel, df_final_ema_f2_Zurich, df_final_ema_f3_Basel, df_final_ema_f3_Zurich]
labels = ['RCP2.6', 'RCP2.6', 'RCP4.5', 'RCP4.5', 'RCP8.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['blue', 'green', 'blue', 'green', 'blue', 'green']
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
axs.set_title('Average number of extreme events per year. Basel-Binningen vs. Zurich-Fluntern \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =15)
axs.set_ylabel('# of events', size=10, labelpad=10)
legend_elements = [Patch(facecolor='blue', label='Basel-Binningen'), Patch(facecolor='green', label='Zurich-Fluntern')]
axs.legend(handles=legend_elements, prop={'size': 10})
fig2.savefig(os.path.join(rootOut, 'Basel_vs_Zurich_NEE_boxplot.png'), format='png', dpi=1000)

'''Create Boxplots length extreme events'''
fig3, axs = plt.subplots(figsize=(8,5), dpi=1000)
data1 = [df_final_ema_f1_L_Basel, df_final_ema_f1_L_Zurich, df_final_ema_f2_L_Basel, df_final_ema_f2_L_Zurich, df_final_ema_f3_L_Basel, df_final_ema_f3_L_Zurich]
labels = ['RCP2.6', 'RCP2.6', 'RCP4.5', 'RCP4.5', 'RCP8.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['blue', 'green', 'blue', 'green', 'blue', 'green']
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
axs.set_title('Average length of extreme events per year. Basel-Binningen vs. Zurich-Fluntern \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =5)
axs.set_ylabel('# of days', size=10, labelpad=10)
legend_elements = [Patch(facecolor='blue', label='Basel-Binningen'), Patch(facecolor='green', label='Zurich-Fluntern')]
axs.legend(handles=legend_elements, prop={'size': 10})
fig3.savefig(os.path.join(rootOut, 'Basel_vs_Zurich_LEE_boxplot.png'), format='png', dpi=1000)

'''Show plots'''
plt.show()


