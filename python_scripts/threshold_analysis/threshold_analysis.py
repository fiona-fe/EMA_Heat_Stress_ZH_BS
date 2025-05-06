'''
script purpose:             visualization of three different WBGT threshold ranges and its impacts on the heat indicators NED, NEE, LEE
outputs:                    '26_30_vs_28_32_vs_33_35_number_of_extreme_days_boxplot.png', 
                            '26_30_vs_28_32_number_of_extreme_events_boxplot.png'
                            '26_30_vs_28_32_length_of_extreme_events_boxplot.png'
                            '33_35_number_of_extreme_days_boxplot.png'
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
import seaborn as sns
from decimal import Decimal, ROUND_DOWN, ROUND_UP
from itertools import groupby
import ema_workbench
import tarfile
from ema_workbench import load_results, ema_logging

'''Import the ema workbench outputs'''
from ema_workbench.analysis.pairs_plotting import (pairs_lines, pairs_scatter,
                                                   pairs_density)
ema_logging.log_to_stderr(level=ema_logging.DEFAULT_LEVEL)

###### Data preparation ######

'''Load the data for CH2018: uncomment either the Basel or Zurich code block, change paths (MODIFY)'''
# Paths Basel
# station = 'Basel-Binningen'
# rootOut_26_30 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/26_30_bernard'
# rootOut_28_32 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/28_32_bernard'
# rootOut_33_35 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/33_35_bernard'

# rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/threshold_analysis_figures'

# Paths Zurich
station = 'Zurich-Fluntern'
rootOut_26_30 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/26_30_bernard'
rootOut_28_32 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/28_32_bernard'
rootOut_33_35 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/33_35_bernard'

rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/threshold_analysis_figures'

out_zip_file_name_26_30 = '20250301_runs.tar.gz'
out_zip_file_name_28_32 = '20250309_runs.tar.gz'
out_zip_file_name_33_35 = '20250309_runs.tar.gz'

fh_26_30 = os.path.join(rootOut_26_30, out_zip_file_name_26_30)
experiments_26_30, outcomes_26_30 = load_results(fh_26_30)
fh_28_32 = os.path.join(rootOut_28_32, out_zip_file_name_28_32)
experiments_28_32, outcomes_28_32 = load_results(fh_28_32)
fh_33_35 = os.path.join(rootOut_33_35, out_zip_file_name_33_35)
experiments_33_35, outcomes_33_35 = load_results(fh_33_35)

'''Create data frames for input data 26_30'''
with tarfile.open(os.path.join(rootOut_26_30, out_zip_file_name_26_30),"r") as zip_ref_26_30:
    zip_ref_26_30.extractall(os.path.join(rootOut_26_30, out_zip_file_name_26_30[:-7]))
	# rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_26_30, out_zip_file_name_26_30[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_26_30, out_zip_file_name_26_30[:-7]), file),
                os.path.join(os.path.join(rootOut_26_30, out_zip_file_name_26_30[:-7]), file[:-4] + ".csv"))

# We have 6 types of outputs
outDaily_26_30 = os.path.join(rootOut_26_30, 'Outputs_py')
#M1
outSeasonTippingPoint_26_30 = os.path.join(rootOut_26_30, 'outSeason')
outSeason_Length_26_30 = os.path.join(rootOut_26_30, 'outSeason_length')
outSeason_StartEvents_26_30 = os.path.join(rootOut_26_30, 'outSeason_StartEvents')
outSeason_EndEvents_26_30 = os.path.join(rootOut_26_30, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_26_30 = os.path.join(rootOut_26_30, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_26_30 = os.path.join(rootOut_26_30, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_26_30 = os.path.join(rootOut_26_30, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_26_30 = os.path.join(rootOut_26_30, 'outSeason_Ave_StartEvents')
out_ema_26_30 = os.path.join(rootOut_26_30, out_zip_file_name_26_30[:-7])

df4_ema_experiment_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'experiments.csv'))
#M1
df4_ema_y0_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_26_30.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'GCM_RCM.csv'), header=None)
df4_ema_y1_26_30.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_26_30.columns = ["Yout2_Ave_Length"]
df4_ema_y3_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_26_30.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_26_30.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_26_30 = pd.read_csv(os.path.join(out_ema_26_30, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_26_30.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_26_30 = pd.concat((df4_ema_experiment_26_30, df4_ema_y0_26_30, df4_ema_y1_26_30, df4_ema_y2_26_30, df4_ema_y3_26_30, df4_ema_y4_26_30, df4_ema_y5_26_30), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_26_30 = df_final_ema_26_30['Xfactor1'].values
xClimateModel_26_30 = df_final_ema_26_30['xClimateModel'].values
xRCP_26_30 = df_final_ema_26_30['xRCP'].values
x7Wbgtthreshold_26_30 = df_final_ema_26_30['x7Wbgtthreshold'].values
#M1
y0_26_30 = df_final_ema_26_30['Yout0_S_Ave_HotDay'].values
y1_26_30 = df_final_ema_26_30['Yout1_GCM_RCM'].values
#M3
y2_26_30 = df_final_ema_26_30['Yout2_Ave_Length'].values
y3_26_30 = df_final_ema_26_30['Yout3_Ave_StartEvent'].values
y4_26_30 = df_final_ema_26_30['Yout4_Ave_EndEvent'].values
#M2
y5_26_30 = df_final_ema_26_30['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_26_30['xRCP']) == 1.0)
filt2 = (round(df_final_ema_26_30['xRCP']) == 2.0)
filt3 = (round(df_final_ema_26_30['xRCP']) == 3.0)
df_final_ema_f1_day_26_30 = df_final_ema_26_30.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_26_30 = df_final_ema_26_30.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_26_30 = df_final_ema_26_30.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_26_30 = df_final_ema_26_30.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_26_30 = df_final_ema_26_30.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_26_30 = df_final_ema_26_30.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_26_30 = df_final_ema_26_30.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_26_30 = df_final_ema_26_30.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_26_30 = df_final_ema_26_30.loc[filt3, 'Yout2_Ave_Length']

'''Create data frames for input data 28_32'''
import tarfile
with tarfile.open(os.path.join(rootOut_28_32, out_zip_file_name_28_32),"r") as zip_ref_28_32:
    zip_ref_28_32.extractall(os.path.join(rootOut_28_32, out_zip_file_name_28_32[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_28_32, out_zip_file_name_28_32[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_28_32, out_zip_file_name_28_32[:-7]), file),
                os.path.join(os.path.join(rootOut_28_32, out_zip_file_name_28_32[:-7]), file[:-4] + ".csv"))
# We have 6 types of outputs
outDaily_28_32 = os.path.join(rootOut_28_32, 'Outputs_py')
#M1
outSeasonTippingPoint_28_32 = os.path.join(rootOut_28_32, 'outSeason')
outSeason_Length_28_32 = os.path.join(rootOut_28_32, 'outSeason_length')
outSeason_StartEvents_28_32 = os.path.join(rootOut_28_32, 'outSeason_StartEvents')
outSeason_EndEvents_28_32 = os.path.join(rootOut_28_32, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_28_32 = os.path.join(rootOut_28_32, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_28_32 = os.path.join(rootOut_28_32, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_28_32 = os.path.join(rootOut_28_32, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_28_32 = os.path.join(rootOut_28_32, 'outSeason_Ave_StartEvents')
out_ema_28_32 = os.path.join(rootOut_28_32, out_zip_file_name_28_32[:-7])

df4_ema_experiment_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'experiments.csv'))
#M1 
df4_ema_y0_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_28_32.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'GCM_RCM.csv'), header=None)
df4_ema_y1_28_32.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_28_32.columns = ["Yout2_Ave_Length"]
df4_ema_y3_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_28_32.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_28_32.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_28_32 = pd.read_csv(os.path.join(out_ema_28_32, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_28_32.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_28_32 = pd.concat((df4_ema_experiment_28_32, df4_ema_y0_28_32, df4_ema_y1_28_32, df4_ema_y2_28_32, df4_ema_y3_28_32, df4_ema_y4_28_32, df4_ema_y5_28_32), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_28_32 = df_final_ema_28_32['Xfactor1'].values
xClimateModel_28_32 = df_final_ema_28_32['xClimateModel'].values
xRCP_28_32 = df_final_ema_28_32['xRCP'].values
x7Wbgtthreshold_28_32 = df_final_ema_28_32['x7Wbgtthreshold'].values
#M1
y0_28_32 = df_final_ema_28_32['Yout0_S_Ave_HotDay'].values
y1_28_32 = df_final_ema_28_32['Yout1_GCM_RCM'].values
#M3
y2_28_32 = df_final_ema_28_32['Yout2_Ave_Length'].values
y3_28_32 = df_final_ema_28_32['Yout3_Ave_StartEvent'].values
y4_28_32 = df_final_ema_28_32['Yout4_Ave_EndEvent'].values
#M2
y5_28_32 = df_final_ema_28_32['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_28_32['xRCP']) == 1.0)
filt2 = (round(df_final_ema_28_32['xRCP']) == 2.0)
filt3 = (round(df_final_ema_28_32['xRCP']) == 3.0)
df_final_ema_f1_day_28_32 = df_final_ema_28_32.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_28_32 = df_final_ema_28_32.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_28_32 = df_final_ema_28_32.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_28_32 = df_final_ema_28_32.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_28_32 = df_final_ema_28_32.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_28_32 = df_final_ema_28_32.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_28_32 = df_final_ema_28_32.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_28_32 = df_final_ema_28_32.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_28_32 = df_final_ema_28_32.loc[filt3, 'Yout2_Ave_Length']

'''Create data frames for input data 33_35'''
import tarfile
with tarfile.open(os.path.join(rootOut_33_35, out_zip_file_name_33_35),"r") as zip_ref_33_35:
    zip_ref_33_35.extractall(os.path.join(rootOut_33_35, out_zip_file_name_33_35[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_33_35, out_zip_file_name_33_35[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_33_35, out_zip_file_name_33_35[:-7]), file),
                os.path.join(os.path.join(rootOut_33_35, out_zip_file_name_33_35[:-7]), file[:-4] + ".csv"))
# We have 6 types of outputs
outDaily_33_35 = os.path.join(rootOut_33_35, 'Outputs_py')
#M1
outSeasonTippingPoint_33_35 = os.path.join(rootOut_33_35, 'outSeason')
outSeason_Length_33_35 = os.path.join(rootOut_33_35, 'outSeason_length')
outSeason_StartEvents_33_35 = os.path.join(rootOut_33_35, 'outSeason_StartEvents')
outSeason_EndEvents_33_35 = os.path.join(rootOut_33_35, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_33_35 = os.path.join(rootOut_33_35, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_33_35 = os.path.join(rootOut_33_35, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_33_35 = os.path.join(rootOut_33_35, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_33_35 = os.path.join(rootOut_33_35, 'outSeason_Ave_StartEvents')
out_ema_33_35 = os.path.join(rootOut_33_35, out_zip_file_name_33_35[:-7])

df4_ema_experiment_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'experiments.csv'))
#M1 
df4_ema_y0_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_33_35.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'GCM_RCM.csv'), header=None)
df4_ema_y1_33_35.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_33_35.columns = ["Yout2_Ave_Length"]
df4_ema_y3_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_33_35.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_33_35.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_33_35 = pd.read_csv(os.path.join(out_ema_33_35, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_33_35.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_33_35 = pd.concat((df4_ema_experiment_33_35, df4_ema_y0_33_35, df4_ema_y1_33_35, df4_ema_y2_33_35, df4_ema_y3_33_35, df4_ema_y4_33_35, df4_ema_y5_33_35), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_33_35 = df_final_ema_33_35['Xfactor1'].values
xClimateModel_33_35 = df_final_ema_33_35['xClimateModel'].values
xRCP_33_35 = df_final_ema_33_35['xRCP'].values
x7Wbgtthreshold_33_35 = df_final_ema_33_35['x7Wbgtthreshold'].values
#M1
y0_33_35 = df_final_ema_33_35['Yout0_S_Ave_HotDay'].values
y1_33_35 = df_final_ema_33_35['Yout1_GCM_RCM'].values
#M3
y2_33_35 = df_final_ema_33_35['Yout2_Ave_Length'].values
y3_33_35 = df_final_ema_33_35['Yout3_Ave_StartEvent'].values
y4_33_35 = df_final_ema_33_35['Yout4_Ave_EndEvent'].values
#M2
y5_33_35 = df_final_ema_33_35['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_33_35['xRCP']) == 1.0)
filt2 = (round(df_final_ema_33_35['xRCP']) == 2.0)
filt3 = (round(df_final_ema_33_35['xRCP']) == 3.0)
df_final_ema_f1_day_33_35 = df_final_ema_33_35.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_33_35 = df_final_ema_33_35.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_33_35 = df_final_ema_33_35.loc[filt3, 'Yout0_S_Ave_HotDay']
df_final_ema_f1_33_35 = df_final_ema_33_35.loc[filt1, 'Yout5_S_Ave_numEvent']
df_final_ema_f2_33_35 = df_final_ema_33_35.loc[filt2, 'Yout5_S_Ave_numEvent']
df_final_ema_f3_33_35 = df_final_ema_33_35.loc[filt3, 'Yout5_S_Ave_numEvent']
df_final_ema_f1_L_33_35 = df_final_ema_33_35.loc[filt1, 'Yout2_Ave_Length']
df_final_ema_f2_L_33_35 = df_final_ema_33_35.loc[filt2, 'Yout2_Ave_Length']
df_final_ema_f3_L_33_35 = df_final_ema_33_35.loc[filt3, 'Yout2_Ave_Length']

'''Making matrices for 26_30'''
outFolder_26_30 = rootOut_26_30
outSeasonFolder_26_30 = os.path.join(outFolder_26_30, 'outSeason')
all_Files_26_30 = []
for filename in os.walk(outSeasonFolder_26_30):
    all_Files_26_30 = filename[2]
    
totalFiles_loc_26_30 = []
for i in range(len(all_Files_26_30)):
    totalFiles_loc_26_30.append(os.path.join(outSeasonFolder_26_30, all_Files_26_30[i]))

emptyMatrixHotDays_26_30 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_26_30))]
    
for i in range(0, len(totalFiles_loc_26_30), 1):
    b_26_30 = []
    with open(totalFiles_loc_26_30[i], 'r') as file:
        outputReaderlines_26_30 = file.readlines()
        for j in range (len(outputReaderlines_26_30)):
            b_26_30.append(outputReaderlines_26_30[j].replace('\n','').split(','))
    emptyMatrixHotDays_26_30[i] = np.array(b_26_30)

matrix_HotDays_26_30 = np.stack(emptyMatrixHotDays_26_30[:len(totalFiles_loc_26_30)], axis=0)

matrix_HotDays26_26_30= np.stack([arr for arr in emptyMatrixHotDays_26_30[:len(totalFiles_loc_26_30)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_26_30= np.stack([arr for arr in emptyMatrixHotDays_26_30[:len(totalFiles_loc_26_30)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_26_30= np.stack([arr for arr in emptyMatrixHotDays_26_30[:len(totalFiles_loc_26_30)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for 28_32'''
outFolder_28_32 = rootOut_28_32
outSeasonFolder_28_32 = os.path.join(outFolder_28_32, 'outSeason')
all_Files_28_32 = []
for filename in os.walk(outSeasonFolder_28_32):
    all_Files_28_32 = filename[2]
    
totalFiles_loc_28_32 = []
for i in range(len(all_Files_28_32)):
    totalFiles_loc_28_32.append(os.path.join(outSeasonFolder_28_32, all_Files_28_32[i]))

emptyMatrixHotDays_28_32 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_28_32))]
    
for i in range(0, len(totalFiles_loc_28_32), 1):
    b_28_32 = []
    with open(totalFiles_loc_28_32[i], 'r') as file:
        outputReaderlines_28_32 = file.readlines()
        for j in range (len(outputReaderlines_28_32)):
            b_28_32.append(outputReaderlines_28_32[j].replace('\n','').split(','))
    emptyMatrixHotDays_28_32[i] = np.array(b_28_32)

matrix_HotDays_28_32 = np.stack(emptyMatrixHotDays_28_32[:len(totalFiles_loc_28_32)], axis=0)

matrix_HotDays26_28_32= np.stack([arr for arr in emptyMatrixHotDays_28_32[:len(totalFiles_loc_28_32)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_28_32= np.stack([arr for arr in emptyMatrixHotDays_28_32[:len(totalFiles_loc_28_32)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_28_32= np.stack([arr for arr in emptyMatrixHotDays_28_32[:len(totalFiles_loc_28_32)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for 33_35'''
outFolder_33_35 = rootOut_33_35
outSeasonFolder_33_35 = os.path.join(outFolder_33_35, 'outSeason')
all_Files_33_35 = []
for filename in os.walk(outSeasonFolder_33_35):
    all_Files_33_35 = filename[2]
    
totalFiles_loc_33_35 = []
for i in range(len(all_Files_33_35)):
    totalFiles_loc_33_35.append(os.path.join(outSeasonFolder_33_35, all_Files_33_35[i]))

emptyMatrixHotDays_33_35 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_33_35))]
    
for i in range(0, len(totalFiles_loc_33_35), 1):
    b_33_35 = []
    with open(totalFiles_loc_33_35[i], 'r') as file:
        outputReaderlines_33_35 = file.readlines()
        for j in range (len(outputReaderlines_33_35)):
            b_33_35.append(outputReaderlines_33_35[j].replace('\n','').split(','))
    emptyMatrixHotDays_33_35[i] = np.array(b_33_35)

matrix_HotDays_33_35 = np.stack(emptyMatrixHotDays_33_35[:len(totalFiles_loc_33_35)], axis=0)

matrix_HotDays26_33_35= np.stack([arr for arr in emptyMatrixHotDays_33_35[:len(totalFiles_loc_33_35)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_33_35= np.stack([arr for arr in emptyMatrixHotDays_33_35[:len(totalFiles_loc_33_35)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_33_35= np.stack([arr for arr in emptyMatrixHotDays_33_35[:len(totalFiles_loc_33_35)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


###### Visualization ######

'''Create Boxplots number of extreme days'''
fig1, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data1 = [df_final_ema_f1_day_26_30, df_final_ema_f1_day_28_32, df_final_ema_f1_day_33_35, df_final_ema_f2_day_26_30, df_final_ema_f2_day_28_32, df_final_ema_f2_day_33_35, df_final_ema_f3_day_26_30, df_final_ema_f3_day_28_32, df_final_ema_f3_day_33_35]
labels = ['26°-30° C', '28°-32° C', '33°-35°C', '26°-30° C', '28°-32° C', '33°-35°C', '26°-30° C', '28°-32° C', '33°-35°C']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['sandybrown', 'mediumpurple', 'lavenderblush', 'chocolate', 'blueviolet', 'palevioletred', 'saddlebrown', 'rebeccapurple', 'mediumvioletred']
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
axs.set_title(f'Average number of extreme days per year in {station} \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =42)
axs.set_ylabel('# of days', size=10, labelpad=10)
os.makedirs(rootOut, exist_ok=True)
fig1.savefig(os.path.join(rootOut, '26_30_vs_28_32_vs_33_35_number_of_extreme_days_boxplot.png'), format='png', dpi=1000)


# '''Create Boxplots number of extreme events'''
fig2, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data2 = [df_final_ema_f1_26_30, df_final_ema_f1_28_32, df_final_ema_f1_33_35, df_final_ema_f2_26_30, df_final_ema_f2_28_32, df_final_ema_f2_33_35, df_final_ema_f3_26_30, df_final_ema_f3_28_32, df_final_ema_f3_33_35]
labels = ['26°-30° C', '28°-32° C', '33°-35°C', '26°-30° C', '28°-32° C', '33°-35°C', '26°-30° C', '28°-32° C', '33°-35°C']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data2, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['sandybrown', 'mediumpurple', 'lavenderblush', 'chocolate', 'blueviolet', 'palevioletred', 'saddlebrown', 'rebeccapurple', 'mediumvioletred']
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
axs.set_title(f'Average number of extreme events per year in {station} \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =12)
axs.set_ylabel('# of events', size=10, labelpad=10)
fig2.savefig(os.path.join(rootOut, '26_30_vs_28_32_number_of_extreme_events_boxplot.png'), format='png', dpi=1000)


'''Create Boxplots length of extreme events'''
fig3, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data3 = [df_final_ema_f1_L_26_30, df_final_ema_f1_L_28_32, df_final_ema_f1_L_33_35, df_final_ema_f2_L_26_30, df_final_ema_f2_L_28_32, df_final_ema_f2_L_33_35, df_final_ema_f3_L_26_30, df_final_ema_f3_L_28_32, df_final_ema_f3_L_33_35]
labels = ['26°-30° C', '28°-32° C', '33°-35°C', '26°-30° C', '28°-32° C', '33°-35°C', '26°-30° C', '28°-32° C', '33°-35°C']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data3, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['sandybrown', 'mediumpurple', 'lavenderblush', 'chocolate', 'blueviolet', 'palevioletred', 'saddlebrown', 'rebeccapurple', 'mediumvioletred']
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
axs.set_title(f'Average length of extreme events per year in {station} \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =4.2)
axs.set_ylabel('# of days', size=10, labelpad=10)
fig3.savefig(os.path.join(rootOut, '26_30_vs_28_32_length_of_extreme_events_boxplot.png'), format='png', dpi=1000)


'''Create Boxplots number of extreme days for very extreme events or for very acclimatized people'''
fig6, axs = plt.subplots(figsize=(8,5), dpi=1000) 
data1 = [df_final_ema_f1_day_33_35, df_final_ema_f2_day_33_35, df_final_ema_f3_day_33_35]
labels = ['RCP2.6', 'RCP4.5', 'RCP8.5']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['lavenderblush', 'palevioletred','mediumvioletred']
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
axs.set_title(f'Average number of extreme days per year in ({station}) \n 1981 - 2099', size=12, pad=10)
axs.set_ylim(bottom=0, top =5)
axs.set_ylabel('# of days', size=10, labelpad=10)
fig6.savefig(os.path.join(rootOut, '33_35_number_of_extreme_days_boxplot.png'), format='png', dpi=1000)

'''Statistical values'''

# Mean values 2.6 NED
print('Mean values 2.6 NED')
print(df_final_ema_f1_day_26_30.mean())
print(df_final_ema_f1_day_28_32.mean())
print(df_final_ema_f1_day_33_35.mean())
# Mean values 4.5 NED
print('Mean values 4.5 NED')
print(df_final_ema_f2_day_26_30.mean())
print(df_final_ema_f2_day_28_32.mean())
print(df_final_ema_f2_day_33_35.mean())
# Mean values 8.5 NED
print('Mean values 8.5 NED')
print(df_final_ema_f3_day_26_30.mean())
print(df_final_ema_f3_day_28_32.mean())
print(df_final_ema_f3_day_33_35.mean())

# Mean values 2.6 NEE
print('Mean values 2.6 NEE')
print(df_final_ema_f1_26_30.mean())
print(df_final_ema_f1_28_32.mean())
print(df_final_ema_f1_33_35.mean())
# Mean values 4.5 NEE
print('Mean values 4.5 NEE')
print(df_final_ema_f2_26_30.mean())
print(df_final_ema_f2_28_32.mean())
print(df_final_ema_f2_33_35.mean())
# Mean values 8.5 NEE
print('Mean values 8.5 NEE')
print(df_final_ema_f3_26_30.mean())
print(df_final_ema_f3_28_32.mean())
print(df_final_ema_f3_33_35.mean())

# Mean values 2.6 LEE
print('Mean values 2.6 LEE')
print(df_final_ema_f1_L_26_30.mean())
print(df_final_ema_f1_L_28_32.mean())
print(df_final_ema_f1_L_33_35.mean())
# Mean values 4.5 LEE
print('Mean values 4.5 LEE')
print(df_final_ema_f2_L_26_30.mean())
print(df_final_ema_f2_L_28_32.mean())
print(df_final_ema_f2_L_33_35.mean())
# Mean values 8.5 LEE
print('Mean values 8.5 LEE')
print(df_final_ema_f3_L_26_30.mean())
print(df_final_ema_f3_L_28_32.mean())
print(df_final_ema_f3_L_33_35.mean())

# Min values 2.6 NED
print('Min values 2.6 NED')
print(df_final_ema_f1_day_26_30.min())
print(df_final_ema_f1_day_28_32.min())
print(df_final_ema_f1_day_33_35.min())
# Min values 4.5 NED
print(df_final_ema_f2_day_26_30.min())
print(df_final_ema_f2_day_28_32.min())
print(df_final_ema_f2_day_33_35.min())
# Min values 8.5 NED
print(df_final_ema_f3_day_26_30.min())
print(df_final_ema_f3_day_28_32.min())
print(df_final_ema_f3_day_33_35.min())

# Min values 2.6 NEE
print('Min values 2.6 NEE')
print(df_final_ema_f1_26_30.min())
print(df_final_ema_f1_28_32.min())
print(df_final_ema_f1_33_35.min())
# Min values 4.5 NEE
print(df_final_ema_f2_26_30.min())
print(df_final_ema_f2_28_32.min())
print(df_final_ema_f2_33_35.min())
# Min values 8.5 NEE
print(df_final_ema_f3_26_30.min())
print(df_final_ema_f3_28_32.min())
print(df_final_ema_f3_33_35.min())

# Min values 2.6 LEE
print('Min values 2.6 LEE')
print(df_final_ema_f1_L_26_30.min())
print(df_final_ema_f1_L_28_32.min())
print(df_final_ema_f1_L_33_35.min())
# Min values 4.5 LEE
print(df_final_ema_f2_L_26_30.min())
print(df_final_ema_f2_L_28_32.min())
print(df_final_ema_f2_L_33_35.min())
# Min values 8.5 LEE
print(df_final_ema_f3_L_26_30.min())
print(df_final_ema_f3_L_28_32.min())
print(df_final_ema_f3_L_33_35.min())

# Max values 2.6 NED
print('Max values 2.6 NED')
print(df_final_ema_f1_day_26_30.max())
print(df_final_ema_f1_day_28_32.max())
print(df_final_ema_f1_day_33_35.max())
# Max values 4.5 NED
print(df_final_ema_f2_day_26_30.max())
print(df_final_ema_f2_day_28_32.max())
print(df_final_ema_f2_day_33_35.max())
# Max values 8.5 NED
print(df_final_ema_f3_day_26_30.max())
print(df_final_ema_f3_day_28_32.max())
print(df_final_ema_f3_day_33_35.max())

# Max values 2.6 NEE
print('Max values 2.6 NEE')
print(df_final_ema_f1_26_30.max())
print(df_final_ema_f1_28_32.max())
print(df_final_ema_f1_33_35.max())
# Max values 4.5 NEE
print(df_final_ema_f2_26_30.max())
print(df_final_ema_f2_28_32.max())
print(df_final_ema_f2_33_35.max())
# Max values 8.5 NEE
print(df_final_ema_f3_26_30.max())
print(df_final_ema_f3_28_32.max())
print(df_final_ema_f3_33_35.max())

# Max values 2.6 LEE
print('Max values 2.6 LEE')
print(df_final_ema_f1_L_26_30.max())
print(df_final_ema_f1_L_28_32.max())
print(df_final_ema_f1_L_33_35.max())
# Max values 4.5 LEE
print(df_final_ema_f2_L_26_30.max())
print(df_final_ema_f2_L_28_32.max())
print(df_final_ema_f2_L_33_35.max())
# Max values 8.5 LEE
print(df_final_ema_f3_L_26_30.max())
print(df_final_ema_f3_L_28_32.max())
print(df_final_ema_f3_L_33_35.max())


# # '''Create time series for RCP2.6'''
# # p_Step = 1
# # alpha_Fade = 0.91
# # fig4, axs = plt.subplots(2,1, figsize=(15,10))
# # plt.style.use('seaborn-pastel')

# # x_axis = np.arange(1981, 2100, step=1)

# # for i in range(0, len(matrix_HotDays26_26_30), p_Step):
# #     a_26_30 = matrix_HotDays26_26_30[i, 1:, 1:2].astype(float)
# #     axs[0].plot(x_axis, a_26_30[:,0], color='chocolate', linestyle = '--', linewidth=0.5, alpha=alpha_Fade)

# # for i in range(0, len(matrix_HotDays26_28_32), p_Step):
# #     a_28_32 = matrix_HotDays26_28_32[i, 1:, 1:2].astype(float)
# #     axs[1].plot(x_axis, a_28_32[:,0], color='blueviolet', linestyle = '--', linewidth=0.5, alpha=alpha_Fade)

# # axs[0].set_title('Number of extreme days per year in {station} for RCP2.6 \n 1981 - 2099', fontsize=25, pad=20)
# # axs[1].set_xlabel('year', fontsize=20, labelpad=15)
# # axs[0].set_ylabel('# of days', fontsize=20, labelpad=15)
# # axs[1].set_ylabel('# of days', fontsize=20, labelpad=15)
# # axs[0].set_ylim(bottom=0, top=110)
# # axs[1].set_ylim(bottom=0, top = 110)
# # axs[0].yaxis.set_ticks(np.arange(0, 110, 20))
# # axs[1].yaxis.set_ticks(np.arange(0, 110, 20))
# # axs[0].xaxis.set_ticks(np.arange(1980, 2101, 10))
# # axs[1].xaxis.set_ticks(np.arange(1980, 2101, 10))
# # axs[0].tick_params(axis='x', labelsize=20)
# # axs[1].tick_params(axis='x', labelsize=20)
# # axs[0].tick_params(axis='y', labelsize=20)
# # axs[1].tick_params(axis='y', labelsize=20)
# # axs[0].annotate('26° - 30° C', fontsize=20, xy=(1985, 100), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# # axs[1].annotate('28° - 32° C', fontsize=20, xy=(1985, 100), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# # fig4.savefig(os.path.join(rootOut, '26_30_vs_28_32_Number_extreme_days_RCP26.png'), format='png', dpi=150)

# # '''Create time series for RCP8.5'''
# # p_Step = 1
# # alpha_Fade = 0.91
# # fig5, axs = plt.subplots(2,1, figsize=(15,10))
# # plt.style.use('seaborn-pastel')

# # x_axis = np.arange(1981, 2100, step=1)

# # for i in range(0, len(matrix_HotDays85_26_30), p_Step):
# #     a_26_30 = matrix_HotDays85_26_30[i, 1:, 1:2].astype(float)
# #     axs[0].plot(x_axis, a_26_30[:,0], color='chocolate', linestyle = '--', linewidth=0.5, alpha=alpha_Fade)

# # for i in range(0, len(matrix_HotDays85_28_32), p_Step):
# #     a_28_32 = matrix_HotDays85_28_32[i, 1:, 1:2].astype(float)
# #     axs[1].plot(x_axis, a_28_32[:,0], color='blueviolet', linestyle = '--', linewidth=0.5, alpha=alpha_Fade)

# # axs[0].set_title(f'Number of extreme days per year in {station} for RCP2.6 \n 1981 - 2099', fontsize=25, pad=20)
# # axs[1].set_xlabel('year', fontsize=20, labelpad=15)
# # axs[0].set_ylabel('# of days', fontsize=20, labelpad=15)
# # axs[1].set_ylabel('# of days', fontsize=20, labelpad=15)
# # axs[0].set_ylim(bottom=0, top=110)
# # axs[1].set_ylim(bottom=0, top = 110)
# # axs[0].yaxis.set_ticks(np.arange(0, 110, 20))
# # axs[1].yaxis.set_ticks(np.arange(0, 110, 20))
# # axs[0].xaxis.set_ticks(np.arange(1980, 2101, 10))
# # axs[1].xaxis.set_ticks(np.arange(1980, 2101, 10))
# # axs[0].tick_params(axis='x', labelsize=20)
# # axs[1].tick_params(axis='x', labelsize=20)
# # axs[0].tick_params(axis='y', labelsize=20)
# # axs[1].tick_params(axis='y', labelsize=20)
# # axs[0].annotate('26° - 30° C', fontsize=20, xy=(1985, 100), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# # axs[1].annotate('28° - 32° C', fontsize=20, xy=(1985, 100), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# # fig5.savefig(os.path.join(rootOut, '26_30_vs_28_32_Number_extreme_days_RCP85.png'), format='png', dpi=150)

plt.show()
