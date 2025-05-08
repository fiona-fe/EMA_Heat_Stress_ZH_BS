'''
script purpose:             visualization of simulated adaptation measures and its impacts on the heat indicators NED, NEE, LEE
outputs:                    norm_vs_Tmp_1_0_vs_Tmp_3_0_vs_Tmp_5_0_NED.png
                            norm_vs_Slr_15_50_vs_Wnd_minus3_vs_Tmp_3_0_NED.png
                            statistics
modifications required:     commented with "MODIFY"

codeauthor: Fabian Weibel (2022)
adapted for publication by: Fiona Federer (2025)
'''

# '''Import libraries'''
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
from ema_workbench import load_results, ema_logging
import tarfile

'''Import the ema workbench outputs'''
from ema_workbench.analysis.pairs_plotting import (pairs_lines, pairs_scatter,
                                                   pairs_density)
ema_logging.log_to_stderr(level=ema_logging.DEFAULT_LEVEL)

###### Data preparation ######

'''Load the data for CH2018: uncomment either the Basel or Zurich code block, change paths (MODIFY)'''
# Paths Basel
# station = 'Basel-Binningen'
# rootOut_norm = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/26_30_bernard'
# rootOut_Slr_15_50 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/Slr_15_50'
# rootOut_Wnd_1 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/Wnd_1'
# rootOut_Tmp_1_0 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/Tmp_1_0'
# rootOut_Tmp_3_0 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/Tmp_3_0'
# rootOut_Tmp_5_0 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/Tmp_5_0'

# rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis/sensitivity_analysis_figures'

# Paths Zurich
station = 'Zurich-Fluntern' 
rootOut_norm = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/26_30_bernard'
rootOut_Slr_15_50 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Slr_15_50'
rootOut_Wnd_1 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Wnd_1'
rootOut_Tmp_1_0 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Tmp_1_0'
rootOut_Tmp_3_0 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Tmp_3_0'
rootOut_Tmp_5_0 = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Tmp_5_0'

rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/sensitivity_analysis_figures'

# EMA output names
out_zip_file_name_norm = '20250301_runs.tar.gz'
out_zip_file_name_Slr_15_50 = '20250325_runs.tar.gz'
out_zip_file_name_Wnd_1 = '20250325_runs.tar.gz'
out_zip_file_name_Tmp_1_0 = '20250325_runs.tar.gz'
out_zip_file_name_Tmp_3_0 = '20250325_runs.tar.gz'
out_zip_file_name_Tmp_5_0 = '20250325_runs.tar.gz'

fh_norm = os.path.join(rootOut_norm, out_zip_file_name_norm)
experiments_norm, outcomes_norm = load_results(fh_norm)
fh_Slr_15_50 = os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50)
experiments_Slr_15_50, outcomes_Slr_15_50 = load_results(fh_Slr_15_50)
fh_Wnd_1 = os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1)
experiments_Wnd_1, outcomes_Wnd_1 = load_results(fh_Wnd_1)
fh_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0)
experiments_Tmp_1_0, outcomes_Tmp_1_0 = load_results(fh_Tmp_1_0)
fh_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0)
experiments_Tmp_3_0, outcomes_Tmp_3_0 = load_results(fh_Tmp_3_0)
fh_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0)
experiments_Tmp_5_0, outcomes_Tmp_5_0 = load_results(fh_Tmp_5_0)

'''Create data frames for input data norm'''
with tarfile.open(os.path.join(rootOut_norm, out_zip_file_name_norm),"r") as zip_ref_norm:
    zip_ref_norm.extractall(os.path.join(rootOut_norm, out_zip_file_name_norm[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_norm, out_zip_file_name_norm[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_norm, out_zip_file_name_norm[:-7]), file),
                os.path.join(os.path.join(rootOut_norm, out_zip_file_name_norm[:-7]), file[:-4] + ".csv"))
# We have 6 types of outputs
outDaily_norm = os.path.join(rootOut_norm, 'Outputs_py')
#M1
outSeasonTippingPoint_norm = os.path.join(rootOut_norm, 'outSeason')
outSeason_Length_norm = os.path.join(rootOut_norm, 'outSeason_length')
outSeason_StartEvents_norm = os.path.join(rootOut_norm, 'outSeason_StartEvents')
outSeason_EndEvents_norm = os.path.join(rootOut_norm, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_norm = os.path.join(rootOut_norm, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_norm = os.path.join(rootOut_norm, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_norm = os.path.join(rootOut_norm, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_norm = os.path.join(rootOut_norm, 'outSeason_Ave_StartEvents')
out_ema_norm = os.path.join(rootOut_norm, out_zip_file_name_norm[:-7])

df4_ema_experiment_norm = pd.read_csv(os.path.join(out_ema_norm, 'experiments.csv'))
#M1
df4_ema_y0_norm = pd.read_csv(os.path.join(out_ema_norm, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_norm.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_norm = pd.read_csv(os.path.join(out_ema_norm, 'GCM_RCM.csv'), header=None)
df4_ema_y1_norm.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_norm = pd.read_csv(os.path.join(out_ema_norm, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_norm.columns = ["Yout2_Ave_Length"]
df4_ema_y3_norm = pd.read_csv(os.path.join(out_ema_norm, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_norm.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_norm = pd.read_csv(os.path.join(out_ema_norm, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_norm.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_norm = pd.read_csv(os.path.join(out_ema_norm, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_norm.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_norm = pd.concat((df4_ema_experiment_norm, df4_ema_y0_norm, df4_ema_y1_norm, df4_ema_y2_norm, df4_ema_y3_norm, df4_ema_y4_norm, df4_ema_y5_norm), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_norm = df_final_ema_norm['Xfactor1'].values
xClimateModel_norm = df_final_ema_norm['xClimateModel'].values
xRCP_norm = df_final_ema_norm['xRCP'].values
x7Wbgtthreshold_norm = df_final_ema_norm['x7Wbgtthreshold'].values
#M1
y0_norm = df_final_ema_norm['Yout0_S_Ave_HotDay'].values
y1_norm = df_final_ema_norm['Yout1_GCM_RCM'].values
#M3
y2_norm = df_final_ema_norm['Yout2_Ave_Length'].values
y3_norm = df_final_ema_norm['Yout3_Ave_StartEvent'].values
y4_norm = df_final_ema_norm['Yout4_Ave_EndEvent'].values
#M2
y5_norm = df_final_ema_norm['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_norm['xRCP']) == 1.0)
filt2 = (round(df_final_ema_norm['xRCP']) == 2.0)
filt3 = (round(df_final_ema_norm['xRCP']) == 3.0)
df_final_ema_f1_day_norm = df_final_ema_norm.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_norm = df_final_ema_norm.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_norm = df_final_ema_norm.loc[filt3, 'Yout0_S_Ave_HotDay']

'''Create data frames for input data Slr_15_50'''
with tarfile.open(os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50),"r") as zip_ref_Slr_15_50:
    zip_ref_Slr_15_50.extractall(os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50[:-7]), file),
                os.path.join(os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50[:-7]), file[:-4] + ".csv"))

# We have 6 types of outputs
outDaily_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'Outputs_py')
#M1
outSeasonTippingPoint_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason')
outSeason_Length_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_length')
outSeason_StartEvents_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_StartEvents')
outSeason_EndEvents_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Slr_15_50 = os.path.join(rootOut_Slr_15_50, 'outSeason_Ave_StartEvents')
out_ema_Slr_15_50 = os.path.join(rootOut_Slr_15_50, out_zip_file_name_Slr_15_50[:-7])

df4_ema_experiment_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'experiments.csv'))
#M1 
df4_ema_y0_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Slr_15_50.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Slr_15_50.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Slr_15_50.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Slr_15_50.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Slr_15_50.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Slr_15_50 = pd.read_csv(os.path.join(out_ema_Slr_15_50, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Slr_15_50.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Slr_15_50 = pd.concat((df4_ema_experiment_Slr_15_50, df4_ema_y0_Slr_15_50, df4_ema_y1_Slr_15_50, df4_ema_y2_Slr_15_50, df4_ema_y3_Slr_15_50, df4_ema_y4_Slr_15_50, df4_ema_y5_Slr_15_50), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Slr_15_50 = df_final_ema_Slr_15_50['Xfactor1'].values
xClimateModel_Slr_15_50 = df_final_ema_Slr_15_50['xClimateModel'].values
xRCP_Slr_15_50 = df_final_ema_Slr_15_50['xRCP'].values
x7Wbgtthreshold_Slr_15_50 = df_final_ema_Slr_15_50['x7Wbgtthreshold'].values
#M1
y0_Slr_15_50 = df_final_ema_Slr_15_50['Yout0_S_Ave_HotDay'].values
y1_Slr_15_50 = df_final_ema_Slr_15_50['Yout1_GCM_RCM'].values
#M3
y2_Slr_15_50 = df_final_ema_Slr_15_50['Yout2_Ave_Length'].values
y3_Slr_15_50 = df_final_ema_Slr_15_50['Yout3_Ave_StartEvent'].values
y4_Slr_15_50 = df_final_ema_Slr_15_50['Yout4_Ave_EndEvent'].values
#M2
y5_Slr_15_50 = df_final_ema_Slr_15_50['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Slr_15_50['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Slr_15_50['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Slr_15_50['xRCP']) == 3.0)
df_final_ema_f1_day_Slr_15_50 = df_final_ema_Slr_15_50.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Slr_15_50 = df_final_ema_Slr_15_50.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Slr_15_50 = df_final_ema_Slr_15_50.loc[filt3, 'Yout0_S_Ave_HotDay']

'''Create data frames for input data Wnd_1'''
with tarfile.open(os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1),"r") as zip_ref_Wnd_1:
    zip_ref_Wnd_1.extractall(os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1[:-7]), file),
                os.path.join(os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1[:-7]), file[:-4] + ".csv"))
            
# We have 6 types of outputs
outDaily_Wnd_1 = os.path.join(rootOut_Wnd_1, 'Outputs_py')
#M1
outSeasonTippingPoint_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason')
outSeason_Length_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_length')
outSeason_StartEvents_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_StartEvents')
outSeason_EndEvents_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Wnd_1 = os.path.join(rootOut_Wnd_1, 'outSeason_Ave_StartEvents')
out_ema_Wnd_1 = os.path.join(rootOut_Wnd_1, out_zip_file_name_Wnd_1[:-7])

df4_ema_experiment_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'experiments.csv'))
#M1 
df4_ema_y0_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Wnd_1.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Wnd_1.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Wnd_1.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Wnd_1.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Wnd_1.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Wnd_1 = pd.read_csv(os.path.join(out_ema_Wnd_1, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Wnd_1.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Wnd_1 = pd.concat((df4_ema_experiment_Wnd_1, df4_ema_y0_Wnd_1, df4_ema_y1_Wnd_1, df4_ema_y2_Wnd_1, df4_ema_y3_Wnd_1, df4_ema_y4_Wnd_1, df4_ema_y5_Wnd_1), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Wnd_1 = df_final_ema_Wnd_1['Xfactor1'].values
xClimateModel_Wnd_1 = df_final_ema_Wnd_1['xClimateModel'].values
xRCP_Wnd_1 = df_final_ema_Wnd_1['xRCP'].values
x7Wbgtthreshold_Wnd_1 = df_final_ema_Wnd_1['x7Wbgtthreshold'].values
#M1
y0_Wnd_1 = df_final_ema_Wnd_1['Yout0_S_Ave_HotDay'].values
y1_Wnd_1 = df_final_ema_Wnd_1['Yout1_GCM_RCM'].values
#M3
y2_Wnd_1 = df_final_ema_Wnd_1['Yout2_Ave_Length'].values
y3_Wnd_1 = df_final_ema_Wnd_1['Yout3_Ave_StartEvent'].values
y4_Wnd_1 = df_final_ema_Wnd_1['Yout4_Ave_EndEvent'].values
#M2
y5_Wnd_1 = df_final_ema_Wnd_1['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Wnd_1['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Wnd_1['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Wnd_1['xRCP']) == 3.0)
df_final_ema_f1_day_Wnd_1 = df_final_ema_Wnd_1.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Wnd_1 = df_final_ema_Wnd_1.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Wnd_1 = df_final_ema_Wnd_1.loc[filt3, 'Yout0_S_Ave_HotDay']

'''Create data frames for input data Tmp_1_0'''
with tarfile.open(os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0),"r") as zip_ref_Tmp_1_0:
    zip_ref_Tmp_1_0.extractall(os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0[:-7]), file),
                os.path.join(os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0[:-7]), file[:-4] + ".csv"))
            
# We have 6 types of outputs
outDaily_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'Outputs_py')
#M1
outSeasonTippingPoint_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason')
outSeason_Length_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_length')
outSeason_StartEvents_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_StartEvents')
outSeason_EndEvents_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, 'outSeason_Ave_StartEvents')
out_ema_Tmp_1_0 = os.path.join(rootOut_Tmp_1_0, out_zip_file_name_Tmp_1_0[:-7])

df4_ema_experiment_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'experiments.csv'))
#M1 
df4_ema_y0_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Tmp_1_0.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Tmp_1_0.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Tmp_1_0.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Tmp_1_0.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Tmp_1_0.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Tmp_1_0 = pd.read_csv(os.path.join(out_ema_Tmp_1_0, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Tmp_1_0.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Tmp_1_0 = pd.concat((df4_ema_experiment_Tmp_1_0, df4_ema_y0_Tmp_1_0, df4_ema_y1_Tmp_1_0, df4_ema_y2_Tmp_1_0, df4_ema_y3_Tmp_1_0, df4_ema_y4_Tmp_1_0, df4_ema_y5_Tmp_1_0), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Tmp_1_0 = df_final_ema_Tmp_1_0['Xfactor1'].values
xClimateModel_Tmp_1_0 = df_final_ema_Tmp_1_0['xClimateModel'].values
xRCP_Tmp_1_0 = df_final_ema_Tmp_1_0['xRCP'].values
x7Wbgtthreshold_Tmp_1_0 = df_final_ema_Tmp_1_0['x7Wbgtthreshold'].values
#M1
y0_Tmp_1_0 = df_final_ema_Tmp_1_0['Yout0_S_Ave_HotDay'].values
y1_Tmp_1_0 = df_final_ema_Tmp_1_0['Yout1_GCM_RCM'].values
#M3
y2_Tmp_1_0 = df_final_ema_Tmp_1_0['Yout2_Ave_Length'].values
y3_Tmp_1_0 = df_final_ema_Tmp_1_0['Yout3_Ave_StartEvent'].values
y4_Tmp_1_0 = df_final_ema_Tmp_1_0['Yout4_Ave_EndEvent'].values
#M2
y5_Tmp_1_0 = df_final_ema_Tmp_1_0['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Tmp_1_0['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Tmp_1_0['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Tmp_1_0['xRCP']) == 3.0)
df_final_ema_f1_day_Tmp_1_0 = df_final_ema_Tmp_1_0.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Tmp_1_0 = df_final_ema_Tmp_1_0.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Tmp_1_0 = df_final_ema_Tmp_1_0.loc[filt3, 'Yout0_S_Ave_HotDay']

'''Create data frames for input data Tmp_3_0'''
with tarfile.open(os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0),"r") as zip_ref_Tmp_3_0:
    zip_ref_Tmp_3_0.extractall(os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0[:-7]), file),
                os.path.join(os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0[:-7]), file[:-4] + ".csv"))
            
# We have 6 types of outputs
outDaily_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'Outputs_py')
#M1
outSeasonTippingPoint_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason')
outSeason_Length_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_length')
outSeason_StartEvents_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_StartEvents')
outSeason_EndEvents_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, 'outSeason_Ave_StartEvents')
out_ema_Tmp_3_0 = os.path.join(rootOut_Tmp_3_0, out_zip_file_name_Tmp_3_0[:-7])

df4_ema_experiment_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'experiments.csv'))
#M1 
df4_ema_y0_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Tmp_3_0.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Tmp_3_0.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Tmp_3_0.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Tmp_3_0.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Tmp_3_0.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Tmp_3_0 = pd.read_csv(os.path.join(out_ema_Tmp_3_0, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Tmp_3_0.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Tmp_3_0 = pd.concat((df4_ema_experiment_Tmp_3_0, df4_ema_y0_Tmp_3_0, df4_ema_y1_Tmp_3_0, df4_ema_y2_Tmp_3_0, df4_ema_y3_Tmp_3_0, df4_ema_y4_Tmp_3_0, df4_ema_y5_Tmp_3_0), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Tmp_3_0 = df_final_ema_Tmp_3_0['Xfactor1'].values
xClimateModel_Tmp_3_0 = df_final_ema_Tmp_3_0['xClimateModel'].values
xRCP_Tmp_3_0 = df_final_ema_Tmp_3_0['xRCP'].values
x7Wbgtthreshold_Tmp_3_0 = df_final_ema_Tmp_3_0['x7Wbgtthreshold'].values
#M1
y0_Tmp_3_0 = df_final_ema_Tmp_3_0['Yout0_S_Ave_HotDay'].values
y1_Tmp_3_0 = df_final_ema_Tmp_3_0['Yout1_GCM_RCM'].values
#M3
y2_Tmp_3_0 = df_final_ema_Tmp_3_0['Yout2_Ave_Length'].values
y3_Tmp_3_0 = df_final_ema_Tmp_3_0['Yout3_Ave_StartEvent'].values
y4_Tmp_3_0 = df_final_ema_Tmp_3_0['Yout4_Ave_EndEvent'].values
#M2
y5_Tmp_3_0 = df_final_ema_Tmp_3_0['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Tmp_3_0['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Tmp_3_0['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Tmp_3_0['xRCP']) == 3.0)
df_final_ema_f1_day_Tmp_3_0 = df_final_ema_Tmp_3_0.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Tmp_3_0 = df_final_ema_Tmp_3_0.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Tmp_3_0 = df_final_ema_Tmp_3_0.loc[filt3, 'Yout0_S_Ave_HotDay']

'''Create data frames for input data Tmp_5_0'''
with tarfile.open(os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0),"r") as zip_ref_Tmp_5_0:
    zip_ref_Tmp_5_0.extractall(os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0[:-7]))
    # rename .cls files to .csv after extraction (reason for .cls = MacOS)
    for file in os.listdir(os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0[:-7])):
        if file.endswith(".cls"):
            os.rename(
                os.path.join(os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0[:-7]), file),
                os.path.join(os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0[:-7]), file[:-4] + ".csv"))
            
# We have 6 types of outputs
outDaily_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'Outputs_py')
#M1
outSeasonTippingPoint_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason')
outSeason_Length_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_length')
outSeason_StartEvents_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_StartEvents')
outSeason_EndEvents_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_EndEvents')
#M2
outSeason_NumEvents_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_NumEvents')
outSeason_Ave_EndEvents_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_Ave_EndEvents')
#M3
outSeason_Ave_LenEvents_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_Ave_LenEvents')
outSeason_Ave_StartEvents_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, 'outSeason_Ave_StartEvents')
out_ema_Tmp_5_0 = os.path.join(rootOut_Tmp_5_0, out_zip_file_name_Tmp_5_0[:-7])

df4_ema_experiment_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'experiments.csv'))
#M1 
df4_ema_y0_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'S_Ave_HotDay.csv'), header=None)
df4_ema_y0_Tmp_5_0.columns = ["Yout0_S_Ave_HotDay"]
df4_ema_y1_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'GCM_RCM.csv'), header=None)
df4_ema_y1_Tmp_5_0.columns = ["Yout1_GCM_RCM"]
#M3
df4_ema_y2_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'S_Ave_Length.csv'), header=None)
df4_ema_y2_Tmp_5_0.columns = ["Yout2_Ave_Length"]
df4_ema_y3_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'S_Ave_StartEvent.csv'), header=None)
df4_ema_y3_Tmp_5_0.columns = ["Yout3_Ave_StartEvent"]
df4_ema_y4_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'S_Ave_EndEvent.csv'), header=None)
df4_ema_y4_Tmp_5_0.columns = ["Yout4_Ave_EndEvent"]
#M2
df4_ema_y5_Tmp_5_0 = pd.read_csv(os.path.join(out_ema_Tmp_5_0, 'S_Ave_numEvent.csv'), header=None)
df4_ema_y5_Tmp_5_0.columns = ["Yout5_S_Ave_numEvent"]

df_final_ema_Tmp_5_0 = pd.concat((df4_ema_experiment_Tmp_5_0, df4_ema_y0_Tmp_5_0, df4_ema_y1_Tmp_5_0, df4_ema_y2_Tmp_5_0, df4_ema_y3_Tmp_5_0, df4_ema_y4_Tmp_5_0, df4_ema_y5_Tmp_5_0), axis = 1)

# Taking the values and storing in separate variables (CH2018)
xRandomness_Tmp_5_0 = df_final_ema_Tmp_5_0['Xfactor1'].values
xClimateModel_Tmp_5_0 = df_final_ema_Tmp_5_0['xClimateModel'].values
xRCP_Tmp_5_0 = df_final_ema_Tmp_5_0['xRCP'].values
x7Wbgtthreshold_Tmp_5_0 = df_final_ema_Tmp_5_0['x7Wbgtthreshold'].values
#M1
y0_Tmp_5_0 = df_final_ema_Tmp_5_0['Yout0_S_Ave_HotDay'].values
y1_Tmp_5_0 = df_final_ema_Tmp_5_0['Yout1_GCM_RCM'].values
#M3
y2_Tmp_5_0 = df_final_ema_Tmp_5_0['Yout2_Ave_Length'].values
y3_Tmp_5_0 = df_final_ema_Tmp_5_0['Yout3_Ave_StartEvent'].values
y4_Tmp_5_0 = df_final_ema_Tmp_5_0['Yout4_Ave_EndEvent'].values
#M2
y5_Tmp_5_0 = df_final_ema_Tmp_5_0['Yout5_S_Ave_numEvent'].values

filt1 = (round(df_final_ema_Tmp_5_0['xRCP']) == 1.0)
filt2 = (round(df_final_ema_Tmp_5_0['xRCP']) == 2.0)
filt3 = (round(df_final_ema_Tmp_5_0['xRCP']) == 3.0)
df_final_ema_f1_day_Tmp_5_0 = df_final_ema_Tmp_5_0.loc[filt1, 'Yout0_S_Ave_HotDay']
df_final_ema_f2_day_Tmp_5_0 = df_final_ema_Tmp_5_0.loc[filt2, 'Yout0_S_Ave_HotDay']
df_final_ema_f3_day_Tmp_5_0 = df_final_ema_Tmp_5_0.loc[filt3, 'Yout0_S_Ave_HotDay']

'''Making matrices for norm'''
outFolder_norm = rootOut_norm
outSeasonFolder_norm = os.path.join(outFolder_norm, 'outSeason')
all_Files_norm = []
for filename in os.walk(outSeasonFolder_norm):
    all_Files_norm = filename[2]
    
totalFiles_loc_norm = []
for i in range(len(all_Files_norm)):
    totalFiles_loc_norm.append(os.path.join(outSeasonFolder_norm, all_Files_norm[i]))

emptyMatrixHotDays_norm = [np.empty([1,1]) for _ in range(len(totalFiles_loc_norm))]
    
for i in range(0, len(totalFiles_loc_norm), 1):
    b_norm = []
    with open(totalFiles_loc_norm[i], 'r') as file:
        outputReaderlines_norm = file.readlines()
        for j in range (len(outputReaderlines_norm)):
            b_norm.append(outputReaderlines_norm[j].replace('\n','').split(','))
    emptyMatrixHotDays_norm[i] = np.array(b_norm)

matrix_HotDays_norm = np.stack(emptyMatrixHotDays_norm[:len(totalFiles_loc_norm)], axis=0) 

matrix_HotDays26_norm= np.stack([arr for arr in emptyMatrixHotDays_norm[:len(totalFiles_loc_norm)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_norm= np.stack([arr for arr in emptyMatrixHotDays_norm[:len(totalFiles_loc_norm)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_norm= np.stack([arr for arr in emptyMatrixHotDays_norm[:len(totalFiles_loc_norm)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for Slr_15_50'''
outFolder_Slr_15_50 = rootOut_Slr_15_50
outSeasonFolder_Slr_15_50 = os.path.join(outFolder_Slr_15_50, 'outSeason')
all_Files_Slr_15_50 = []
for filename in os.walk(outSeasonFolder_Slr_15_50):
    all_Files_Slr_15_50 = filename[2]
    
totalFiles_loc_Slr_15_50 = []
for i in range(len(all_Files_Slr_15_50)):
    totalFiles_loc_Slr_15_50.append(os.path.join(outSeasonFolder_Slr_15_50, all_Files_Slr_15_50[i]))

emptyMatrixHotDays_Slr_15_50 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_Slr_15_50))]
    
for i in range(0, len(totalFiles_loc_Slr_15_50), 1):
    b_Slr_15_50 = []
    with open(totalFiles_loc_Slr_15_50[i], 'r') as file:
        outputReaderlines_Slr_15_50 = file.readlines()
        for j in range (len(outputReaderlines_Slr_15_50)):
            b_Slr_15_50.append(outputReaderlines_Slr_15_50[j].replace('\n','').split(','))
    emptyMatrixHotDays_Slr_15_50[i] = np.array(b_Slr_15_50)

matrix_HotDays_Slr_15_50 = np.stack(emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)], axis=0) 

matrix_HotDays26_Slr_15_50= np.stack([arr for arr in emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_Slr_15_50= np.stack([arr for arr in emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_Slr_15_50= np.stack([arr for arr in emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for Wnd_1'''
outFolder_Wnd_1 = rootOut_Wnd_1
outSeasonFolder_Wnd_1 = os.path.join(outFolder_Wnd_1, 'outSeason')
all_Files_Wnd_1 = []
for filename in os.walk(outSeasonFolder_Wnd_1):
    all_Files_Wnd_1 = filename[2]
    
totalFiles_loc_Wnd_1 = []
for i in range(len(all_Files_Wnd_1)):
    totalFiles_loc_Wnd_1.append(os.path.join(outSeasonFolder_Wnd_1, all_Files_Wnd_1[i]))

emptyMatrixHotDays_Wnd_1 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_Wnd_1))]
    
for i in range(0, len(totalFiles_loc_Wnd_1), 1):
    b_Wnd_1 = []
    with open(totalFiles_loc_Wnd_1[i], 'r') as file:
        outputReaderlines_Wnd_1 = file.readlines()
        for j in range (len(outputReaderlines_Wnd_1)):
            b_Wnd_1.append(outputReaderlines_Wnd_1[j].replace('\n','').split(','))
    emptyMatrixHotDays_Wnd_1[i] = np.array(b_Wnd_1)

matrix_HotDays_Wnd_1 = np.stack(emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)], axis=0) 

matrix_HotDays26_Wnd_1= np.stack([arr for arr in emptyMatrixHotDays_Wnd_1[:len(totalFiles_loc_Wnd_1)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_Wnd_1= np.stack([arr for arr in emptyMatrixHotDays_Wnd_1[:len(totalFiles_loc_Wnd_1)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_Wnd_1= np.stack([arr for arr in emptyMatrixHotDays_Wnd_1[:len(totalFiles_loc_Wnd_1)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for Tmp_1_0'''
outFolder_Tmp_1_0 = rootOut_Tmp_1_0
outSeasonFolder_Tmp_1_0 = os.path.join(outFolder_Tmp_1_0, 'outSeason')
all_Files_Tmp_1_0 = []
for filename in os.walk(outSeasonFolder_Tmp_1_0):
    all_Files_Tmp_1_0 = filename[2]
    
totalFiles_loc_Tmp_1_0 = []
for i in range(len(all_Files_Tmp_1_0)):
    totalFiles_loc_Tmp_1_0.append(os.path.join(outSeasonFolder_Tmp_1_0, all_Files_Tmp_1_0[i]))

emptyMatrixHotDays_Tmp_1_0 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_Tmp_1_0))]
    
for i in range(0, len(totalFiles_loc_Tmp_1_0), 1):
    b_Tmp_1_0 = []
    with open(totalFiles_loc_Tmp_1_0[i], 'r') as file:
        outputReaderlines_Tmp_1_0 = file.readlines()
        for j in range (len(outputReaderlines_Tmp_1_0)):
            b_Tmp_1_0.append(outputReaderlines_Tmp_1_0[j].replace('\n','').split(','))
    emptyMatrixHotDays_Tmp_1_0[i] = np.array(b_Tmp_1_0)

matrix_HotDays_Tmp_1_0 = np.stack(emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)], axis=0) 

matrix_HotDays26_Tmp_1_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_1_0[:len(totalFiles_loc_Tmp_1_0)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_Tmp_1_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_1_0[:len(totalFiles_loc_Tmp_1_0)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_Tmp_1_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_1_0[:len(totalFiles_loc_Tmp_1_0)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for Tmp_3_0'''
outFolder_Tmp_3_0 = rootOut_Tmp_3_0
outSeasonFolder_Tmp_3_0 = os.path.join(outFolder_Tmp_3_0, 'outSeason')
all_Files_Tmp_3_0 = []
for filename in os.walk(outSeasonFolder_Tmp_3_0):
    all_Files_Tmp_3_0 = filename[2]
    
totalFiles_loc_Tmp_3_0 = []
for i in range(len(all_Files_Tmp_3_0)):
    totalFiles_loc_Tmp_3_0.append(os.path.join(outSeasonFolder_Tmp_3_0, all_Files_Tmp_3_0[i]))

emptyMatrixHotDays_Tmp_3_0 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_Tmp_3_0))]
    
for i in range(0, len(totalFiles_loc_Tmp_3_0), 1):
    b_Tmp_3_0 = []
    with open(totalFiles_loc_Tmp_3_0[i], 'r') as file:
        outputReaderlines_Tmp_3_0 = file.readlines()
        for j in range (len(outputReaderlines_Tmp_3_0)):
            b_Tmp_3_0.append(outputReaderlines_Tmp_3_0[j].replace('\n','').split(','))
    emptyMatrixHotDays_Tmp_3_0[i] = np.array(b_Tmp_3_0)

matrix_HotDays_Tmp_3_0 = np.stack(emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)], axis=0) 

matrix_HotDays26_Tmp_3_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_3_0[:len(totalFiles_loc_Tmp_3_0)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_Tmp_3_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_3_0[:len(totalFiles_loc_Tmp_3_0)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_Tmp_3_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_3_0[:len(totalFiles_loc_Tmp_3_0)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)


'''Making matrices for Tmp_5_0'''
outFolder_Tmp_5_0 = rootOut_Tmp_5_0
outSeasonFolder_Tmp_5_0 = os.path.join(outFolder_Tmp_5_0, 'outSeason')
all_Files_Tmp_5_0 = []
for filename in os.walk(outSeasonFolder_Tmp_5_0):
    all_Files_Tmp_5_0 = filename[2]
    
totalFiles_loc_Tmp_5_0 = []
for i in range(len(all_Files_Tmp_5_0)):
    totalFiles_loc_Tmp_5_0.append(os.path.join(outSeasonFolder_Tmp_5_0, all_Files_Tmp_5_0[i]))

emptyMatrixHotDays_Tmp_5_0 = [np.empty([1,1]) for _ in range(len(totalFiles_loc_Tmp_5_0))]
    
for i in range(0, len(totalFiles_loc_Tmp_5_0), 1):
    b_Tmp_5_0 = []
    with open(totalFiles_loc_Tmp_5_0[i], 'r') as file:
        outputReaderlines_Tmp_5_0 = file.readlines()
        for j in range (len(outputReaderlines_Tmp_5_0)):
            b_Tmp_5_0.append(outputReaderlines_Tmp_5_0[j].replace('\n','').split(','))
    emptyMatrixHotDays_Tmp_5_0[i] = np.array(b_Tmp_5_0)

matrix_HotDays_Tmp_5_0 = np.stack(emptyMatrixHotDays_Slr_15_50[:len(totalFiles_loc_Slr_15_50)], axis=0) 

matrix_HotDays26_Tmp_5_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_5_0[:len(totalFiles_loc_Tmp_5_0)] if 'RCP26' in arr[0, 1] or '_26_' in arr[0, 1]], axis=0)
matrix_HotDays45_Tmp_5_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_5_0[:len(totalFiles_loc_Tmp_5_0)] if 'RCP45' in arr[0, 1] or '_45_' in arr[0, 1]], axis=0)
matrix_HotDays85_Tmp_5_0= np.stack([arr for arr in emptyMatrixHotDays_Tmp_5_0[:len(totalFiles_loc_Tmp_5_0)] if 'RCP85' in arr[0, 1] or '_85_' in arr[0, 1]], axis=0)

###### Visualization ######

'''Create Boxplots Norm, Tmp_1_0, Tmp_3_0, Tmp_5_0'''
fig1, axs = plt.subplots(figsize=(15,5), dpi=1000)

data1 = [df_final_ema_f1_day_norm, df_final_ema_f1_day_Tmp_1_0, df_final_ema_f1_day_Tmp_3_0, df_final_ema_f1_day_Tmp_5_0, df_final_ema_f2_day_norm, df_final_ema_f2_day_Tmp_1_0, df_final_ema_f2_day_Tmp_3_0, df_final_ema_f2_day_Tmp_5_0, df_final_ema_f3_day_norm, df_final_ema_f3_day_Tmp_1_0, df_final_ema_f3_day_Tmp_3_0, df_final_ema_f3_day_Tmp_5_0]
labels = ['Baseline', 'Tmp: -1°C', 'Tmp: -3°C', 'Tmp: -5°C', 'Baseline', 'Tmp: -1°C', 'Tmp: -3°C', 'Tmp: -5°C', 'Baseline', 'Tmp: -1°C', 'Tmp: -3°C', 'Tmp: -5°C',]
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['grey', 'cornflowerblue', 'blue', 'midnightblue', 'grey', 'cornflowerblue', 'blue', 'midnightblue', 'grey', 'cornflowerblue', 'blue', 'midnightblue']
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
axs.set_ylim(bottom=0, top =25)
axs.set_xlabel('climate scenario', size=10, labelpad=10)
axs.set_ylabel('# of days', size=10, labelpad=10)
os.makedirs(rootOut, exist_ok=True)
fig1.savefig(os.path.join(rootOut, 'norm_vs_Tmp_1_0_vs_Tmp_3_0_vs_Tmp_5_0_ave_NED.png'), format='png', dpi=1000)



'''Create Boxplots Norm, Slr_15_50, Wnd_1, Tmp_1_0'''
fig3, axs = plt.subplots(figsize=(15,5))

data1 = [df_final_ema_f1_day_norm, df_final_ema_f1_day_Slr_15_50, df_final_ema_f1_day_Wnd_1, df_final_ema_f1_day_Tmp_1_0, df_final_ema_f2_day_norm, df_final_ema_f2_day_Slr_15_50, df_final_ema_f2_day_Wnd_1, df_final_ema_f2_day_Tmp_1_0, df_final_ema_f3_day_norm, df_final_ema_f3_day_Slr_15_50, df_final_ema_f3_day_Wnd_1, df_final_ema_f3_day_Tmp_1_0]
labels = ['Baseline', 'Glr: 15_50', 'Wnd: +1', 'Tmp: -3°C', 'Baseline', 'Glr: 15_50', 'Wnd: +1', 'Tmp: -3°C', 'Baseline', 'Glr: 15_50', 'Wnd +1', 'Tmp: -3°C']
flierprops = {'marker': 'o', 'markersize': 3, 'markeredgewidth': 0.5}
box = axs.boxplot(data1, vert=True, patch_artist=True, labels=labels, flierprops=flierprops)
plt.xticks(fontsize=8, rotation=45)
plt.yticks(fontsize=8)
colors = ['grey', 'orange', 'darkgreen', 'cornflowerblue', 'grey', 'orange', 'darkgreen', 'cornflowerblue', 'grey', 'orange', 'darkgreen', 'cornflowerblue']
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
axs.set_ylim(bottom=0, top =25)
axs.set_xlabel('climate scenario', size=10, labelpad=10)
axs.set_ylabel('# of days', size=10, labelpad=10)
if not os.path.exists(rootOut):
    os.makedirs(rootOut)
fig3.savefig(os.path.join(rootOut, 'norm_vs_Slr_15_50_vs_Wnd_minus3_vs_Tmp_3_0_ave_NED.png'), format='png', dpi=1000)


'''Create time series for RCP4.5''' # löschen: not sure what this part is for
# p_Step = 1
# alpha_Fade = 0.91
# fig2, axs = plt.subplots(2,1, figsize=(15,10))
# plt.style.use('seaborn-pastel')

# x_axis = np.arange(1981, 2100, step=1)

# for i in range(0, len(matrix_HotDays45_Slr_15_50), p_Step):
#     a_Slr_15_50 = matrix_HotDays45_Slr_15_50[i, 1:, 1:2].astype(float)
#     axs[0].plot(x_axis, a_Slr_15_50[:,0], color='orange', linestyle = '--', alpha=alpha_Fade)

# for i in range(0, len(matrix_HotDays45_Tmp_3_0), p_Step):
#     a_Tmp_3_0 = matrix_HotDays45_Tmp_3_0[i, 1:, 1:2].astype(float)
#     axs[1].plot(x_axis, a_Tmp_3_0[:,0], color='orange', linestyle = '--', alpha=alpha_Fade)

# axs[0].set_title('Number of extreme days per year in Basel-Binningen for RCP4.5 \n 1981 - 2099', fontsize=25, pad=20)
# axs[0].set_xlabel('year', fontsize=20, labelpad=20)
# axs[0].set_ylabel('# of days', fontsize=20, labelpad=20)
# axs[1].set_ylabel('# of days', fontsize=20, labelpad=20)
# axs[0].set_ylim(bottom=0, top=110)
# axs[1].set_ylim(bottom=0, top = 110)
# axs[0].yaxis.set_ticks(np.arange(0, 110, 20))
# axs[1].yaxis.set_ticks(np.arange(0, 110, 20))
# axs[0].xaxis.set_ticks(np.arange(1980, 2101, 10))
# axs[1].xaxis.set_ticks(np.arange(1980, 2101, 10))
# axs[0].tick_params(axis='x', labelsize=20)
# axs[1].tick_params(axis='x', labelsize=20)
# axs[0].tick_params(axis='y', labelsize=20)
# axs[1].tick_params(axis='y', labelsize=20)
# axs[0].annotate('26° - 30° C', fontsize=20, xy=(1985, 100), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# axs[1].annotate('28° - 32° C', fontsize=20, xy=(1985, 100), xytext=(0, 0), textcoords='offset points', ha='center', va='center')
# fig2.savefig(os.path.join(rootOut, 'Number_extreme_days_26_45_85.png'), format='png', dpi=150)


'''Statistical values'''

# Min values for each scenario 
print("min values 2.6")
print(df_final_ema_f1_day_norm.min())
print(df_final_ema_f1_day_Slr_15_50.min())
print(df_final_ema_f1_day_Wnd_1.min())
print(df_final_ema_f1_day_Tmp_1_0.min())
print(df_final_ema_f1_day_Tmp_3_0.min())
print(df_final_ema_f1_day_Tmp_5_0.min())
print("min values 4.5")
print(df_final_ema_f2_day_norm.min())
print(df_final_ema_f2_day_Slr_15_50.min())
print(df_final_ema_f2_day_Wnd_1.min())
print(df_final_ema_f2_day_Tmp_1_0.min())
print(df_final_ema_f2_day_Tmp_3_0.min())
print(df_final_ema_f2_day_Tmp_5_0.min())
print("min values 8.5")
print(df_final_ema_f3_day_norm.min())
print(df_final_ema_f3_day_Slr_15_50.min())
print(df_final_ema_f3_day_Wnd_1.min())
print(df_final_ema_f3_day_Tmp_1_0.min())
print(df_final_ema_f3_day_Tmp_3_0.min())
print(df_final_ema_f3_day_Tmp_5_0.min())

# Max values for each scenario
print("max values 2.6")
print(df_final_ema_f1_day_norm.max())
print(df_final_ema_f1_day_Slr_15_50.max())
print(df_final_ema_f1_day_Wnd_1.max())
print(df_final_ema_f1_day_Tmp_1_0.max())
print(df_final_ema_f1_day_Tmp_3_0.max())
print(df_final_ema_f1_day_Tmp_5_0.max())
print("max values 4.5")
print(df_final_ema_f2_day_norm.max())
print(df_final_ema_f2_day_Slr_15_50.max())
print(df_final_ema_f2_day_Wnd_1.max())
print(df_final_ema_f2_day_Tmp_1_0.max())
print(df_final_ema_f2_day_Tmp_3_0.max())
print(df_final_ema_f2_day_Tmp_5_0.max())
print("max values 8.5")
print(df_final_ema_f3_day_norm.max())
print(df_final_ema_f3_day_Slr_15_50.max())
print(df_final_ema_f3_day_Wnd_1.max())
print(df_final_ema_f3_day_Tmp_1_0.max())
print(df_final_ema_f3_day_Tmp_3_0.max())
print(df_final_ema_f3_day_Tmp_5_0.max())

# Mean values for each scenario
print("mean values 2.6")
print(df_final_ema_f1_day_norm.mean())
print(df_final_ema_f1_day_Slr_15_50.mean())
print(df_final_ema_f1_day_Wnd_1.mean())
print(df_final_ema_f1_day_Tmp_1_0.mean())
print(df_final_ema_f1_day_Tmp_3_0.mean())
print(df_final_ema_f1_day_Tmp_5_0.mean())
print("mean values 4.5")
print(df_final_ema_f2_day_norm.mean())
print(df_final_ema_f2_day_Slr_15_50.mean())
print(df_final_ema_f2_day_Wnd_1.mean())
print(df_final_ema_f2_day_Tmp_1_0.mean())
print(df_final_ema_f2_day_Tmp_3_0.mean())
print(df_final_ema_f2_day_Tmp_5_0.mean())
print("mean values 8.5")
print(df_final_ema_f3_day_norm.mean())
print(df_final_ema_f3_day_Slr_15_50.mean())
print(df_final_ema_f3_day_Wnd_1.mean())
print(df_final_ema_f3_day_Tmp_1_0.mean())
print(df_final_ema_f3_day_Tmp_3_0.mean())
print(df_final_ema_f3_day_Tmp_5_0.mean())

plt.show()
