'''
script purpose:             inputdata preparation for adaptation measures sensitivity analysis (transformation)
outputs:                    modifies the values in the 'input' folder (Climate_ref, Climate_ref_no_randomness_0. Climate_ref_randomness_1/2)
best practice:              create a folder 'setup3_sensitivity_analysis' > subfolders for the modelled adaptation measures
                            > copy 'input' into the folders
                            > careful: in each run the csv get overwritten!                                     
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
import seaborn as sns
import time
import os


'''Load the data for CH2018: uncomment either the Basel or Zurich code block, change paths (MODIFY)'''
# base_path = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_sensitivity_analysis'
base_path = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis'

'''Read csv data and store it as pandas dataframe'''
# Define foldernames and corresponding wildcard (*) (csv-name endings) for different modelled adaptation measures
options = {
    "Wnd_1": "w",
    "Tmp_5_0": "t",
    "Tmp_3_0": "t",
    "Tmp_1_0": "t",
    "Corrfactor_Tmp_1_3_5": "t",
    "Slr_15_50": "s",
    "Corrfactor_Wnd_3": "w",
}

# Define datasets to be modified
subfolders = [
    "Climate_ref",
    "Climate_ref_no_randomness_0",
    "Climate_ref_randomness_1",
    "Climate_ref_randomness_2"
]

'''Functions'''
# Modification function
def process_csv(file_path, option):
    df = pd.read_csv(file_path, sep=',', header=0)
    print(df)

    if option == "w":  # Wind-related modifications (Wnd_1, Corrfactor_Wnd_3)
        if "Wnd_1" in file_path:
            df.loc[1:,:] += 1
        elif "Corrfactor_Wnd_3" in file_path:
              for col in df.columns:
                df.loc[1:, col] = df.loc[1:, col].apply(lambda x: x - 3 if x >= 3 else 0)

    elif option == "t":  # Temperature-related modifications (Tmp_* and station_vs_center_comparison)
        # Add month column
        df['Date'] = pd.NaT  # Initialize with NaT (Not a Time value)
        df.loc[1:, 'Date'] = pd.date_range(start="1981-01-01", periods=len(df) - 1, freq="D")
        df['year'] = df['Date'].dt.year
        df['month'] = df['Date'].dt.month
        df['day'] = df['Date'].dt.day

        # Transform temperatures (simulates adaptation measures)
        # Apply transformation to all cells of the rows where "month" is between April and September
        mask = (df['month'] >= 4) & (df['month'] <= 9)
        columns_to_transform = df.columns.difference(['Date', 'year', 'month', 'day'])
        
        if "Corrfactor_Tmp_1_3_5" in file_path:
            df.loc[1:, columns_to_transform] = df.loc[1:, columns_to_transform].applymap(lambda x: 
                x if x <= 20 else
                x + 1 if 20 < x <= 25 else
                x + 3 if 25 < x <= 30 else
                x + 5)
            df.drop(columns=['Date', 'year', 'month', 'day'], inplace=True)


        if "Tmp_5_0" in file_path:
            df.loc[mask, columns_to_transform] -= 5
            df.drop(columns=['Date', 'year', 'month', 'day'], inplace=True)

        elif "Tmp_3_0" in file_path:
            df.loc[mask, columns_to_transform] -= 3
            df.drop(columns=['Date', 'year', 'month', 'day'], inplace=True)

        elif "Tmp_1_0" in file_path:
            df.loc[mask, columns_to_transform] -= 1
            df.drop(columns=['Date', 'year', 'month', 'day'], inplace=True)


    elif option == "s": # Radiation-related modifications (Slr_15_50)
        # Add month column
        df['Date'] = pd.NaT  # Initialize with NaT (Not a Time value)
        df.loc[1:, 'Date'] = pd.date_range(start="1981-01-01", periods=len(df) - 1, freq="D")
        df['year'] = df['Date'].dt.year
        df['month'] = df['Date'].dt.month
        df['day'] = df['Date'].dt.day

        # Transform radiation (simulates adaptation measures)  
        columns_to_transform = df.columns.difference(['Date', 'year', 'month', 'day'])

        if "Slr_15_50" in file_path:
            # Apply transformations based on the month
            df.loc[df['month'].between(4, 9), columns_to_transform] *= 0.15  # April to September
            df.loc[df['month'].isin([1, 2, 3, 10, 11, 12]), columns_to_transform] *= 0.5  # October to March
            df.drop(columns=['Date', 'year', 'month', 'day'], inplace=True)



    df.to_csv(file_path, sep=',', index=False)
    print("end")

# Main function to modify all datasets and all modelled adapation measures
def process_all():
    for option, wildcard in options.items():
        for subfolder in subfolders:
            folder_path = os.path.join(base_path, option, subfolder)
            if os.path.exists(folder_path):
                file_pattern = f"47-5410000_7-5836000{wildcard}.csv"
                for file_name in os.listdir(folder_path):
                    if file_name.endswith(file_pattern):
                        file_path = os.path.join(folder_path, file_name)
                        process_csv(file_path, wildcard)
                        print(f"Processed: {file_path}")
            else:
                print(f"Skipping missing folder: {folder_path}")


process_all()

