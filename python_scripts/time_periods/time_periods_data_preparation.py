'''
script purpose:             csv-preparation for time_periods_....py
outputs:                    e.g. time_periods_NED_RCP26_{city}_2.csv
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
from numpy.polynomial.polynomial import polyfit
import os
import csv as csv


'''Load the data for CH2018: uncomment either the Basel or Zurich code block, change paths (MODIFY)'''
# Define the paths for CH2018 Basel
# city = 'Basel' # name used for directories
# rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Basel/setup3_bernard_136_30/26_30_bernard' 

# Define the paths for CH2018 Basel
city = 'Zurich'
rootOut = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_bernard_136_30/26_30_bernard' 

# List of CSV file paths for all the matrices 
files = [
    os.path.join(rootOut, f'time_periods_NED_RCP26_{city}.csv'), os.path.join(rootOut, f'time_periods_NED_RCP45_{city}.csv'), os.path.join(rootOut, f'time_periods_NED_RCP85_{city}.csv'),
    os.path.join(rootOut, f'time_periods_NEE_RCP26_{city}.csv'), os.path.join(rootOut, f'time_periods_NEE_RCP45_{city}.csv'), os.path.join(rootOut, f'time_periods_NEE_RCP85_{city}.csv'),
    os.path.join(rootOut, f'time_periods_LEE_RCP26_{city}.csv'), os.path.join(rootOut, f'time_periods_LEE_RCP45_{city}.csv'), os.path.join(rootOut, f'time_periods_LEE_RCP85_{city}.csv')
]

''' Modify data'''
# Loop through each file path and process the DataFrame
for file in files:
    # Read the CSV into a DataFrame
    df = pd.read_csv(file, sep=',', header=None)
    
    # Slice the elements in each cell (excluding the first 9 characters and last 2 characters)
    for i in range(df.shape[1]):  # Number of columns
        for j in range(df.shape[0]):  # Number of rows
            # Apply the slicing operation on each cell
            df.iloc[j, i] = df.iloc[j, i][9:-2]  # Modify the cell content
    
    
    # Save the modified DataFrame to a new CSV
    df.to_csv(f'{file[:-4]}_2.csv', index=False)

