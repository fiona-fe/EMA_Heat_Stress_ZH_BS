'''
script purpose:             EMA with WBGT calculations using simulated climate data (1981-2099) for RCP2.6/4.5/8.5    
outputs:                    'wbgt.csv', 
                            '**name**_runs.tar.gz',
                            'Outputs_py'
                            'outSeason',
                            'outSeason_Ave_EndEvents',
                            'outSeason_Ave_LenEvents',
                            'outSeason_Ave_StartEvents',
                            'outSeason_EndEvents',
                            'outSeason_length',
                            'outSeason_NumEvents',
                            'outSeason_StartEvents
modifications required:     commented with "MODIFY"

codeauthor: Dr. Seyed Saeid Ashraf Vaghefi (2022)
adapted for Master's thesis by: Fabian Weibel (2022)
adapted for publication by: Fiona Federer (2025)
'''

# import libraries
import os
import os.path
import random
from operator import add
from datetime import datetime, date, timedelta
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import shutil
import ema_workbench
import time
import platform

# Detect operating system (apply changes if MacOS)
system = platform.system()

# Function for initiating the main dictionary of climate stations
def create_dic(a):
    '''Function: creating a dictionary for each climate station'''
    a = {}
    keys = ['DayTmin', 'TempTmax', 'TempTmin', 'DayTmax', 'DayHmd', 'ValHmd', 'DayPcp', 'ValPcp',
    'DaySlr', 'ValSlr', 'DayWnd', 'ValWnd', 'Tmax_weight', 'Tmin_weight', 'elev', 'lat', 'long', 'fileName']
    a = {key: None for key in keys}
    return a

# Function reading the metadata of the input files
def initialize_input_dict (mainFolderHDNs):
    ''' This function returns a dictionary , and addresses of all folders'''

    '''Step 1: Defining paths of the data folders''' 
    rootFolder = mainFolderHDNs
    inputFolder = os.path.join(rootFolder, 'input')
    TmaxFolder = os.path.join(inputFolder, 'Tmax')
    TminFolder = os.path.join(inputFolder, 'Tmin')
    HmdFolder = os.path.join(inputFolder, 'Hmd')
    PcpFolder = os.path.join(inputFolder, 'Pcp')
    SlrFolder = os.path.join(inputFolder, 'Slr')
    WndFolder = os.path.join(inputFolder, 'Wnd')
    climate_ref_Folder = os.path.join(inputFolder, 'Climate_ref')
    climate_Ref_Folder_org = os.path.join(inputFolder, 'Climate_ref_no_randomness_0')
    climate_ref_Folder_rand_1 = os.path.join(inputFolder, 'Climate_ref_randomness_1')
    climate_ref_Folder_rand_2 = os.path.join(inputFolder, 'Climate_ref_randomness_2')

    '''Step 2: Reading all files' names inside the Tmax, Tmin, ..., and climate folders''' 
    for filename in os.walk(TmaxFolder):
        files = filename[2]  # extract file list
        if system == "Darwin":  # if macOS
            files = [f for f in files if not f.startswith(".DS_")] # filter out .DS_Store
            files.sort() # avoid sorting problem induced by different file system handling
        TmaxFiles = files

    TminFiles = list()
    for filename in os.walk(TminFolder):
        files = filename[2]  
        if system == "Darwin": # if macOS
            files = [f for f in files if not f.startswith(".DS_")] 
            files.sort() 
        TminFiles = files

    HmdFiles = []
    for filename in os.walk(HmdFolder):
        files = filename[2]  
        if system == "Darwin":  
            files = [f for f in files if not f.startswith(".DS_")] 
            files.sort() 
        HmdFiles = files

    PcpFiles = list()
    for filename in os.walk(PcpFolder):
        files = filename[2]  
        if system == "Darwin": # if macOS
            files = [f for f in files if not f.startswith(".DS_")] 
            files.sort() 
        PcpFiles = files

    SlrFiles = list()
    for filename in os.walk(SlrFolder):
        files = filename[2]  
        if system == "Darwin": # if macOS
            files = [f for f in files if not f.startswith(".DS_")] 
            files.sort() 
        SlrFiles = files

    WndFiles = list()
    for filename in os.walk(WndFolder):
        files = filename[2]  
        if system == "Darwin": # if macOS
            files = [f for f in files if not f.startswith(".DS_")] 
            files.sort() 
        WndFiles = files

    climate_ref_Files = list()
    for filename in os.walk(climate_ref_Folder):
       climate_ref_Files = filename[2]

    '''Step 3_1: Reading files inside Tmax folder '''
    os.chdir(TmaxFolder)
    with open(TmaxFiles[0], 'r') as file: 
        weights_Inputs = file.read()
    with open(TmaxFiles[1], 'r') as file: 
        day_Tmax = file.read()
    with open(TmaxFiles[2], 'r') as file: 
        x1Tmaxthreshold = file.read()

    '''Step 3_2: Reading the lines of files inside Tmax folder'''
    weights_Inputs = weights_Inputs.replace('\n', '\t')
    weights_Inputs = weights_Inputs.split('\t')
    day_Tmax = day_Tmax.replace('\n', '\t').split('\t')
    x1Tmaxthreshold = x1Tmaxthreshold.replace('\n', '\t').split('\t')

    '''Step 3_3: Reading files and lines of files inside Tmin folder'''
    os.chdir(TminFolder)

    with open(TminFiles[0], 'r') as file:
        day_Tmin = file.read()
    with open(TminFiles[1], 'r') as file:
        x2Tminthreshold = file.read()
    day_Tmin = day_Tmin.replace('\n', '\t')
    day_Tmin = day_Tmin.split('\t')
    x2Tminthreshold = x2Tminthreshold.replace('\n', '\t').split('\t')

    '''Step 3_4: Reading files and lines of files inside Pcp folder'''
    os.chdir(PcpFolder)
    with open(PcpFiles[0], 'r') as file:
        day_Pcp = file.read()
    with open(PcpFiles[1], 'r') as file:
        x3Pcpthreshold = file.read()
    day_Pcp = day_Pcp.replace('\n', '\t')
    day_Pcp = day_Pcp.split('\t')
    x3Pcpthreshold = x3Pcpthreshold.replace('\n', '\t').split('\t')

    '''Step 3_5: Reading files and lines of files inside Hmd folder'''
    os.chdir(HmdFolder)
    with open(HmdFiles[0], 'r') as file:
        day_Hmd = file.read()
    with open(HmdFiles[1], 'r') as file:
        x4Hmdthreshold = file.read()
    day_Hmd = day_Hmd.replace('\n', '\t')
    day_Hmd = day_Hmd.split('\t')
    x4Hmdthreshold = x4Hmdthreshold.replace('\n', '\t').split('\t')

    '''Step 3_6: Reading files and lines of files inside Slr folder'''
    os.chdir(SlrFolder)
    with open(SlrFiles[0], 'r') as file:
        day_Slr = file.read()
    with open(SlrFiles[1], 'r') as file:
        x5Slrthreshold = file.read()
    day_Slr = day_Slr.replace('\n', '\t')
    day_Slr = day_Slr.split('\t')
    x5Slrthreshold = x5Slrthreshold.replace('\n', '\t').split('\t')

    '''Step 3_7: Reading files and lines of files inside Wnd folder'''
    os.chdir(WndFolder)
    with open(WndFiles[0], 'r') as file:
        day_Wnd = file.read()
    with open(WndFiles[1], 'r') as file:
        x6Wndthreshold = file.read()
    day_Wnd = day_Wnd.replace('\n', '\t')
    day_Wnd = day_Wnd.split('\t')
    x6Wndthreshold = x6Wndthreshold.replace('\n', '\t').split('\t')

    '''Step 4: Reading the lines of each file inside climate_ref folder'''
    os.chdir(climate_ref_Folder)
    with open('pcp.txt', 'r') as file:
        pcpData = file.read()
    with open('tmp.txt', 'r') as file:
        tmpData = file.read()
    with open('hmd.txt', 'r') as file:
        hmdData = file.read()
    with open('slr.txt', 'r') as file:
        srdData = file.read()
    with open('wnd.txt', 'r') as file:
        wndData = file.read()

    # Since the information of pcp.txt, tmp.txt, wnd.txt, hmd.txt, slr.txt are the same we only read pcp
    pcpData = pcpData.split('\n')
    for i in range(len(pcpData)):
        pcpData[i] = pcpData[i].split(',')

    '''Step 5: Initialazing the input dictionary of climate stations which holds the information of the stations'''
    nameStn = []
    for file in climate_ref_Files:
        if 'p.csv' in file:
            nameStn.append(file[-25: -5])
    stnDicts = []
    for i in range(len(nameStn)):
        stnDicts.append(create_dic(nameStn[i]))

    '''Step 6: Assigning the file names to the dictionary'''
    for i in range (len(nameStn)):
        stnDicts[i]['fileName'] = nameStn[i]
    
    '''Step 7: Assigning the Tamx, Tmin, Hmd, Pcp, Slr and Wnd values'''
    for stnDict in stnDicts:
        # Tmax
        for i, element in enumerate(day_Tmax):
            if element == stnDict['fileName'][:]:
                stnDict['DayTmax'] = day_Tmax[i+1]
        for i, element in enumerate(x1Tmaxthreshold):
            if element == stnDict['fileName'][:]:
                stnDict['TempTmax'] = x1Tmaxthreshold[i+1]
        for i, element in enumerate(weights_Inputs):
            stnDict['Tmax_weight'] = weights_Inputs[1]
            stnDict['Tmin_weight'] = weights_Inputs[3]

        # Tmin
        for i, element in enumerate(day_Tmin):
            if element == stnDict['fileName'][:]:
                stnDict['DayTmin'] = day_Tmin[i+1]
        for i, element in enumerate(x2Tminthreshold):
            if element == stnDict['fileName'][:]:
                stnDict['TempTmin'] = x2Tminthreshold[i+1]

        # Hmd
        for i, element in enumerate(day_Hmd):
            if element == stnDict['fileName'][:]:
                stnDict['DayHmd'] = day_Hmd[i+1]
        for i, element in enumerate(x4Hmdthreshold):
            if element == stnDict['fileName'][:]:
                stnDict['ValHmd'] = x4Hmdthreshold[i+1]

        # Pcp
        for i, element in enumerate(day_Pcp):
            if element == stnDict['fileName'][:]:
                stnDict['DayPcp'] = day_Pcp[i+1]
        for i, element in enumerate(x3Pcpthreshold):
            if element == stnDict['fileName'][:]:
                stnDict['ValPcp'] = x3Pcpthreshold[i+1]

        # Slr
        for i, element in enumerate(day_Slr):
            if element == stnDict['fileName'][:]:
                stnDict['DaySlr'] = day_Slr[i+1]
        for i, element in enumerate(x5Slrthreshold):
            if element == stnDict['fileName'][:]:
                stnDict['ValSlr'] = x5Slrthreshold[i+1]

        # Wnd
        for i, element in enumerate(day_Wnd):
            if element == stnDict['fileName'][:]:
                stnDict['DayWnd'] = day_Wnd[i+1]
        for i, element in enumerate(x6Wndthreshold):
            if element == stnDict['fileName'][:]:
                stnDict['ValWnd'] = x6Wndthreshold[i+1]

    '''Step 8: Assigning the elevation, Lat and long to the dictionaries'''
    for i in range(len(stnDicts)):
        for j in range(1, len(pcpData)):
            if pcpData[j][1][:-1] == stnDicts[i]['fileName'][:]:
                stnDicts[i]['lat']= pcpData[j][2]
                stnDicts[i]['long']= pcpData[j][3]
                stnDicts[i]['elev']= pcpData[j][4]
    return stnDicts, inputFolder, TmaxFolder, TminFolder, HmdFolder, PcpFolder, SlrFolder, WndFolder, climate_ref_Folder, \
    climate_Ref_Folder_org, climate_ref_Folder_rand_1, climate_ref_Folder_rand_2

# HDNs Model
# Initializing the main dictionary for a case study
caseStudyStns = {}
inputFolder = ''
TmaxFolder = ''
TminFolder = ''
HmdFolder = ''
PcpFolder = ''
SlrFolder = ''
WndFolder = ''
climateFolder = ''
climateFolder_org = ''
climateFolder1 = ''
climateFolder2 = ''

# Defining the root of the climate data and for the results to be stored
root = r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/Corrfactor_Tmp_1_3_5'


# calling the function with multiple return values
caseStudyStns, inputFolder, TmaxFolder, TminFolder, HmdFolder, PcpFolder, SlrFolder, WndFolder, \
climateFolder, climateFolder_org, climateFolder1, climateFolder2 = initialize_input_dict(root)

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

# 1st column as index: creating date from 01 01 1981 to 2099 12 31
from datetime import timedelta, date

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date ).days + 1)):
        yield start_date + timedelta(n)

# class Health_Policy_HDNs
class Policy_Health:
    def __init__(self, x1Tmaxthreshold, x2Tminthreshold):
        self.x1Tmaxthreshold = x1Tmaxthreshold
        self.x2Tminthreshold = x2Tminthreshold

    def policy_release2(self):
        return(self.x1Tmaxthreshold, self.x2Tminthreshold)


# class creating the RCP-ClimateModel combination
class RCP_Model:
    def __init__(self, xRCP, xClimateModel):
        self.input1 = round(xRCP)
        self.input2 = xClimateModel

    def rcpGenerator(self):
        if self.input1 == 1:
            RCP = str(2.6)
            rcpInt = 1
        if self.input1 == 2:
            RCP = str(4.5)
            rcpInt = 2
        if self.input1 == 3:
            RCP = str(8.5)
            rcpInt = 3
        return(RCP, rcpInt)

    def climateModel(self):
        a, b = RCP_Model.rcpGenerator(self)
        if b == 1:
            climateModel = round(self.input2*7)
        elif b == 2:
            climateModel = 7 + max(1,round(self.input2*17))
        else:
            climateModel = 24 + max(1, round(self.input2*22))
        return (int(climateModel))

# main function to calculate the WBGT
def HDNs_Model (xRCP=None, xClimateModel=None, Xfactor1 = None,
                 x1Tmaxthreshold = None, x2Tminthreshold = None, x7Wbgtthreshold = None):

    '''' This function controls the HDNs model in the XLR framework'''

    ''' VERY IMPORTANT --- Controling the randomness --- VERY IMPORTANT'''
    xClimateRandomness = round(Xfactor1)
    if (xClimateRandomness == 1):
        os.chdir(climateFolder_org)
        src = os.getcwd()
        os.chdir(climateFolder)
        dst = os.getcwd()
        print('Original CH2018 is being used')
    elif (xClimateRandomness == 2) :
        os.chdir(climateFolder1)
        src = os.getcwd()
        os.chdir(climateFolder)
        dst = os.getcwd()
        print('Random Climate realization version 1 is being used')
    else:
        os.chdir(climateFolder2)
        src = os.getcwd()
        os.chdir(climateFolder)
        dst = os.getcwd()
        print('Random Climate realization version 2 is being used')

    os.chdir(climateFolder)
    fnames = os.listdir()

    print('HDNs_DMDU: Matching the station names values with CSV files!')

    # Matching the station names values in the dictionary of stations with CSV files in Climate folder of the case Study
    pcpCaseStudy = []
    tmpCaseStudy = []
    hmdCaseStudy = []
    slrCaseStudy = []
    wndCaseStudy = []

    if (xClimateRandomness == 1):
        for i in range(len(caseStudyStns)):
            pcpCaseStudy.append(os.path.join(climateFolder, caseStudyStns[i]['fileName'] + 'p.csv'))
            tmpCaseStudy.append(os.path.join(climateFolder, caseStudyStns[i]['fileName'] + 't.csv'))
            hmdCaseStudy.append(os.path.join(climateFolder, caseStudyStns[i]['fileName'] + 'h.csv'))
            slrCaseStudy.append(os.path.join(climateFolder, caseStudyStns[i]['fileName'] + 's.csv'))
            wndCaseStudy.append(os.path.join(climateFolder, caseStudyStns[i]['fileName'] + 'w.csv'))
    elif (xClimateRandomness == 2) :
        for i in range(len(caseStudyStns)):
            pcpCaseStudy.append(os.path.join(climateFolder1, caseStudyStns[i]['fileName'] + 'p.csv'))
            tmpCaseStudy.append(os.path.join(climateFolder1, caseStudyStns[i]['fileName'] + 't.csv'))
            hmdCaseStudy.append(os.path.join(climateFolder1, caseStudyStns[i]['fileName'] + 'h.csv'))
            slrCaseStudy.append(os.path.join(climateFolder1, caseStudyStns[i]['fileName'] + 's.csv'))
            wndCaseStudy.append(os.path.join(climateFolder1, caseStudyStns[i]['fileName'] + 'w.csv'))
    else:
        for i in range(len(caseStudyStns)):
            pcpCaseStudy.append(os.path.join(climateFolder2, caseStudyStns[i]['fileName'] + 'p.csv'))
            tmpCaseStudy.append(os.path.join(climateFolder2, caseStudyStns[i]['fileName'] + 't.csv'))
            hmdCaseStudy.append(os.path.join(climateFolder2, caseStudyStns[i]['fileName'] + 'h.csv'))
            slrCaseStudy.append(os.path.join(climateFolder2, caseStudyStns[i]['fileName'] + 's.csv'))
            wndCaseStudy.append(os.path.join(climateFolder2, caseStudyStns[i]['fileName'] + 'w.csv'))

    print('HDNs_DMDU: Building a database for each csv file (tmp, pcp, hmd, slr, wnd)!')

    '''Step 11: building a database for each precipitation/temperature/humidity/solar radiation/wind file in Climate folder and saving them in a list'''
    '''11.1 reading the csv files as databases'''
    dfpcp = [None for _ in range(len(pcpCaseStudy))]
    dftmp = [None for _ in range(len(tmpCaseStudy))]
    dfhmd = [None for _ in range(len(hmdCaseStudy))]
    dfslr = [None for _ in range(len(slrCaseStudy))]
    dfwnd = [None for _ in range(len(wndCaseStudy))]

    for i in range(len(pcpCaseStudy)):
        dfpcp[i] = pd.read_csv(pcpCaseStudy[i])
        dftmp[i] = pd.read_csv(tmpCaseStudy[i])
        dfhmd[i] = pd.read_csv(hmdCaseStudy[i])
        dfslr[i] = pd.read_csv(slrCaseStudy[i])
        dfwnd[i] = pd.read_csv(wndCaseStudy[i])

    '''11.2 making a header for output files'''
    dfpcpCol = dfpcp[0].columns
    dftmpCol = dftmp[0].columns
    dfhmdCol = dfhmd[0].columns
    dfslrCol = dfslr[0].columns
    dfwndCol = dfwnd[0].columns

    '''11.3 defining the length of simulations and scenarios'''
    scenariosLength = len(dfpcpCol)
    simulationLength = len(dftmp[0][dftmpCol[0]]) - 1

    '''11.4 Reading the beginning and end of the simulation'''
    start_date = date(1981, 1, 1)
    end_date = date(2099, 12, 31)
    dateList = []
    for single_date in daterange(start_date, end_date):
        dateList.append(single_date.strftime("%m/%d/%Y"))
    seasonList = []
    for n in range (1981, 2100, 1):
        seasonList.append(str(n))

    print('HDNs_DMDU: Part 1 Running the model, looking for extreme events and printing the output!')
    '''################################ PART1 ################################'''
    '''Running the model for each climate station:'''

    for k in range(len(caseStudyStns)):

        '''making a header for output files'''
        dfpcpCol = dfpcp[k].columns
        dftmpCol = dftmp[k].columns
        dfhmdCol = dfhmd[k].columns
        dfslrCol = dfslr[k].columns
        dfwndCol = dfwnd[k].columns

        '''defining the length of simulations and scenarios'''
        scenariosLength = 1
        simulationLength = len(dftmp[0][dftmpCol[0]]) - 1

        '''declaring the initial arrays'''
        wbgt = [0 for _ in range(simulationLength)]
        is_extreme_wbgt = [0 for _ in range(simulationLength)]
        total = np.zeros([simulationLength, 3*scenariosLength])

        '''RCP and Climate Model Controler'''
        rcp_Model = RCP_Model(xRCP, xClimateModel)
        RCP, intRCP = rcp_Model.rcpGenerator()
        climateModel = rcp_Model.climateModel()

        '''Running the model for each climate scenario:'''
        for j in range(climateModel, climateModel + 1, 1):

            # Reading the information and inputs of the FIRST DAY of simulation
            todayPCP = dfpcp[k][dfpcpCol[j]].iloc[1] if (dfpcp[k][dfpcpCol[j]].iloc[1] != -99) else 0
            todayTMPMAX = round(dftmp[k][dftmpCol[2*j]].iloc[1],2) if(dftmp[k][dftmpCol[2*j]].iloc[1] != -99) else 0
            todayTMPMIN = round(dftmp[k][dftmpCol[2*j+1]].iloc[1],2) if(dftmp[k][dftmpCol[2*j+1]].iloc[1] != -99) else 0
            todayTMPAVE = round((todayTMPMAX+todayTMPMIN)/2,2) if((todayTMPMAX+todayTMPMIN)/2 != -99) else 0
            todayHMD = round(dfhmd[k][dfhmdCol[j]].iloc[1],2) if(dfhmd[k][dfhmdCol[j]].iloc[1] != -99) else 0
            todayWND = dfwnd[k][dfwndCol[j]].iloc[1] if(dfwnd[k][dfwndCol[j]].iloc[1] != -99) else 0
            todaySLR = round(dfslr[k][dfslrCol[j]].iloc[1],2) if(dfslr[k][dfslrCol[j]].iloc[1] != -99) else 0

            '''threshold Temperatures C
            EMA_workbench_controler for the threshold of hot days and hot nights'''

            policyCities = Policy_Health(x1Tmaxthreshold, x2Tminthreshold) 
            x1Tmaxthreshold, x2Tminthreshold = policyCities.policy_release2()

            '''Equations to calculate WBGT ### Input data units: todayTMPMAX --> °C, todayHMD --> %, todayWND --> m/s, todaySLR --> W/m2'''

            '''1.) http://www.bom.gov.au/info/thermal_stress/#approximation, WBGT'''
            '''wbgt = 0.567 × Ta + 0.393 × e + 3.94, '''
            '''e = rh / 100 × 6.105 × exp ( 17.27 × Ta / ( 237.7 + Ta ) )'''

            # wbgt[0] =  0.567 * todayTMPMAX + 0.393 * (todayHMD / 100 * 6.105 * np.exp(17.27 * todayTMPMAX / (237.7 + todayTMPMAX))) + 3.94

            '''2.) http://www.bom.gov.au/info/thermal_stress/#approximation, AT'''
            '''AT = Ta + 0.348×e - 0.70×ws + 0.70×Q/(ws + 10) - 4.25'''
            '''e = rh / 100 × 6.105 × exp ( 17.27 × Ta / ( 237.7 + Ta ) )'''

            # wbgt[0] =  todayTMPMAX + 0.348 * (todayHMD / 100 * 6.105 * np.exp(17.27 * todayTMPMAX / (237.7 + todayTMPMAX))) - 0.7 * float(todayWND) + 0.7 * (todaySLR/((float(todayWND))+10)) - 4.25

            # e = todayHMD/100 * 6.105 * np.exp(17.27*todayTMPMAX/(237.7+todayTMPMAX))

            # wbgt[0] =  todayTMPMAX + 0.348 * e - 0.7 * float(todayWND) + 0.7 * todaySLR/((float(todayWND))+10) - 4.25

            # if wbgt[0] >= 45:
            #     wbgt[0] = 45
            # elif wbgt[0] <= -10:
            #     wbgt[0] = -10

            '''3.) Bernard, 1999 and Carter et al., 2020'''
            '''WBGT = 0.7*Tnwb + 0.2*Tg + 0.1*Ta'''

            Tg = [0 for _ in range(simulationLength)]
            vapre = [0 for _ in range(simulationLength)]
            Tpwb = [0 for _ in range(simulationLength)]
            Tnwb = [0 for _ in range(simulationLength)]
            C = [0 for _ in range(simulationLength)]

            Tg[0] = 0.009624*todaySLR + 1.102*todayTMPMAX - 0.00404*todayHMD - 2.2776
            vapre[0] = (todayHMD/100)*(0.6107*np.exp((17.27*todayTMPMAX)/(todayTMPMAX + 237.3)))
            Tpwb[0] = 0.376 + 5.79*vapre[0] + (0.388 - 0.0465*vapre[0])*todayTMPMAX

            if (Tg[0] - todayTMPMAX) < 4:
                if todayWND < 0.03:
                    C[0] = 0.85
                elif todayWND > 3.0:
                    C[0] = 1.0
                else:
                    C[0] = 0.96 + 0.069*np.log10(todayWND)
                Tnwb[0] = todayTMPMAX - C[0]*(todayTMPMAX - Tpwb[0])

            if (Tg[0] - todayTMPMAX) >= 4:
                if todayWND < 0.1:
                    C[0] = 1.1
                elif todayWND > 1.0:
                    C[0] = -0.1
                else:
                    C[0] = 0.10/(todayWND**1.1) - 0.2
                Tnwb[0] = Tpwb[0] + 0.25*(Tg[0] - todayTMPMAX) + C[0]

            wbgt[0] = 0.7*Tnwb[0] + 0.2*Tg[0] + 0.1*todayTMPMAX

            '''4.) Ono and Tonouchi, 2014; Heo, Bell and Lee, 2019'''
            '''WBGT=0.735×Ta + 0.0374×RH + 0.00292×Ta×RH + 7.619×SR - 4.557×SR^2 - 0.0572×WS - 4.064'''

            # wbgt[0] = 0.735*todayTMPMAX + 0.0374*todayHMD + 0.00292*todayTMPMAX*todayHMD + 7.619*(todaySLR/1000) - 4.557*((todaySLR/1000)*(todaySLR/1000)) - 0.0572*todayWND - 4.064

            '''Checking whether WBGT is above a given threshold for first day:'''
            if wbgt[0] >= x7Wbgtthreshold:
                is_extreme_wbgt[0] = 1
            else:
                is_extreme_wbgt[0] = 0

            '''storing three values in a list for the first day ready fo printing in the csv file:'''
            total[0,2] = round(is_extreme_wbgt[0], 2)

            # For the SECOND DAY to the END DAY of Simulation
            i = 0
            for i in range(2, simulationLength + 1, 1):
                '''# precipitation and temperature missing values were handled'''
                todayPCP = dfpcp[k][dfpcpCol[j]].iloc[i] if (dfpcp[k][dfpcpCol[j]].iloc[i] != -99) else 0
                todayTMPMAX = round(dftmp[k][dftmpCol[2*j]].iloc[i],2) if(dftmp[k][dftmpCol[2*j]].iloc[i] != -99) else 0
                todayTMPMIN = round(dftmp[k][dftmpCol[2*j+1]].iloc[i],2) if(dftmp[k][dftmpCol[2*j+1]].iloc[i] != -99) else 0
                todayTMPAVE = round((todayTMPMAX+todayTMPMIN)/2,2) if((todayTMPMAX+todayTMPMIN)/2 != -99) else 0
                todayHMD = round(dfhmd[k][dfhmdCol[j]].iloc[i],2) if(dfhmd[k][dfhmdCol[j]].iloc[i] != -99) else 0
                todayWND = dfwnd[k][dfwndCol[j]].iloc[i] if(dfwnd[k][dfwndCol[j]].iloc[i] != -99) else 0
                todaySLR = round(dfslr[k][dfslrCol[j]].iloc[i],2) if(dfslr[k][dfslrCol[j]].iloc[i] != -99) else 0

                '''Equations to calculate WBGT ### Input data units: todayTMPMAX --> °C, todayHMD --> %, todayWND --> m/s, todaySLR --> W/m2'''

                '''1.) http://www.bom.gov.au/info/thermal_stress/#approximation, WBGT'''
                '''wbgt = 0.567 × Ta + 0.393 × e + 3.94, '''
                '''e = rh / 100 × 6.105 × exp ( 17.27 × Ta / ( 237.7 + Ta ) )'''

                # wbgt[i-1] =  0.567 * todayTMPMAX + 0.393 * (todayHMD / 100 * 6.105 * np.exp(17.27 * todayTMPMAX / (237.7 + todayTMPMAX))) + 3.94

                '''2.) http://www.bom.gov.au/info/thermal_stress/#approximation, AT'''
                '''AT = Ta + 0.348×e - 0.70×ws + 0.70×Q/(ws + 10) - 4.25'''
                '''e = rh / 100 × 6.105 × exp ( 17.27 × Ta / ( 237.7 + Ta ) )'''

                # wbgt[i-1] =  todayTMPMAX + 0.348 * (todayHMD / 100 * 6.105 * np.exp(17.27 * todayTMPMAX / (237.7 + todayTMPMAX))) - 0.7 * float(todayWND) + 0.7 * (todaySLR/((float(todayWND))+10)) - 4.25

                # e = todayHMD/100 * 6.105 * np.exp(17.27*todayTMPMAX/(237.7+todayTMPMAX))

                # wbgt[i-1] =  todayTMPMAX + 0.348 * e - 0.7 * float(todayWND) + 0.7 * todaySLR/((float(todayWND))+10) - 4.25

                # if wbgt[i-1] >= 45:
                #     wbgt[i-1] = 45
                # elif wbgt[i-1] <= -10:
                #     wbgt[i-1] = -10

                '''3.) Carter et al., 2020'''
                '''WBGT = 0.7*Tnwb + 0.2*Tg + 0.1*Ta'''

                Tg[i-1] = 0.009624*todaySLR + 1.102*todayTMPMAX - 0.00404*todayHMD - 2.2776
                vapre[i-1] = (todayHMD/100)*(0.6107*np.exp((17.27*todayTMPMAX)/(todayTMPMAX + 237.3)))
                Tpwb[i-1] = 0.376 + 5.79*vapre[i-1] + (0.388 - 0.0465*vapre[i-1])*todayTMPMAX

                if (Tg[i-1] - todayTMPMAX) < 4:
                    if todayWND < 0.03:
                        C[i-1] = 0.85
                    elif todayWND > 3.0:
                        C[i-1] = 1.0
                    else:
                        C[i-1] = 0.96 + 0.069*np.log10(todayWND)
                    Tnwb[i-1] = todayTMPMAX - C[i-1]*(todayTMPMAX - Tpwb[i-1])

                if (Tg[i-1] - todayTMPMAX) >= 4:
                    if todayWND < 0.1:
                        C[i-1] = 1.1
                    elif todayWND > 1.0:
                        C[i-1] = -0.1
                    else:
                        C[i-1] = 0.10/(todayWND**1.1) - 0.2
                    Tnwb[i-1] = Tpwb[i-1] + 0.25*(Tg[i-1] - todayTMPMAX) + C[i-1]

                wbgt[i-1] = 0.7*Tnwb[i-1] + 0.2*Tg[i-1] + 0.1*todayTMPMAX

                '''4.) Ono and Tonouchi, 2014; Heo, Bell and Lee, 2019'''
                '''WBGT=0.735×Ta + 0.0374×RH + 0.00292×Ta×RH + 7.619×SR - 4.557×SR^2 - 0.0572×WS - 4.064'''

                # wbgt[i-1] = 0.735*todayTMPMAX + 0.0374*todayHMD + 0.00292*todayTMPMAX*todayHMD + 7.619*(todaySLR/1000) - 4.557*((todaySLR/1000)*(todaySLR/1000)) - 0.0572*todayWND - 4.064

                ''' Checking whether WBGT is above a given threshold'''
                if wbgt[i-1] >= x7Wbgtthreshold:
                    is_extreme_wbgt[i-1] = 1

                else:
                    is_extreme_wbgt[i-1] = 0

                '''storing three values in a list for the first day ready fo printing in the csv file:'''
                total[i-1,2] = round(is_extreme_wbgt[i-1], 2)

        '''Saving the Outputs of total list in a CSV file in a specific path'''
        # 1st row as the column names:
        columnsDF = []
        nameHeader = dfpcpCol[climateModel]
        columnsDF.append('is_Tmax_exEve' + nameHeader)
        columnsDF.append('is_Tmin_exEve' + nameHeader)
        columnsDF.append('is_Wbgt_exEve' + nameHeader)

        '''Extreme analyses daily'''
        columnsDF0 = ['DATE']
        dfnew0 = pd.DataFrame(dateList, columns = columnsDF0)
        dfnew1 = pd.DataFrame(total, columns = columnsDF)
        df1 = pd.concat([dfnew0, dfnew1], axis=1, sort=False)

        if os.path.isdir(os.path.join(root, 'Outputs_py')):
            pass
        else: os.mkdir(os.path.join(root, 'Outputs_py'))

        '''Make CSvs for daily extreme Outputs'''
        outfolder =os.path.join(root, 'Outputs_py')
        outfileName = 'Total_daily_' + caseStudyStns[k]['fileName'] + '.csv'
        outputFile = os.path.join(outfolder, outfileName )
        df1.to_csv(outputFile, index = False)

        '''Saving raw wbgt values in a csv file'''
        wbgt_df = pd.DataFrame(wbgt)
        wbgt_df.to_csv(r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/wbgt.csv', index=False)



        print('##############################End of Part 1 Calculations!##############################')

        '''################################ PART2 ################################'''
        '''##### PART 2 Annual outputs #####'''

        print('####HDNS_Model: Starting Part 2 Running the model, Annual (seasonal) outputs, reading files!####')

        total_Daily_FilesAll = list()
        total_Daily_Files = []

        for filename in os.walk(outfolder):
            total_Daily_FilesAll = filename[2]

        for bIndex in range (len(total_Daily_FilesAll)):
            if 'Total_daily_' in total_Daily_FilesAll[bIndex]:
                total_Daily_Files.append(total_Daily_FilesAll[bIndex])
            else: continue

        # Adding the whole address of directory to the name of total daily files
        totalFiles = []
        for i in range(len(total_Daily_Files)):
            totalFiles.append(os.path.join(outfolder, total_Daily_Files[i]))

        print('#### HDNS Model: Continuing of Part 2 Annual Outputs, Performing XLRM Analyses!####')

        # databases are read here:
        dfSeason = [ None for _ in range(len(totalFiles))]

        # Here we calcluate seasonal tipping points
        for i in range(len(totalFiles)):
            dfSeason[i] = pd.read_csv(totalFiles[i], low_memory=False)

            start_date = date(1981, 1, 1)
            end_date = date(2099, 12, 31)

            dateList = []
            for single_date in daterange(start_date, end_date):
                dateList.append(single_date.strftime("%m/%d/%Y"))

            start_season = []
            end_season = []

            for pp in range (1981, 2100, 1):
                start_season.append(date(pp, 1, 1))
                end_season.append(date(pp, 12, 31))

            df2 = dfSeason[i]
            df2.set_index('DATE', inplace = True)
            df2Col = df2.columns

            df2ColCal = []

            '''We perform our analysis for one climate scenario at a time inseated of 68'''
            '''m+0 --> Tmin, m+1 --> Tmax, m+2 --> WBGT'''
            for m in range(1):
                df2ColCal.append(df2Col[3*m+2])

            sumHotCondition = np.zeros([len(start_season), len(df2ColCal)])
            sumRows = np.zeros(len(start_season))
            Average_indexLength =  np.zeros(len(start_season))
            Average_StartEvents =  np.zeros(len(start_season))
            Average_EndEvents =  np.zeros(len(start_season))
            num_Events = np.zeros(len(start_season))

            for j in range(len(df2ColCal)):
                # matrix with 118 years but variable length
                indexLength = [[] for i in range(len(start_season))]
                indexStartEvents = [[] for i in range(len(start_season))]
                indexEndEvents = [[] for i in range(len(start_season))]

                for k in range(len(start_season)):
                    start_date = start_season[k]
                    end_date = end_season[k]
                    'How to calculate different length of concurrent HDNs'
                    PP = 0
                    qqDayYear = -1
                    for single_date in daterange(start_date, end_date):
                        sumHotCondition[k,j] += df2[df2ColCal[j]].loc[single_date.strftime("%m/%d/%Y")]
                        qqDayYear += 1
                        if df2[df2ColCal[j]].loc[single_date.strftime("%m/%d/%Y")] == 1:
                            PP += 1
                        else:
                            if PP > 0:
                                indexLength[k].append(PP)
                                indexStartEvents[k].append(qqDayYear-PP)
                                indexEndEvents[k].append(qqDayYear)
                                PP = 0
                    # we need to consider the last day and say if the single_date == end_date >> indexLength[k].append(PP)

                    'Here we calculates the outputs'
                    sumRows[k] +=  sumHotCondition[k,j]  ### aa = 11 if len(b) != 0 else 0
                    Average_indexLength[k] = np.average(indexLength[k]) if len(indexLength[k]) !=0 else 0
                    Average_StartEvents[k] = np.average(indexStartEvents[k]) if len(indexStartEvents[k]) !=0 else 0
                    Average_EndEvents[k] = np.average(indexEndEvents[k]) if len(indexEndEvents[k]) !=0 else 0
                    num_Events[k]= len(indexLength[k])

            # First average all events in one year and then average all years (1981-2100) to get one number
            # These Averages of averages ar not usefull in my opinion !!
            AveragesumRows = np.average(sumRows) ## Saeed 2020/08/17
            Ave_Average_indexLength = np.average(Average_indexLength)
            Ave_Average_StartEvents = np.average(Average_StartEvents)
            Ave_Average_EndEvents = np.average(Average_EndEvents)
            Ave_num_Events = np.average(num_Events)

            df3 = pd.DataFrame(sumHotCondition, columns = df2ColCal)
            df4 = pd.DataFrame(indexLength)
            df4_Start_Date_Event = pd.DataFrame(indexStartEvents)
            df4_End_Date_Event = pd.DataFrame(indexEndEvents)
            df4_num_Events = pd.DataFrame(num_Events, columns = df2ColCal)
            df5_Ave_length =  pd.DataFrame(Average_indexLength, columns = df2ColCal)
            df5_Ave_StartEvents =  pd.DataFrame(Average_StartEvents, columns = df2ColCal)
            df5_Ave_EndEvents =  pd.DataFrame(Average_EndEvents, columns = df2ColCal)

            firstCol = []
            for o in range (len(seasonList)):
                firstCol.append(seasonList[o])

            columnsDF1 = ['Season']
            dfnew3 = pd.DataFrame(firstCol, columns = columnsDF1)

            '''Building 4 seasonal output'''
            dfFinalSeason = pd.concat([dfnew3, df3], axis=1, sort=False)
            dfFinalSeason_length = pd.concat([dfnew3, df4], axis=1, sort=False)
            dfFinalSeason_StartEvents = pd.concat([dfnew3, df4_Start_Date_Event], axis=1, sort=False)
            dfFinalSeason_EndEvents = pd.concat([dfnew3, df4_End_Date_Event], axis=1, sort=False)
            dfFinalSeason_NumEvents = pd.concat([dfnew3, df4_num_Events], axis=1, sort=False)

            '''Building 3 average seasonal output'''
            dfFinalSeason_Ave_length = pd.concat([dfnew3, df5_Ave_length], axis=1, sort=False)
            dfFinalSeason_Ave_StartEvents = pd.concat([dfnew3, df5_Ave_StartEvents], axis=1, sort=False)
            dfFinalSeason_Ave_EndEvents = pd.concat([dfnew3, df5_Ave_EndEvents], axis=1, sort=False)

            '''OutPut number of individual extreme HDNs days '''
            if os.path.isdir(os.path.join(root, 'outSeason')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason'))

            outfileNameSeason = 'season_' + total_Daily_Files[i]
            outFolderSeason = os.path.join(root, 'outSeason')
            outputFileSeason = os.path.join(outFolderSeason, outfileNameSeason)

            outFilesFinal = []
            for filename in os.walk(outFolderSeason):
                outFilesFinal = filename[2]
                iii = len(outFilesFinal)
                if os.path.isfile(outputFileSeason):
                    newOutFileNameSeason = outputFileSeason[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason.to_csv(newOutFileNameSeason, index = False)
                else:
                    dfFinalSeason.to_csv(outputFileSeason, index = False)

            '''OutPut Lengths of individual events'''
            if os.path.isdir(os.path.join(root, 'outSeason_length')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_length'))

            outfileNameSeason_length = 'season_' + total_Daily_Files[i]
            outFolderSeason_length = os.path.join(root, 'outSeason_length')
            outputFileSeason_length = os.path.join(outFolderSeason_length, outfileNameSeason_length)

            outFilesFinal_length = []
            for filename_length in os.walk(outFolderSeason_length):
                outFilesFinal_length = filename_length[2]
                iii = len(outFilesFinal_length)
                if os.path.isfile(outputFileSeason_length):
                    newOutFileNameSeason_length = outputFileSeason_length[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_length.to_csv(newOutFileNameSeason_length, index = False)
                else:
                    dfFinalSeason_length.to_csv(outputFileSeason_length, index = False)

            '''OutPut Start Dates of individual events'''
            if os.path.isdir(os.path.join(root, 'outSeason_StartEvents')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_StartEvents'))

            outfileNameSeason_StartEvents = 'season_' + total_Daily_Files[i]
            outFolderSeason_StartEvents = os.path.join(root, 'outSeason_StartEvents')
            outputFileSeason_StartEvents = os.path.join(outFolderSeason_StartEvents, outfileNameSeason_StartEvents)

            outFilesFinal_StartEvents = []
            for filename_StartEvents in os.walk(outFolderSeason_StartEvents):
                outFilesFinal_StartEvents = filename_StartEvents[2]
                iii = len(outFilesFinal_StartEvents)
                if os.path.isfile(outputFileSeason_StartEvents):
                    newOutFileNameSeason_StartEvents = outputFileSeason_StartEvents[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_StartEvents.to_csv(newOutFileNameSeason_StartEvents, index = False)
                else:
                    dfFinalSeason_StartEvents.to_csv(outputFileSeason_StartEvents, index = False)

            '''OutPut End Dates of individual events'''
            if os.path.isdir(os.path.join(root, 'outSeason_EndEvents')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_EndEvents'))

            outfileNameSeason_EndEvents = 'season_' + total_Daily_Files[i]
            outFolderSeason_EndEvents = os.path.join(root, 'outSeason_EndEvents')
            outputFileSeason_EndEvents = os.path.join(outFolderSeason_EndEvents, outfileNameSeason_EndEvents)

            outFilesFinal_EndEvents = []
            for filename_EndEvents in os.walk(outFolderSeason_EndEvents):
                outFilesFinal_EndEvents = filename_EndEvents[2]
                iii = len(outFilesFinal_EndEvents)
                if os.path.isfile(outputFileSeason_EndEvents):
                    newOutFileNameSeason_EndEvents = outputFileSeason_EndEvents[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_EndEvents.to_csv(newOutFileNameSeason_EndEvents, index = False)
                else:
                    dfFinalSeason_EndEvents.to_csv(outputFileSeason_EndEvents, index = False)

            '''OutPut Number of individual events'''
            if os.path.isdir(os.path.join(root, 'outSeason_NumEvents')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_NumEvents'))

            outfileNameSeason_NumEvents = 'season_' + total_Daily_Files[i]
            outFolderSeason_NumEvents = os.path.join(root, 'outSeason_NumEvents')
            outputFileSeason_NumEvents = os.path.join(outFolderSeason_NumEvents, outfileNameSeason_NumEvents)

            outFilesFinal_NumEvents = []
            for filename_NumEvents in os.walk(outFolderSeason_NumEvents):
                outFilesFinal_NumEvents = filename_NumEvents[2]
                iii = len(outFilesFinal_NumEvents)
                if os.path.isfile(outputFileSeason_NumEvents):
                    newOutFileNameSeason_NumEvents = outputFileSeason_NumEvents[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_NumEvents.to_csv(newOutFileNameSeason_NumEvents, index = False)
                else:
                    dfFinalSeason_NumEvents.to_csv(outputFileSeason_NumEvents, index = False)

            '''OutPut Ave Length of events'''
            if os.path.isdir(os.path.join(root, 'outSeason_Ave_LenEvents')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_Ave_LenEvents'))

            outfileNameSeason_Ave_LenEvents = 'season_' + total_Daily_Files[i]
            outFolderSeason_Ave_LenEvents = os.path.join(root, 'outSeason_Ave_LenEvents')
            outputFileSeason_Ave_LenEvents = os.path.join(outFolderSeason_Ave_LenEvents, outfileNameSeason_Ave_LenEvents)

            outFilesFinal_Ave_LenEvents = []
            for filename_Ave_LenEvents in os.walk(outFolderSeason_Ave_LenEvents):
                outFilesFinal_Ave_LenEvents = filename_Ave_LenEvents[2]
                iii = len(outFilesFinal_Ave_LenEvents)
                if os.path.isfile(outputFileSeason_Ave_LenEvents):
                    newOutFileNameSeason_Ave_LenEvents = outputFileSeason_Ave_LenEvents[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_Ave_length.to_csv(newOutFileNameSeason_Ave_LenEvents, index = False)
                else:
                    dfFinalSeason_Ave_length.to_csv(outputFileSeason_Ave_LenEvents, index = False)

            '''OutPut Ave Start of events'''
            if os.path.isdir(os.path.join(root, 'outSeason_Ave_StartEvents')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_Ave_StartEvents'))

            outfileNameSeason_Ave_StartEvents = 'season_' + total_Daily_Files[i]
            outFolderSeason_Ave_StartEvents = os.path.join(root, 'outSeason_Ave_StartEvents')
            outputFileSeason_Ave_StartEvents = os.path.join(outFolderSeason_Ave_StartEvents, outfileNameSeason_Ave_StartEvents)

            outFilesFinal_Ave_StartEvents = []
            for filename_Ave_StartEvents in os.walk(outFolderSeason_Ave_StartEvents):
                outFilesFinal_Ave_StartEvents = filename_Ave_StartEvents[2]
                iii = len(outFilesFinal_Ave_StartEvents)
                if os.path.isfile(outputFileSeason_Ave_StartEvents):
                    newOutFileNameSeason_Ave_StartEvents = outputFileSeason_Ave_StartEvents[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_Ave_StartEvents.to_csv(newOutFileNameSeason_Ave_StartEvents, index = False)
                else:
                    dfFinalSeason_Ave_StartEvents.to_csv(outputFileSeason_Ave_StartEvents, index = False)

            '''OutPut Ave End of events'''
            if os.path.isdir(os.path.join(root, 'outSeason_Ave_EndEvents')):
                pass
            else:
                os.mkdir(os.path.join(root, 'outSeason_Ave_EndEvents'))

            outfileNameSeason_Ave_EndEvents = 'season_' + total_Daily_Files[i]
            outFolderSeason_Ave_EndEvents = os.path.join(root, 'outSeason_Ave_EndEvents')
            outputFileSeason_Ave_EndEvents = os.path.join(outFolderSeason_Ave_EndEvents, outfileNameSeason_Ave_EndEvents)

            outFilesFinal_Ave_EndEvents = []
            for filename_Ave_EndEvents in os.walk(outFolderSeason_Ave_EndEvents):
                outFilesFinal_Ave_EndEvents = filename_Ave_EndEvents[2]
                iii = len(outFilesFinal_Ave_EndEvents)
                if os.path.isfile(outputFileSeason_Ave_EndEvents):
                    newOutFileNameSeason_Ave_EndEvents = outputFileSeason_Ave_EndEvents[0 : -4] + '_' + str(iii) + '.csv'
                    dfFinalSeason_Ave_EndEvents.to_csv(newOutFileNameSeason_Ave_EndEvents, index = False)
                else:
                    dfFinalSeason_Ave_EndEvents.to_csv(outputFileSeason_Ave_EndEvents, index = False)

            print('############################## End of all calculations ##############################')

        return{'S_Ave_HotDay' : AveragesumRows, 'GCM_RCM' : climateModel, 'y2' : dfpcpCol[climateModel], 'S_HotDay' : sumRows,
        'S_Ave_Length' : Ave_Average_indexLength, 'S_Length' : Average_indexLength,
        'S_Ave_StartEvent' : Ave_Average_StartEvents, 'S_StartEvent' : Average_StartEvents,
        'S_Ave_EndEvent' : Ave_Average_EndEvents, 'S_EndEvent' : Average_EndEvents,
        'S_Ave_numEvent' : Ave_num_Events, 'S_numEvent' : num_Events}

# EMA_Workbench connector
'''
Created on 20 dec. 2010

This file illustrated the use the EMA classes for a contrived example
It's main purpose has been to test the parallel processing functionality

.. codeauthor:: jhkwakkel <j.h.kwakkel (at) tudelft (dot) nl>
'''
from ema_workbench import (Model, RealParameter, Constant, ScalarOutcome, ema_logging, IntegerParameter,
                          CategoricalParameter, perform_experiments, TimeSeriesOutcome, ArrayOutcome)
from ema_workbench import (MultiprocessingEvaluator)

from ema_workbench.util import ema_logging

''' Try logging'''
# Enable logging to stderr
ema_logging.log_to_stderr()
# Get the logger
logger = ema_logging.get_rootlogger()
''' End logging'''


# import time
start_time = time.time()

if __name__ == '__main__':
    ema_logging.LOG_FORMAT = '[%(name)s/%(levelname)s/%(processName)s] %(message)s'
    ema_logging.log_to_stderr(ema_logging.INFO)

    model = Model('UZHModel', function = HDNs_Model)  # instantiate the model

    # specify process model parameters  xRCP=None, xClimateModel=None
    model.uncertainties = [RealParameter("Xfactor1",  0.51, 3.49), # default: ("xfactor1", 0.51, 3.49), defines climate randomness
                           RealParameter("xRCP", 0.51, 3.49), # default: ("xRCP", 0.51, 3.49) cf. rcpGenerator
                           RealParameter("xClimateModel", 0, 1), # default: ("xClimateModel", 0, 1) cf. climateModel
                           ]

    # specify polices IntegerParameter
    model.levers = [
                    RealParameter("x7Wbgtthreshold", 26, 30) # default: ("x7Wbgtthreshold", 26, 30), defines threshold range indicating extreme events
                    ]

    # specify outcomes
    model.outcomes = [ScalarOutcome('S_Ave_HotDay'),
                      ScalarOutcome('GCM_RCM'),
                      ScalarOutcome('S_Ave_Length'),
                      ScalarOutcome('S_Ave_StartEvent'),
                      ScalarOutcome('S_Ave_EndEvent'),
                      ScalarOutcome('S_Ave_numEvent'),
                      ArrayOutcome('S_HotDay'),
                      ArrayOutcome('S_Length'),
                      ArrayOutcome('S_StartEvent'),
                      ArrayOutcome('S_EndEvent'),
                      ArrayOutcome('S_numEvent')
                      ]

    try:
        results = perform_experiments(model, 136, 30)
    except Exception as e:
        logger.error(f"An error occurred during the model run: {e}")


print('end!')
training_time = time.time() - start_time

print("--- %s seconds ---" % (training_time))
print('training time : {} mins and {} seconds'.format((training_time // 60) , round((training_time % 60), 1)))
print('training time : {} hours {} mins and {} seconds '.format(training_time // 3600 , round((training_time % 3600 // 60), 1), round((training_time % 3600) % 60 ,1)))

# Save the outputs
from ema_workbench import save_results

save_results(results, r'/Users/erdtiger/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Studium/EClim/Publication_Thesis_fweibel/ema_case_Basel_HDNs/case2_Zurich/setup3_sensitivity_analysis/20250325_runs.tar.gz')

print('end!')

# Sources for WBGT equations:

# - Steadman, R. G. (1994) ‘Norms of apparent temperature in Australia’, Australian Meteorological Magazine, 43, pp. 1–16.

# - Bernard, T. E. (1999) ‘Prediction of Workplace Wet Bulb Global Temperature’, Applied Occupational and Environmental Hygiene, 14(2), pp. 126–134. doi: https://doi.org/10.1080/104732299303296.

# - Ono, M. and Tonouchi, M. (2014) ‘Estimation of wet-bulb globe temperature using generally measured meteorological indices’, Japanese Journal of Biometeorology, 50(4), pp. 147–157. doi: 10.11227/seikisho.50.147.
