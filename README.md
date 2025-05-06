


# comment: The structure is copied from the switchdrive folders and not changed. Therefore, the some of the old folder names might be confusing.

Folder structure:
	•	ema_case_Basel_HDNs: data folder with input data, simulations, figures
	⁃	case2_Base = folder for projections 
	⁃	setup3 > input = input for EMA
	⁃	setup3_bernard_136_30 = results from the runs using 136X, 30L and the equation by Bernard
	⁃	setup3_sensitivity_analysis = consists of all input data for adaptation measures / city-centre to measurement station correction, results (always with 136_30_bernard, threshold range = 26-30)
	⁃	case2_Basel_obs: folder for observed data
	⁃	setup3 > input = input for EMA
	⁃	setup3_bernard_136_30 = results from the runs using 136X, 30L and the equation by Bernard
	⁃	case2_Zurich: folder for projections
	⁃	setup3 > input = input for EMA
	⁃	setup3_bernard_136_30 = results from the runs using 136X, 30L and the equation by Bernard
	⁃	setup3_sensitivity_analysis = consists of all input data for adaptation measures / city-centre to measurement station correction, results (always with 136_30_bernard, threshold range = 26-30)
	•	python_scripts: data folder with scripts
	⁃	WF_FF_EMA_Workbench_future = EMA-workbench script for projections
	⁃	WF_FF_EMA_Workbench_obs = EMA-workbench script for observed period
	⁃	WF_FF_Publication_WBGT_136_30_Bernard = visualisation script for time series, feature scoring
	⁃	folders: self-explanatory for visualisation and preprocessing, cf. comment in scripts


Order of running scripts:
1. WF_FF_EMA_Workbench_future/obs
2. WF_FF_Publication_WBGT_136_30_Bernard or other visualisation scripts in the folders
