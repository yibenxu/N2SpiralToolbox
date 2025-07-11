# N2SpiralToolbox_SleepEEG
MATLAB toolbox to automatically detect and analyse spiral wave patterns in sleep EEG data, developed by Prof. Pulin Gong's group at University of Sydney.
## Instructions for use
System dependencies: N/A <br />
Software version: MATLAB 2016b and above (has been tested on MATLAB 2016b and 2020a) <br />
Hardware requirement: N/A

Data format: Sleep EEG data (256-channel HdEEG) in .fdt and .set formats (EEGLAB) <br />
Data tested: 9 patients with obstructive sleep apnea (OSA), recorded during stage 2 NREM (N2) sleep <br />

### Launch: <br />
The high-density EEG source data used in this study are not publicly available but can be accessed upon request. 

Please download EEGLAB2024 (https://sccn.ucsd.edu/eeglab/) before launch.

Run 'HdEEG_Processing_main_GitHub.m' in matlab. 

### Main function: 

HdEEG_Processing_main_GitHub.m (includes the entire pipeline, from preprocessing to spiral detection to sprial-based analysis)

### Subfunctions:
anglesubtract.m <br />
pattDetection_v4.m <br />
opticalFlow2.m <br />

## Expected output <br />

The following outputs can be expected along the pipeline:<br />

EEG signal time series in N2 epochs.<br />
Sigma (11-15Hz) band EEG signal time series in N2 epochs.<br />
Sigma (11-15Hz) band EEG signal time series in spindle epochs.<br />
Sigma (11-15Hz) band EEG signal time series in spindle epochs, interpolated across a 20*21 grid.<br />
Phase velocity field calculated from the sigma (11-15Hz) band EEG time series in spindle epochs.<br />
Spiral centre locations & distribution map calculated from the phase velocity field.<br />
Spiral short-term (overnight) and long-term (3-month) consistencies and their correlations with age and memory retention.<br />
