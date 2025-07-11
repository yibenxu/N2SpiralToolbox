# N2SpiralToolbox
MATLAB toolbox to automatically detect and analyse spiral wave patterns in sleep EEG data, developed by Prof. Pulin Gong's group at University of Sydney.
## Instructions for use
System dependencies: N/A <br />
Software version: MATLAB 2016b and above (has been tested on MATLAB 2016b and 2020a) <br />
Hardware requirement: N/A

Data format: Sleep EEG data (256-channel HdEEG) in .fdt and .set formats (EEGLAB) <br />
Data tested: 9 patients with obstructive sleep apnea (OSA), recorded during stage 2 NREM (N2) sleep <br />

### Launch: <br />
Download the HdEEG data (not publicly available, but can be accessed upon request) and save to the main folder. <br />
Download EEGLAB2024 (https://sccn.ucsd.edu/eeglab/). <br />
Run 'HdEEG_Processing_main_GitHub.m' in matlab. <br />

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


## Authors

* **Yiben Xu** - yixu4976@uni.sydney.edu.au
* **Alex McInnes** - amci2773@uni.sydney.edu.au
* **Pulin Gong** - pulin.gong@sydney.edu.au
  
## Cite
Xu, Y., McInnes, A., Kao, CH. et al. Spatiotemporal dynamics of sleep spindles form spiral waves that predict overnight memory consolidation and age-related memory decline. Commun Biol 8, 1014 (2025). https://doi.org/10.1038/s42003-025-08447-4
