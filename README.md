This code is specific to data collected for the VIPA Study at the University of East Anglia, led by Dr Jordan Tsigarides (MD).
Written by Vanessa Grove (MSc)

Purpose of code: To run spectral analysis on resting state EEG data and compare data during/following use of Virtual Reality technology to a baseline dataset. Main analysis occurs using data averaged across the whole scalp in seven frequency bands. If the power in a frequency band is found to be significantly different to the baseline dataset, a supplementary analysis will run which will analyse each of the 64 channels individually to examine which channels/groups of channels are demonstrating differences. The results of this code can be visualised in the visualisation code, which is in a separate file. 
Purpose of code: To compare sprectral EEG data between groups. Just as in the static spectral analysis, the main analysis occurs using data averaged across the whole scalp in seven frequency bands. If the power in a frequency band is found to be significantly different between the groups, a supplementary analysis will run which will analyse each of the 64 channels individually to examine which channels/groups of channels are demonstrating differences. The results of this code will be visualised in a figure that will be produced when code is run completely. 

IMPORTANT- This code can only be run on a device that is connected to the S: drive on a UEA-based server.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Setup: 
To use this code on new data, user will update the first few lines of code with the appropriate information. This must be done manually before running the code.

IMPORTANT- This code will run on fft_data variables that have been created through running the staticspec.m code NOT on EEG files. Be sure that the data you would like to run has an associated fft_data.mat structure.

Code Structure:
Specify Type of Analysis
* analysis_type- will either be '1' for a UEA analysis or '2' for UoL
* visualisation_only- if you would like to rerun statistics with different parameters, make sure this is set to '0'. If you would just like to produce figures of previously run statistics, set this variable to '1' (default)

SECTION 1a- Data Parameters- UEA Analysis Only
* fig_title- Give analysis a title to be displayed on figures (e.g. 'FMS v HC baseline'). Be sure to use ''
* dataset1- ID number of first participant or appended dataset. IMPORTANT: this variable must match the processed EEG data file of interest exactly (i.e. if EEG file is called APPENDED_EC_pre, then 'ppt_id' MUST be 'APPENDED').
* dataset2- as above.
* existing_var- If only looking to produce figures, this is the title of the saved gp_compare variable (see Variables file for options)
* Exp- Number of experiments being analysed (either 2 or 3). IMPORTANT: HC data is located under 3
* stats_title- name of statistical analysis (e.g. 'FMS:HC baseline stats', similar to fig_title but will be stored in the gp_compare variable_
* Comparison- Type of comparison of interest (1 = sequential order (default), 2 = condition (i.e. by game or by environment).
Other information in this section regarding statistical analysis details will get stored in the output variable for future reference. Only update this section if rerunning statistical analysis with new parameters. 

 SECTION 1b- Data Setup- UEA Analysis Only
Loads in appropriate spectral data from the fft_data files entered into the first section of code:
* 'gp_compare.powvals' = spectral power in dB for each channel at 1 Hz resolution from 2-45 Hz
* 'gp_compare.scalpavg' = average spectral power in each frequency band across the scalp in microvolts
* 'gp_compare.scalppow' = average scalp spectral power at 1 Hz resolution from 2-45 Hz

SECTION 2- Data Parameters- UoL Anlaysis Only
This section will load existing UoL gp_compare variable with UEA and UoL data stored for analysis and/or visualisation

SECTION 3- Statistics Setup
This section of code prepares the output variable ('gp_compare') to store the results of the statistical analysis.

SECTION 4- Whole-Scalp Analysis
* ‘For’ loop cycles through each of the seven frequency bands of interest
* Data is selected from fft_data(1).allpower and fft_data(n).allpower for each of the frequencies included within the frequency range (output: 64 x m matrices, where m = number of frequencies in specified band at a 1Hz resolution).
* Data is selected from 'gp_compare' for each of the frequencies included within the frequency range (output: 64 x m matrices, where m = number of frequencies in specified band at a 1Hz resolution).
* The average across channels is computed for each of the datasets (output: 1 x m matrices)
* The ‘true’ test statistic is calculated via a Wilcoxon Signed Rank test. The p-value of this test is stored as the test statistic and is what is used to build the data-based distribution for permutation testing.
* The data is permuted by shuffling the data points randomly between the two datasets 5000 times and performing a Wilcoxon Signed Rank test for each permute. The output of this step is a 1 x 5000 matrix in which each entry corresponds to a p-value for the shuffled data.
* A final p-value for the permutation testing is calculated by summing the number of times permuted test statistic values were less than the true test statistic, divided by the total number of permutes: (∑(permuted statistic < test statistic)/N). This value is stored in decimal form in 'fft_data.scalppvals'. Additionally, if the p-value is less than the Bonferroni-corrected critical value, this information is stored as a ‘1’ in 'fft_data.scalpsig' and if the value is greater than the critical p-value, it is stored as a ‘0’.
* Optional visualization: The null hypotheses distribution can be visualized as histograms with the true test statistic plotted on top.
* The ‘true’ test statistic is calculated via a Mann Whitney U/Rank Sum test. The p-value of this test is stored as the test statistic and is what is used to build the data-based distribution for permutation testing.
* The data is permuted by shuffling the data points randomly between the two datasets 5000 times and performing a Mann Whitney U/Rank Sum test for each permute. The output of this step is a 1 x 5000 matrix in which each entry corresponds to a p-value for the shuffled data.
* A final p-value for the permutation testing is calculated by summing the number of times permuted test statistic values were less than the true test statistic, divided by the total number of permutes: (∑(permuted statistic < test statistic)/N). This value is stored in decimal form in 'gp_compare.scalppvals'. Additionally, if the p-value is less than the Bonferroni-corrected critical value, this information is stored as a ‘1’ in 'gp_compare.scalpsig' and if the value is greater than the critical p-value, it is stored as a ‘0’.

SECTION 5- Single Channel Analysis
Note: This section is only performed if the result of the whole-scalp data is significant.
*	The relevant value from 'fft_data.scalpsig' is evaluated. If the value is ‘0’, this section is skipped.
*	The same procedure as above is repeated, however, this time it is computed for each of the 64 channels individually, therefore the ‘true’ test statistic output is a 1 x 64 matrix of p-values corresponding to each channel. The result of permutation testing is a 5000 x 64 matrix. In order to calculate the final p-value for each of the channels and to account for multiple comparisons, the minimum value of each of the permutes is selected (across channels), resulting in a final distribution of 5000 p-values. The true test statistic of each of the 64 channels is then compared to this single distribution using the formula above. The values are stored in 'fft_data.chan_pvals', and the logical values representing significance are stored in 'fft_data.sigdif'. 

*	The relevant value from 'gp_compare.scalpsig' is evaluated. If the value is ‘0’, this section is skipped.
*	The same procedure as above is repeated, however, this time it is computed for each of the 64 channels individually, therefore the ‘true’ test statistic output is a 1 x 64 matrix of p-values corresponding to each channel. The result of permutation testing is a 5000 x 64 matrix. In order to calculate the final p-value for each of the channels and to account for multiple comparisons, the minimum value of each of the permutes is selected (across channels), resulting in a final distribution of 5000 p-values. The true test statistic of each of the 64 channels is then compared to this single distribution using the formula above. The values are stored in 'gp_compare.chan_pvals', and the logical values representing significance are stored in 'gp_compare.sigdif'.

SECTION 6- Visualisation
Produces a figure to visualise results.
Top Panel: Shows a scalp map of the average power difference distribution between groups (dataset2 - dataset1) in each frequency band. Significant channels will be displayed on these scalp maps if any were found. 
Middle Panel: Shows a bar graph of average power in each frequency band between groups. If a significant difference was demonstrated, it will be displayed on this panel.
Bottom Panel: Shows the full power spectrum for 2-45 Hz in dB for each of the groups. The frequency bands are marked with vertical lines. 
