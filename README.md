# hippocampus-prefrontal-social-information
This software package contains code for population decoding, synergy decoding and pseudo-population decoding of single-unit responses from freely behaving mice exposed to various social stimuli. Example data files are provided to demonstrate the execution of each function. All scripts were written as MATLAB .m files and tested in MATLAB R2023a under Windows 10 v22H2.

1. "mPFC_vHPC_decoding.m"
This script evaluates decoding performances using neurons from the medial prefrontal cortex (mPFC) and the ventral hippocampus (vHPC). Population decoding was conducted within each session using binned population vectors from the during stimulus period. A naive Bayesian decoder was trained to classify firing patterns associated with four distinct social stimuli. The effects of bin length and population dimensions on classification accuracy were compared between mPFC and vHPC.

By defalt, the script displays pre-computed multi-session accuracy results for mPFC and vHPC, loaded from '\data\decode_results_Eular_dist_4_class_all.mat'. To retrain and evaluate the model using the provided example firing data, set the execution condition statement in "Line 36" to "true".

Example data includes single-unit firing from two sessions each for mPFC and vHPC located in
"\data\mPFC\produce_single_cell_data_all_mPFC"
"\data\vHPC\produce_single_cell_data_all_vHPC"

It should be noted that, the pre-computed file will be overwritten after re-sunning the decoding.

2. "synergy_decoding.m"
This script assesses synergy-based decoding performance. Synergy decoding was applied per session using balanced sets of eight neurons. A naive Bayesian decoder was trained to discriminate among four social stimuli. Decoding accuracy was compared between using both spike+LFP datasets versus spike-only datasets in mPFC and vHPC.

Example data includes all sessions containing at least eight neurons for mPFC and vHPC, stored in:
"\data\synergy\mPFC\"
"\data\synergy\vHPC\"

3. "pseudo_decoding.m"
This script evaluates decoding performance for discriminateing between two social stimulus identities (S and N) across the time course of social investigation. Pseudo-population decoding was implemented via single-trial population decoding, constructing pseudo-populations from neurons or LFPs across multiple recording sessions. A naive Bayesian decoder was trained to differentiate between two social stimuli. Decoding accuracy was compared between mPFC and vHPC.

Example data includes five sessions each for mPFC and vHPC located in:
"\data\pseudo\mPFC\"
"\data\pseudo\vHPC\"

All file paths was specified relative to the root directory of the project. Ensure MATLAB's current working directory is set to the project root when executing these scripts.

