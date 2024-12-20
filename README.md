# NCD_Global_2024
Files related to the paper 'Particle-associated N2 fixation by heterotrophic bacteria in the global ocean'

# Description of files
Python codes related to the paper 'Particle-associated N2 fixation by heterotrophic bacteria in the global ocean'. The provided code will help to (1) generate the dynamics of N2 fixation by non-cyanobacterial diazotrophs (NCDs) inside a particle of radius 0.25 cm at the temperature 17°C (Supplementary Figure S1) and (2) calculate cellular rates  at a distance of 0.05 cm from the center of a particle (Fig. 3). One can also find the cell specific N2 fixation rates (mug N/cell/d) and total amount of fixed N2 inside the particle (mug N/particle) by looking at the files 'N_fix_dat' and the output file 'J_N2_tot', respectively.

# How to run the code
Running the code requires two scripts in the 'NCD_Python_Code' folder, the 'main_script.py' and the 'parameters.py'. The first script is for decribing and running the model, while the second script contains all the parameter values.

Running the code also need three '.csv' data files: 'Size_vs_halfsaturation_NO3_uptake', 'Size_vs_NO3_uptake', and 'viscocity_temperature_Jumars_1993'. These '.csv' files are also included in the folder 'NCD_Python_Code'.
