# ECA_Multireactor_Incubations

## Workflow for processing data from optode results through creating Respiration file for publishing

# Do in this order:

1. Pull out saturated DO values for each site from the calibration file. Make this into a study.code_fast_rates.xlsx file

2. Run 01_Optode_Raw_DO script

   a. **Script:** This pulls in raw data from optode results, mapping sheets, saturated 0 minute value, and flags if the sample was put on at the same time a picture was taken by adding a RATE_006 Method Deviation to the 2 minute point

   b. **You:** Add other deviations to the Raw DO .csv after it is exported

4. Run 02a_ECA_MOI_Tins

   a. **Script:** Calculates gravimetric moisture from individual moisture tins taken right before optode incubations. Checks CV of replicates. If CV is high, run 02b_ECA_MOI_outliers.
       Note: 02b_ECA_MOI_outliers was written for ECA_EC samples. You may need to manually change this script to fit your data.
   
   b. **You:** Check for deviations manually (e.g. spills) to make sure that you are calculating correct gravimetric moisture. Add other deviations to the moisture tin .csv after it is exported 
