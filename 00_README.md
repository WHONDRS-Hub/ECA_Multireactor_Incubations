# ECA_Multireactor_Incubations

## Workflow for processing data from optode results through creating Respiration file for publishing

# Do in this order:

1. Pull out saturated DO values for each site from the calibration file. Make this into a study.code_fast_rates.xlsx file.\
     Note: This is done internally at PNNL

2. Run 01_Optode_Raw_DO script\
   Note: This is done internally at PNNL

   a. **Script:** This pulls in raw data from optode results, mapping sheets, saturated 0 minute value, and flags if the sample was put on at the same time a picture was taken by adding a RATE_006 Method Deviation to the 2 minute point

   b. **You:** Add other deviations to the Raw DO .csv after it is exported

3. Run 01a_ECA_MOI_Tins\
   Note: This is done internally at PNNL. Script is located in the river-corridors-sfa/QAQC_scripts folder. Alternatively, this data may have already been published and you can pull from ESS-DIVE.
   
   a. **Script:** Calculates gravimetric moisture from individual moisture tins taken right before optode incubations. Checks CV of replicates. If CV is high, run 01b_ECA_MOI_outliers.
       Note: 01b_ECA_MOI_outliers was written for ECA_EC samples. You may need to manually change this script to fit your data.
   
   b. **You:** Check for deviations manually (e.g. spills) to make sure that you are calculating correct gravimetric moisture. Add other deviations to the moisture tin .csv after it is exported

4. Run 03_ECA_drying\
   Note: This is done internally at PNNL. Moisture tin data will need to be located in Data folder on Github.
   
   a. **Script:** Calculates wet sediment mass, dry sediment mass, added water, total water, and gravimetric moisture for each day when samples were weighed over 21-day incubation. Also calculates summary file for initial water mass, final water mass, dry sediment mass, initial gravimetric moisture, and final gravimetric moisture for each sample.

   b. **You:** Add other deviations to the drying .csv after it is exported

5. Run 04_Wet_Sediment_Model\
   Note: This is done internally at PNNL
   
   a. **Script:** Gets amount of water from eca incubation and amends to ECA drying summary file. For ECA, this is an actual model built using grain size data and samples for which we could calculate actual incubation water masses. For EV, we have all weights, so no model needed. Just scroll to the bottom of the script. 

6. Run 05_ECA_SpC_pH_Temp\
   Note: This is done internally at PNNL
   a. **Script:** Gets SpC, pH, Temp values from ECA Mapping files.

   b. **You:** Add other deviations to the .csv after it is exported

7. Run 06a/06b_Sensitivity Analyses to Publish\
   Note: This can be run by public and internally at PNNL
   
   a. **Script:** 06a for ECA. 06b for EV data. Includes descriptions of parameters used in rate fitting, and options used for CV cut offs for rates.

   b. **You:** Check that Method Deviations are correct. 
