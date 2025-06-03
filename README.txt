Readme file for softwares, codes, and data used in the manuscript entitled "Increasing constraint of aridity on tree intrinsic water use efficiency"

Software used for data processing: 
Rstudio 2021.09.1 Build 372
R version 4.1.2 (2021-11-01) -- "Bird Hippie"

Platform: aarch64-apple-darwin20 (64-bit)

Programs: 
(1) 1_iWUE_TR_iwue_recalculate.R
Calculate the intrinsic water use efficiency according to 13C. The variable "iWUE_new" is used in the study. "ca_recom" and "d13Catm_recom" represent atmospheric CO2 concentration and stable carbon isotope data, respectively. 
(2) 2_iWUE_TR_iwue_iav_splinefit.R
Detrend the estimated iWUE and get the annual anomalies of iWUE. 
(3) 3_iWUE_TR_data_iav_prepare.R
Connect the anomalies of iWUE with the annual SPEI. Prepare for sensitivity estimations
(4) 4_iWUE_TR_iav_figure.R
Plot the Figture 1 in the manuscript. 
(5) 5_iWUE_TR_iwue_win_response_spei
Estimate the sensitivity of iWUE to SPEI at each period window from 1951 to 2010. 
(6) 6_iWUE_TR_winsens_diff_analysis_spei.R
Estimate the temporal change in the sensitivity between the first and last periods. Also, include the differences of potential influencing factors. 
(7) 7_iWUE_TR_winsens_diff_figure.R
Plot the Figure 2 in the manuscript.
(8) 8_iWUE_TR_winsens_diff_factor_figure.R
Plot the Figure 3 in the manuscript. 
(9) 9_1_iWUE_model_simulate_figure.R
Simulations of the first-principle based model for the sensitivity changes under different scenarios. Plot the Figure 4 in the manuscript. 

Data: 
(1) TR_iwue_newco2_processed_iwue_full_iav_20nyrs_all.csv
The file contains the new estimated iWUE and its annual anomalies at each tree site. 
(2) TR_spei6mon_gs3.csv
6-month scale SPEI averaged over the growing season (i.e, from April to October). 
(3) TR_cru_climate_gs3.csv
Climate data averaged over the growing season. 
(4) TR_iwueiav_v3.csv
Combined anomaly data of both iWUE and SPEI. 
(5) TR_iwue_sensdiff_win22sens_1951_2010.csv
Estimated difference of the sensitivity at each tree between the first and last periods. 
(6) df_sim_sens_historical_meanco2.csv, df_sim_sens_ssp126_meanco2.csv, df_sim_sens_ssp245_meanco2.csv, df_sim_sens_ssp370_meanco2.csv, df_sim_sens_ssp585_meanco2.csv
Multiple files store simulated sensitivities of iWUE under different scenarios as well as input variables. 

* Source_Data: raw data for generating each figure.

