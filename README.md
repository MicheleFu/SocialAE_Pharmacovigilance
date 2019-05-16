# SocialAE_Pharmacovigilance

Repository for MD thesis on the use of PharmacoVigilance-Pharmacodynamic approach to assess social adverse events.

The repository contains:
- Data preprocessing script
- Data analysis script

N.B. the actual data is not shared to avoid redistribution rights issues.

Materials and preprocessing:
* *FAERS ID/*: 1 xlsx with all the observations for each Adverse Event considered (42 events [AE_list.Rda]). All the observations were merged, and duplicates were removed (175590 raw ICSR), and then the resulting ICSR_df.Rda was filtered for higher quality data (only ICSR with known reporter type and sex: 151612 ICSR [56898 HCP + 94714 Consumer])
* *ATC.csv*: csv that links the ATC-code to the drug name used by WHO and to different search terms that can be found in the FAERS. The ATC code has been integrated with drugs reported in the FAERS but not in the ATC code. (3103 Search Terms, 2157 drugs [D_list.Rda])
* *pKi_mean.csv*: csv with drug and target interacting, human pCHEMBL (-pLog of molar IC50, XC50, EC50, AC50, Ki, Kd or Potency; gathered from CHEMBL and Guidetopharmacology.org databases) and action (primarily agonists, antagonists and partial agonists) (9185 raw-pCHEMBL). PDSP was not considered because of the absence of action informations. Where multiple data for the same interaction was found, the geometrical mean was calculated (1887 pCHEMBL: 497 Targets x 20 actions x 2157 drugs)

Method:
* For each combination (drug x event) a contingency table was built, and when the reports with both the drug and the event were at least three the Rerporting Odds Ratio was calculated ((drug_AE x no_drug_no_AE)/(no-drug_AE x drug_no_AE)) together with the confidence interval (s <- sqrt(1/F_EA + 1/F_nEA + 1/nF_EA + 1/nF_nEA);  ROR_m <- exp(log(ROR) - 1.96 x s); ROR_M <- exp(log(ROR) + 1.96 x s)). The results were reported in the dataframe *ROR_df.Rda*.
* With the data from ROR_df.Rda a matrix was printed, with drugs as rows (gathered in drug families by ATC code) and events as column, and with the specific ROR in every cell. Every cell was then coloured using a gradient from white (ROR = 0 or not available), to red (ROR = 10), to brown(ROR = 50), to violet (ROR = 100 or above). The heatmap thus built has been saved as *HM_(ROR_df).png*
* The passages above were conducted also for HCP and Consumer (PZ) reports separated.
* For each combination of event, target and action a linear regression model was produced, if possible, for the entire database, and only lrm with p-value less then 0.1 were reported (*LRM.pdf*). The same process was conducted on some of the individual drug families (antiParkinson [LRMA_N04.pdf], antiPsychotic [LRMA_N05A.pdf]). Slope, intercept and p-value were reported inside the *LRM_csv* folder. The *LRM_2plots.pdf* shows side by side the linear model with all the drugs and the linear models of individual drug families with at least two points in the plot.