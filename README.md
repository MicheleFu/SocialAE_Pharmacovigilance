# SocialAE_Pharmacovigilance

Repository for MD thesis on the use of PharmacoVigilance-Pharmacodynamic approach to assess social adverse events.

The repository contains:
- Data preprocessing script
- Data analysis script

N.B. the actual data is not shared to avoid redistribution rights issues.

Materials and preprocessing:
* FAERS ID/: 1 xlsx with all the observations for each Adverse Event considered (42 events [AE_list.Rda]). All the observations were merge, and duplicates were removed (175590 raw ICSR), and then the resulting ICSR_df.Rda was filtered for higher quality data (only ICSR with known reporter type and sex: 151612 ICSR [56898 HCP + 94714 Consumer])
* ATC.csv: csv that links the ATC-code to the drug name used by WHO and to different search terms that can be found in the FAERS. The ATC code has been integrated with drugs reported in the FAERS but not in the ATC code. (3103 Search Terms, 2157 drugs [D_list.Rda])
* pKi_mean.csv: csv with drug and target interacting, pCHEMBL (-pLog of molar IC50, XC50, EC50, AC50, Ki, Kd or Potency; gathered from CHEMBL and Guidetopharmacology.org databases) and action (primarily agonists, antagonists and partial agonists) (9185 raw-pCHEMBL). PDSP was not considered because of the absence of action informations. Where multiple data for the same interaction was found, the geometrical mean was calculated (1887 pCHEMBL: 497 Targets x 20 actions x 2157 drugs)