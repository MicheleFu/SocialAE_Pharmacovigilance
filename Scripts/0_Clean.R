# Introduction ----------------------------------------------------------------
# This script creates the environment used in the study
# Fusaroli Michele, 2019, SDR of Social ADR.

# Libraries -------------------------------------------------------------------
library(tidyverse)
library(readxl)

# Files needed ----------------------------------------------------------------
# ATC code integrated
ATC <- read_delim("ATC.csv", ";", escape_double = FALSE,
                  trim_ws = TRUE)

# fileList --------------------------------------------------------------------
# List of xlsx with observations from FAERS public dashboard
fileList <- list.files(path="FAERS ID",pattern="xlsx",full.names = T)

# AE_list ---------------------------------------------------------------------
# List of AE created from fileList
AE_list <- list()
for (f in fileList){
  s <- gsub("FAERS ID/","",f)
  s <- gsub(".xlsx","",s)
  s <- gsub("_", " ",s)
  AE_list <- c(AE_list, s)
}
AE_list <- tolower(AE_list)
save(AE_list, file = "Rda/AE_list.Rda")

# df_ICSR----------------------------------------------------------------------
# Dataframe containing all the observations obtained from FAERS-pd,
# removing duplicates and low quality data (no reporter or sex specified).
# Then the Suspected Product Active Ingredients column is rewritten
# with ATC code
ICSR_df <- data.frame(matrix(ncol = 24, nrow = 0))
colnames(ICSR_df) <- c("Case ID",
        "Suspect Product Names",
        "Suspect Product Active Ingredients",
        "Reason for Use", "Reactions", "Serious",
        "Outcomes",
        "Sex",
        "Event Date",
        "Latest FDA Received Date",
        "Case Priority",
        "Patient Age",
        "Patient Weight",
        "Sender",
        "Reporter Type",
        "Report Source",
        "Concomitant Product Names",
        "Latest Manufacturer Received Date",
        "Initial FDA Received Date",
        "Country where Event occurred",
        "Reported to Manufacturer?",
        "Manufacturer Control Number",
        "Literature Reference",
        "Compounded Flag")
for (f in fileList) {
  d  <- read_xlsx(f)
  ICSR_df <- rbind(ICSR_df, d)
}

ICSR_df <- ICSR_df %>%
  distinct()
ICSR_df <- ICSR_df %>%
  filter(`Reporter Type` != "Not Specified") %>%
  filter(`Sex` != "Not Specified")
ICSR_df$`Suspect Product Active Ingredients` <- tolower(ICSR_df$`Suspect Product Active Ingredients`)
ICSR_df$`Reactions` <- tolower(ICSR_df$`Reactions`)

ICSR_df <- ICSR_df %>% 
  separate_rows(`Suspect Product Active Ingredients`, sep = ";") %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = "\\\\")
ICSR_df$`Suspect Product Active Ingredients` <- trimws(ICSR_df$`Suspect Product Active Ingredients`)
for (i in seq(nrow(ATC))){
  print(i)
  ICSR_df$`Suspect Product Active Ingredients`[ICSR_df$`Suspect Product Active Ingredients`== ATC$`Search term`[i]] = ATC$Code[i]
}
ICSR_df <- ICSR_df %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = " & ") %>%
  distinct() %>%
  group_by(`Case ID`) %>%
  summarise(`Suspect Product Active Ingredients` = paste(`Suspect Product Active Ingredients`, collapse = ";"),
            `Reactions` = first(`Reactions`),
            `Reason for Use` = first(`Reason for Use`),
            `Reporter Type` = first(`Reporter Type`),
            `Serious` = first(`Serious`),
            `Outcomes` = first(`Outcomes`),
            `Sex` = first(`Sex`),
            `Patient Age`= first(`Patient Age`),
            `Patient Weight`= first(`Patient Weight`),
            `Concomitant Product Names` = first(`Concomitant Product Names`))
save(ICSR_df, file = "RDA/ICSR_df.RDA")

# D_list ----------------------------------------------------------------------
# List of suspect product active ingredients reported in FAERS,
# reported to an ATC code integrated.
D_list <- list(ATC$`Code`) %>%
  flatten() %>%
  unlist %>% 
  unique
D_list <- as.data.frame(D_list)
D_list <- list(D_list[!grepl("&", D_list$D_list),]) %>%
  unique() %>% 
  unlist() %>%
  trimws()
save(D_list, file = "Rda/D_list.Rda")

# HCP_df and PZ_df-------------------------------------------------------------
HCP_df <- ICSR_df %>%
  filter(`Reporter Type` == "Healthcare Professional")
save(HCP_df, file = "RDA/HCP_df.RDA")
PZ_df <- ICSR_df %>%
  filter(`Reporter Type` == "Consumer")
save(PZ_df, file = "RDA/PZ_df.RDA")

# Pharmacodynamics Database ---------------------------------------------------
# CHEMBL
CHEMBL_Mechanisms <- read_delim("Ki Databases/Chembl/Mechanisms.csv", 
                                ";", escape_double = FALSE, trim_ws = TRUE)
Chembl_pKi <- read_delim("Ki Databases/Chembl/pKi.csv", 
                         ";", escape_double = FALSE, trim_ws = TRUE)
Chembl_pKi <- Chembl_pKi %>%
  select(Molecule, pCHEMBL = `pChEMBL Value`,Target = `Target Pref. Name`)
CHEMBL_Mechanisms <- CHEMBL_Mechanisms %>%
  select(Molecule = `Molecule ChEMBL ID`, `Molecule Name`,
         `ATC Classifications`, Target =`Target Name`, Mechanism = `Mechanism Type`) %>%
  mutate(`ATC Classifications` = sub(".*level5': '", "",`ATC Classifications`)) %>%
  mutate(`ATC Classifications` = sub("',.*", "",`ATC Classifications`))
Chembl <- left_join(Chembl_pKi,CHEMBL_Mechanisms, by = c("Molecule","Target"))
Chembl <- Chembl %>%
  filter(is.na(`Molecule Name`)==FALSE) %>%
  filter(is.na(Mechanism)==FALSE)
Chembl$`Molecule Name` <- tolower(Chembl$`Molecule Name`)
Chembl$Mechanism <- tolower(Chembl$Mechanism)
Chembl <- Chembl %>%
  select(Molecule = `Molecule Name`, ATC = `ATC Classifications`, Target, Action = Mechanism, pCHEMBL)
Chembl$Target <- tolower(Chembl$Target)
write_csv2(Chembl, "Ki Databases/Cleaned/CHEMBL.csv")

# PDSP
PDSP <- KiDatabase <- read_delim("Ki Databases/PDSP/KiDatabase.csv", 
                                 ",", escape_double = FALSE, trim_ws = TRUE)
PDSP <- PDSP %>%
  select(Target = Name, Molecule = `Ligand Name`, species, Ki = `ki Val`) %>%
  filter(species == "HUMAN") %>%
  select(-species)
PDSP$Molecule <- tolower(PDSP$Molecule)
PDSP <- PDSP %>%
  mutate(pCHEMBL = -log10(Ki/1000000000))
PDSP$Target <- tolower(PDSP$Target)
write_csv2(PDSP, "Ki Databases/Cleaned/PDSP.csv")

#Guidetopharmacology
GtP <- read_csv("Ki Databases/Guidetopharmacology/interactions.csv")
GtP <- GtP %>%
  select(Target = target, species = `target_species`, Molecule =ligand,
         Action = action, affinity_high, affinity_median, affinity_low) %>%
  filter(species == "Human") %>%
  select(-species) %>% 
  unite(affinity_high, affinity_median, affinity_low, col=pCHEMBL, sep = ";") %>%
  separate_rows(pCHEMBL)
GtP$Molecule <- tolower(GtP$Molecule)
GtP$Target <- tolower(GtP$Target)
GtP$Action <- tolower(GtP$Action)
GtP$ATC <- NA
GtP <- GtP %>%
  select(Molecule, ATC, Target, Action, pCHEMBL)
GtP <- GtP %>%
  filter(pCHEMBL != "NA")
write_csv2(GtP, "Ki Databases/Cleaned/GtP.csv")

# Merge databases in pKi
pKi <- pKi %>%
  mutate(pCHEMBL =pCHEMBL/100)
pKi <- pKi %>%
  filter(Action != "none")
write_csv2(pKi, "pKi.csv")
pKi <- pKi %>%
  group_by(Molecule, Target) %>%
  summarize(ATC = first(ATC), Action = first(Action), pCHEMBL = exp(mean(log(pCHEMBL))))
write_csv2(pKi, "pKi_mean.csv")
