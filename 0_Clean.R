# Introduction ----------------------------------------------------------------
# This script creates the environment used in the study
# Fusaroli Michele, 2019, SDR of Social ADR.

# Libraries -------------------------------------------------------------------
library(tidyverse)
library(readxl)

# fileList ----------------------------------------------------------------------
# List of xlsx with observations from FAERS public dashboard
fileList <- list.files(path="FAERS ID",pattern="xlsx",full.names = T)

# AE_list ----------------------------------------------------------------------
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

# df_ICSR---------------------------------------------------------------------------
# Dataframe containing all the observations obtained from FAERS-pd,
# removing duplicates
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
ICSR_df <- distinct(ICSR_df)
ICSR_df$`Suspect Product Active Ingredients` <- tolower(ICSR_df$`Suspect Product Active Ingredients`)
ICSR_df$`Reactions` <- tolower(ICSR_df$`Reactions`)
save(ICSR_df, file = "Rda/ICSR_df.Rda")

# D_list ----------------------------------------------------------------------
# List of suspect product active ingredients reported in FAERS, needed to
# confront with problems like different spelling or name,
# or absence from ATC
ATC <- read_delim("ATC.csv", ";", escape_double = FALSE,
                  trim_ws = TRUE)
D_df <- list(ICSR_df$`Suspect Product Active Ingredients`) %>%
  unlist() %>%
  strsplit(";") %>%
  unlist() %>%
  trimws() %>%
  unique() %>%
  strsplit("\\\\") %>% 
  unlist() %>% 
  trimws() %>%
  unique() %>% 
  as.data.frame() %>%
  rename("Search term" = ".") %>%
  left_join(ATC, by = "Search term") %>% 
  arrange(Code)
D_df$Substance <- tolower(D_df$Substance)
### Manual integration of missing ATCs in the file ATC.csv
D_list <- list(ATC$`Code`) %>%
  flatten() %>%
  unlist %>% 
  unique() %>%
  save(D_list, file = "Rda/D_list.Rda")

# Tidy ICSR_df-----------------------------------------------------------------
# Suspect Product Active Ingredients referred to with different names
# are converted to unique ATCs, and the df loose unused data
ATC <- read_delim("ATC.csv", ";",
                  escape_double = FALSE, trim_ws = TRUE)
load("Rda/ICSR_df.Rda")

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

# HCP_df and PZ_df-------------------------------------------------------------
load("Rda/ICSR_df.Rda")
HCP_df <- ICSR_df %>%
  filter(`Reporter Type` == "Healthcare Professional")
save(HCP_df, file = "RDA/HCP_df.RDA")
PZ_df <- ICSR_df %>%
  filter(`Reporter Type` == "Consumer")
save(PZ_df, file = "RDA/PZ_df.RDA")
