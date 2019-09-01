# Introduction ----------------------------------------------------------------
# This script creates the environment used in the study
# Fusaroli Michele, 2019, SDR of Social ADR.

# Libraries -------------------------------------------------------------------
library(tidyverse)
library(readxl)

# Files needed ----------------------------------------------------------------
# ATC code integrated and adapted
ATC <- read_delim("Databases Used/ATC.csv", ";", escape_double = FALSE,
                  trim_ws = TRUE)

# fileList --------------------------------------------------------------------
# List of xlsx with observations from FAERS public dashboard
fileList <- list.files(path="FAERS ID_AE",pattern="xlsx",full.names = T)

# AE_list.Rda ---------------------------------------------------------------------
# List of Adverse Events created from fileList
AE_list <- list()
for (f in fileList){
  s <- gsub("FAERS ID_AE/","",f)
  s <- gsub(".xlsx","",s)
  s <- gsub("_", " ",s)
  AE_list <- c(AE_list, s)
}
AE_list <- tolower(AE_list)
save(AE_list, file = "Rda/AE_list.Rda")

# ICSR_df.Rda----------------------------------------------------------------------
# Dataframe containing all the observations obtained from FAERS-pd,
# removing duplicates.
# Then the "Suspected Product Active Ingredients" column is rewritten
# with standardized ATC code
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
ICSR_df$`Suspect Product Active Ingredients` <- tolower(ICSR_df$`Suspect Product Active Ingredients`)
ICSR_df$`Reactions` <- tolower(ICSR_df$`Reactions`)

ICSR_df <- ICSR_df %>% 
  separate_rows(`Suspect Product Active Ingredients`, sep = ";") %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = "\\\\")
ICSR_df$`Suspect Product Active Ingredients` <- trimws(ICSR_df$`Suspect Product Active Ingredients`)
for (i in seq(nrow(ATC))){
  print(i)
  ICSR_df$`Suspect Product Active Ingredients Code`[ICSR_df$`Suspect Product Active Ingredients`== ATC$`Search term`[i]] = ATC$Code[i]
  ICSR_df$`Suspect Product Active Ingredients`[ICSR_df$`Suspect Product Active Ingredients`== ATC$`Search term`[i]] = ATC$Substance[i]
  }
ICSR_df <- ICSR_df %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = " & ") %>%
  distinct() %>%
  group_by(`Case ID`) %>%
  summarise(`Suspect Product Active Ingredients` = paste(`Suspect Product Active Ingredients`, collapse = ";"),
            `Suspect Product Active Ingredients Code` = paste(`Suspect Product Active Ingredients Code`, collapse = ";"),
            `Suspect Product Names` = first(`Suspect Product Names`),
            `Reason for Use` = first(`Reason for Use`),
            `Reactions` = first(`Reactions`),
            `Reporter Type` = first(`Reporter Type`),
            `Serious` = first(`Serious`),
            `Outcomes` = first(`Outcomes`),
            `Sex` = first(`Sex`),
            `Patient Age`= first(`Patient Age`),
            `Patient Weight`= first(`Patient Weight`),
            `Concomitant Product Names` = first(`Concomitant Product Names`),
            `Event Date` = first(`Event Date`),
            `Initial FDA Event Date` = first(`Initial FDA Received Date`),
            `Country where Event occurred` = first(`Country where Event occurred`),
            `Literature Reference` = first(`Literature Reference`))
ICSR_df$`Patient Age` <-  na_if(ICSR_df$`Patient Age`, "Not Specified")
ICSR_df <- ICSR_df %>%
  separate(`Patient Age`, into = c("Age (YR)","Age_Unit"), sep = " ", convert = TRUE)
for (i in 1:nrow(ICSR_df)){
  if (ICSR_df$Age_Unit[i] == "DEC" & is.na(ICSR_df$`Age (YR)`[i]) == FALSE){
    ICSR_df$`Age (YR)`[i] <- round(ICSR_df$`Age (YR)`[i]*10,1)
  }
  else if (ICSR_df$Age_Unit[i] == "MTH" & is.na(ICSR_df$`Age (YR)`[i]) == FALSE){
    ICSR_df$`Age (YR)`[i] <- round(ICSR_df$`Age (YR)`[i]/12,1)
  }
  else if (ICSR_df$Age_Unit[i] == "DAY" & is.na(ICSR_df$`Age (YR)`[i]) == FALSE){
    ICSR_df$`Age (YR)`[i] <- round(ICSR_df$`Age (YR)`[i]/365,1)
  }
  else if (ICSR_df$Age_Unit[i] == "HR" & is.na(ICSR_df$`Age (YR)`[i]) == FALSE){
    ICSR_df$`Age (YR)`[i] <- round(ICSR_df$`Age (YR)`[i]/8760,1)
  }
  else if (ICSR_df$Age_Unit[i] == "WEEK" & is.na(ICSR_df$`Age (YR)`[i]) == FALSE){
    ICSR_df$`Age (YR)`[i] <- round((ICSR_df$`Age (YR)`[i]/365)*7,1)
  }
}
ICSR_df$`Patient Weight` <-  na_if(ICSR_df$`Patient Weight`, "Not Specified")
ICSR_df <- ICSR_df %>%
  separate(`Patient Weight`, into = c("Weight (KG)","Weight_Unit"), sep = " ", convert = TRUE)
ICSR_df$`Weight (KG)` <- parse_number(ICSR_df$`Weight (KG)`)
for (i in 1:nrow(ICSR_df)){
  print(i)
  if (ICSR_df$Weight_Unit[i] == "LB" & is.na(ICSR_df$`Weight (KG)`[i]) == FALSE & is.na(ICSR_df$Weight_Unit[i]) == FALSE){
    ICSR_df$`Weight (KG)`[i] <- round(ICSR_df$`Weight (KG)`[i]*0.4536,1)
  }
  else if (ICSR_df$Weight_Unit[i] == "GMS" & is.na(ICSR_df$`Weight (KG)`[i]) == FALSE & is.na(ICSR_df$Weight_Unit[i]) == FALSE){
    ICSR_df$`Weight (KG)`[i] <- round(ICSR_df$`Weight (KG)`[i]/1000,1)
  }
}
ICSR_df <- ICSR_df %>%
  select(-Weight_Unit, -Age_Unit)
ICSR_df$`Weight (KG)` <- round(ICSR_df$`Weight (KG)`,1)
# when not specified the date of the event,
# the first report to FDA was considered as the event date
for (i in 1:nrow(ICSR_df)){
  if (ICSR_df$`Event Date`[i] == "-"){
    ICSR_df$`Event Date`[i] <- ICSR_df$`Initial FDA Event Date`[i]
  }
  x <- ICSR_df$`Event Date`[i]
  ICSR_df$`Event Date`[i] <- substr(x, nchar(x)-2+1, nchar(x))
  }
ICSR_df <- ICSR_df %>%
  select(-`Initial FDA Event Date`)
save(ICSR_df, file = "RDA/ICSR_df.RDA")

# D_list ----------------------------------------------------------------------
# List of suspect product active ingredients reported in FAERS
ICSR_df <- ICSR_df %>% 
  separate_rows(`Suspect Product Active Ingredients`, sep = ";") %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = "\\\\")
ICSR_df$`Suspect Product Active Ingredients` <- trimws(ICSR_df$`Suspect Product Active Ingredients`)
D_list <- list(unique(ICSR_df$`Suspect Product Active Ingredients`)) %>%
  flatten() %>%
  unlist %>% 
  unique %>%
  sort()
D_list <- D_list[D_list != "altro"]
save(D_list, file = "Rda/D_list.Rda")

# HCP_df and C_df-------------------------------------------------------------
# divide the dataframe by reporter type
load("Rda/ICSR_df.Rda")

HCP_df <- ICSR_df %>%
  filter(`Reporter Type` == "Healthcare Professional")
save(HCP_df, file = "RDA/HCP_df.RDA")
C_df <- ICSR_df %>%
  filter(`Reporter Type` == "Consumer")
save(PZ_df, file = "RDA/C_df.RDA")
NS_df <- ICSR_df %>%
  filter(`Reporter Type` == "Not Specified")
save(NS_df, file = "RDA/NS_df.RDA")

#Population description -------------------------------------------------------
Pop_char_df <- as.data.frame(matrix(ncol = 8, nrow = 1), stringsAsFactors = FALSE)
colnames(Pop_char_df) <- c("df","Tot","F(%)","M(%)","HCP(%)","PZ(%)",
                           "mean age","mean weight")
Descr <- function(df){
  Tot <- nrow(df)
  F_num <-sum(df$Sex == "Female")
  Fem <- paste(round((F_num/Tot)*100,2))
  M_num <- sum(df$Sex == "Male")
  Mal <- paste(round((M_num/Tot)*100,2))
  HCP_num <- sum(df$`Reporter Type` == "Healthcare Professional")
  HCP <- paste(round((HCP_num/Tot)*100,2))
  PZ_num <- sum(df$`Reporter Type` == "Consumer")
  PZ <- paste(round((PZ_num/Tot)*100,2))
  Age <- round(mean(df$`Age (YR)`, na.rm = TRUE),0)
  Weight <- round(mean(df$`Weight (KG)`, na.rm = TRUE),1)
  dataf <- deparse(substitute(df))
  new_row <- c(dataf,Tot,Fem,Mal,HCP,PZ,Age,Weight)
  rbind(Pop_char_df, new_row)
  }
Pop_char_df <- Descr(ICSR_df)
Pop_char_df <- Descr(HCP_df)
Pop_char_df <- Descr(PZ_df)
Pop_char_df <- Descr(NS_df)
write_csv2(Pop_char_df,"Population characteristics/Population_characteristics.csv")





