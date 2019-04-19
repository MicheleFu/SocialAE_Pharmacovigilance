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
# Manual integration of missing ATCs
D_list <- list(D_df$`Search term`) %>%
  flatten()
save(D_list, file = "Rda/D_list.Rda")

#Crea lista di farmaci effettivamente usata in df_ATC
ATC_used <- read_delim("ATC-used.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE)
Flist <- list(ATC_used$`Code`) %>%
  flatten() %>%
  unique()
save(Flist, file = "Rda/F_ATC.Rda")

#Creazione df-sep -------------------------------------------------------------
ATC_used <- read_delim("ATC-used.csv", ";",
                       escape_double = FALSE, trim_ws = TRUE)
ATC_used$Substance <- tolower(ATC_used$Substance)

x <- df %>% 
  separate_rows(`Suspect Product Active Ingredients`, sep = ";") %>% 
  separate_rows(`Reactions`, sep = ";")

for (i in seq(nrow(ATC_used))){
  print(i)
  x$Substance[str_detect(x$`Suspect Product Active Ingredients`,ATC_used$`Search term`[i])]=ATC_used$Substance[i]
}

# elimino search term e duplicati
x <- x %>% select(-`Suspect Product Active Ingredients`) %>% distinct()

EAlist <- tolower(EAlist)
x$Reactions <- tolower(x$Reactions)
x$Reactions <- trimws(x$Reactions)
t <- x
for (i in seq(length(EAlist))){
  print(i)
  x$AE[str_detect(x$`Reactions`, EAlist[[i]])] = EAlist[[i]]
}
x <- x %>% select(-`Reactions`) %>% distinct()
df_sep <- x

save(df_sep, file = "RDA/df_separated.RDA")

# Eventualmente riguardare gli altri NA di Substance per altri eventuali farmaci per cui impostare ATC
## unique(x$`Suspect Product Active Ingredients`[is.na(x$Substance)])
## Se son misspellings possibile anche gsub all'inizio dello script

#Creazione df_ATC--------------------------------------------------------------
ATC_used <- read_delim("ATC-used.csv", ";",
                       escape_double = FALSE, trim_ws = TRUE)
ATC_used$Substance <- tolower(ATC_used$Substance)

x <- df %>% 
  separate_rows(`Suspect Product Active Ingredients`, sep = ";")

for (i in seq(nrow(ATC_used))){
  print(i)
  x$`Suspect Product Active Ingredients`[str_detect(x$`Suspect Product Active Ingredients`,ATC_used$`Search term`[i])]=ATC_used$Code[i]
}
x <- x %>% distinct()
t <- x %>%
  group_by(`Case ID`) %>%
  summarise(`Suspect Product Active Ingredients` = paste(`Suspect Product Active Ingredients`, collapse = ";"),
            `Reactions` = first(`Reactions`),
            `Reason for Use` = first(`Reason for Use`),
            `Reporter Type` = first(`Reporter Type`),
            `Serious` = first(`Serious`),
            `Outcomes` = first(`Outcomes`),
            `Sex` = first(`Sex`))
x <- t
x <- x %>% separate_rows(`Reactions`, sep = ";")
EAlist <- tolower(EAlist)
x$Reactions <- tolower(x$Reactions)
x$Reactions <- trimws(x$Reactions)

for (i in seq(length(EAlist))){
  print(i)
  x$AE[str_detect(x$`Reactions`, EAlist[[i]])] = EAlist[[i]]
}
x <- x %>% select(-`Reactions`) %>% distinct()
x <- x %>%
  subset(!is.na(`AE`))

t <- x %>%
  group_by(`Case ID`) %>%
  summarise(`Suspect Product Active Ingredients` = first(`Suspect Product Active Ingredients`),
            `AE` = paste(`AE`, collapse = ";"),
            `Reason for Use` = first(`Reason for Use`),
            `Reporter Type` = first(`Reporter Type`),
            `Serious` = first(`Serious`),
            `Outcomes` = first(`Outcomes`),
            `Sex` = first(`Sex`))
df_ATC <- t
save(df_ATC, file = "RDA/df_ATC.RDA")

# df_ATC_HCP/PZ-------------------------------------------------------------------
df_ATC_HCP <- df_ATC %>%
  filter(`Reporter Type` == "Healthcare Professional")
df_ATC_PZ <- df_ATC %>%
  filter(`Reporter Type` == "Consumer")
df_ATC_NS <- df_ATC %>%
  filter(`Reporter Type` == "Not Specified")
save(df_ATC_HCP, file = "RDA/df_ATC_HCP.RDA")
save(df_ATC_PZ, file = "RDA/df_ATC_PZ.RDA")
