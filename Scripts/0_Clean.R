# Introduction ----------------------------------------------------------------
# This script creates the environment used in the study
# Fusaroli Michele, 2019, SDR of Social ADR.

# Libraries -------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(scales)

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
ICSR_df$`Age (YR)` <- abs(ICSR_df$`Age (YR)`)
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
for (i in 1:nrow(ICSR_df)){
  if(isTRUE((ICSR_df$`Age (YR)`[i] <= 0) & (ICSR_df$`Weight (KG)`[i] >= 12))){
    ICSR_df$`Age (YR)`[i] <- NA
  }
}
for (i in 1:nrow(ICSR_df)){
  if(isTRUE((ICSR_df$`Weight (KG)`[i] >= 200))){
    ICSR_df$`Weight (KG)`[i] <- NA
  }
}
for (i in 1:nrow(ICSR_df)){
  if(isTRUE((ICSR_df$`Weight (KG)`[i] <= 20) & (ICSR_df$`Age (YR)`[i]) >= 20)){
    ICSR_df$`Weight (KG)`[i] <- NA
  }
}
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
ICSR_df$`Event Date` <- as.numeric(ICSR_df$`Event Date`)
for (i in 1:nrow(ICSR_df)){
  print(i)
  if(ICSR_df$`Event Date`[i] > 19){
    ICSR_df$`Event Date`[i] <- ICSR_df$`Event Date`[i] + 1900
  } else {
    ICSR_df$`Event Date`[i] <- ICSR_df$`Event Date`[i] + 2000
  }
}
for (i in 1:nrow(ICSR_df)){
  print(i)
  if(ICSR_df$`Event Date`[i] <1970){
    ICSR_df$`Event Date`[i] <- NA
  }
}
save(ICSR_df, file = "RDA/ICSR_df.RDA")

# D_list ----------------------------------------------------------------------
# List of suspect product active ingredients reported in FAERS
load("Rda/ICSR_df.Rda")
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
D_Code_list <- D_list
for (i in 1:length(D_Code_list)){
  D_Code_list[i] <- ATC$Code[ATC$Substance == D_Code_list[i]]
  }
save(D_Code_list, file = "Rda/D_Code_list.Rda")

# HCP_df and C_df-------------------------------------------------------------
# divide the dataframe by reporter type
load("Rda/ICSR_df.Rda")

HCP_df <- ICSR_df %>%
  filter(`Reporter Type` == "Healthcare Professional")
save(HCP_df, file = "RDA/HCP_df.RDA")
C_df <- ICSR_df %>%
  filter(`Reporter Type` == "Consumer")
save(C_df, file = "RDA/C_df.RDA")
NS_df <- ICSR_df %>%
  filter(`Reporter Type` == "Not Specified")
save(NS_df, file = "RDA/NS_df.RDA")

# Population description -------------------------------------------------------
# descriptive analysis of the dataframe
Pop_char_df <- as.data.frame(matrix(ncol = 12, nrow = 1), stringsAsFactors = FALSE)
colnames(Pop_char_df) <- c("Reporter","Tot","F(%)","M(%)","HCP(%)","C(%)",
                           "età media","età mediana", "età DS", "peso media","peso mediana", "peso DS")
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
  Med_A <- round(median(df$`Age (YR)`, na.rm = TRUE),0)
  Age_sd <- round(sd(df$`Age (YR)`, na.rm = TRUE),0)
  Weight <- round(mean(df$`Weight (KG)`, na.rm = TRUE),1)
  Med_W <- round(median(df$`Weight (KG)`, na.rm = TRUE),1)
  Weight_sd <- round(sd(df$`Weight (KG)`, na.rm = TRUE),1)
  dataf <- deparse(substitute(df))
  dataf <- gsub("_.*$","", dataf)
  new_row <- c(dataf,Tot,Fem,Mal,HCP,PZ,Age, Med_A, Age_sd,Weight,Med_W, Weight_sd)
  rbind(Pop_char_df, new_row)
}
Pop_char_df <- Descr(ICSR_df)
Pop_char_df <- Descr(HCP_df)
Pop_char_df <- Descr(C_df)
Pop_char_df <- Descr(NS_df)
Pop_char_df <- Pop_char_df[-1,]
write_csv(Pop_char_df,"Population characteristics/Population_characteristics.csv")

load("Rda/ICSR_df.Rda")
draw_sex <- function(df){
  x <- df
  x <- x %>%
    group_by(Sex) %>%
    count() %>%
    ungroup() %>%
    mutate(per = `n`/sum(`n`)) %>%
    arrange(desc(Sex)) %>%
    mutate(label = scales::percent(per)) %>%
    mutate(r = rank(desc(n)))
  p <- ggplot(data = subset(x, r < 7)) +
    geom_bar(aes(x="", y = per, fill = Sex), stat = "identity", width = 0.5) +
    geom_text(aes(x=1, y = cumsum(per) - per/2, label = label), size = 5, color = "white") +
    ggtitle("Genere dei pazienti") +
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("#E06B6F","#1282A2","#93905B")) +
    labs(y = "Frazione") +
    scale_fill_discrete(name = "Genere", labels = c("Femmine", "Maschi", "Non specificato")
                        )
  return(p)
}
draw_reporter <- function(df){
  x <- df
  x <- x %>%
    group_by(`Reporter Type`) %>%
    count() %>%
    ungroup() %>%
    mutate(per = `n`/sum(`n`)) %>%
    arrange(desc(`Reporter Type`)) %>%
    mutate(label = scales::percent(per)) %>%
    mutate(r = rank(desc(n)))
  p <- ggplot(data = subset(x, r < 7)) +
    geom_bar(aes(x="", y = per, fill = `Reporter Type`), stat = "identity", width = 0.5) +
    geom_text(aes(x=1, y = cumsum(per) - per/2, label = label), size = 5, color = "white") +
    ggtitle("Tipologia dei segnalatori") +
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("#EFC000FF","#0073C2FF","#605E3C")) +
    labs(y = "Frazione") +
    scale_fill_discrete(name = "Segnalatore", labels = c("Consumatore", "Professionista sanitario", "Non specificato")
    )
  return(p)
}
draw_serious <- function(df){
  x <- df
  x <- x %>%
    group_by(Serious) %>%
    count() %>%
    ungroup() %>%
    mutate(per = `n`/sum(`n`)) %>%
    arrange(desc(Serious)) %>%
    mutate(label = scales::percent(per)) %>%
    mutate(r = rank(desc(n)))
  p <- ggplot(data = subset(x, r < 7)) +
    geom_bar(aes(x="", y = per, fill = Serious), stat = "identity", width = 0.5) +
    geom_text(aes(x=1, y = cumsum(per) - per/2, label = label), size = 5, color = "white") +
    ggtitle("Gravità riportata") +
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c("#95C623","#E55812","")) +
    labs(y = "Frazione") +
    scale_fill_discrete(name = "Gravità"#, labels = c("Consumatore", "Professionista sanitario", "Non specificato")
    )
  return(p)
}
draw_date <- function(df){
  x <- df
  x <- x %>%
    group_by(`Event Date`) %>%
    count() %>%
    ungroup() %>%
    mutate(per = `n`/sum(`n`)) %>%
    arrange(desc(`Event Date`)) %>%
    mutate(label = scales::percent(per)) %>%
    mutate(r = rank(desc(n)))
  p <- ggplot(data = x, aes(x= `Event Date`)) +
    geom_histogram(aes(y = n), fill = "blue", stat = "identity", width = 1) +
    ggtitle("Data di occorrenza dell'evento") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Anno di occorrenza", y = "Numero segnalazioni")
  return(p)
}
draw_age <- function(df){
  x <- df
  x <- x %>%
    group_by(`Age (YR)`) %>%
    count() %>%
    ungroup() %>%
    arrange(desc(`Age (YR)`))
  ggplot(data = x, aes(x= `Age (YR)`)) +
    ggtitle("Distribuzione dell'età") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_density(alpha=.2, fill="blue") +
    labs(x = "Età in anni", y = "Densità")
}
draw_weight <- function(df){
  x <- df
  x <- x %>%
    group_by(`Weight (KG)`) %>%
    count() %>%
    ungroup() %>%
    arrange(desc(`Weight (KG)`))
  ggplot(data = x, aes(x= `Weight (KG)`)) +
    ggtitle("Distribuzione del peso in Kg") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_density(alpha=.2, fill="blue") +
    labs(x = "Peso in Kg", y = "Densità")
}

pdf("Population characteristics/Tot/Sex.pdf")
draw_sex(ICSR_df)
dev.off()

pdf("Population characteristics/Tot/Reporter.pdf")
draw_reporter(ICSR_df)
dev.off()

pdf("Population characteristics/Tot/Serious.pdf")
draw_serious(ICSR_df)
dev.off()

pdf("Population characteristics/Tot/Dates.pdf")
draw_date(ICSR_df)
dev.off()

pdf("Population characteristics/Tot/Age.pdf")
draw_age(ICSR_df)
dev.off()

pdf("Population characteristics/Tot/Weight.pdf")
draw_weight(ICSR_df)
dev.off()
