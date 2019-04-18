# Introduzione ----------------------------------------------------------------
# Questo codice contiene script utilizzati per creare df e liste utili alle DA
# Fusaroli Michele, 2019, Progetto SDR di ADR sociali.

# Librerie necessarie ---------------------------------------------------------
library(tidyverse)
library(readxl)

# fileList ----------------------------------------------------------------------
fileList <- list.files(path="FAERS ID",pattern="xlsx",full.names = T)

# EAlist ----------------------------------------------------------------------
# si crea la lista dei nomi degli EA a partire dai nomi dei file
EAlist <- list()
for (f in fileList){
  s <- gsub("FAERS ID/","",f)
  s <- gsub(".xlsx","",s)
  s <- gsub("_", " ",s)
  EAlist <- c(EAlist, s)
}
# si toglie gambling disorder (duplicato)
EAlist[[13]] <- NULL
# si salva
EAlist <- tolower(EAlist)
save(EAlist, file = "Rda/EA.Rda")

# df---------------------------------------------------------------------------
# si crea una matrice vuota con 24 colonne e senza osservazioni,
# e si impostano i nomi delle variabili rispettando la struttura del FAERS
df <- data.frame(matrix(ncol = 24, nrow = 0))
x  <- c("Case ID",
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
colnames(df) <- x
# si riempie la matrice con gli ICSR dei df EA
for (f in fileList) {
  d  <- read_xlsx(f)
  df <- rbind(df, d)
}
# si rimuovono i duplicati (gli ICSR possono essere presenti in più df)
df <- distinct(df)
# si trasforma tutta la colonna dei farmaci in minuscolo in modo da confrontarla
# poi con la lista ATC
df$`Suspect Product Active Ingredients` <- tolower(df$`Suspect Product Active Ingredients`)
# si salva
save(df, file = "Rda/df_unico.Rda")

# F_list ----------------------------------------------------------------------
# crea la lista che viene poi confrontata con l'ATC per affrontare problemi di
# ortografia, definizione molteplice, e mancata comparsa nel sistema atc
# si crea una lista con tutti i valori della variabile farmaci sospettati
F_list <- list(df$`Suspect Product Active Ingredients`)
# si sistema la lista in modo che tutti i farmaci siano separati,
# che siano sullo stesso livello
# e che non vi siano duplicati
F_list <- F_list %>%
  # si appiattisce la lista (i farmaci sospettati insieme sono in cluster: "a;b""c"
  unlist() %>%
  # si separano i farmaci segnalati insieme (a;b)
  strsplit(";") %>%
  # si appiattisce nuovamente ("a""b")
  unlist() %>%
  # si rimuovono i duplicati
  unique() %>%
  # si separano le combinazioni di farmaci (a1\a2)
  strsplit("\\\\") %>%
  # si appiattisce nuovamente ("a1""a2")
  unlist() %>%
  # si rimuovono i duplicati
  unique() %>% {
    # si rimuove la punteggiatura finale
    gsub("[[:punct:]]*$", "", .) %>%
      # si rimuove tutto ciò che sta dopo frontslash
      gsub("/.*$", "", .) %>%
      # si rimuovono slash e parentesi
      gsub("\\(.*$", "", .) %>%
      # si rimuove lo spazio finale
      gsub(" *$", "", .) %>%
      # si rimuove la punteggiatura finale
      gsub("[[:punct:]]*$", "", .)
  } %>%
  # si rimuovono i duplicati
  unique(F_list)
# Si ritrasforma in una lista
Flist <- list()
for (f in F_list) {
  Flist <- c(Flist, f)
}
# si salva
save(Flist, file = "Rda/F.Rda")

#Crea lista di farmaci effettivamente usata in df unico
ATC_used <- read_delim("ATC-used.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE)
Flist <- list(ATC_used$`Search term`) %>%
  flatten()
save(Flist, file = "Rda/F.Rda")

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
