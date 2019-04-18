# Introduzione ----------------------------------------------------------------
# Questo codice esegue esegue analisi di disproporzionalità con ROR + IC
# e crea una matrice con per righe i farmaci e per colonne gli EA
# e in ogni cella, se significativo (i.e. IC-min > 1), una lista
# contenente (IC_min,ROR,IC_max), F_EA, F_nEA
# Fusaroli Michele, 2019, Progetto SDR di ADR sociali.

# Libraries Needed ------------------------------------------------------------

library(tidyverse)

# Rda Needed ------------------------------------------------------------------
load("Rda/df_ATC.Rda")
load("Rda/EA.Rda")
load("Rda/F_ATC.Rda")
ATC_used <- read_delim("ATC-used.csv", 
                       ";", escape_double = FALSE, trim_ws = TRUE)

# Disproportionality Analysis (tabelle EA)-------------------------------------

DA <- function(dataframe, farmaci) {
  # Esegue le analisi di disproporzionalità e riporta i ROR significativi
  # in una matrice.
  #
  # Args:
  #   dataframe : dataframe unico contenente tutti i farmaci e tutti gli EA
  #   farmaci   : lista dei farmaci da ricercare
  #
  # Returns:
  #   Significance_Matrix : una matrice con per righe i farmaci e per colonne gli EA
  #                         e in ogni cella, se significativo (i.e. IC-min > 1), una lista
  #                         contenente (IC_min,ROR,IC_max), F_EA, F_nEA
  # Si crea la matrice vuota
  Significance_Matrix <- matrix(ncol = length(EAlist), nrow = length(farmaci))
  rownames(Significance_Matrix) <- farmaci
  colnames(Significance_Matrix) <- EAlist
  for (e in EAlist) {
    # si procede solo se l'EA è riportato almeno una volta
    df_EA <- filter(dataframe, str_detect(`AE`, e))
    v_temp <- nrow(df_EA)
    if (v_temp > 0) {
      print(e)
      # si procede per ogni F, solo se segnalato almeno una volta per l'EA
      for (f in farmaci) {
        df_EA_F <- filter(df_EA, str_detect(`Suspect Product Active Ingredients`, f))
        v_temp <- nrow(df_EA_F)
        if (v_temp>0) {
          # si crea una tabella di frequenza F*EA
          mytable <- table(str_detect(dataframe$`Suspect Product Active Ingredients`, f),
                           str_detect(dataframe$AE, e))
          colnames(mytable)[1] <- "nEA" #ICSR in cui non è segnalato l'EA
          colnames(mytable)[2] <- "EA" #ICSR in cui è segnalato l'EA
          rownames(mytable)[1] <- "nF" #ICSR in cui non è segnalato il F
          rownames(mytable)[2] <- "F" #ICSR in cui è segnalato il F
          d <- as.data.frame.matrix(mytable)
          # si procede all'analisi di disproporzionalità se
          # vi sono almeno 3 casi F_EA segnalati
          if (d["F","EA"]>=3) {
            r <- Calculate_IC(d)
            # se ROR infinito si riporta senza IC
            if (is.infinite(r[[2]])) {
              IC <- "[|Inf|]"
              Significance_Matrix[[f,e]]   <-   IC
            } else if (r[[1]] > 1) {
              # se ROR non è infinito e sta con tutto il suo IC sopra 1 si riporta col suo IC
              # si approssimano i valori a 1 cifra decimale
              r[[1]] <- round(r[[1]], digits = 1)
              r[[2]] <- round(r[[2]], digits = 1)
              r[[3]] <- round(r[[3]], digits = 1)
              # si crea l'etichetta con tali valori
              IC <- paste("[", r[[1]],
                          "|", r[[2]], "|",
                          r[[3]], "]",
                          sep = "")
              # si riporta nella matrice
              Significance_Matrix[[f,e]]   <-   IC
            }
          }
        }
      }
    }
  }
  save(Significance_Matrix, file = "Rda/Matrix_ATC.Rda")
}


Calculate_IC <- function(m) {
  #
  # Calcola i ROR con i rispettivi IC
  #
  # Args:
  #   d : matrice di contingenza di un singolo F con un singolo EA
  #
  # Returns:
  #   r  : lista che contiene IC_min (estremo inf), ROR e IC_max (estremo sup)
  #
  # si impostano le combinazioni di F e EA
  F_EA <- m["F","EA"]
  nF_nEA <- m["nF","nEA"]
  nF_EA <- m["nF","EA"]
  F_nEA <- m["F","nEA"]
  # si calcola il ROR
  ROR <- F_EA * nF_nEA / nF_EA / F_nEA
  # si calcola la deviazione standard
  s <- sqrt(1/F_EA + 1/F_nEA + 1/nF_EA + 1/nF_nEA)
  # si calcola l'estremo inferiore dell'intervallo di confidenza
  IC_min <- exp(log(ROR) - 1.96*s)
  # si calcola anche l'estremo superiore
  IC_max <- exp(log(ROR) + 1.96*s)
  # si crea la lista che risulta dalla funzione
  r <- list("IC_min" <- IC_min, "ROR" <- ROR, "IC_max" <- IC_max)
}

# Execute ---------------------------------------------------------------------
DA(df_ATC,Flist)

# Clean------------------------------------------------------------------------
load("Rda/Matrix_ATC.Rda")
# si rimuovono le colonne senza ROR significativi
onlyNAcolumns_idx <- Significance_Matrix %>%
  is.na() %>%
  apply(MARGIN = 2, FUN = all)
Significance_Matrix <- Significance_Matrix[,!onlyNAcolumns_idx]
# si rimuovono le righe senza ROR significativi
Significance_Matrix <- Significance_Matrix[rowSums(is.na(Significance_Matrix)) != ncol(Significance_Matrix), ]
save(Significance_Matrix, file = "Rda/Matrix_ATC.Rda")
