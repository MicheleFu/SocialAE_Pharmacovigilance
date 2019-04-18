# Introduzione ----------------------------------------------------------------
# Questo codice esegue esegue analisi di disproporzionalità con ROR + IC
# e crea il dataframe usato poi per ogni visualizzazione, con tutt i dati su
# per ogni F e EA (F_EA, F_nEA, ROR, IC; tutto per HCP; tutto per PZ; differenza).
# Fusaroli Michele, 2019, Progetto SDR di ADR sociali.

# Libraries Needed ------------------------------------------------------------

library(tidyverse)

# Rda Needed ------------------------------------------------------------------
load("~/Desktop/PhV_of_Social_Impact/Rda/df_ATC.Rda")
load("~/Desktop/PhV_of_Social_Impact/Rda/EA.Rda")
load("~/Desktop/PhV_of_Social_Impact/Rda/F_ATC.Rda")

# Disproportionality Analysis (tabelle EA)--------------------------------------------------

Wrangle <- function(dataframe, farmaci) {
  # Crea un database con tutti i dati che servono.
  #
  # Args:
  #   dataframe : dataframe unico con tutte le osservazioni
  #   farmaci   : lista dei farmaci da ricercare
  #
  # Returns:
  #   dataframe con i dati elaborati
  #
  Wrangle <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_EA_sign) <- c("Drug", "AE", "F_EA", "F_nEA", "ROR", "ROR_m", "ROR_M", "IC", "")
  for (e in EAlist) {
    # si procede solo se l'EA è riportato almeno una volta
    df_EA <- filter(dataframe, str_detect(`Reactions`, e))
    v_temp <- nrow(df_EA)
    if (v_temp > 0) {
      # si crea una matrice senza osservazioni per riportare i risultati
      df_EA_sign <- data.frame(matrix(ncol = 4, nrow = 0))
      colnames(df_EA_sign) <- c("Drug", "IC", "F_EA", "F_nEA")
      # si procede per ogni F, solo se segnalato almeno una volta per l'EA
      for (f in farmaci) {
        df_EA_F <- filter(df_EA, str_detect(`Suspect Product Active Ingredients`, f))
        v_temp <- nrow(df_EA_F)
        if (v_temp>0) {
          # si crea una tabella di frequenza F*EA
          mytable <- table(str_detect(dataframe$`Suspect Product Active Ingredients`, f),
                           str_detect(dataframe$Reactions, e))
          colnames(mytable)[1] <- "nEA" #ICSR in cui non è segnalato l'EA
          colnames(mytable)[2] <- "EA" #ICSR in cui è segnalato l'EA
          rownames(mytable)[1] <- "nF" #ICSR in cui non è segnalato il F
          rownames(mytable)[2] <- "F" #ICSR in cui è segnalato il F
          d <- as.data.frame.matrix(mytable)
          # si procede all'analisi di disproporzionalità se
          # vi sono almeno 3 casi F_EA segnalati
          if (d["F","EA"]>=3) {
            r <- Calculate_IC(d)
            # se ROR infinito si riporta con IC infinito
            if (is.infinite(r[[2]])) {
              IC <- "[Inf]"
              df_EA_sign [nrow(df_EA_sign)+1,] <- c(f, IC, d["F","EA"], d["F","nEA"])
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
              # si riporta nella tabella
              df_EA_sign [nrow(df_EA_sign)+1,] <- c(f, IC, d["F","EA"], d["F","nEA"])
            }
          }
        }
      }
      # se la tabella ha osservazioni si stampa, con F in ordine per numero casi
      if (nrow(df_EA_sign) > 0) {
        # si crea l'etichetta per salvare la tabella (METTERE IN ORDINE PER IC_min)
        l <- paste("Tables/", e,".csv", sep ="")
        # si salva la tabella
        write.csv2(df_EA_sign, l)
      }
    }
  }
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
DA(df,Flist)
