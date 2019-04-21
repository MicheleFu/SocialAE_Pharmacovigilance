# Introduction ----------------------------------------------------------------
# Questo codice esegue esegue analisi di disproporzionalità con ROR + IC
# e crea una matrice con per righe i farmaci e per colonne gli EA
# e in ogni cella, se significativo (i.e. IC-min > 1), una lista
# contenente (IC_min,ROR,IC_max), F_EA, F_nEA
# Fusaroli Michele, 2019, SDR of Social ADR.

# Libraries Needed ------------------------------------------------------------

library(tidyverse)
library(superheat)
library(RColorBrewer)

# Rda Needed ------------------------------------------------------------------
load("Rda/ICSR_df.Rda")
load("Rda/AE_list.Rda")
load("Rda/D_list.Rda")
ATC <- read_delim("ATC.csv", ";", escape_double = FALSE, 
                  trim_ws = TRUE)
ATC <- select(ATC, -"Search term") %>%
  unique()
ATC <- ATC[!grepl("&", ATC$Code),]
D_list <- as.data.frame(D_list)
D_list <- list(D_list[!grepl("&", D_list$D_list),]) %>%
  unique() %>% 
  unlist() %>%
  trimws()

# Wrangle_df-------------------------------------------------------------------
# Df with all the data needed already summarized
Wrangle_df <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(Wrangle_df) <- c("Drug_Code", "Drug_Name", "AE", "F_EA", "F_nEA", "ROR", "s", "ROR_m", "ROR_M", "IC")
Wrangle <- function(df, D) {
  #
  # Args:
  #   df  : ICSR_df
  #   D   : D_list
  #
  # Returns:
  #   Wrangle_df : Disproportionality Analysis and summarized data
  for (i in 1:length(AE_list)) {
    print(i)
    for (d in D) {
      print(d)
      tab <- table(str_detect(df$`Suspect Product Active Ingredients`, d),
                   str_detect(df$Reactions, AE_list[i]))
      colnames(tab)[1] <- "nEA" # other reactions
      colnames(tab)[2] <- "EA"  # reaction e
      rownames(tab)[1] <- "nF"  # other drugs
      rownames(tab)[2] <- "F"   # drug d
      t <- as.data.frame.matrix(tab)
      F_EA <- t["F","EA"]
      F_nEA <- t["F","nEA"]
      ROR_m <- NA
      ROR_M <- NA
      ROR <- NA
      IC <- NA
      s <- NA
      if (F_EA>=3) {
        nF_nEA <- t["nF","nEA"]
        nF_EA <- t["nF","EA"]
        ROR <- F_EA * nF_nEA / nF_EA / F_nEA
        s <- sqrt(1/F_EA + 1/F_nEA + 1/nF_EA + 1/nF_nEA)
        ROR_m <- exp(log(ROR) - 1.96*s)
        ROR_M <- exp(log(ROR) + 1.96*s)
        if (is.infinite(ROR)) {IC <- "[|Inf|]"} else {
          ROR_m <- round(ROR_m, digits = 1)
          ROR <- round(ROR, digits = 1)
          ROR_M <- round(ROR_M, digits = 1)
          IC <- paste("[", ROR_m,
                      "|", ROR, "|",
                      ROR_M, "]",
                      sep = "")}
      }
      new_row <- list(d, ATC$Substance[ATC$Code == d], AE_list[i], F_EA, F_nEA, ROR, s, ROR_m, ROR_M, IC)
      Wrangle_df[nrow(Wrangle_df)+1,] <-  new_row
    }
  }
  return(Wrangle_df)
}
Wrangle_df <- Wrangle(ICSR_df,D_list)
save(Wrangle_df, file = "Rda/Wrangle_df.Rda")

# Heatmatrix ------------------------------------------------------------------
load("Rda/Wrangle_df.Rda")
Heatmatrix <- matrix(ncol = length(AE_list), nrow = length(D_list))
rownames(Heatmatrix) <- D_list
colnames(Heatmatrix) <- AE_list
for (e in AE_list) {
  for (d in D_list){
    Heatmatrix[d,e] = Wrangle_df$ROR[Wrangle_df$Drug_Code == d, Wrangle_df$AE == e]
  }
}
save(Heatmatrix, file = "Rda/Heatmatrix.Rda")

## remove not significant columns and rows?
onlyNAcolumns_idx <- Heatmatrix %>%
  is.na() %>%
  apply(MARGIN = 2, FUN = all)
Heatmatrix <- Heatmatrix[,!onlyNAcolumns_idx]
Heatmatrix <- Heatmatrix[rowSums(is.na(Heatmatrix)) != ncol(Heatmatrix), ]

#Si sostituisce, nella matrice, al search term "Code: Substance"
MAT <- Significance_Matrix
MAT <- MAT %>%
  as.data.frame()
MAT <- setNames(cbind(rownames(MAT), MAT, row.names = NULL), c("Code",colnames(MAT)))
MAT <- MAT %>%
  left_join(ATC_unique, by = "Code") %>%
  select(Code, "Substance", everything()) %>%
  unite(Drug, c(Code, "Substance"), sep = ": ")
MAT <- unique(MAT)
# Si crea il codice ATC di 4 caratteri che viene usato nella heatmap per clusterizzare
Code_member <- substr(MAT$Drug, start = 1, stop = 4)

# si ritrasforma in una matrice
Heatmatrix <- MAT %>%
  remove_rownames() %>%
  column_to_rownames(var = "Drug") %>%
  as.matrix()

#Creo la matrice alternativa con le etichette da apporre all'heatmap
MAT1 <- MAT %>%
  remove_rownames() %>%
  column_to_rownames(var = "Drug") %>%
  as.matrix()
MAT1[is.na(MAT1)] <- ""

# Heatmatrix[!is.na(Heatmatrix)] <- 1
# Heatmatrix[is.na(Heatmatrix)] <- 0
# Si prende il solo ROR senza IC e si converte gli Inf e tutto ciò che sta
# sopra a 100 a 100 (solo per determinare il colore)
Heatmatrix[,] <- sub(".*?\\|", "", Heatmatrix[,])
Heatmatrix[,] <- sub("\\|.*$", "", Heatmatrix[,])
Heatmatrix[Heatmatrix[,] == "Inf"] <- 100
Heatmatrix <- as.data.frame(Heatmatrix)
Heatmatrix <- type.convert(Heatmatrix)
Heatmatrix[Heatmatrix[,] > 100] <- 100
Heatmatrix <- as.matrix(Heatmatrix)

png("Color_ATC.png", height = 15000, width = 8000) #if desired 5lvl remember to turn active left.label
superheat(Heatmatrix,
          heat.pal = c("white", "red", "#b35806","#542788"),
          heat.pal.values = c(0, 0.1, 0.5, 1),
          heat.col.scheme = "red",
          heat.lim = c(1,100),
          bottom.label.text.angle = 90,
          bottom.label.text.size = 10,
          bottom.label.size = 0.1,
          force.left.label = TRUE,
          left.label.text.size = 4,
          left.label.size = 0.4,
          force.grid.hline = TRUE,
          grid.hline.col = "gray",
          grid.vline.col = "gray",
          heat.na.col= "white",
          X.text = MAT1,
          X.text.size = 3,
          left.label.text.alignment = "right",
          pretty.order.rows = FALSE,
          membership.rows = Code_member,
          left.label = "variable",
          row.title = "Substance",
          row.title.size = 6,
          column.title = "Adverse Event",
          column.title.size = 6)
dev.off()
