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

# Data Needed ------------------------------------------------------------------
load("Rda/ICSR_df.Rda")
load("Rda/AE_list.Rda")
load("Rda/HCP_df.Rda")
load("Rda/PZ_df.Rda")
load("Rda/D_list.Rda")
D_list <- as.data.frame(D_list)
D_list <- list(D_list[!grepl("&", D_list$D_list),]) %>%
  unique() %>% 
  unlist() %>%
  trimws()

ATC <- read_delim("ATC.csv", ";", escape_double = FALSE, 
                  trim_ws = TRUE) %>%
  select(-"Search term") %>%
  unique()
ATC <- ATC[!grepl("&", ATC$Code),]

# Functions -------------------------------------------------------------------

Wrangle <- function(df, D) {
  #
  # Args:
  #   df  : Dataframe of observations
  #   D   : List of drugs
  #
  # Returns:
  #   Wrangle_df : Disproportionality Analysis and summarized data
  Wrangle_df <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(Wrangle_df) <- c("Drug_Code", "Drug_Name", "AE", "F_EA", "F_nEA", "nF_EA","nF_nEA", "ROR", "s", "ROR_m", "ROR_M", "IC")
  for (i in 1:length(AE_list)) {
    EA_Name <- AE_list[i]
    print(EA_Name)
    x <- subset(df, str_detect(Reactions, AE_list[i]))
    for (d in D) {
      print(d)
      D_Name <- ATC$Substance[ATC$Code == d]
      y <- subset(df,!str_detect(`Suspect Product Active Ingredients`,d))
      if (sum(str_detect(x$`Suspect Product Active Ingredients`, d))==0) {
        F_nEA <- sum(str_detect(df$`Suspect Product Active Ingredients`, d))
        nF_EA <- sum(str_detect(y$Reactions,EA_Name))
        nF_nEA <- sum(!str_detect(y$Reactions,EA_Name))
        new_row <- list(d, D_Name, EA_Name, 0, F_nEA, nF_EA, nF_nEA, NA, NA, NA, NA, NA)
      } else {
        tab <- table(str_detect(df$`Suspect Product Active Ingredients`, d),
                     str_detect(df$Reactions, AE_list[i]))
        colnames(tab)[1] <- "nEA" # other reactions
        colnames(tab)[2] <- "EA"  # reaction e
        rownames(tab)[1] <- "nF"  # other drugs
        rownames(tab)[2] <- "F"   # drug d
        t <- as.data.frame.matrix(tab)
        F_EA <- t["F","EA"]
        F_nEA <- t["F","nEA"]
        nF_EA <- t["nF","EA"]
        nF_nEA <- t["nF","nEA"]
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
                        sep = "")
            }
          }
        new_row <- list(d, D_Name, EA_Name, F_EA, F_nEA, nF_EA, nF_nEA, ROR, s, ROR_m, ROR_M, IC)
      }
      Wrangle_df[nrow(Wrangle_df)+1,] <-  new_row
    }
  }
  return(Wrangle_df)
}

# Wrangled_dfs-------------------------------------------------------------------
# Df with all the data needed already summarized
TOT_df <- Wrangle(ICSR_df,D_list)
save(TOT_df, file = "Rda/Tot_df.Rda")

HCP_Wrangle_df <- Wrangle(HCP_df,D_list) %>%
  mutate(Reporter_Type = "HCP")
PZ_Wrangle_df <- Wrangle(PZ_df,D_list) %>%
  mutate(Reporter_Type = "PZ")
HCP_PZ_Wrangled_df <- rbind(HCP_Wrangle_df,PZ_Wrangle_df) %>%
  arrange(Drug_Code, AE)
save(HCP_PZ_Wrangled_df, file = "Rda/HCP_PZ_Wrangled_df.Rda")
## HCP_PZ_Wrangled <- HCP_PZ_Wrangled_df %>%
##  filter (!is.na(ROR))

# Print_Heatmap ------------------------------------------------------------------

Create_Matrix <- function(df, index) {
  m <- matrix(ncol = length(AE_list), nrow = length(D_list))
  rownames(m) <- D_list
  colnames(m) <- AE_list
  for (e in AE_list) {
    x <- subset(df, AE == e)
    for (d in D_list){
      print(d)
      if (!is.na(x$ROR_m[x$Drug_Code == d])) {
        if(x$ROR_m[x$Drug_Code == d] > 1) {
          m[d,e] = x[[index]][x$Drug_Code == d]
        }
      }
    }
  }
  onlyNAcolumns_idx <- m %>%
    is.na() %>%
    apply(MARGIN = 2, FUN = all)
  m <- m[,!onlyNAcolumns_idx]
  m <- m[rowSums(is.na(m)) != ncol(m), ]
  m <- m %>%
    as.data.frame()
  m <- setNames(cbind(rownames(m), m, row.names = NULL), c("Code",colnames(m)))
  m <- m %>%
    left_join(ATC, by = "Code") %>%
    select(Code, "Substance", everything()) %>%
    unite(Drug, c(Code, "Substance"), sep = ": ")
  return(m)
}

Print_Heatmap <- function(df) {
  Heatmatrix <- Create_Matrix(df,"ROR")
  Code_member <- substr(Heatmatrix$Drug, start = 1, stop = 4)
  Heatmatrix <- Heatmatrix %>%
    remove_rownames() %>%
    column_to_rownames(var = "Drug") %>%
    as.matrix()
  IC_matrix <- Create_Matrix(df, "IC")
  IC_matrix <- IC_matrix %>%
    remove_rownames() %>%
    column_to_rownames(var = "Drug") %>%
    as.matrix()
  IC_matrix[is.na(IC_matrix)] <- ""
  Heatmatrix[Heatmatrix[,] == "Inf"] <- 100
  Heatmatrix[Heatmatrix[,] > 100] <- 100
  l <- paste("Hm_",deparse(substitute((df))),".png", sep = "")
  png(l, height = 15000, width = 8000)
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
            X.text = IC_matrix,
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
}

Print_Heatmap(Wrangle_df)
HCP <- filter(HCP_PZ_Wrangled_df, HCP_PZ_Wrangled_df$Reporter_Type == "HCP")
Print_Heatmap(HCP)
PZ <- filter(HCP_PZ_Wrangled_df, HCP_PZ_Wrangled_df$Reporter_Type == "PZ")
Print_Heatmap(PZ)

# Comparation HCP/PZ ----------------------------------------------------------
Comparation_df <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(Comparation_df) <- c("Drug_Code", "Drug_Name", "AE", "F_EA_HCP",
                              "F_nEA_HCP", "nF_EA_HCP","nF_nEA_HCP",
                              "F_EA_PZ","F_nEA_PZ", "nF_EA_PZ", "nF_nEA_PZ",
                              "Log(odds)_EA_HCP(nF)", "Log(odds)_EA_PZ/HCP(nF)",
                              "Log(Odds)_F_HCP","Comparation_Index", "Pr")
for (e in AE_list) {
  print(e)
  for (d in D_list) {
    print(d)
    x <- subset(HCP_PZ_Wrangled_df, Drug_Code == d & AE == e)
    y <- matrix(ncol=4, nrow=4)
    colnames(y)  <-  c("EA","nEA", "Reporter_Type", "Drug")
    rownames(y) <- c("F_HCP", "nF_HCP", "F_PZ", "nF_PZ")
    y[1,1] <- x$F_EA[x$Reporter_Type == "HCP"]
    y[1,2] <- x$F_nEA[x$Reporter_Type == "HCP"]
    y[1,3] <- "HCP"
    y[1,4] <- "1"
    y[2,1] <- x$nF_EA[x$Reporter_Type == "HCP"]
    y[2,2] <- x$nF_nEA[x$Reporter_Type == "HCP"]
    y[2,3] <- "HCP"
    y[2,4] <- "0"
    y[3,1] <- x$F_EA[x$Reporter_Type == "PZ"]
    y[3,2] <- x$F_nEA[x$Reporter_Type == "PZ"]
    y[3,3] <- "PZ"
    y[3,4] <- "1"
    y[4,1] <- x$nF_EA[x$Reporter_Type == "PZ"]
    y[4,2] <- x$nF_nEA[x$Reporter_Type == "PZ"]
    y[4,3] <- "PZ"
    y[4,4] <- "0"
    y <- as.data.frame(y)
    m <- glm(cbind(EA,nEA) ~ 1 + Reporter_Type*Drug, family=binomial, data=y)
    y <- as.matrix(y)
    new_row <- list(d, ATC$Substance[ATC$Code == d], e, y[1,1],y[1,2], y[2,1],y[2,2],
                    y[3,1],y[3,2], y[4,1], y[4,2], m[[1]][[1]], m[[1]][[2]],
                    m[[1]][[3]], m[[1]][[4]], coef(summary(m))[4,4])
    Comparation_df[nrow(Comparation_df)+1,] <-  new_row
  }
}
save(Comparation_df, file="Rda/Comparation_df.Rda")

##Intercetta = valore in HCP (Reporter=0) e nF (F=0)
##Reporter_TypePZ = ruolo del Reporter Type (F=0)
##F1 = ruolo del farmaco (Reporter=0)
##Reporter_TypePZ : F = ruolo combinato di F e Reporter

mP <- read_delim("mP.csv", ";", escape_double = FALSE, 
                  trim_ws = TRUE)
png(filename = "Parkinson_Gambling")
ggplot(data=mP) +
  geom_boxplot(mapping=aes(x=reorder(Drug_Name, -ROR), middle = ROR, lower = ROR_m, ymin=ROR_m, upper = ROR_M, ymax=ROR_M), stat = "identity") +
  labs(title    = "Gambling RORs of anti-Parkinson Drugs") +
  xlab("Drug") +
  geom_label(aes(x = reorder(`Drug_Name`, -ROR), ROR, label = ROR), size = 3, colour = "red", fill="white") +
  coord_flip()
dev.off()
