library(tidyverse)
library(superheat)
library(RColorBrewer)

# si caricano i file utilizzati
load("Rda/Matrix_ATC.Rda")
ATC_unique <- read_delim("ATC-unique.csv", 
                         +     ";", escape_double = FALSE, trim_ws = TRUE)

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

# extracting clusters----------------------------------------------------------
t <- MAT %>%
  filter(str_detect(Drug,"N04B")) %>%
  select(Drug, gambling)

# si ritrasforma in una matrice
HM_t <- t %>%
  remove_rownames() %>%
  column_to_rownames(var = "Drug") %>%
  as.matrix()

#Creo la matrice alternativa con le etichette da apporre all'heatmap
HM_t1 <- t%>%
  remove_rownames() %>%
  column_to_rownames(var = "Drug") %>%
  as.matrix()
HM_t1[is.na(HM_t1)] <- ""

# Si prende il solo ROR senza IC e si converte gli Inf e tutto ciò che sta
# sopra a 100 a 100 (solo per determinare il colore)
HM_t[,] <- sub(".*?\\|", "", HM_t[,])
HM_t[,] <- sub("\\|.*$", "", HM_t[,])
HM_t[HM_t[,] == "Inf"] <- 100
HM_t <- as.data.frame(HM_t)
HM_t <- type.convert(HM_t)
HM_t$gambling[HM_t$gambling > 100] <- 100
HM_t <- as.matrix(HM_t)
png("HM_1.png", height = 15000, width = 8000) #if desired 5lvl remember to turn active left.label
superheat(HM_t,
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
          X.text = HM_t1,
          X.text.size = 3,
          left.label.text.alignment = "right",
          pretty.order.rows = FALSE,
          left.label = "variable",
          row.title = "Substance",
          row.title.size = 6,
          column.title = "Adverse Event",
          column.title.size = 6)
dev.off()

#HCP_PZ------------------------------------------------------------------------
x <- HCP_PZ
y <- HCP_PZ
x[x[,] == "Inf"] <- 100
x[x[,] > 100] <- 100
x[x[,] < (-100)] <- (-100)
x[is.nan(x)] <- -100
y[is.nan(y)] <- "-100"
Code_member <- substr(rownames(x), start = 1, stop = 4)


png("HCP_PZ.png", height = 20000, width = 10000) #if desired 5lvl remember to turn active left.label
superheat(x,
          heat.pal = c("brown", "red", "white", "green", "blue"),
          heat.pal.values = c(0,0.4, 0.5, 0.6,1),
#          heat.col.scheme = "red",
          heat.lim = c(-100,100),
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
          X.text = y,
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
