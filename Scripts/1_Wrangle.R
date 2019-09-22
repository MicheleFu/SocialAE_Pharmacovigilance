# Introduction ----------------------------------------------------------------
# This code execute a DA using ROR with CI 95%.
# It prints a heatmap with the results.
# Fusaroli Michele, 2019, SDR of Social ADR.

# Libraries Needed ------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(superheat)
library(RColorBrewer)

# Data Needed ------------------------------------------------------------------
load("Rda/ICSR_df.Rda")
load("Rda/AE_list.Rda")
load("Rda/HCP_df.Rda")
load("Rda/C_df.Rda")
load("Rda/D_list.Rda")
load("Rda/D_Code_list.Rda")

ATC <- read_delim("Databases Used/ATC.csv", ";", escape_double = FALSE, 
                  trim_ws = TRUE) %>%
  select(-"Search term") %>%
  unique()
ATC <- ATC[!grepl("&", ATC$Code),]

# Calculate ROR ---------------------------------------------------------------

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
      D_Code <- ATC$Code[ATC$Substance == d]
      y <- subset(df,!str_detect(`Suspect Product Active Ingredients Code`,D_Code))
      if (sum(str_detect(x$`Suspect Product Active Ingredients Code`, D_Code))==0) {
        F_nEA <- sum(str_detect(df$`Suspect Product Active Ingredients Code`, D_Code))
        nF_EA <- sum(str_detect(y$Reactions,EA_Name))
        nF_nEA <- sum(!str_detect(y$Reactions,EA_Name))
        new_row <- list(D_Code, d, EA_Name, 0, F_nEA, nF_EA, nF_nEA, NA, NA, NA, NA, NA)
      } else {
        tab <- table(str_detect(df$`Suspect Product Active Ingredients Code`, D_Code),
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
        new_row <- list(D_Code, d, EA_Name, F_EA, F_nEA, nF_EA, nF_nEA, ROR, s, ROR_m, ROR_M, IC)
      }
      Wrangle_df[nrow(Wrangle_df)+1,] <-  new_row
    }
  }
  return(Wrangle_df)
}

ROR_df <- Wrangle(ICSR_df,D_list)
save(ROR_df, file = "Rda/ROR_df.Rda")

# Heatmap ---------------------------------------------------------------------
D_Code_list <- D_Code_list %>%
  sort()
Create_Matrix <- function(df, index) {
  m <- matrix(ncol = length(AE_list), nrow = length(D_Code_list))
  rownames(m) <- D_Code_list
  colnames(m) <- AE_list
  for (e in AE_list) {
    x <- subset(df, AE == e)
    for (d in D_Code_list){
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
  l <- paste("Visualization/Heatmap_",deparse(substitute((df))),".png", sep = "")
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

Print_Heatmap(ROR_df)

# Specific Heatmaps -----------------------------------------------------------

N04_code <- subset(D_Code_list, str_detect(D_Code_list,"N04B"))
AE_ICD_list <- c("compulsive shopping","gambling","gambling disorder","hypersexuality","impulse-control disorder", "impulsive behaviour")
Parkinson <- ROR_df %>%
  filter(Drug_Code %in% N04_code) %>%
  filter(AE %in% AE_ICD_list)

opioids_code <- subset(D_Code_list, str_detect(D_Code_list,"N02A"))
Diversion <- c("drug diversion")
Opioids <- ROR_df %>%
  filter(Drug_Code %in% opioids_code) %>%
  filter(AE %in% Diversion)

bisfosfonati_code <- subset(D_Code_list, str_detect(D_Code_list,"M05B"))
Diversion <- c("poor personal hygiene")
Bisfosfonati <- ROR_df %>%
  filter(Drug_Code %in% bisfosfonati_code) %>%
  filter(AE %in% Diversion)

lassativi_code <- subset(D_Code_list, str_detect(D_Code_list,"A06A"))
Diversion <- c("social avoidant behaviour")
Lassativi <- ROR_df %>%
  filter(Drug_Code %in% lassativi_code) %>%
  filter(AE %in% Diversion)

psicostimolanti_code <- subset(D_Code_list, str_detect(D_Code_list,"N06B"))
Diversion <- c("educational problems","impulsive behaviour","fight in school","drug diversion","crime")
ADHD <- ROR_df %>%
  filter(Drug_Code %in% psicostimolanti_code) %>%
  filter(AE %in% Diversion)


draw_tile <- function(df,h,w){
  l <- paste("Visualization/", (deparse(substitute(df))), ".png", sep="")
  png(filename = l,height = h, width = w)
  x <- ggplot(data = df, aes(x = AE, y = Drug_Name)) +
    geom_tile(aes(fill = ROR)) +
    geom_text(aes(label = IC)) +
    scale_fill_gradient(low = "yellow", high = "red") +
    labs(x = "Evento Avverso", y = "Farmaco")
  print(x)
  dev.off()
}
draw_tile(Parkinson,800,1200)
draw_tile(Opioids,800,400)
draw_tile(Bisfosfonati,800,400)
draw_tile(Lassativi,800,400)
draw_tile(ADHD,800,1200)

# HCP-PZ ROR ------------------------------------------------------------------
HCP_ROR_df <- Wrangle(HCP_df,D_list) %>%
  mutate(Reporter_Type = "HCP")
save(HCP_ROR_df,file="RDA/HCP_ROR_df")
C_ROR_df <- Wrangle(C_df,D_list) %>%
  mutate(Reporter_Type = "C")
save(C_ROR_df,file="RDA/C_ROR_df")
HCP_C_ROR_df <- rbind(HCP_ROR_df,C_ROR_df) %>%
  arrange(Drug_Code, AE)
save(HCP_C_ROR_df, file = "Rda/HCP_C_ROR_df.Rda")

HCP <- filter(HCP_C_ROR_df, HCP_C_ROR_df$Reporter_Type == "HCP")
Print_Heatmap(HCP)
C <- filter(HCP_C_ROR_df, HCP_C_ROR_df$Reporter_Type == "C")
Print_Heatmap(C)
Comparation_df <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(Comparation_df) <- c("Drug_Code", "Drug_Name", "AE", "F_EA_HCP",
                              "F_nEA_HCP", "nF_EA_HCP","nF_nEA_HCP",
                              "F_EA_C","F_nEA_C", "nF_EA_C", "nF_nEA_C",
                              "Log(odds)_EA_HCP(nF)", "Log(odds)_EA_PZ/HCP(nF)",
                              "Log(Odds)_F_HCP","Comparation_Index", "Pr")
for (e in AE_list) {
  print(e)
  for (d in D_Code_list) {
    print(d)
    x <- subset(HCP_C_ROR_df, Drug_Code == d & AE == e)
    y <- matrix(ncol=4, nrow=4)
    colnames(y)  <-  c("EA","nEA", "Reporter_Type", "Drug")
    rownames(y) <- c("F_HCP", "nF_HCP", "F_C", "nF_C")
    y[1,1] <- x$F_EA[x$Reporter_Type == "HCP"]
    y[1,2] <- x$F_nEA[x$Reporter_Type == "HCP"]
    y[1,3] <- "HCP"
    y[1,4] <- "1"
    y[2,1] <- x$nF_EA[x$Reporter_Type == "HCP"]
    y[2,2] <- x$nF_nEA[x$Reporter_Type == "HCP"]
    y[2,3] <- "HCP"
    y[2,4] <- "0"
    y[3,1] <- x$F_EA[x$Reporter_Type == "C"]
    y[3,2] <- x$F_nEA[x$Reporter_Type == "C"]
    y[3,3] <- "C"
    y[3,4] <- "1"
    y[4,1] <- x$nF_EA[x$Reporter_Type == "C"]
    y[4,2] <- x$nF_nEA[x$Reporter_Type == "C"]
    y[4,3] <- "C"
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



# Boxplots --------------------------------------------------------------------
pdf("Visualization/Boxplots.pdf")
for (i in AE_list){
  x <- ROR_df %>%
    filter(AE == i) %>%
    filter(ROR_m > 1) %>%
    mutate(r = rank(ROR)) %>%
    filter(r <= 20) %>% 
    mutate(Drug_Family = substr(Drug_Code, start = 1, stop = 4))
  print(ggplot(data=x) +
          geom_boxplot(mapping=aes(x=reorder(Drug_Name, -ROR), middle = ROR, lower = ROR_m, ymin=ROR_m, upper = ROR_M, ymax=ROR_M, fill = Drug_Family), stat = "identity") +
          labs(title    = i) +
          xlab("Drug") +
          geom_label(aes(x = reorder(`Drug_Name`, -ROR), ROR, label = ROR), size = 1, colour = "red", fill="white") +
          coord_flip())
}
dev.off()

# LRM -------------------------------------------------------------------------
load("Rda/ROR_df.Rda")
PhD <- read_delim("Databases Used/PhD.csv", ";", escape_double = FALSE, trim_ws = TRUE)
PhD$pCHEMBL <- PhD$pCHEMBL/1000000
LRM <- function(df, PhD_df){
  Targets_list <- as.list(unique(PhD_df$`Target`))
  Action_list <- as.list(unique(PhD_df$Action))
  df <- df %>%
    mutate(Drug_Family = substr(df$Drug_Code, start = 1, stop = 4))
  LRM_df <- as.data.frame(matrix(nrow=0, ncol=8))
  colnames(LRM_df) <- c("AE", "Target","Action", "Intercept", "Slope", "SE", "p_value", "Pearson")
  pdf("Visualization/LRM.pdf")
  for (e in AE_list){
    x <- subset(df, df$AE == e)
    x <- subset(x, is.na(x$ROR) == FALSE)
    x <- subset(x, is.infinite(x$ROR) == FALSE)
    for (t in Targets_list) {
      y <- subset(PhD_df, `Target` == t) %>%
        rename(Drug_Name = Substance)
      z1 <- left_join(x,y, by ="Drug_Name")
      z1 <- subset(z1, is.na(z1$pCHEMBL) == FALSE)
      for (m in Action_list){
        z <- subset(z1, z1$Action == m)
        if (dim(z)[1] >= 4){
          Intercept <- round(coefficients(summary(lm(z$ROR~z$pCHEMBL)))[1],2)
          Slope <- round(coefficients(summary(lm(z$ROR~z$pCHEMBL)))[2],2)
          SE <- round(coefficients(summary(lm(z$ROR~z$pCHEMBL)))[4],2)
          p_value <- round(coefficients(summary(lm(z$ROR~z$pCHEMBL)))[8],6)
          P <- cor.test(z$ROR, z$pCHEMBL, method="pearson")
          Pearson <- round(P$estimate,2)
          LRM_df[nrow(LRM_df)+1,] <- c(e, t, m, Intercept, Slope, SE, p_value, Pearson)
          if (p_value <= 1){
            plot1 <- ggplot(data = z, aes(x=z$pCHEMBL, y=z$ROR, main= paste("ROR ~ ", t))) +
              geom_smooth(method ="lm") +
              geom_point(aes(color = z$Drug_Family)) +
              xlab("pCHEMBL") +
              ylab("ROR")
            title <- ggdraw() + draw_label(paste(e," ~ ", t, "[",m,"]"), fontface = "bold")
            results <- ggdraw() + draw_label(paste("Intercept: ", Intercept,
                                                   "     Slope: ", Slope,
                                                   "     SE: ", SE,
                                                   "     p-value: ", p_value,
                                                   "     Pearson: ", Pearson), , size = 10)
            plot2 <- ggplot(data = z, aes(x=z$pCHEMBL, y=z$ROR, main= paste("ROR ~ ", t))) +
              geom_smooth(aes(color = Drug_Family), method ="lm") +
              geom_point(aes(color = Drug_Family)) +
              xlab("pCHEMBL") +
              ylab(NULL)
            legenda <- get_legend(plot2)
            p <- plot_grid(plot_grid(plot1 + theme(legend.position = "none"),
                                     plot2 + theme(legend.position = "none"),
                                     labels = "AUTO",
                                     rel_widths = c(1,1), align = "h"),
                           legenda,
                           rel_widths = c(2,.3))
            print(plot_grid(title, p, results, ncol=1, rel_heights=c(0.1, 1, 0.1)))
          }
        }
      }
    }
  }
  dev.off()
  LRM_df <- LRM_df %>%
    arrange(p_value) %>%
    mutate(Rank = rank(p_value))
  LRM_df <- LRM_df %>%
    mutate(BH20 = (Rank/nrow(LRM_df))*0.20) %>%
    mutate(Sign20 = (p_value <= BH20))
  write_csv2(LRM_df, "Results/LRM.csv")
}

LRM(ROR_df, PhD)
