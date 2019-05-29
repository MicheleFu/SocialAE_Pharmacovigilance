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
load("Rda/PZ_df.Rda")
load("Rda/D_list.Rda")

ATC <- read_delim("ATC.csv", ";", escape_double = FALSE, 
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

ROR_df <- Wrangle(ICSR_df,D_list)
save(ROR_df, file = "Rda/ROR_df.Rda")

# Heatmap ---------------------------------------------------------------------
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

Print_Heatmap(ROR_df)

# HCP-PZ ROR ------------------------------------------------------------------
HCP_ROR_df <- Wrangle(HCP_df,D_list) %>%
  mutate(Reporter_Type = "HCP")
PZ_ROR_df <- Wrangle(PZ_df,D_list) %>%
  mutate(Reporter_Type = "PZ")
HCP_PZ_ROR_df <- rbind(HCP_ROR_df,PZ_ROR_df) %>%
  arrange(Drug_Code, AE)
save(HCP_PZ_ROR_df, file = "Rda/HCP_PZ_ROR_df.Rda")

HCP <- filter(HCP_PZ_ROR_df, HCP_PZ_ROR_df$Reporter_Type == "HCP")
Print_Heatmap(HCP)
PZ <- filter(HCP_PZ_ROR_df, HCP_PZ_ROR_df$Reporter_Type == "PZ")
Print_Heatmap(PZ)

# LRM -------------------------------------------------------------------------

PhD <- read_delim("Universe.csv", ";", escape_double = FALSE, trim_ws = TRUE)
PhD$pCHEMBL <- PhD$pCHEMBL/1000000
LRM <- function(df, PhD_df){
  Targets_list <- as.list(unique(PhD_df$`Target`))
  Action_list <- as.list(unique(PhD_df$Action))
  df <- df %>%
    mutate(Drug_Family = substr(df$Drug_Code, start = 1, stop = 4))
  LRM_df <- as.data.frame(matrix(nrow=0, ncol=7))
  colnames(LRM_df) <- c("AE", "Target","Action", "Intercept", "Slope", "SE", "p_value")
  pdf("LRM.pdf")
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
          LRM_df[nrow(LRM_df)+1,] <- c(e, t, m, Intercept, Slope, SE, p_value)
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
                                                   "     p-value: ", p_value))
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
  write_csv2(LRM_df, "LRM_csv/LRM.csv")
}

LRM(ROR_df, PhD)

