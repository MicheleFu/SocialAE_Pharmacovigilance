# Introduction ----------------------------------------------------------------
# This code explores the putative mechanisms by which
# anti-Parkinson drugs can cause gambling, compulsive shopping and
# hypersexuality
# Fusaroli Michele, 2019, SDR of Social ADR.

# LRM -------------------------------------------------------------------------

DopamineR <- read_delim("Databases Used/DopamineR.csv", 
                        +     ";", escape_double = FALSE, trim_ws = TRUE)
AE_list <-  c("gambling", "gambling disorder", "compulsive shopping",
              "hypersexuality", "compulsive sexual behaviour",
              "impulsive behaviour")

LRM <- function(df, PhD_df){
  Targets_list <- as.list(unique(PhD_df$Target))
  Action_list <- as.list(unique(PhD_df$Action))
  LRM_df <- as.data.frame(matrix(nrow=0, ncol=8))
  colnames(LRM_df) <- c("AE", "Target","Action", "Intercept", "Slope", "SE", "p_value", "Pearson")
  pdf("Visualization/ICD.pdf")
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
              geom_point(aes(color = z$Drug_Name)) +
              xlab("pCHEMBL") +
              ylab("ROR")
            title <- ggdraw() + draw_label(paste(e," ~ ", t, "[",m,"]"), fontface = "bold")
            results <- ggdraw() + draw_label(paste("Intercept: ", Intercept,
                                                   "     Slope: ", Slope,
                                                   "     SE: ", SE,
                                                   "     p-value: ", p_value,
                                                   "     Pearson: ", Pearson), size = 10)
            print(plot_grid(title, plot1, results, ncol=1, rel_heights=c(0.1, 1, 0.1)))
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
  ICD_df <- LRM_df
  write_csv2(ICD_df, "Results/ICD.csv")
}

LRM(ROR_df, DopamineR)
