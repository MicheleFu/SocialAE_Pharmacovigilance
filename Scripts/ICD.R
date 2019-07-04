# Introduction ----------------------------------------------------------------
# This code explores the putative mechanisms by which
# anti-Parkinson drugs can cause gambling, compulsive shopping and
# hypersexuality
# Fusaroli Michele, 2019, SDR of Social ADR.

fileList <- list.files(path="FAERS ID_ICD",pattern="xlsx",full.names = T)

ICD_df <- data.frame(matrix(ncol = 24, nrow = 0))
colnames(ICD_df) <- c("Case ID",
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
  ICD_df <- rbind(ICD_df, d)
}
ICD_df <- ICD_df %>%
  distinct()
ICD_df$`Suspect Product Active Ingredients` <- tolower(ICD_df$`Suspect Product Active Ingredients`)
ICD_df$`Reactions` <- tolower(ICD_df$`Reactions`)

ICD_df <- ICD_df %>% 
  separate_rows(`Suspect Product Active Ingredients`, sep = ";") %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = "\\\\")
ICD_df$`Suspect Product Active Ingredients` <- trimws(ICD_df$`Suspect Product Active Ingredients`)
for (i in seq(nrow(ATC))){
  print(i)
  ICD_df$`Suspect Product Active Ingredients`[ICD_df$`Suspect Product Active Ingredients`== ATC$`Search term`[i]] = ATC$Code[i]
}
ICD_df <- ICD_df %>%
  separate_rows(`Suspect Product Active Ingredients`, sep = " & ") %>%
  distinct() %>%
  group_by(`Case ID`) %>%
  summarise(`Suspect Product Active Ingredients` = paste(`Suspect Product Active Ingredients`, collapse = ";"),
            `Reactions` = first(`Reactions`),
            `Reason for Use` = first(`Reason for Use`),
            `Reporter Type` = first(`Reporter Type`),
            `Serious` = first(`Serious`),
            `Outcomes` = first(`Outcomes`),
            `Sex` = first(`Sex`),
            `Patient Age`= first(`Patient Age`),
            `Patient Weight`= first(`Patient Weight`),
            `Concomitant Product Names` = first(`Concomitant Product Names`),
            `Date` = first(`Latest FDA Received Date`))
ICD_df$`Patient Age` <-  na_if(ICD_df$`Patient Age`, "Not Specified")
ICD_df <- ICD_df %>%
  separate(`Patient Age`, into = c("Age (YR)","Age_Unit"), sep = " ", convert = TRUE)
for (i in 1:nrow(ICD_df)){
  if (ICD_df$Age_Unit[i] == "DEC" & is.na(ICD_df$`Age (YR)`[i]) == FALSE){
    ICD_df$`Age (YR)`[i] <- round(ICD_df$`Age (YR)`[i]*10,1)
  }
  else if (ICD_df$Age_Unit[i] == "MTH" & is.na(ICD_df$`Age (YR)`[i]) == FALSE){
    ICD_df$`Age (YR)`[i] <- round(ICD_df$`Age (YR)`[i]/12,1)
  }
  else if (ICD_df$Age_Unit[i] == "DAY" & is.na(ICD_df$`Age (YR)`[i]) == FALSE){
    ICD_df$`Age (YR)`[i] <- round(ICD_df$`Age (YR)`[i]/365,1)
  }
  else if (ICD_df$Age_Unit[i] == "HR" & is.na(ICD_df$`Age (YR)`[i]) == FALSE){
    ICD_df$`Age (YR)`[i] <- round(ICD_df$`Age (YR)`[i]/8760,1)
  }
  else if (ICD_df$Age_Unit[i] == "WEEK" & is.na(ICD_df$`Age (YR)`[i]) == FALSE){
    ICD_df$`Age (YR)`[i] <- round((ICD_df$`Age (YR)`[i]/365)*7,1)
  }
}
ICD_df$`Patient Weight` <-  na_if(ICD_df$`Patient Weight`, "Not Specified")
ICD_df <- ICD_df %>%
  separate(`Patient Weight`, into = c("Weight (KG)","Weight_Unit"), sep = " ", convert = TRUE)
ICD_df$`Weight (KG)` <- parse_number(as.character(ICD_df$`Weight (KG)`))
for (i in 1:nrow(ICD_df)){
  if (ICD_df$Weight_Unit[i] == "LB" & is.na(ICD_df$`Weight (KG)`[i]) == FALSE & is.na(ICD_df$Weight_Unit[i]) == FALSE){
    ICD_df$`Weight (KG)`[i] <- round(ICD_df$`Weight (KG)`[i]*0.4536,1)
  }
  else if (ICD_df$Weight_Unit[i] == "GMS" & is.na(ICD_df$`Weight (KG)`[i]) == FALSE & is.na(ICD_df$Weight_Unit[i]) == FALSE){
    ICD_df$`Weight (KG)`[i] <- round(ICD_df$`Weight (KG)`[i]/1000,1)
  }
}
ICD_df <- ICD_df %>%
  select(-Weight_Unit, -Age_Unit) %>%
  mutate(Date = substr(Date, nchar(Date)-3, nchar(Date))) %>%
  mutate(gambling = str_detect(string = Reactions, pattern = "gambling")) %>%
  mutate(compulsive_shopping = str_detect(string = Reactions, pattern = "compulsive shopping")) %>%
  mutate(hypersexuality = str_detect(string = Reactions, pattern = "hypersexuality"))
ICD_df$`Weight (KG)` <- round(ICD_df$`Weight (KG)`,1)
ICD_df$Date <- parse_number(ICD_df$Date)
save(ICD_df, file = "RDA/ICD_df.RDA")



#Population description -------------------------------------------------------
Char_df <- as.data.frame(matrix(ncol = 45, nrow = 1), stringsAsFactors = FALSE)
variables <- c("Set","Number","N%", "Mean_Age", "<2y","%","2-12y","%","12-18y","%",
               "18-30y","%","30-45y","%","45-65y","%","65-75y","%", ">75y","%",
               "Female","F%", "Male","M%", "Mean_Weight", "HCP", "HCP%", "PZ","PZ%",
               "<=2004","%", "2005-2006","%", "2007-2008","%","2009-2010","%",
               "2011-2012","%","2013-2014","%","2015-2016","%","2017-2018","%")
colnames(Char_df) <- c(variables)
Counting <- function(set,x1){
  number <- nrow(x1)
  Np <- round(number/nrow(ICD_df)*100,2)
  Mean <- round(mean(x1$`Age (YR)`, na.rm = TRUE),0)
  less2 <- sum(x1$`Age (YR)` < 2, na.rm = TRUE)
  less2p <- round((less2/number)*100,2)
  less11 <- sum(x1$`Age (YR)`< 12 & x1$`Age (YR)`>= 2, na.rm = TRUE)
  less11p <- round((less11/number)*100,2)
  less17 <- sum(x1$`Age (YR)`< 18 & x1$`Age (YR)`>= 12, na.rm = TRUE)
  less17p <- round((less17/number)*100,2)
  less30 <- sum(x1$`Age (YR)`< 30 & x1$`Age (YR)`>= 18, na.rm = TRUE)
  less30p <- round((less30/number)*100,2)
  less45 <- sum(x1$`Age (YR)`< 45 & x1$`Age (YR)`>= 30, na.rm = TRUE)
  less45p <- round((less45/number)*100,2)
  less65 <- sum(x1$`Age (YR)`< 65 & x1$`Age (YR)`>= 45, na.rm = TRUE)
  less65p <- round((less65/number)*100,2)
  less74 <- sum(x1$`Age (YR)`< 75 & x1$`Age (YR)`>= 65, na.rm = TRUE)
  less74p <- round((less74/number)*100,2)
  more74 <- sum(x1$`Age (YR)`>= 75, na.rm = TRUE)
  more74p <- round((more74/number)*100,2)
  F_num <-sum(x1$Sex == "Female")
  Fp <- round((F_num/number)*100,2)
  M_num <- sum(x1$Sex == "Male")
  Mp <- round((M_num/number)*100,2)
  Weight <- round(mean(x1$`Weight (KG)`, na.rm = TRUE),0)
  HCP_num <- sum(x1$`Reporter Type` == "Healthcare Professional")
  HCPp <- round((HCP_num/number)*100,2)
  PZ_num <- sum(x1$`Reporter Type` == "Consumer")
  PZp <- round((PZ_num/number)*100,2)
  less2004 <- sum(x1$`Date` <= 2004, na.rm = TRUE)
  p2004 <- round((less2004/number)*100,2)
  less2006 <- sum(x1$`Date` > 2004 & x1$`Date` <= 2006, na.rm = TRUE)
  p2006 <- round((less2006/number)*100,2)
  less2008 <- sum(x1$`Date` > 2006 & x1$`Date` <= 2008, na.rm = TRUE)
  p2008 <- round((less2008/number)*100,2)
  less2010 <- sum(x1$`Date` > 2008 & x1$`Date` <= 2010, na.rm = TRUE)
  p2010 <- round((less2010/number)*100,2)
  less2012 <- sum(x1$`Date` > 2010 & x1$`Date` <= 2012, na.rm = TRUE)
  p2012 <- round((less2012/number)*100,2)
  less2014 <- sum(x1$`Date` > 2012 & x1$`Date` <= 2014, na.rm = TRUE)
  p2014 <- round((less2014/number)*100,2)
  less2016 <- sum(x1$`Date` > 2014 & x1$`Date` <= 2016, na.rm = TRUE)
  p2016 <- round((less2016/number)*100,2)
  more2016 <- sum(x1$`Date` > 2016, na.rm = TRUE)
  p2018 <- round((more2016/number)*100,2)
  new_row <- c(set, number, Np, Mean, less2, less2p, less11, less11p, 
               less17, less17p,less30,less30p,less45,less45p,less65,less65p,
               less74,less74p,more74,more74p,F_num,Fp,M_num,Mp,Weight,
               HCP_num,HCPp,PZ_num,PZp,less2004,p2004,less2006,p2006,
               less2008,p2008,less2010,p2010,less2012,p2012,less2014,
               p2014,less2016,p2016,more2016,p2018)
  return(new_row)
}

Descr <- function(AE){
  set <- "Case"
  x <-subset(ICD_df, ICD_df[,AE] == TRUE)
  row_Case <- Counting(set,x)
  y <-subset(ICD_df, ICD_df[,AE] == FALSE)
  set <- "Non_Case"
  row_NonCase <- Counting(set,y)
  Char_df <- rbind(Char_df, row_NonCase)
  Char_df <- rbind(Char_df, row_Case)
}
  
Char_df <- Descr("gambling")
Char_df <- Char_df[-1,]
Char_df <- as.data.frame(sapply(Char_df, as.numeric))
Char_df$Set <- c("Non_Case","Case")
summary(glm(cbind(HCP,PZ)~Set, data = Char_df, family="binomial"))
summary(results)
chisq.test(cbind(Char_df$HCP,Char_df$PZ), Char_df$Set, correct = FALSE)
exp(coefficients(results))


Char_df <- Descr("compulsive_shopping")
Char_df <- Descr("hypersexuality")
save(Char_df, file = "RDA/Char_df.RDa")

D_list <- c("G04BE07, N04BC07", "N04BC01, G02CB01","N04BC06, G02CB03",
            "G02CB02, N02CA07","N04BC02", "N04BC08", "N04BC05",
            "N04BC04", "N04BC09", "N04BX-a")
AE_list <- c("gambling", "compulsive shopping", "hypersexuality")
Calculate_ROR <- function(df) {
  #
  # Args:
  #   df  : Dataframe of observations
  #
  # Returns:
  #   Results_df : Disproportionality Analysis and summarized data
  ICD_Res_df <- data.frame(matrix(ncol = 20, nrow = 0))
  colnames(ICD_Res_df) <- c("D_Name","AE", "F_EA", "F_nEA","nF_EA","nF_nEA",
                            "Intercept", "Intercept_min", "Intercept_max",
                            "Exposure1", "Exposure1_min", "Exposure1_max")
  for (i in 1:length(AE_list)) {
    AE_Name <- AE_list[[i]]
    print(AE_Name)
    y <- df
    y <- y %>%
      mutate(AE = str_detect(string = Reactions, pattern=AE_Name))
    EA <- sum()
    nEA <- sum(y$AE == FALSE)
    for (d in D_list){
      y <- y %>%
        mutate(D = str_detect(string = `Suspect Product Active Ingredients`, pattern=d))
      F_EA <- sum(y$D == TRUE & y$AE == TRUE)
      nF_EA <- sum(y$D == FALSE & y$AE == TRUE)
      F_nEA <- sum(y$D == TRUE & y$AE == FALSE)
      nF_nEA <- sum(y$D == FALSE & y$AE == FALSE)
      D_Name <- ATC$Substance[ATC$Code == d]
    RORraw=glm(AE ~ factor(D), data = y, family = binomial) 
    Intercept <- exp(RORraw$coefficients[1])
    Exposure1 <- exp(RORraw$coefficients[2])
    Intercept_min <- exp(confint(RORraw, level = 0.95)[1]) 
    Exposure1_min <- exp(confint(RORraw, level = 0.95)[2])
    Intercept_max <- exp(confint(RORraw, level = 0.95)[3])
    Exposure1_max <- exp(confint(RORraw, level = 0.95)[4])
    new_row <- list(D_Name,AE_Name, F_EA, F_nEA, nF_EA, nF_nEA,
                    Intercept, Intercept_min, Intercept_max,
                    Exposure1, Exposure1_min, Exposure1_max)
    ICD_Res_df[nrow(ICD_Res_df)+1,] <-  new_row
    }
  }
  return(ICD_Res_df)
}

ICD_Res_df <- Calculate_ROR(ICD_df)
write_csv2(ICD_Res_df, "ICD_Res_df")

pdf("Visualization/Boxplots_ICD.pdf")
for (i in AE_list){
  x <- ROR_df %>%
    filter(AE == i) %>%
    filter(ROR_m > 1) %>%
    filter(Drug_Name %in% c("apomorphine", "bromocriptine","cabergoline",
                            "lisuride","pergolide","piribedil","pramipexole",
                            "ropinirole","rotigotine"))
  print(ggplot(data=x) +
          geom_boxplot(mapping=aes(x=reorder(Drug_Name, -ROR), middle = ROR, lower = ROR_m, ymin=ROR_m, upper = ROR_M, ymax=ROR_M), stat = "identity") +
          labs(title    = i) +
          xlab("Drug") +
          geom_label(aes(x = reorder(`Drug_Name`, -ROR), ROR, label = ROR), size = 1, colour = "red", fill="white") +
          coord_flip())
}
dev.off()

ICD_PhD <- read_delim("Databases Used/ICD_PhD.csv", 
                      ";", escape_double = FALSE, trim_ws = TRUE)
h <- matrix(nrow = 9, ncol = 17)
rownames(h) <- c("apomorphine", "bromocriptine","cabergoline",
                 "lisuride","pergolide","piribedil","pramipexole",
                 "ropinirole","rotigotine")
colnames(h) <- c("ADRA1a","ADRA1b","ADRA1d","ADRA2a","ADRA2b","ADRA2c",
                 "D1","D2","D3","D4","D5","5HT1a","5HT1b","5HT1d",
                 "5HT2a","5HT2b","5HT2c")
for (x in rownames(h)){
  for (y in colnames(h)){
    print(paste(x,y))
    h[[x,y]] <- ICD_PhD$Occupancy[ICD_PhD$Substance==x & ICD_PhD$Target==y]
    if(is.na(h[x,y])){h[[x,y]] <- 0}
  }
}
png("Visualization/Occupancy.png")
superheat(h,
          title = "Occupancy",
          heat.pal = c("white","yellow", "green", "blue"),
          heat.pal.values = c(0, 0.25, 0.5, 1),
          heat.col.scheme = "red",
          heat.lim = c(0,100),
          bottom.label.text.angle = 90,
          bottom.label.text.size = 2,
          bottom.label.size = 0.1,
          force.left.label = TRUE,
          left.label.text.size = 4,
          left.label.size = 0.4,
          force.grid.hline = TRUE,
          grid.hline.col = "gray",
          grid.vline.col = "gray",
          heat.na.col= "white",
          X.text = h,
          X.text.size = 2,
          left.label.text.alignment = "right",
          pretty.order.rows = TRUE,
          left.label = "variable",
          row.title = "Substance",
          row.title.size = 6,
          column.title = "Target",
          column.title.size = 6)
dev.off()
 

LRM <- function(df, PhD_df){
  Targets_list <- as.list(unique(PhD_df$Target))
  Action_list <- as.list(unique(PhD_df$Action))
  LRM_df <- as.data.frame(matrix(nrow=0, ncol=8))
  colnames(LRM_df) <- c("AE", "Target","Action", "Intercept", "Slope", "SE", "p_value", "Pearson")
  pdf("Visualization/ICD2_occ.pdf")
  for (e in AE_list){
    x <- subset(df, df$AE == e)
    x <- subset(x, is.na(x$ROR) == FALSE)
    x <- subset(x, is.infinite(x$ROR) == FALSE)
    for (t in Targets_list) {
      y <- subset(PhD_df, `Target` == t) %>%
        rename(Drug_Name = Substance)
      z1 <- left_join(x,y, by ="Drug_Name")
      z1 <- subset(z1, is.na(z1$Occupancy) == FALSE)
      for (m in Action_list){
        z <- subset(z1, z1$Action == m)
        if (dim(z)[1] >= 4){
          Intercept <- round(coefficients(summary(lm(z$ROR~z$Occupancy)))[1],2)
          Slope <- round(coefficients(summary(lm(z$ROR~z$Occupancy)))[2],2)
          SE <- round(coefficients(summary(lm(z$ROR~z$Occupancy)))[4],2)
          p_value <- round(coefficients(summary(lm(z$ROR~z$Occupancy)))[8],6)
          P <- cor.test(z$ROR, z$Occupancy, method="pearson")
          Pearson <- round(P$estimate,2)
          LRM_df[nrow(LRM_df)+1,] <- c(e, t, m, Intercept, Slope, SE, p_value, Pearson)
          if (p_value <= 1){
            plot1 <- ggplot(data = z, aes(x=z$Occupancy, y=z$ROR, main= paste("ROR ~ ", t))) +
              geom_smooth(method ="lm") +
              geom_point(aes(color = z$Drug_Name)) +
              xlab("Occupancy") +
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
  ICD2_occ_df <- LRM_df
  write_csv2(ICD2_occ_df, "Results/ICD2_occ.csv")
}

LRM(ROR_df, ICD_PhD)

plot_pearson <- function(event){
  x <- ICD2_occ %>% 
    arrange(Pearson) %>%
    filter(AE == event)
  ggplot(data = x, mapping = aes(x = Pearson, y = reorder(Target,Pearson), color = Sign20)) +
    geom_point() +
    geom_text(aes(label= p_value),hjust=0, vjust=0) +
    coord_cartesian(xlim = c(-1, 1)) +
    ggtitle(event)
}

png("Visualization/gambling.png")
print(plot_pearson("gambling"))
dev.off()

png("Visualization/compulsive shopping.png")
print(plot_pearson("compulsive shopping"))
dev.off()

png("Visualization/hypersexuality.png")
print(plot_pearson("hypersexuality"))
dev.off()
