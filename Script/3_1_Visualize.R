# Introduzione ----------------------------------------------------------------
# Questo codice crea grafici per facilitare la visualizzazione dei dati di PhV.
# Fusaroli Michele, 2019, Progetto SDR di ADR sociali.

# Librerie necessarie ---------------------------------------------------------

library(tidyverse)

# Oggett necessari ------------------------------------------------------------
load("Rda/df_separated.Rda")
load("Rda/Matrix.Rda")

# Functions -------------------------------------------------------------------
Segnalazioni <- function(dataframe) {
  # Crea un grafico a barre che confronta n dei 25 farmaci più segnalati
  #
  # Args:
  #   dataframe: dataframe contenente tutte le segnalazioni di un singolo EA,
  #              pulito
  #
  # Returns:
  #   Grafico a barre orizzontali con lunghezza dipendente da n
  #
  #
  # Seleziona i 25 farmaci con più segnalazioni per lo specifico EA,
  # per dare grafici più chiari (meno affollati)
  subset(dataframe, r<25) %>%
    # crea il grafico con stacking per tipologia di reporter
    ggplot(mapping=aes(x    = reorder(`Substance`, n),
                       fill = `Reporter Type`)) +
    geom_bar() +
    theme(axis.text = element_text(size = 5)) +
    # imposta layout (stringa s è generata dal comando che richiama questa ƒ)
    labs(title    = "Segnalazioni per Farmaco",
         subtitle =  s) +
    # ribalta le coordinate, con farmaci sulle y e conteggio sulle x
    coord_flip() +
    xlab("Substance") +
    geom_text(aes(x    = reorder(`Substance`, n), n, label = ROR), size = 2)
}

x <- df_sep
x <- x %>%
  subset(!is.na(`AE`)) %>%
  subset(!is.na(Substance))
##



##
w <- as.data.frame(Significance_Matrix)
advev <- colnames(Significance_Matrix)


pdf("Plots.pdf")
for (e in advev) {
  dr <- rownames(w[!is.na(w[,e]),])
  dr <- list(dr)
  dr <- flatten(dr)
  s <- e
  v <- select(w, e) %>%
    rownames_to_column(var = "Substance") %>% 
    rename("ROR" = e)
  y <- subset(x, AE == e & Substance %in% dr)
  y <- merge(y,v,all.x = TRUE)
  y <- y %>% 
    add_count(Substance) %>%
    # crea una variabile r che ordina per n decrescente i diversi farmaci,
    # per permettere poi di creare grafici chiari (con i farmaci più segnalati)
    mutate(r=dense_rank(desc(n))) %>%
    # ordina i farmaci in base a n decrescente
    arrange(desc(n))
  print(Segnalazioni(y))
}
while (!is.null(dev.list()))  dev.off()
