setwd("C:/Users/scej4/OneDrive/Desktop/diplomska_koda/cbdasB_regulationPredicition")
library(tidyverse)


plus_strand2 <- function(x) {
  read.table(file = x, sep = "\t") |>
    filter(V5 == '+')
}

minus_strand2 <- function(x) {
  read.table(file = x, sep = "\t") |>
    filter(V5 == '-')
}

#Regulation prediction - vedno samo s promotorjem

#Plus strand----------------------------------------------------------------------------------------------------

abacus_reg <- plus_strand2("abacus_cbdasB_reg.txt")

finola_reg <- plus_strand2("finola_cbdasB_reg.txt")

JL5_reg <- plus_strand2("JL5_cbdasB_reg.txt")

JLdash_reg <- plus_strand2("JLdash_cbdasB_reg.txt")

pBanana_reg <- plus_strand2("pBanana_cbdasB_reg.txt")

PK_reg <- plus_strand2("PK_cbdasB_reg.txt")

SL2_reg <- plus_strand2("SL2_cbdasB_reg.txt")

#Minus strand--------------------------------------------------------------------------------------------------

anc1_reg <- minus_strand2('anc1_cbdasB_reg.txt')

anc2_reg <- minus_strand2('anc2_cbdasB_reg.txt')

cannatonic_reg <- minus_strand2('cannatonic_cbdasB_reg.txt')

cannbio_reg <- minus_strand2('cannbio2_cbdasB_reg.txt')

CP_reg <- minus_strand2('CP_cbdasB_reg.txt')

hOG1_reg <- minus_strand2('hOG1_cbdasB_reg.txt') 

hOG2_reg <- minus_strand2('hOG2_cbdasB_reg.txt')

JL_reg <- minus_strand2('JL_cbdasB_reg.txt') 

JLfather_reg <- minus_strand2('JLfather_cbdasB_reg.txt')

JLmother_reg <- minus_strand2('JLmother_cbdasB_reg.txt')

PR_reg <- minus_strand2('PR_cbdasB_reg.txt')

SL1_reg <- minus_strand2('SL1_cbdasB_reg.txt')

travisCBD1_reg <- minus_strand2('travisCBD1_cbdasB_reg.txt')

travisCBD2_reg <- minus_strand2('travisCBD2_cbdasB_reg.txt')

Ppepper <- minus_strand2('Ppepper_cbdasB_reg.txt') #dodaj (plantregmap not working)




#Plus strand: Heatmap-----------------------------------------------------------------------------------------




plus_reg_presence_matrix <- 
  bind_rows(
  abacus_reg %>% mutate(Cultivar = "abacus"),
  SL2_reg %>% mutate(Cultivar = "SL2"),
  finola_reg %>% mutate(Cultivar = "finola"),
  PK_reg %>% mutate(Cultivar = "PKush"),
  JLdash_reg %>% mutate(Cultivar = "JLdash"),
  pBanana_reg %>% mutate(Cultivar = "Pbanana")) %>%
  distinct(Cultivar, V1) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  pivot_longer(-V1, names_to = "Cultivar", values_to = "Present")


heatmap_plus_reg <- 
  ggplot(plus_reg_presence_matrix, aes(x = Cultivar, y = V1, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkgreen")) +
  theme_minimal() +
  labs(fill = "Present", 
       title = "Predicted Regulatory sites presence across different Cannabis S. cultivars (cbdas) - plus st.")
  +theme(axis.text.y = element_text(size = 6, angle = 45, lineheight = 1.2),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("heatmap_plus_reg.png", plot = heatmap_plus_reg, width = 8, 
       height = max(6, length(unique(plus_reg_presence_matrix$V1)) * 0.15), 
       dpi = 300, bg = "white")




#Minus strand: Heatmap-----------------------------------------------------------------------------------------
#Manjka pPepper!  

minus_reg_presence_matrix <- 
  bind_rows(
  anc1_reg %>% mutate(Cultivar = "anc1"),
  anc2_reg %>% mutate(Cultivar = "anc2"),
  cannatonic_reg %>% mutate(Cultivar = "cannatonic"),
  cannbio_reg %>% mutate(Cultivar = "cannbio"),
  CP_reg %>% mutate(Cultivar = "CP"),
  hOG1_reg %>% mutate(Cultivar = "hOG1"),
  hOG2_reg %>% mutate(Cultivar = "hOG2"),
  JL_reg %>% mutate(Cultivar = "JL"),
  PR_reg %>% mutate(Cultivar = "PR"),
  SL1_reg %>% mutate(Cultivar = "SL1"),
  travisCBD1_reg %>% mutate(Cultivar = "t1"),
  travisCBD2_reg %>% mutate(Cultivar = "t2")) %>% 
  distinct(Cultivar, V1) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  pivot_longer(-V1, names_to = "Cultivar", values_to = "Present")



heatmap_minus_reg <- 
  ggplot(minus_reg_presence_matrix, aes(x = Cultivar, y = V1, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkgreen")) +
  theme_minimal() +
  labs(fill = "Present", 
       title = "Predicted Regulatory sites presence across different Cannabis S. cultivars (cbdas) - minus st.") + 
  theme(axis.text.y = element_text(size = 6, angle = 45, lineheight = 1.2),
       axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("heatmap_minus_reg.png", plot = heatmap_minus_reg, width = 8, 
       height = max(6, length(unique(minus_reg_presence_matrix$V1)) * 0.15), 
       dpi = 300, bg = "white")





