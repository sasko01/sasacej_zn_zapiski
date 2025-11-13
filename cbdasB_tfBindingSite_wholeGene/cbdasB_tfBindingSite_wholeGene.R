setwd("C:/Users/scej4/OneDrive/Desktop/diplomska_koda/cbdasB_tfBindingSite_wholeGene")
library(tidyverse)

#TF binding sites prediction dobljeni s celotnim genom (files: ..._cbdasB.fasta)



#Plus strand---------------------------------------------------------------------------------------------------

abacus_tfb <- plus_strand('abacus_cbdasB_tfb.txt')

finola_tfb <- plus_strand('finola_cbdasB_tfb.txt')

JL5_tfb <- plus_strand('JL5_cbdasB_tfb.txt')

JLdash_tfb <- plus_strand('JLdash_cbdasB_tfb.txt')

Pbanana_tfb <- plus_strand('Pbanana_cbdasB_tfb.txt')

PK_tfb <- plus_strand('PK_cbdasB_tfb.txt')

SL2_tfb <- plus_strand('SL2_cbdasB_tfb.txt')

#Minus strand--------------------------------------------------------------------------------------------------

anc1_tfb <- minus_strand('anc1_cbdasB_tfb.txt')

anc2_tfb <- minus_strand('anc2_cbdasB_tfb.txt')

cannatonic_tfb <- minus_strand('cannatonic_cbdasB_tfb.txt')

cannbio_tfb <- minus_strand('cannbio2_cbdasB_tfb.txt')

CP_tfb <- minus_strand('CP_cbdasB_tfb.txt')

hOG1_tfb <- minus_strand('hOG1_cbdasB_tfb.txt') 

hOG2_tfb <- minus_strand('hOG2_cbdasB_tfb.txt')

JL_tfb <- minus_strand('JL_cbdasB_tfb.txt') 

JLfather_tfb <- minus_strand('JLfather_cbdasB_tfb.txt')

JLmother_tfb <- minus_strand('JLmother_cbdasB_tfb.txt')

PR_tfb <- minus_strand('PR_cbdasB_tfb.txt')

SL1_tfb <- minus_strand('SL1_cbdasB_tfb.txt')

travisCBD1_tfb <- minus_strand('travisCBD1_cbdasB_tfb.txt')

travisCBD2_tfb <- minus_strand('travisCBD2_cbdasB_tfb.txt')

pPepper_tfb <- minus_strand('pPepper_cbdasB_tfb.txt')




#Plus strand: Heatmap-----------------------------------------------------------------------------------------


plus_tfb_presence_matrix <- 
  bind_rows(
    abacus_tfb %>% mutate(Cultivar = "abacus"),
    SL2_tfb %>% mutate(Cultivar = "SL2"),
    finola_tfb %>% mutate(Cultivar = "finola"),
    PK_tfb %>% mutate(Cultivar = "PKush"),
    JLdash_tfb %>% mutate(Cultivar = "JLdash"),
    JL5_tfb %>% mutate(Cultivar = "JL5"),
    Pbanana_tfb %>% mutate(Cultivar = "Pbanana")) %>%
  distinct(Cultivar, V1) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  pivot_longer(-V1, names_to = "Cultivar", values_to = "Present") %>% 
  add_column(WholeGene = 1)


heatmap_plus_tfb <- 
  ggplot(plus_tfb_presence_matrix, aes(x = Cultivar, y = V1, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkgreen")) +
  theme_minimal() +
  labs(fill = "Present", 
       title = "Predicted TF binding sites presence across different Cannabis S. cultivars (cbdas, whole gene) - plus st.") + 
  theme(axis.text.y = element_text(size = 6, angle = 45, lineheight = 1.2),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("heatmap_plus_tfb.png", plot = heatmap_plus_tfb, width = 10, 
       height = max(6, length(unique(plus_tfb_presence_matrix$V1)) * 0.15), 
       dpi = 300, bg = "white")




#Minus strand: Heatmap-----------------------------------------------------------------------------------------


minus_tfb_presence_matrix <- 
  bind_rows(
    anc1_tfb %>% mutate(Cultivar = "anc1"),
    anc2_tfb %>% mutate(Cultivar = "anc2"),
    cannatonic_tfb %>% mutate(Cultivar = "cannatonic"),
    cannbio_tfb %>% mutate(Cultivar = "cannbio"),
    CP_tfb %>% mutate(Cultivar = "CP"),
    hOG1_tfb %>% mutate(Cultivar = "hOG1"),
    hOG2_tfb %>% mutate(Cultivar = "hOG2"),
    JL_tfb %>% mutate(Cultivar = "JL"),
    JLfather_tfb %>% mutate(Cultivar = "JLfather"),
    JLmother_tfb %>% mutate(Cultivar = "JLmother"),
    PR_tfb %>% mutate(Cultivar = "PR"),
    SL1_tfb %>% mutate(Cultivar = "SL1"),
    travisCBD1_tfb %>% mutate(Cultivar = "t1"),
    travisCBD2_tfb %>% mutate(Cultivar = "t2"),
    pPepper_tfb %>% mutate(Cultivar = "pPepper")) %>% 
  distinct(Cultivar, V1) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  pivot_longer(-V1, names_to = "Cultivar", values_to = "Present") %>% 
  add_column(WholeGene = 1)



heatmap_minus_tfb <- 
  ggplot(minus_tfb_presence_matrix, aes(x = Cultivar, y = V1, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkgreen")) +
  theme_minimal() +
  labs(fill = "Present", 
       title = "Predicted TF binding sites presence across different Cannabis S. cultivars (cbdas, whole gene) - minus st.") + 
  theme(axis.text.y = element_text(size = 6, angle = 45, lineheight = 1.2),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("heatmap_minus_tfb.png", plot = heatmap_minus_tfb, width = 10, 
       height = max(6, length(unique(minus_tfb_presence_matrix$V1)) * 0.15), 
       dpi = 300, bg = "white")




#Plus strand: Heatmap Gene Family-----------------------------------------------------------------------------


plus_family_presence_wG <- 
  bind_rows(
    abacus_tfb %>% mutate(Cultivar = "abacus"),
    SL2_tfb %>% mutate(Cultivar = "SL2"),
    finola_tfb %>% mutate(Cultivar = "finola"),
    PK_tfb %>% mutate(Cultivar = "PKush"),
    JLdash_tfb %>% mutate(Cultivar = "JLdash"),
    JL5_tfb %>% mutate(Cultivar = "JL5"),
    Pbanana_tfb %>% mutate(Cultivar = "Pbanana")) %>%
  distinct(Cultivar, V2) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  pivot_longer(-V2, names_to = "Cultivar", values_to = "Present") %>% 
  add_column(WholeGene = 1)


heatmap_plus_family_wG <- 
  ggplot(plus_family_presence_wG, aes(x = Cultivar, y = V2, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkgreen")) +
  theme_minimal() +
  labs(fill = "Present", 
       title = "Predicted TF families presence across different Cannabis S. cultivars (cbdas, whole gene) - plus st.") + 
  theme(axis.text.y = element_text(size = 6, angle = 45, lineheight = 1.2),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("heatmap_plus_family_wG.png", plot = heatmap_plus_family_wG, width = 10, 
       height = max(6, length(unique(plus_family_presence_wG$V2)) * 0.15), 
       dpi = 300, bg = "white")


#Minus strand: Heatmap Gene Family-----------------------------------------------------------------------------


minus_family_presence_wG <- 
  bind_rows(
    anc1_tfb %>% mutate(Cultivar = "anc1"),
    anc2_tfb %>% mutate(Cultivar = "anc2"),
    cannatonic_tfb %>% mutate(Cultivar = "cannatonic"),
    cannbio_tfb %>% mutate(Cultivar = "cannbio"),
    CP_tfb %>% mutate(Cultivar = "CP"),
    hOG1_tfb %>% mutate(Cultivar = "hOG1"),
    hOG2_tfb %>% mutate(Cultivar = "hOG2"),
    JL_tfb %>% mutate(Cultivar = "JL"),
    JLfather_tfb %>% mutate(Cultivar = "JLfather"),
    JLmother_tfb %>% mutate(Cultivar = "JLmother"),
    PR_tfb %>% mutate(Cultivar = "PR"),
    SL1_tfb %>% mutate(Cultivar = "SL1"),
    travisCBD1_tfb %>% mutate(Cultivar = "t1"),
    travisCBD2_tfb %>% mutate(Cultivar = "t2"),
    pPepper_tfb %>% mutate(Cultivar = "pPepper")) %>% 
  distinct(Cultivar, V2) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  pivot_longer(-V2, names_to = "Cultivar", values_to = "Present") %>% 
  add_column(WholeGene = 1)

heatmap_minus_family_wG <- 
  ggplot(minus_family_presence_wG, aes(x = Cultivar, y = V2, fill = factor(Present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "darkgreen")) +
  theme_minimal() +
  labs(fill = "Present", 
       title = "Predicted TF families presence across different Cannabis S. cultivars (cbdas, whole gene) - minus st.") + 
  theme(axis.text.y = element_text(size = 6, angle = 45, lineheight = 1.2),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("heatmap_minus_family_wG.png", plot = heatmap_minus_family_wG, width = 10, 
       height = max(6, length(unique(minus_family_presence_wG$V2)) * 0.15), 
       dpi = 300, bg = "white")




#Plus strand: Heatmap - Gene Family + TFB sites------------------------------------------------------------------


plus_tfb_family_matrix <- 
  bind_rows(
    abacus_tfb %>% mutate(Cultivar = "abacus"),
    SL2_tfb %>% mutate(Cultivar = "SL2"),
    finola_tfb %>% mutate(Cultivar = "finola"),
    PK_tfb %>% mutate(Cultivar = "PKush"),
    JLdash_tfb %>% mutate(Cultivar = "JLdash"),
    JL5_tfb %>% mutate(Cultivar = "JL5"),
    Pbanana_tfb %>% mutate(Cultivar = "Pbanana")) %>%
  distinct(Cultivar, V1, V2) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  add_column(WholeGene = 1) %>% 
  pivot_longer(cols = abacus:Pbanana, names_to = "Cultivar", values_to = "Presence")

ggplot(plus_tfb_family_matrix, aes(x = Cultivar, y = V1, fill = V2, alpha = Presence)) +
  geom_tile(color = "white") +
  scale_alpha_manual(values = c("0" = 0.1, "1" = 1)) +  
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Cultivar",
    y = "Gene ID",
    fill = "Family",
    alpha = "Presence"
  )

