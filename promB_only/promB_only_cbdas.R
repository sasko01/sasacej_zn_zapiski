setwd("C:/Users/scej4/OneDrive/Desktop/diplomska_koda/promB_only/TF_binding_sites_promB")
library(tidyverse)

plus_strand <- function(x) {
  read.table(file = x, sep = "\t") |>
    filter(V6 == '+')
}

minus_strand <- function(x) {
  read.table(file = x, sep = "\t") |>
    filter(V6 == '-')
}


#Tukaj so uporabljeni le promotorji tj. 2028bp pred genom CBDAS




#PLUS str._______________________________________________________________________________________________________


abacus_promB <- plus_strand('abacus_promB_tfb.txt')

finola_promB <- plus_strand('finola_promB_tfb.txt')

JL5_promB <- plus_strand('JL5_promB_tfb.txt')

JLdash_promB <- plus_strand('JLdash_promB_tfb.txt')

Pbanana_promB <- plus_strand('Pbanana_promB_tfb.txt')

PK_promB <- plus_strand('PK_promB_tfb.txt')

SL2_promB <- plus_strand('SL2_promB_tfb.txt')

anc1_promB <- plus_strand('anc1_promB_tfb.txt')

anc2_promB <- plus_strand('anc2_promB_tfb.txt')

cannatonic_promB <- plus_strand('cannatonic_promB_tfb.txt')

cannbio_promB <- plus_strand('cannbio_promB_tfb.txt')

CP_promB <- plus_strand('CP_promB_tfb.txt')

hOG1_promB <- plus_strand('hOG1_promB_tfb.txt') 

hOG2_promB <- plus_strand('hOG2_promB_tfb.txt')

JLfather_promB <- plus_strand('JLfather_promB_tfb.txt')

JLmother_promB <- plus_strand('JLmother_promB_tfb.txt')

PR_promB <- plus_strand('PR_promB_tfb.txt')

SL1_promB <- plus_strand('SL1_promB_tfb.txt')

travis1_promB <- plus_strand('travis1_promB_tfb.txt')

travis2_promB <- plus_strand('travis2_promB_tfb.txt')

pPepper_promB <- plus_strand('pPepper_promB_tfb.txt')



#Heatmap - Plus st. TF binding sites_______________________________________________________________________


plus_promB_presence_matrix <- 
  bind_rows(
    abacus_promB  %>% mutate(Cultivar = "abacus"),
    
    finola_promB %>% mutate(Cultivar = "finola"),
    
    JL5_promB %>% mutate(Cultivar = "JL5"),
    
    JLdash_promB %>% mutate(Cultivar = "JLdash"),
    
    Pbanana_promB %>% mutate(Cultivar = "Pbanana"),
    
    PK_promB %>% mutate(Cultivar = "PK"),
    
    SL2_promB %>% mutate(Cultivar = "SL2"),
    
    anc1_promB %>% mutate(Cultivar = "anc1"),
    
    anc2_promB %>% mutate(Cultivar = "anc2"),
    
    cannatonic_promB %>% mutate(Cultivar = "cannatonic"),
    
    cannbio_promB %>% mutate(Cultivar = "cannbio"),
    
    CP_promB %>% mutate(Cultivar = "CP"),
    
    hOG1_promB %>% mutate(Cultivar = "hOG1"), 
    
    hOG2_promB %>% mutate(Cultivar = "hOG2"),
    
    JLfather_promB %>% mutate(Cultivar = "JLfather"),
    
    JLmother_promB %>% mutate(Cultivar = "JLmother"),
    
    PR_promB %>% mutate(Cultivar = "PR"),
    
    SL1_promB %>% mutate(Cultivar = "SL1"), 
    
    travis1_promB %>% mutate(Cultivar = "travis1"),
    
    travis2_promB %>% mutate(Cultivar = "travis2"),
    
    pPepper_promB %>% mutate(Cultivar = "pPepper")
    ) %>%
  distinct(Cultivar, V1, V2) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cultivar, values_from = value, values_fill = 0) %>% 
  add_column(WholeGene = 1) %>% 
  pivot_longer(cols = anc1:pPepper, names_to = "Cultivar", values_to = "Presence")


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
