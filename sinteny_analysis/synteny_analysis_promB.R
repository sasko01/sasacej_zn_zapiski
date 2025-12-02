setwd("C:/Users/scej4/OneDrive/Desktop/diplomska_koda/sinteny_analysis")
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(dplyr)



all_cultivars <- bind_rows(
  read_tsv('abacus_promB_tfb.txt') %>% mutate(Cultivar = "abacus"),
  
  read_tsv('finola_promB_tfb.txt') %>% mutate(Cultivar = "finola"),
  
  read_tsv('JL5_promB_tfb.txt') %>% mutate(Cultivar = "JL5"),
  
  read_tsv('JLdash_promB_tfb.txt') %>% mutate(Cultivar = "JLdash"),
  
  read_tsv('Pbanana_promB_tfb.txt') %>% mutate(Cultivar = "Pbanana"),
  
  read_tsv('PK_promB_tfb.txt') %>% mutate(Cultivar = "PK"),
  
  read_tsv('anc1_promB_tfb.txt') %>% mutate(Cultivar = "anc1"),
  
  read_tsv('anc2_promB_tfb.txt') %>% mutate(Cultivar = "anc2"),
  
  read_tsv('cannatonic_promB_tfb.txt') %>% mutate(Cultivar = "cannatonic"),
  
  read_tsv('cannbio_promB_tfb.txt') %>% mutate(Cultivar = "cannbio"),
  
  read_tsv('CP_promB_tfb.txt') %>% mutate(Cultivar = "CP"),
  
  read_tsv('hOG1_promB_tfb.txt') %>% mutate(Cultivar = "hOG1"),
  
  read_tsv('hOG2_promB_tfb.txt') %>% mutate(Cultivar = "hOG2"),
  
  read_tsv('JLfather_promB_tfb.txt') %>% mutate(Cultivar = "JLfather"),
  
  read_tsv('JLmother_promB_tfb.txt') %>% mutate(Cultivar = "JLmother"),
  
  read_tsv('PR_promB_tfb.txt') %>% mutate(Cultivar = "PR"),
  
  read_tsv('SL1_promB_tfb.txt') %>% mutate(Cultivar = "SL1"),
  
  read_tsv('travis1_promB_tfb.txt') %>% mutate(Cultivar = "travis1"),
  
  read_tsv('travis2_promB_tfb.txt') %>% mutate(Cultivar = "travis2"),
  
  read_tsv('pPepper_promB_tfb.txt') %>% mutate(Cultivar = "pPepper")
) %>% mutate(direction = ifelse(strand == "+", 1, -1))


names(all_cultivars)[names(all_cultivars) == "matched sequence"] <- "matched_sequence"
names(all_cultivars)[names(all_cultivars) == "sequence name"] <- "sequence_name"

#Synteny plot poskus___________________________________________________________________________________________

(
ggplot(all_cultivars,
       aes(x = start, xend = stop,
           y = Cultivar, yend = Cultivar,
           color = Family)) +
  geom_segment(
    size = 5,
    arrow = arrow(
      length = unit(0.15, "cm"),
      type = "closed",
      ends = "last"
    )
  ) +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Promoter TFBS Synteny Map — CBDA-synthase",
    x = "Promoter Position (bp)",
    y = "Cultivar",
    color = "TF Family"
  )
) %>%  ggsave("synteny.png", plot = .,
                width = 12, height = 10, dpi = 300,
                bg = "white")


(
  ggplot(all_cultivars,
         aes(x = start, xend = stop,
             y = Cultivar, yend = Cultivar,
             color = motif)) +
    geom_segment(
      size = 5,
      arrow = arrow(
        length = unit(0.15, "cm"),
        type = "closed",
        ends = "last"
      )
    ) +
    scale_color_viridis_d(option = "turbo") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Promoter TFBS Synteny Map — CBDA-synthase",
      x = "Promoter Position (bp)",
      y = "Cultivar",
      color = "TF Motif"
    )
) %>%  ggsave("synteny_motifs.png", plot = .,
              width = 25, height = 15, dpi = 300,
              bg = "white")


#With GenomicRanges______________________________________________________________________________________


gr_cultivars <- with(all_cultivars,
                     GRanges(
                       seqnames = Cultivar,       
                       ranges   = IRanges(start, stop),
                       strand   = ifelse(direction == 1, "+", "-"),
                       motif    = motif,
                       family   = Family,
                       matched  = matched_sequence
                     ))


gr <- GRanges(
  seqnames = all_cultivars$Cultivar,       
  ranges = IRanges(start = all_cultivars$start, end = all_cultivars$stop),
  strand = all_cultivars$strand,
  Family = all_cultivars$Family,
  motif = all_cultivars$motif
)





#_______________________________________________________________________________________________________

library(rtracklayer)

export.bed(gr_cultivars, "cultivar_TFBS.bed")
