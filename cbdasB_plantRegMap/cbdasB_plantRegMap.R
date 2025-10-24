setwd("C:/Users/scej4/OneDrive/Desktop/diplomska_koda/cbdasB_plantRegMap")
library(tidyverse)
 
plus_strand <- function(x) {
  read.table(file = x, sep = "\t") |>
  filter(V6 == '+')
}

minus_strand <- function(x) {
  read.table(file = x, sep = "\t") |>
  filter(V6 == '-')
}

#Uporabljeni fasta files SAMO s promotorsko regijo (2000bp pred genom)


#Plus strand---------------------------------------------------------------------------------------------------

abacus_cbdasPromB_prm <- plus_strand('abacus_cbadsPromB_plantRegMap.txt')

finola_cbdasPromB_rpm.txt <- plus_strand('finola_cbdasPromB_rpm.txt')

JL5_cbdasPromB_prm <- plus_strand('JL5_cbdasPromB_prm.txt')

JLdash_cbdasPromB_prm <- plus_strand('JLdash_cbdasPromB_prm.txt')

Pbanana_cbdasPromB_prm <- plus_strand('Pbanana_cbdasPromB_prm.txt')

PK_cbdasPromB_prm <- plus_strand('PK_cbdasPromB_prm.txt')

SL2_cbdasPromB_prm <- plus_strand('SL2_cbdasPromB_prm.txt')

#Minus strand--------------------------------------------------------------------------------------------------

anc1_cbdasPromB_prm <- minus_strand('anc1_cbdasPromB_prm.txt')

anc2_cbdasPromB_prm <- minus_strand('anc2_cbdasPromB_prm.txt')

cannatonic_cbdasPromB_prm <- minus_strand('cannatonic_cbdasPromB_prm.txt')

cannbio_cbdasPromB_prm <- minus_strand('cannbio_cbdasPromB_prm.txt')

CP_cbdasPromB_rpm <- minus_strand('CP_cbdasPromB_rpm.txt')

hOG1_cbdasPromB_prm <- minus_strand('hOG1_cbdasPromB_prm.txt') 

hOG2_cbdasPromB_prm <- minus_strand('hOG2_cbdasPromB_rpm.txt')

JL_cbdasPromB_prm.txt <- minus_strand('JL_cbdasPromB_prm.txt') 

JLfather_cbdasPromB_prm <- minus_strand('JLfather_cbdasPromB_prm.txt')

JLmother_cbdasPromB_prm <- minus_strand('JLmother_cbdasPromB_prm.txt')

PR_cbdasPromB_prm <- minus_strand('PR_cbdasPromB_prm.txt')

SL1_cbdasPromB_prm <- minus_strand('SL1_cbdasPromB_prm.txt')

travisCBD1_cbdasPromB_prm <- minus_strand('travisCBD1_cbdasPromB_prm.txt')

travisCBD2_cbdasPromB_prm <- minus_strand('travisCBD2_cbdasPromB_prm.txt')
