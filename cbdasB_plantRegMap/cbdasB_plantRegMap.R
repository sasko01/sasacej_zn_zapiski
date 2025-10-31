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

abacus <- plus_strand('abacus_cbadsPromB_plantRegMap.txt')

finola <- plus_strand('finola_cbdasPromB_rpm.txt')

JL5 <- plus_strand('JL5_cbdasPromB_prm.txt')

JLdash <- plus_strand('JLdash_cbdasPromB_prm.txt')

Pbanana <- plus_strand('Pbanana_cbdasPromB_prm.txt')

PKush <- plus_strand('PK_cbdasPromB_prm.txt')

SL2 <- plus_strand('SL2_cbdasPromB_prm.txt')

#Minus strand--------------------------------------------------------------------------------------------------

anc1 <- minus_strand('anc1_cbdasPromB_prm.txt')

anc2 <- minus_strand('anc2_cbdasPromB_prm.txt')

cannatonic <- minus_strand('cannatonic_cbdasPromB_prm.txt')

cannbio <- minus_strand('cannbio_cbdasPromB_prm.txt')

CP <- minus_strand('CP_cbdasPromB_rpm.txt')

hOG1 <- minus_strand('hOG1_cbdasPromB_prm.txt') 

hOG2 <- minus_strand('hOG2_cbdasPromB_rpm.txt')

JL <- minus_strand('JL_cbdasPromB_prm.txt') 

JLfather <- minus_strand('JLfather_cbdasPromB_prm.txt')

JLmother <- minus_strand('JLmother_cbdasPromB_prm.txt')

PR <- minus_strand('PR_cbdasPromB_prm.txt')

SL1 <- minus_strand('SL1_cbdasPromB_prm.txt')

travisCBD1 <- minus_strand('travisCBD1_cbdasPromB_prm.txt')

travisCBD2 <- minus_strand('travisCBD2_cbdasPromB_prm.txt')

Ppepper <- minus_strand('Ppepper_cbdasPromB_prm.txt')


#PLOT----------------------------------------------------------------------------------------------------------





