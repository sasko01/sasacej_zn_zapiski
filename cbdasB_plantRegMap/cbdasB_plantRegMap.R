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

#Uporabljeni fasta files SAMO s promotorsko regijo (2000bp pred genom) - files: eg. abacus_prom


#Plus strand---------------------------------------------------------------------------------------------------

abacus_prom <- plus_strand('abacus_cbadsPromB_plantRegMap.txt')

finola_prom <- plus_strand('finola_cbdasPromB_rpm.txt')

JL5_prom <- plus_strand('JL5_cbdasPromB_prm.txt')

JLdash_prom <- plus_strand('JLdash_cbdasPromB_prm.txt')

Pbanana_prom <- plus_strand('Pbanana_cbdasPromB_prm.txt')

PKush_prom <- plus_strand('PK_cbdasPromB_prm.txt')

SL2_prom <- plus_strand('SL2_cbdasPromB_prm.txt')

#Minus strand--------------------------------------------------------------------------------------------------

anc1_prom <- minus_strand('anc1_cbdasPromB_prm.txt')

anc2_prom <- minus_strand('anc2_cbdasPromB_prm.txt')

cannatonic_prom <- minus_strand('cannatonic_cbdasPromB_prm.txt')

cannbio_prom <- minus_strand('cannbio_cbdasPromB_prm.txt')

CP_prom <- minus_strand('CP_cbdasPromB_rpm.txt')

hOG1_prom <- minus_strand('hOG1_cbdasPromB_prm.txt') 

hOG2_prom <- minus_strand('hOG2_cbdasPromB_rpm.txt')

JL_prom <- minus_strand('JL_cbdasPromB_prm.txt') 

JLfather_prom <- minus_strand('JLfather_cbdasPromB_prm.txt')

JLmother_prom <- minus_strand('JLmother_cbdasPromB_prm.txt')

PR_prom <- minus_strand('PR_cbdasPromB_prm.txt')

SL1_prom <- minus_strand('SL1_cbdasPromB_prm.txt')

travisCBD1_prom <- minus_strand('travisCBD1_cbdasPromB_prm.txt')

travisCBD2_prom <- minus_strand('travisCBD2_cbdasPromB_prm.txt')

Ppepper_prom <- minus_strand('Ppepper_cbdasPromB_prm.txt')


#---------------------------------------------------------------------------------------------------------------




