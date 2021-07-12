#By: Patrick O'Connell
#06052020

#Identify important SNPs in ERAP1 as pathogenic or benign
#from output of mpileup SNP call file
#just change txt file that is read in each time.

setwd("~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/scRNA_SNP_project/real_snp_call_files")
library(tidyverse)

#read in txt file w/ calls
snp_calls <- read.table("ERAP1_snps_132.txt", header = FALSE)
names(snp_calls) <- c("chr", "position", "ref_allele", "end_position", "alt_allele")


#Q730 SNP. If line is printed out it is present
{snp_calls %>%
  filter(., position == "96783148") %>%
  filter(., alt_allele == "C") %>%
  print.table(.)

#K528R SNP. If line is printed out it is present
snp_calls %>%
  filter(., position == "96788627") %>%
  filter(., alt_allele == "C") %>%
  print.table(.)

#R725Q SNP. If line is printed out it is present
snp_calls %>%
  filter(., position == "96783162") %>%
  filter(., alt_allele == "T") %>%
  print.table(.)

#M349V SNP. If line is printed out it is present
snp_calls %>%
  filter(., position == "96793832") %>%
  filter(., alt_allele == "C") %>%
  print.table(.) 

#D575N SNP. If line is printed out it is present
snp_calls %>%
  filter(., position == "96786506") %>%
  filter(., alt_allele == "T") %>%
  print.table(.)
}
