if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("GenomicAlignments") 
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)
install.packages("phangorn")
library(phangorn)
install.packages("UniprotR")
install.packages("protti", dependencies = TRUE)

#setwd()

#Set Working directory for Data and Homework folder ####

getwd()

#Get DNA Sequence
mySequences01 <- readDNAStringSet("sequence01.fasta")

view(dna_sequences)

#DNA to AA
dna_sequences <- readDNAStringSet("Sequence01.fasta")
aa<-Biostrings::translate(dna_sequences, genetic.code=GENETIC_CODE, no.init.codon=FALSE,
          if.fuzzy.codon="error")
aa

#Make AA file

writeXStringSet(aa,"aa.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

#If two packages that do same thing (match/vector arguments) specify package before
#function, then put two colons like above

#AccessionNumbers

df <- read.csv("HW6.csv", header = TRUE)
df[,1]
view(df)



#df=df%>%
 # mutate(O78681=as.factor(O78681))%>%
  #mutate(I3PAK4=as.factor(I3PAK4))%>%
 # mutate(A0A0F7H0C5=as.factor(A0A0F7H0C5))%>%
 # mutate(A0A0F7H0D2=as.factor(A0A0F7H0D2))%>%
 # mutate(Q85DY2=as.factor(Q85DY2))

df=df%>%
  mutate(O78681=as.factor(O78681))

df=df%>%
  mutate(I3PAK4=as.factor(I3PAK4))

df=df%>%
  mutate(A0A0F7H0C5=as.factor(A0A0F7H0C5))

df=df%>%
  mutate(A0A0F7H0D2=as.factor(A0A0F7H0D2))

df=df%>%
  mutate(Q85DY2=as.factor(Q85DY2))


GetProteinGOInfo(df)

# Load the package 
library(ProteinGO)



