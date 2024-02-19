if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings", force = TRUE)
BiocManager::install("GenomicAlignments", force = TRUE)
library(Biostrings)
library(msa)
library(dplyr)
library(tidyverse)
library(genepop)
library(tidyr)
install.packages("seqinr")
library(seqinr)
library(ape)
install.packages("UniprotR")
install.packages("protti")
install.packages("r3dmol")
library(UniprotR)
library(protti)
library(r3dmol)
library(snakecase)

remotes::install_github("Tazinho/snakecase")
#load Packages####
#pacman::p_load(pacman, tidyverse, BiocManager, Biostrings, GenomicAlignments, UniprotR, protti, r3dmol, seqinr,phangorn, genepop, tidyr, ape, dplyr)

install.packages("protti", dependencies = TRUE)

#Setting the Working Directory####
#setwd("Bioinformatics/")
#setwd("Data/")
#setwd("Homework06/")

#Directory Check####
getwd()

mySequences01 <- readDNAStringSet("sequence01.fasta")
mySequences01

#Translating DNA into AA ####
aa_sequence <- Biostrings::translate(mySequences01)
aa_sequence
as.character(aa_sequence)

#AA Sequence to a FASTA file
output_file <- "amino_acid_sequence.fasta"
writeXStringSet(aa_sequence, file = output_file,
                format = "fasta", width = 60)

#4.Putting the Accession Numbers into R ####
accession_numbers<- read.table("Accession5.txt")

#5.List of Accession numbers ####
accession_numbers <- c("O78681", "I3PAK4", "A0A0F7H0C5", "A0A0F7H0D2", "Q85DY2")

# Converting it into a character string
accession_string <- paste(accession_numbers, collapse = ",")

# Print/look at the string
print(accession_string)

#6. Reading accession numbers into the GO format ####
AccessionNumbersGO <- GetProteinGOInfo(accession_numbers)
str(AccessionNumbersGO)
#write.csv(AccessionNumbersGO, "AccessionNumbersGO.csv", row.names = FALSE)

#7. Extract GO from AccessionNumbersGO ####
View(AccessionNumbersGO)
df <-- read.csv(AccessionNumbersGO.csv)
#PlotGoInfo(AccessionNumbersGO) did not work for me
go_terms <- unlist(strsplit(AccessionNumbersGO$Gene.Ontology..GO., ";"))
go_terms <- gsub("\\[.*?\\]", "", go_terms)  # Remove GO IDs from GO terms

# Creating a DF with GO terms and counts
go_counts <- data.frame(GoTerm = go_terms, Count = rep(1, length(go_terms)))

# Count summary from each GO term
go_counts <- aggregate(Count ~ GoTerm, go_counts, sum)

# Plot GO
barplot(go_counts$Count, names.arg = go_counts$GoTerm,
        xlab = "GO Terms", ylab = "Count", main = "GO Term Distribution")

#Finding associated diseases or pathologies of the gene####
GetPathology_Biotech(accession_numbers)
#NA on all counts
Get.diseases(accession_numbers)
#Error


#Structure using Protti####


#This point onwards did not work for me
#Protti would not load saying it is using a previous version of R


viewtibble <-protti::fetch_uniprot(accession_numbers)
View(viewtibble)

#Looking for available structure info from Protein DataBase
#fetch_pdb("Accession Number")


#Get Proposed 3D Structure of protti
fetch_alphafold_prediction(accession_numbers)


