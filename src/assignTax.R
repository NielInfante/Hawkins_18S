# script to classify 28S amplicon data


#library(DECIPHER)
library(dada2)

setwd('~/depot/projects/Hawkins/Metagenomics_Brzostek/MyGo/18S_fromGit/')


taxa <- assignTaxonomy(seqtab, 'Data/fromJen/taxa.fa',multithread = T)




# Try using https://zenodo.org/record/835855#.XLSFbFNKi5M
seqtab <- readRDS('Data/seqtab.rds')


taxa <- assignTaxonomy(seqtab, 'Data/RDP_LSU/fw/fungiLSU_train_012014.fa',multithread = T)


path<-"Data/RDP_LSU/RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata"

dada2:::makeTaxonomyFasta_RDP(file.path(path, "fungiLSU_taxid_012014.fa"),
							file.path(path, "fungiLSU_taxid_012014.txt"),
							"Data.RDP_LSU/RDP_LSU_fixed_train_set_v2.fa",compress=FALSE)

dada2:::makeSpeciesFasta_RDP("/media/lauren/96BA-19E6/RDPClassifierLSU/current_Fungi_unaligned.fa", "/media/lauren/96BA-19E6/Upload/rdp_species_assignment_LSU_v2.fa", compress=FALSE)





taxa <- assignTaxonomy(seqtab, "Data/RDP_LSU/current_Fungi_unaligned.fa", multithread=TRUE)







# From Dada2 pipeline
seqtab <- readRDS('Data/seqtab.rds')

dna <- DNAStringSet(getSequences(seqtab)) # Create a DNAStringSet from the ASVs

load("Data/SILVA_SSU_r132_March2018.RData") 
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
	m <- match(ranks, x$rank)
	taxa <- x$taxon[m]
	taxa[startsWith(taxa, "unclassified_")] <- NA
	taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)


