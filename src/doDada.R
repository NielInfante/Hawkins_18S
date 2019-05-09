library(dada2); packageVersion("dada2")
library(tidyverse)
library(phyloseq)



baseDir <- "~/depot/projects/Hawkins/Metagenomics_Brzostek/MyGo/18S_fromGit/"
setwd(baseDir)

# List of where the files are

fnFs <- sort(list.files('~/depot/projects/Hawkins/Metagenomics_Brzostek/18Sqiime/18S_forward', pattern=".R1.fq", full.names = T))
fnRs <- sort(list.files('~/depot/projects/Hawkins/Metagenomics_Brzostek/18Sqiime/18S_reverse', pattern=".R2.fq", full.names = T))

# Get Sample names
sample.names <- sort(list.files('~/depot/projects/Hawkins/Metagenomics_Brzostek/18Sqiime/18S_forward', pattern=".R1.fq", full.names = F))
# <- sort(list.files('~/depot/projects/Hawkins/Metagenomics_Brzostek/16Sqiime/16S_reverse', pattern=".R2.fq", full.names = T))

sample.names <- gsub(x = sample.names, pattern="_S.*", replacement = '', perl=T)
sample.names <- gsub(x = sample.names, pattern="-", replacement = '')


# Take a look at some quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Create filtered files
filt_path <- paste0(baseDir,"/filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230), 
										 trimLeft = 10,     
										 maxN=0,  truncQ=2, rm.phix=TRUE,
										 compress=TRUE, multithread=TRUE,
										 #                     maxEE=c(3,4)
										 maxEE=2
)
head(out)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Learn Error Rates   - Might consider pooling for small sample sets (pool=T)
errF <- learnErrors(filtFs, multithread=T, nreads=10000000)
png('pics/error_F.png')
plotErrors(errF, nominalQ = T)
dev.off()

errR <- learnErrors(filtRs, multithread=T, nreads=10000000)
png('pics/error_R.png')
plotErrors(errR, nominalQ = T)
dev.off()


# Sample Inference
# Maybe run again with pool=T when I've got a lot of extra time
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=F)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=F)
dadaFs[[1]]
dadaRs[[1]]

# Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Skip this step for now
# Remove sequences that are inappropriatly sized
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,400)]
dim(seqtab2)
#table(nchar(getSequences(seqtab2)))
seqtab <- seqtab2

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
100 * sum(seqtab.nochim)/sum(seqtab)

# 95% of the reads are still there, even though 96% of identified sequences were removed as chimeas

# Overview of where things were filtered - see if you need to go back and tweek some parameters
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
track$Sample <- sample.names
#track$Survived <- round(100 * track$nonchim / track$input)
head(track)

tidyTrack <- track %>% gather(Step, Count, input:nonchim)
tidyTrack$Step <- factor(tidyTrack$Step, levels = c('input', 'filtered','denoised','merged','tabled','nonchim')) 
head(tidyTrack)

png('pics/filtering.png', width=1400, height=1000)
ggplot(tidyTrack) + geom_line(aes(x=Step, y=Count, group=Sample, color=Sample))
dev.off()

seqtab <- seqtab.nochim
saveRDS(seqtab, file='Data/seqtab.rds')



# Assign Taxonomy   - If you want species, need to use RDP or Silva, not green genes

taxa <- assignTaxonomy(seqtab, "Data/silva_132.18s.99_rep_set.dada2.fa.gz", multithread = T)

saveRDS(taxa, 'Data/taxa.rds')




#write.table(seqtab.nochim, file="test", sep="\t", col.names = NA, quote=F)

# Make a phylogenetic tree Using phangorn
library(phangorn)
library(DECIPHER)
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point.

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



# Some Extra
library(phyloseq)
#packageVersion("phyloseq")
library(ggplot2)

meta <- read.table('meta', header=T, sep="\t", stringsAsFactors = F)
rownames(meta) <- meta$IndexPair

# Issue making tree, try it without for the moment
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(meta), tax_table(taxa),phy_tree(fitGTR$tree))

#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(meta), tax_table(taxa))



saveRDS(ps, file="Data/PhyloseqObject.rds")
save(ps, file="Data/PS.Rdata")
