# Export to otu table


library(phyloseq)
library(reshape2)


# ps - phyloseq object

psm <- psmelt(psm)

otu <- reshape2::dcast(psm,  Kingdom + Phylum + Class + Order + Family + Genus + Species
											 ~ Sample, value.var = "Abundance", fun.aggregate = sum)

write_tsv(otu, 'Data/otu.txt')


getSrcDirectory(function(x) {x})


library(rstudioapi)

getActiveDocumentContext()$path 





