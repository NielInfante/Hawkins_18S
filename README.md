# 18S Analysis

### Start with doDada.R

* Point to wwhere the raw files are
* Created meta file
* Filtering:
  * Trim 10 bases from the start of every read
  * Truncate forward reads to 240, reverse to 230
  * Remove all Ns
  * Expected error of 2
* Filter merged sequences to between 250 and 400 bases
* Remove chimeras
* Assign taxonomy using silva 132 18S


### Do Filtering with doPhyloseq.R

* Set working directory
* Step through to see how you want to filter data

Output a filtered phyloseq object

### Do differential abundance with doDeseq.R

* Set working directory
* Set output directory
* Filter the samples you want
* Set the output prefix, intGroup, contrast
* If there are too many zeros, you may need to use the alternate code to load the dds object

Outputs many pictures, and a list of significantly different taxa

### Do Beta diversity with doVeganBeta.R

### Do Alpha analysis with doAlpha.R

### Make pretty pictures with drawPics.R

