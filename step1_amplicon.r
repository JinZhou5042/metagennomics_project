# This is inspired by https://benjjneb.github.io/dada2/tutorial_1_8.html
# The script requires paired end data
# Script can be run with Rscript dada2_refdb_PE_scriptMode.R --input inputfolder --output --outputfolder ... etc
# Comment out lines for arguments you don't need
# Script takes extra input of refdbfile (in addition to refdir folder)
# conda activate /afs/crc/group/Bioinformatics/SP24_BIO60132/software/conda_envs/dada2

library(dada2); 
packageVersion("dada2")
library(argparser, quietly=TRUE)

p <- arg_parser("DADA2 in R")

p = add_argument(p, "--input", help="Enter input folder", type='character', default="/scratch365/zhuang8/R16S_data/workspace/input_data")
p = add_argument(p, "--output", help="Enter output folder", type='character', default="/scratch365/zhuang8/R16S_data/workspace/output_data")
p = add_argument(p, "--plots", help="Enter plot folder input_path", type='character', default="/scratch365/zhuang8/R16S_data/workspace/plot_data")
p = add_argument(p, "--refdir", help="Enter ref folder", type='character', default="/scratch365/zhuang8/R16S_data/ref_DBs")
p = add_argument(p, "--refdbfile", help="Enter nr99_v138 database", type='character', default="silva_nr99_v138.1_train_set.fa.gz")
p = add_argument(p, "--spfile", help="Enter species database", type='character', default="silva_species_assignment_v138.1.fa.gz")

argv <- parse_args(p)
input_path = argv$input
output_data_path = argv$output
output_plot_path = argv$plots
ref_path = argv$refdir
db_nr99_v138 = argv$refdbfile
db_species_assignment = argv$spfile

clear_folder <- function(folder) {
  files <- list.files(folder, recursive = TRUE, full.names = TRUE)
  
  if (length(files) > 0) {
    file.remove(files)
  }
  subdirs <- list.dirs(folder, recursive = TRUE, full.names = TRUE)
  if (length(subdirs) > 1) {  
    subdirs <- subdirs[subdirs != folder]  
    sapply(subdirs, function(x) {
      if (length(list.files(x)) == 0) {
        dir.remove(x)
      }
    })
  }
}

initial_folder <- function(folder) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  } else {
    clear_folder(folder)
  }
}

initial_folder(output_data_path)
initial_folder(output_plot_path)

######################################################################################################
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
######################################################################################################
fnFs <- sort(list.files(input_path, pattern="_1.*\\.fq$", full.names = TRUE))
fnRs <- sort(list.files(input_path, pattern="_2.*\\.fq$", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####################################################################
# Visualizing the quality profiles of the forward and reverse reads
####################################################################
pdf(file.path(output_plot_path, "dada2_plotQualityProfile.pdf"), onefile=T)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
dev.off()

#############################################################
# Perform filtering and trimming on forward and reverse reads
##############################################################
# place under filtered subdirectory
filtFs <- file.path(input_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(input_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Modify and add new parameter based on your data.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280, 220),
                    maxN=0, maxEE=c(2, 3), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE, minLen=40)
print(out)

## Examine quality profiles of filtered reads
pdf(file.path(output_plot_path, "QualityProfile.filt_plot.pdf"), onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()

#####################################################
# Learn the Error Rates for forward and reverse reads
#####################################################
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## Plot estimated error as sanity check
pdf(file.path(output_plot_path, "ErrorsRates_F.pdf"), onefile=T)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(file.path(output_plot_path, "ErrorsRates_R.pdf"), onefile=T)
plotErrors(errR, nominalQ=TRUE)
dev.off()


################
# Dereplication
################
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

###################
# Sample Inference
###################
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

####################
# Merge paired reads
####################
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###########################
# Construct sequence table (ASV)
###########################
seqtab <- makeSequenceTable(mergers)
## Get dimensions
dim(seqtab)

## Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##################
# Remove chimeras
##################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)

#################
#Save into Rds
#################
saveRDS(seqtab, file.path(output_data_path, "seqtab.Rds"))
saveRDS(seqtab.nochim, file.path(output_data_path, "seqtab.nochim.Rds"))


#Save seqtab as otu table
otutab <- t(seqtab.nochim)
write.table(otutab, file=file.path(output_data_path, "otutab.tsv"), quote=FALSE)


#################################
# Track reads through the pipeline
#################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file.path(output_data_path, "track_reads.txt"), sep=" ", quote=FALSE)

#################
# Assign taxonomy
#################
# IMP : download taxonomy file wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1
taxa <- assignTaxonomy(seqtab.nochim, file.path(
  ref_path, db_nr99_v138), multithread=TRUE, tryRC=TRUE)
taxa <- addSpecies(taxa, file.path(ref_path, db_species_assignment))
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, file.path(output_data_path, "taxa.Rds"))
write.table(taxa, file=file.path(output_data_path, "taxa.tsv"), quote=FALSE)
