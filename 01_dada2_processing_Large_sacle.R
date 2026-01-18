############################################################
# Script 01: DADA2-based ASV inference and preprocessing
#
# Purpose:
#   This script processes raw paired-end 16S rRNA amplicon
#   sequencing reads to generate high-quality ASVs
#   (Amplicon Sequence Variants), taxonomic assignments,
#   and a phylogenetic tree for downstream microbiome analyses.
#
# Major steps:
#   1. Quality inspection of raw FASTQ files
#   2. Primer identification and removal (cutadapt)
#   3. Quality filtering and trimming
#   4. ASV inference using DADA2
#   5. Chimera removal
#   6. Taxonomic assignment using IDTAXA (SILVA)
#   7. Construction of ASV tables and sample metadata
#   8. Phylogenetic tree construction
#   9. Generation of phyloseq objects and coverage estimation
#
# Input:
#   - Paired-end FASTQ files (16S rRNA V4 region)
#
# Output:
#   - ASV count table (ASV_counts.tsv)
#   - Taxonomy table (ASV_identification.tsv)
#   - ASV sequences (ASV.fasta)
#   - Sample metadata (Information.tsv)
#   - Phylogenetic tree (Newick format)
#   - Read processing statistics (track.tsv)
#
# This script should be executed prior to downstream
# ecological analyses (e.g., microeco, SourceTracker).
#
# Author: [Your Name]
# Affiliation: [Your Institution]
# Date: [YYYY-MM-DD]
############################################################


############################
# 1. Load required packages
############################
library(dada2)
library(DECIPHER)
library(Biostrings)
library(phyloseq)
library(tidyverse)
library(ape)
library(ShortRead)
library(doParallel)
library(entropart)


############################
# 2. Define raw data paths
############################
# Each directory corresponds to a biological stage or tissue
path <- c(
  "~/raw_60/Bud",
  "~/raw_60/F2",
  "~/raw_60/R3",
  "~/raw_60/R5"
)

fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))


############################
# 3. Quality inspection
############################
# Visual inspection of raw read quality profiles
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)


############################
# 4. Primer identification
############################
# 16S rRNA V4 primers (515F–805R)
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GACTACHVGGGTATCTAATCC"

# Generate all primer orientations
allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(
    Forward = dna,
    Complement = complement(dna),
    Reverse = reverse(dna),
    RevComp = reverseComplement(dna)
  )
  sapply(orients, toString)
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


############################
# 5. Primer removal (cutadapt)
############################
setwd("~/raw_60")

fnFs.filtN <- file.path("sequence_files", "N-filtered", basename(fnFs))
fnRs.filtN <- file.path("sequence_files", "N-filtered", basename(fnRs))

# Remove reads containing ambiguous bases
filterAndTrim(
  fnFs, fnFs.filtN,
  fnRs, fnRs.filtN,
  maxN = 0,
  compress = TRUE,
  multithread = TRUE
)

# Remove primer sequences using cutadapt
cutadapt <- "/custom/miniconda3/bin/cutadapt"

path.cut <- file.path("sequence_files", "cutadapt")
if (!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

for (i in seq_along(fnFs)) {
  system2(
    cutadapt,
    args = c(
      R1.flags, R2.flags,
      "-n", 2,
      "-o", fnFs.cut[i],
      "-p", fnRs.cut[i],
      fnFs.filtN[i],
      fnRs.filtN[i]
    )
  )
}


############################
# 6. Quality filtering and trimming
############################
filtFs <- file.path("sequence_files", "filtered", basename(fnFs.cut))
filtRs <- file.path("sequence_files", "filtered", basename(fnRs.cut))

out <- filterAndTrim(
  fnFs.cut, filtFs,
  fnRs.cut, filtRs,
  maxN = 0,
  maxEE = c(5, 5),
  truncLen = c(240, 200),
  truncQ = 2,
  minLen = 50,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)


############################
# 7. ASV inference (DADA2)
############################
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE
)


############################
# 8. Taxonomic assignment
############################
dna <- DNAStringSet(getSequences(seqtab.nochim))
load("~/16s_rRNA_microbial_community_data/SILVA_SSU_r138_2019.RData")

ids <- IdTaxa(dna, trainingSet, strand = "top")
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)


############################
# 9. Generate ASV tables and metadata
############################
taxid_bacteria <- cbind(
  data.frame(taxid),
  as.data.frame(t(seqtab.nochim))
)

# Remove non-bacterial sequences
taxid_bacteria <- subset(taxid_bacteria, domain == "Bacteria")
taxid_bacteria <- subset(taxid_bacteria, order != "Chloroplast" | is.na(order))
taxid_bacteria <- subset(taxid_bacteria, family != "Mitochondria" | is.na(family))


############################
# 10. Phylogenetic tree construction
############################
seqs <- DNAStringSet(rownames(taxid_bacteria))
aligned_seqs <- AlignSeqs(seqs)

d <- DistanceMatrix(aligned_seqs, correction = "Jukes-Cantor")
tree <- TreeLine(d, method = "complete", cutoff = 0.03)

WriteDendrogram(
  tree,
  file = "~/raw_60/DADA/phylogenetic_tree.newick",
  quoteLabels = FALSE
)


############################
# 11. Phyloseq object & coverage
############################
OTU <- otu_table(t(seqtab.nochim), taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(taxid))
SAM <- sample_data(read.table("~/raw_60/Information.tsv", header = TRUE, row.names = 1))

physeq <- phyloseq(OTU, TAX, SAM, phy_tree(tree))


############################################################
# End of Script 01
############################################################
