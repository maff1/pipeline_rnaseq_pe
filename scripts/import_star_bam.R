#!/usr/bin/env Rscript
################################################################################
## RNAseq workflow - PART A:: read BAM and summarise overlaps ##################
################################################################################

rm(list = ls())

library("R.utils")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
library("BiocParallel")
register(BPPARAM=SnowParam(24))

setwd("/mnt/marcodata/PROJECTS/snakemake_uzh")
metadata <- read.csv("./config/metadata.csv")

fls <- list.files("./results/star", pattern=".bam", full.names = TRUE)
gtffile <- file.path("./reference","Homo_sapiens.GRCh38.105.gtf")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtffile, format="gtf")
ebg <- exonsBy(txdb, by="gene")
bamLst <- Rsamtools::BamFileList(fls, yieldSize=2000000)
################################################################################
sequ <- GenomicAlignments::summarizeOverlaps(ebg,
                          bamLst,
                          mode="Union",
                          singleEnd=FALSE,
                          ignore.strand=FALSE,
                          fragments=TRUE,
                          BPPARAM=SnowParam())
colnames(sequ) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(sequ))
  indexColdata <- match(colnames(sequ), metadata$sample)
    colData(sequ) <- cbind(colData(sequ), metadata[indexColdata,])
saveRDS(sequ, "./results/counts_from_star.rds")
# ---------------------------------------------------------------------------- #

#BiocManager::install("DESeq2")
library("DESeq2")
sequ <- readRDS("./results/counts_from_star.rds")
dds <- DESeqDataSet(sequ, design = as.formula("~mutant"))
dds <- DESeq(dds, betaPrior=FALSE, parallel=TRUE)
resultsNames(dds)

res = results(dds, 
        contrast=c("mutant", "DM", "BRAF"),
        independentFiltering = TRUE,
        alpha = 0.1,
        pAdjustMethod = "BH",
        format = "DataFrame", tidy = TRUE)