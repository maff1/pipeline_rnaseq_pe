#!/usr/bin/env Rscript

################################################################################
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
gse = readRDS(args[1])
lsFiles <- list.files(path = "/mnt/marcodata/PROJECTS/snakemake_uzh/results/star",
                      pattern = "_ReadsPerGene.out.tab",
                      full.names = TRUE); names(lsFiles) <- gsub("_ReadsPerGene.out.tab", "", basename(lsFiles))
lsReads = lapply(lsFiles, function(r) read.delim(file = r, skip = 4, header = FALSE))
star_reads = do.call(cbind,lsReads)
  star_reads2 = star_reads[,-seq(1, 21, by=4)[-1]]
    unstrand = star_reads2[, c(1, grep("\\.V2$", colnames(star_reads2)))]
    rownames(unstrand) <- unstrand[, 1]; unstrand[, 1] <- NULL
      forward = star_reads2[, c(1, grep("\\.V3$", colnames(star_reads2)))]
      rownames(forward) <- forward[, 1]; forward[, 1] <- NULL
       reverse = star_reads2[, c(1, grep("\\.V4$", colnames(star_reads2)))]
       rownames(reverse) <- reverse[, 1]; reverse[, 1] <- NULL

       lapply(list(unstrand,forward,reverse), function(x) colSums(x))       
