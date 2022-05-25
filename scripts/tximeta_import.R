#!/usr/bin/env Rscript

################################################################################
# 1. import salmon quantification files --> tximeta

rm(list = ls())

library("tximeta")

args = commandArgs(trailingOnly=TRUE)

# parses salmon quantification files and combine
tximeta::setTximetaBFC(args[1], quiet = TRUE)

makeLinkedTxome(
  indexDir = args[5],
  source = "Ensembl",
  organism = "Homo sapiens",
  release = "105",
  genome = "GRCh38",
  fasta = args[6],
  gtf = args[7],
  write = TRUE,
  jsonFile = args[8]
)

# "/mnt/marcodata/PROJECTS/snakemake/config/metadata.csv"
# "/mnt/marcodata/PROJECTS/snakemake/results/salmon_quant"
# "/mnt/marcodata/PROJECTS/snakemake/results/deseq2/tximeta.rds"

samples_tbl <- read.csv(args[2])
    coldata <- data.frame(
      files = unlist(lapply(samples_tbl["sample"], 
                            function(x) file.path(args[3], x, "quant.sf"))),
      names = samples_tbl["sample"], stringsAsFactors = FALSE, row.names = NULL
      ); colnames(coldata) <- c("files", "names")
tse <- tximeta::tximeta(
          coldata = coldata, 
            type = "salmon", 
              txOut = TRUE, 
          skipMeta = FALSE,
        skipSeqinfo = FALSE,
      useHub = FALSE,
    markDuplicateTxps = FALSE,
  cleanDuplicateTxps = FALSE,
customMetaInfo = NULL
)
# Summarize abundances, counts, lengths from transcripts to gene-level
gse <- tximeta::summarizeToGene(tse)
saveRDS(gse, file = args[4])