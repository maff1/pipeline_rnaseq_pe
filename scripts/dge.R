#!/usr/bin/env Rscript

################################################################################
# 1. add annotations
# 2. differential gene expression analysis --> DESeq2

rm(list = ls())

library("ensembldb")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("tidyverse")
library("DESeq2")
#library("limma")

args = commandArgs(trailingOnly=TRUE)
gse = readRDS(args[1])
gtfDir = args[2]
confDir = args[3]
destDir = args[2]
designFormula = args[4]
out_1 = args[5]
out_2 = args[6]
out_3 = args[7]


metadata = read.csv(confDir)
targetPheno = unlist(strsplit(designFormula, "\\ \\+\\ "))[2]
gtffile <- file.path(gtfDir,"Homo_sapiens.GRCh38.105.gtf.gz")
DB <- ensDbFromGtf(gtf = gtffile)
EDB <- EnsDb(DB)
# makeEnsembldbPackage(ensdb=EDB, version="0.0.1",
#                      maintainer="MAFF",
#                      author="M Fernandes", destDir=destDir)
# ---------------------------------------------------------------------------- #
annot <- data.frame(ensembl_id = sapply( strsplit( rownames(gse), split="\\." ), "[", 1 ))
annot$symbol <- mapIds(x = EDB,
                       keys = annot$ensembl_id,
                       column = "SYMBOL",
                       keytype = "GENEID",
                       multiVals = "first")
annot$biotype <- mapIds(x = EDB,
                        keys = annot$ensembl_id,
                        column = "GENEBIOTYPE",
                        keytype = "GENEID",
                        multiVals = "first")
annot$entrez_1 <- mapIds(x = org.Hs.eg.db,
                         keys = annot$ensembl_id,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")
annot$entrez_2 <- mapIds(x = org.Hs.eg.db,
                         keys = annot$symbol,
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")
# bind annotations with the deseq2 object

apply(annot[, c("entrez_1", "entrez_2")], 2, is.character)
apply(annot[, c("entrez_1", "entrez_2")], 2, is.numeric)
# annotations <- annot %>%
#   as.data.frame() %>%
#   mutate(entrez = coalesce(entrez_1, entrez_2)) %>%
#   select(-entrez_1, -entrez_2)

all(rownames(gse) == annot$ensembl_id)
mcols(gse) <- cbind(mcols(gse), annot)
colnames(gse) <- metadata["sample"]
colData(gse) <- merge(colData(gse), metadata, by.x = "names", by.y = "sample")
# ---------------------------------------------------------------------------- #

dds <- DESeqDataSet(gse, design = as.formula(paste0(designFormula)))
dds <- DESeq(dds, betaPrior=FALSE, parallel=TRUE)
saveRDS(dds, out_1)

contr_pairs = lapply(as.data.frame(combn(unique(metadata[[targetPheno]]), 2)),
                     function(x) c(targetPheno, x)
); names(contr_pairs) <- sapply(as.data.frame(combn(unique(metadata[[targetPheno]]), 2)), 
                function(x) paste(x, collapse = "_vs_"))

lsRes = lapply(contr_pairs, function(x) results(dds, 
              contrast=paste0(x),
              independentFiltering = TRUE,
              alpha = 0.1,
              pAdjustMethod = "BH",
              format = "DataFrame")
)

lsResAnno <- lapply(lsRes, function(res) {
  df = cbind(mcols(gse), res)
    df["symbol"] = apply(df["symbol"], 1, function(x) paste(na.omit(x), na.omit(x), sep = "|"))
    }
  )
openxlsx::write.xlsx(lsResAnno, 
                     file = out_2)
# ---------------------------------------------------------------------------- #

# counts batch corrected
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
#mat <- limma::removeBatchEffect(mat, vsd$batch)
assay(vsd) <- mat
counts_batch_corrected <- assay(vsd)
write.table(counts_batch_corrected, file = out_3, sep = "\t")