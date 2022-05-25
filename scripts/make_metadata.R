################################################################################
# METADATA :: SAMPLES TABLE

pathMeta <- "/mnt/marcodata/PROJECTS/UZH/DATA"
metadata <- read.delim(file = file.path(pathMeta, "metadata.tab"), sep = "\t")

fq1 = list.files(path = pathMeta, pattern = "_1.fq.gz", full.names = TRUE)
names(fq1) = gsub("_[0-9].fq.gz|_", "", basename(fq1))
fq2 = list.files(path = pathMeta, pattern = "_2.fq.gz", full.names = TRUE)
names(fq2) = gsub("_[0-9].fq.gz|_", "", basename(fq2))

metadata_tbl = data.frame(
  metadata,
  sample = paste0(gsub("C", "C_", metadata$names)),
  fq1 = fq1[match(metadata$names, names(fq1))],
  fq2 = fq2[match(metadata$names, names(fq2))]
)

write.csv(metadata_tbl, 
          "/mnt/marcodata/PROJECTS/snakemake_uzh/config/metadata.csv",
          row.names = FALSE)