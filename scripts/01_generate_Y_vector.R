# 01_generate_Y_vector.R
# Generates Y vectors (VISTA positives mm10 + negatives via bedtools shuffle)
# and applies mm10 blacklist filtering.

library(Biostrings)
library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(IRanges)

stopifnot("Run this script from the repository root" = dir.exists("scripts"))

ml_dir    <- "data/raw"
bed_dir   <- "data/bed_files"
out_dir   <- "data/Y_vectors_generated"

dir.create(bed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

blacklist_su <- "tutto"

fasta_path  <- file.path(ml_dir, "all_enhancer.fasta")
mm10_path   <- file.path(ml_dir, "ENHANCERS_mm10.bed")
bl_path     <- file.path(bed_dir, "mm10-blacklist.v2.bed.gz")
chrom_sizes <- file.path(bed_dir, "mm10.chrom.sizes")

if (!file.exists(bl_path)) {
  message("Downloading official mm10 blacklist (Boyle Lab)...")
  try(download.file("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
                    bl_path, mode = "wb"), silent = TRUE)
}
if (!file.exists(chrom_sizes)) {
  message("Downloading mm10.chrom.sizes from UCSC...")
  download.file("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes", chrom_sizes)
}

fastaFile3 <- readDNAStringSet(fasta_path)
seq_name3 <- as.data.frame(fastaFile3@ranges@NAMES)
seq_name_all <- data.frame("region" = seq_name3, "value" = 1, "seq" = paste(fastaFile3))
sep3 <- separate(seq_name_all, fastaFile3.ranges.NAMES, into = c("a", "chr", "start", "end", "element", "number", "tipo"))
all_df <- sep3[, -c(1, 5, 6)]
remove(seq_name3, seq_name_all, sep3)

all_df$bed <- paste(all_df$chr, all_df$start, all_df$end, sep = "-")
all_df <- all_df[, -5]
all_df <- all_df[-114, ]

all_mm10 <- fread(mm10_path)[, 1:3]
colnames(all_mm10) <- c("chr", "start", "end")
all_mm10$bed <- paste(all_mm10$chr, all_mm10$start, all_mm10$end, sep = "-")

new <- cbind(all_mm10, all_df)
new <- new[, -c(5, 6, 7, 9)]
new$width <- new$end - new$start

new_positive <- new[new$tipo == "positive", ]
bed2 <- transform(new_positive, 
                  chromStart = ifelse(start > end, end, start),
                  chromEnd   = ifelse(end < start, start, end))

bed3 <- GenomicRanges::makeGRangesFromDataFrame(bed2)
bed4 <- rtracklayer::export(bed3, format = "bed", ignore.strand = FALSE)

pos_full <- file.path(out_dir, "enhancers_only_positive_mm10.bed")
write.table(bed4, pos_full, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

pos3 <- file.path(out_dir, "positive_enhancers_mm10.bed")
all3 <- file.path(out_dir, "all_enhancers_mm10.bed")
fwrite(fread(pos_full, select = 1:3), pos3, sep = "\t", col.names = FALSE)
fwrite(fread(mm10_path, select = 1:3), all3, sep = "\t", col.names = FALSE)

enhancers <- as.data.frame(new[new$tipo == "positive", -c(5, 6)])
enhancers$value <- 1

bl <- fread(bl_path, header = FALSE, select = 1:3)
setnames(bl, c("chr", "start", "end"))
bl_gr <- makeGRangesFromDataFrame(bl)

rimuovi_bl <- function(df) {
  ov <- findOverlaps(makeGRangesFromDataFrame(df), bl_gr, type = "any")
  if (length(ov) > 0) df <- df[-unique(queryHits(ov)), ]
  df
}

shuffle_neg <- function(seed) {
  out <- tempfile(fileext = ".bed")
  system(sprintf("bedtools shuffle -seed %d -excl %s -i %s -g %s > %s",
                 seed, all3, pos3, chrom_sizes, out))
  n <- fread(out, select = 1:3)
  setnames(n, c("chr", "start", "end"))
  n
}

costruisci_Y <- function(neg) {
  neg <- as.data.frame(neg)
  neg$bed <- paste(neg$chr, neg$start, neg$end, sep = "-")
  neg$value <- 0
  neg <- neg[, c("chr", "start", "end", "bed", "value")]
  
  if (blacklist_su == "negativi") neg <- rimuovi_bl(neg)
  Y <- rbind(enhancers[, c("chr", "start", "end", "bed", "value")], neg)
  if (blacklist_su == "tutto") Y <- rimuovi_bl(Y)
  
  Y
}

for (i in 0:9) {
  message("Generating Y_vector", i, "_BL.csv ...")
  seed_val <- if (i == 0) 42 else i
  Yi <- costruisci_Y(shuffle_neg(seed_val))
  write.csv(Yi, file.path(out_dir, paste0("Y_vector", i, "_BL.csv")), row.names = FALSE)
}

message("Completed successfully.")
