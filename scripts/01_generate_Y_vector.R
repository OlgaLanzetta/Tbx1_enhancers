# 01_generate_Y_vector.R
# Builds the labelled region sets (Y vectors): VISTA positives (mm10) + shuffled negatives,
# with ENCODE blacklist regions removed from the final Y vector.
# Output: data/Y_vectors_generated/Y_vector{0..9}_BL.csv
#
# Run from the repository root:  setwd("path/to/enhancer-lasso-tbx1")

library(Biostrings)
library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(IRanges)

stopifnot("Run this script from the repository root" = dir.exists("scripts"))

# ---- PATHS (relative to the repository root) -------------------------------------
raw_dir     <- "data/raw"                    # see data/sources.xlsx for where each file comes from
data_dir    <- "data/Y_vectors_generated"    # output
bed_dir     <- raw_dir
old_bed_dir <- raw_dir                       # optional: opposite_enhancers_senza_negativi.bed (original negatives)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Blacklist is removed from the FINAL Y vector ("tutto") or only from negatives ("negativi")
blacklist_su <- "tutto"

# ---- INPUT FILES -------------------------------------------------------------------
fasta_path  <- file.path(raw_dir, "all_enhancer.fasta")     # VISTA elements (mouse, positive + negative)
mm10_path   <- file.path(raw_dir, "ENHANCERS_mm10.bed")     # same elements after UCSC liftOver mm9 -> mm10
bl_path     <- file.path(raw_dir, "mm10-blacklist.v2.bed.gz")
chrom_sizes <- file.path(raw_dir, "mm10.chrom.sizes")

stopifnot("Missing data/raw/all_enhancer.fasta" = file.exists(fasta_path),
          "Missing data/raw/ENHANCERS_mm10.bed" = file.exists(mm10_path))

if (!file.exists(bl_path)) {
  download.file("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
                bl_path, mode = "wb")
}
if (!file.exists(chrom_sizes)) {
  download.file("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes", chrom_sizes)
}

# ---- COSTRUZIONE POSITIVI (identica al tuo script originale) -------------------
fastaFile3 <- readDNAStringSet(fasta_path)
seq_name3 <- as.data.frame(fastaFile3@ranges@NAMES)
seq_name_all <- data.frame("region" = seq_name3, "value" = 1, "seq" = paste(fastaFile3))
sep3 <- separate(seq_name_all, fastaFile3.ranges.NAMES, into = c("a", "chr", "start", "end", "element", "number", "tipo"))
all_df <- sep3[, -c(1, 5, 6)]
remove(seq_name3, seq_name_all, sep3)
all_df$bed <- paste(all_df$chr, all_df$start, all_df$end, sep = "-")
all_df <- all_df[, -5]
all_df <- all_df[-114, ]

all_mm10 <- fread(mm10_path)
all_mm10 <- all_mm10[, -c(4, 5, 6)]
colnames(all_mm10) <- c("chr", "start", "end")
all_mm10$bed <- paste(all_mm10$chr, all_mm10$start, all_mm10$end, sep = "-")
all_df <- all_df[, -5]

new <- cbind(all_mm10, all_df)
new <- new[, -c(5, 6, 7, 9)]
new$width <- new$end - new$start

new_positive <- new[new$tipo == "positive"]
bed2 <- new_positive |>
  transform(chromStart = ifelse(start > end, end, start),
            chromEnd   = ifelse(end < start, start, end))
bed3 <- GenomicRanges::makeGRangesFromDataFrame(bed2)
bed4 <- rtracklayer::export(bed3, format = "bed", ignore.strand = FALSE)

# come nell'originale, ma SENZA append = TRUE (rieseguendo lo script il file non si duplica)
pos_full <- file.path(data_dir, "enhancers_only_positive_mm10.bed")
write.table(bed4, pos_full, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# equivalente di: cut -f 1,2,3 ...
pos3 <- file.path(data_dir, "positive_enhancers_mm10.bed")
all3 <- file.path(data_dir, "all_enhancers_mm10.bed")
fwrite(fread(pos_full, select = 1:3), pos3, sep = "\t", col.names = FALSE)
fwrite(fread(mm10_path, select = 1:3), all3, sep = "\t", col.names = FALSE)

enhancers <- new[new$tipo == "positive", ]
enhancers <- enhancers[, -c(5, 6)]
enhancers$value <- 1
enhancers <- as.data.frame(enhancers)
cat("Positivi:", nrow(enhancers), "\n")

# ---- NEGATIVI ------------------------------------------------------------------
# Se esiste ancora il file dei negativi del tuo vecchio script lo uso per il Y_vector base.
opp_orig <- file.path(old_bed_dir, "opposite_enhancers_senza_negativi.bed")

shuffle_neg <- function(seed) {
  out <- tempfile(fileext = ".bed")
  system(sprintf("bedtools shuffle -seed %d -excl %s -i %s -g %s > %s",
                 seed, all3, pos3, chrom_sizes, out))
  n <- fread(out, select = 1:3)
  setnames(n, c("chr", "start", "end"))
  n
}

carica_neg_base <- function() {
  if (file.exists(opp_orig)) {
    message("Uso il file originale dei negativi: ", opp_orig)
    n <- fread(opp_orig, select = 1:3)
    setnames(n, c("chr", "start", "end"))
    n
  } else {
    message("File originale dei negativi non trovato: rigenero con seed 42")
    shuffle_neg(42)
  }
}

# ---- BLACKLIST (equivalente di bedtools intersect -v) ---------------------------
bl <- fread(bl_path, header = FALSE, select = 1:3)
setnames(bl, c("chr", "start", "end"))
bl_gr <- makeGRangesFromDataFrame(bl)

rimuovi_bl <- function(df) {
  ov <- findOverlaps(makeGRangesFromDataFrame(df), bl_gr, type = "any")
  if (length(ov) > 0) df <- df[-unique(queryHits(ov)), ]
  df
}

costruisci_Y <- function(neg) {
  neg <- as.data.frame(neg)
  neg$bed <- paste(neg$chr, neg$start, neg$end, sep = "-")
  neg$value <- 0
  neg <- neg[, c("chr", "start", "end", "bed", "value")]
  n_neg0 <- nrow(neg)

  if (blacklist_su == "negativi") neg <- rimuovi_bl(neg)
  Y <- rbind(enhancers[, c("chr", "start", "end", "bed", "value")], neg)
  if (blacklist_su == "tutto") Y <- rimuovi_bl(Y)

  cat(sprintf("Y_vector finale: prima della blacklist %d righe (%d positivi + %d negativi) -> dopo %d righe (positivi %d, negativi %d)\n",
              nrow(enhancers) + n_neg0, nrow(enhancers), n_neg0,
              nrow(Y), sum(Y$value == 1), sum(Y$value == 0)))
  cat("chiavi in comune positivi/negativi:", length(intersect(Y$bed[Y$value == 1], Y$bed[Y$value == 0])), "\n")
  Y
}

# ---- SCRITTURA ------------------------------------------------------------------
# Y_vector = base; Y_vector_neg1..9 = stessi positivi, negativi con seed diverso.
# (pipeline_percorsi_nuovi.R legge questi file se esistono)
Y0 <- costruisci_Y(carica_neg_base())
write.csv(Y0, file.path(data_dir, "Y_vector0_BL.csv"), row.names = FALSE)   # iterazione 0 = base

for (i in 1:9) {
  message("Y_vector", i, "_BL.csv")
  Yi <- costruisci_Y(shuffle_neg(i))
  write.csv(Yi, file.path(data_dir, paste0("Y_vector", i, "_BL.csv")), row.names = FALSE)
}
