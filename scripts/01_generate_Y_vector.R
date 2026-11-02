# LEGACY version: histone columns computed with the original ifelse(keys %in% agg$bed, agg$x, 0),
# which does not align the values with the rows of the Y vector. Kept only to compare with older results.
# The main script is scripts/02_predict_tbx1_enhancers.R.

# 02_predict_tbx1_enhancers.R
# Trains a LASSO logistic regression (glmnet) on each Y vector set and predicts the enhancer
# probability of 14 regions at the Tbx1 locus (mm10, chr16). Final score = mean over iterations.
#
# Run from the repository root:  setwd("path/to/enhancer-lasso-tbx1")

library(dplyr)
library(GenomicRanges)
library(data.table)
library(Biostrings)
library(glmnet)
library(caret)
library(ROCR)

stopifnot("Run this script from the repository root" = dir.exists("scripts"))

# ---- PATHS (relative to the repository root) -------------------------------------
# "original"  -> data/Y_vectors_original  (the Y vectors used for the published results)
# "generated" -> data/Y_vectors_generated (produced by 01_generate_Y_vector.R)
Y_source  <- "original"

ml_dir    <- if (Y_source == "original") "data/Y_vectors_original" else "data/Y_vectors_generated"
peaks_dir <- "data/peaks"     # peak / ChIP-seq files (one column of X each); see data/sources.xlsx
bins_dir  <- "data/bins"      # H3K27me3 and H3K4me1 2000-bp bin files (score in column 5)
out_dir   <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 14 TBX1 regions defined directly (no external RData file)
enhancers <- data.frame(
  seqnames = "chr16",
  start = c(18576567, 18577384, 18578371, 18580095, 18581868, 18586388, 18587866, 18590110, 18591397, 18601232, 18603959, 18607904, 18610483, 18617351),
  end   = c(18577247, 18578222, 18579005, 18580974, 18582636, 18587203, 18588576, 18590935, 18592314, 18602126, 18604876, 18608777, 18611385, 18618243)
)
enhancers$bed <- paste(enhancers$seqnames, enhancers$start, enhancers$end, sep = "-")
y.tbx1_base <- makeGRangesFromDataFrame(enhancers, keep.extra.columns = TRUE)

predictions_list <- list()

for (i in 0:9) {
  y_path <- file.path(ml_dir, paste0("Y_vector", i, "_BL.csv"))
  if (!file.exists(y_path)) next
  
  message("Processing: Y_vector", i, "_BL.csv")
  
  Y_vector <- read.csv(y_path)
  if ("X" %in% colnames(Y_vector) || colnames(Y_vector)[1] == "") Y_vector <- Y_vector[, -1]
  if (!"bed" %in% colnames(Y_vector)) {
    chr_col <- intersect(c("seqnames", "chr"), colnames(Y_vector))[1]
    Y_vector$bed <- paste(Y_vector[[chr_col]], Y_vector$start, Y_vector$end, sep = "-")
  }
  y <- makeGRangesFromDataFrame(Y_vector)
  
  old_wd <- getwd()
  setwd(peaks_dir)
  filenames <- list.files(path = peaks_dir)
  
  peak_counts <- sapply(filenames, simplify = FALSE, function(file_name) {
    df <- fread(file_name)[, 1:3]
    setnames(df, c("seqnames", "chromStart", "chromEnd"))
    df1 <- makeGRangesFromDataFrame(df)
    res <- data.frame(countOverlaps(y, df1))
    colnames(res) <- gsub(".gz", "", file_name)
    res
  })
  setwd(old_wd)
  
  X <- do.call("cbind", peak_counts)
  
  file_h1 <- fread(file.path(bins_dir, "H3K27me3_Mast_wt_bin_2000_mod.bed.gz"))[, c(1, 2, 3, 5)]
  setnames(file_h1, c("seqnames", "chromStart", "chromEnd", "score"))
  file_h1$bed <- paste(file_h1$seqnames, file_h1$chromStart, file_h1$chromEnd, sep = "-")
  
  y_extra <- makeGRangesFromDataFrame(Y_vector, keep.extra.columns = TRUE)
  f1 <- makeGRangesFromDataFrame(file_h1, keep.extra.columns = TRUE)
  
  intersection1 <- findOverlaps(y_extra, f1, type = "any")
  df_int1 <- data.frame(y_extra[queryHits(intersection1), ], f1[subjectHits(intersection1), ])
  df_filter1 <- df_int1[, c(6, 14, 13)]
  match1 <- aggregate(df_filter1$score, by = list(bed = df_filter1$bed), FUN = sum)
  
  col_h1 <- ifelse(Y_vector$bed %in% match1$bed, match1$x, 0)
  col_h1 <- as.numeric(as.data.frame(col_h1)$col_h1)
  
  file_h2 <- fread(file.path(bins_dir, "H3K4me1_Mast_wt_bin_2000_mod.bed.gz"))[, c(1, 2, 3, 5)]
  setnames(file_h2, c("seqnames", "chromStart", "chromEnd", "score"))
  file_h2$bed <- paste(file_h2$seqnames, file_h2$chromStart, file_h2$chromEnd, sep = "-")
  
  f2 <- makeGRangesFromDataFrame(file_h2, keep.extra.columns = TRUE)
  intersection2 <- findOverlaps(y_extra, f2, type = "any")
  df_int2 <- data.frame(y_extra[queryHits(intersection2), ], f2[subjectHits(intersection2), ])
  df_filter2 <- df_int2[, c(6, 14, 13)]
  match2 <- aggregate(df_filter2$score, by = list(bed = df_filter2$bed), FUN = sum)
  
  col_h2 <- ifelse(Y_vector$bed %in% match2$bed, match2$x, 0)
  col_h2 <- as.numeric(as.data.frame(col_h2)$col_h2)
  
  X <- cbind(X, H3K27me3 = col_h1, H3K4me1 = col_h2)
  X <- as.matrix(X)
  
  Y <- Y_vector$value
  
  set.seed(1)
  inTrain <- createDataPartition(Y, p = 0.80, list = FALSE)
  X.train <- X[inTrain, ]
  Y.train <- Y[inTrain]
  
  cv.out <- cv.glmnet(X.train, Y.train, alpha = 1, family = "binomial")
  lasso.mod.train <- glmnet(X.train, Y.train, lambda = cv.out$lambda.min, alpha = 1, family = "binomial")
  
  old_wd <- getwd()
  setwd(peaks_dir)
  peak_counts_tbx1 <- sapply(filenames, simplify = FALSE, function(file_name) {
    df <- fread(file_name)[, 1:3]
    setnames(df, c("seqnames", "chromStart", "chromEnd"))
    df1 <- makeGRangesFromDataFrame(df)
    res <- data.frame(countOverlaps(y.tbx1_base, df1))
    colnames(res) <- gsub(".gz", "", file_name)
    res
  })
  setwd(old_wd)
  
  X.tbx1 <- do.call("cbind", peak_counts_tbx1)
  
  intersection_t1 <- findOverlaps(y.tbx1_base, f1, type = "any")
  df_t1 <- data.frame(y.tbx1_base[queryHits(intersection_t1), ], f1[subjectHits(intersection_t1), ])
  df_filter_t1 <- df_t1[, c(6, 12, 13)]
  match_t1 <- aggregate(df_filter_t1$score, by = list(bed = df_filter_t1$bed), FUN = sum)
  
  col_tbx1_h1 <- ifelse(enhancers$bed %in% match_t1$bed, match_t1$x, 0)
  col_tbx1_h1 <- as.numeric(as.data.frame(col_tbx1_h1)$col_tbx1_h1)
  
  intersection_t2 <- findOverlaps(y.tbx1_base, f2, type = "any")
  df_t2 <- data.frame(y.tbx1_base[queryHits(intersection_t2), ], f2[subjectHits(intersection_t2), ])
  df_filter_t2 <- df_t2[, c(6, 12, 13)]
  match_t2 <- aggregate(df_filter_t2$score, by = list(bed = df_filter_t2$bed), FUN = sum)
  
  col_tbx1_h2 <- ifelse(enhancers$bed %in% match_t2$bed, match_t2$x, 0)
  col_tbx1_h2 <- as.numeric(as.data.frame(col_tbx1_h2)$col_tbx1_h2)
  
  X.tbx1 <- cbind(X.tbx1, H3K27me3 = col_tbx1_h1, H3K4me1 = col_tbx1_h2)
  X.tbx1 <- as.matrix(X.tbx1)
  
  # Only change to the model: if the single-lambda refit does not converge and returns an empty
  # model, use the cv.glmnet fit at lambda.min. Otherwise this is the original code.
  n_coef <- sum(as.numeric(coef(lasso.mod.train))[-1] != 0)
  if (n_coef == 0) {
    warning("Iteration ", i, ": empty refit model, using cv.glmnet at lambda.min")
    Y.lasso.pred.tbx1 <- predict(cv.out, newx = X.tbx1, s = "lambda.min", type = "response")
  } else {
    Y.lasso.pred.tbx1 <- predict(lasso.mod.train, newx = X.tbx1, type = "response")
  }
  pred_df <- as.data.frame(Y.lasso.pred.tbx1)
  rownames(pred_df) <- enhancers$bed
  
  predictions_list[[paste0("iter_", i)]] <- pred_df
}

predictions_all <- do.call(cbind, predictions_list)
predictions_mean <- rowMeans(predictions_all)
print(predictions_mean)

saveRDS(predictions_all, file.path(out_dir, "predictions_all_iterations.rds"))
saveRDS(predictions_mean, file.path(out_dir, "predictions_mean.rds"))
