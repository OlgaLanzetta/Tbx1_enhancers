# 02_predict_tbx1_enhancers.R
# ------------------------------------------------------------------------------------------
# For each Y vector set (Y_vector{i}_BL.csv):
#   1. builds the feature matrix X for every labelled region
#   2. trains a LASSO logistic regression (glmnet, alpha = 1) with an 80/20 train/test split
#   3. predicts the enhancer probability of 14 regions at the Tbx1 locus (mm10, chr16)
# Final score of each region = mean of the predictions over the available iterations.
#
# Features (columns of X):
#   - one column per peak / ChIP-seq file in data/peaks/ : number of peaks overlapping the region
#   - two histone columns (H3K27me3, H3K4me1)             : sum of the 2000-bp bin scores overlapping the region
#
# Run from the repository root:  setwd("path/to/enhancer-lasso-tbx1")
# ------------------------------------------------------------------------------------------

library(GenomicRanges)
library(data.table)
library(glmnet)
library(caret)

stopifnot("Run this script from the repository root" = dir.exists("scripts"))

# ---- PATHS (relative to the repository root) ---------------------------------------------
# Y_source = "original"  -> data/Y_vectors_original  (Y vectors already prepared and used for the results)
#            "generated" -> data/Y_vectors_generated (produced by scripts/01_generate_Y_vector.R)
Y_source <- "original"

y_dir     <- if (Y_source == "original") "data/Y_vectors_original" else "data/Y_vectors_generated"
peaks_dir <- "data/peaks"     # peak / ChIP-seq files, one column of X each (see data/sources.xlsx)
bins_dir  <- "data/bins"      # H3K27me3 and H3K4me1 2000-bp bin files (score in column 5)
out_dir   <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- REGIONS TO PREDICT: 14 regions at the Tbx1 locus (mm10) -----------------------------
tbx1_regions <- data.frame(
  seqnames = "chr16",
  start = c(18576567, 18577384, 18578371, 18580095, 18581868, 18586388, 18587866, 18590110, 18591397, 18601232, 18603959, 18607904, 18610483, 18617351),
  end   = c(18577247, 18578222, 18579005, 18580974, 18582636, 18587203, 18588576, 18590935, 18592314, 18602126, 18604876, 18608777, 18611385, 18618243)
)
tbx1_regions$bed <- paste(tbx1_regions$seqnames, tbx1_regions$start, tbx1_regions$end, sep = "-")
tbx1_gr <- makeGRangesFromDataFrame(tbx1_regions, keep.extra.columns = TRUE)

# ---- FEATURE FILES (each file is read only once) -------------------------------------------
# Peak files: only the first three columns (chr, start, end) are used.
peak_files <- list.files(path = peaks_dir, pattern = "\\.(bed|narrowPeak|broadPeak|txt)(\\.gz)?$", full.names = TRUE)
if (length(peak_files) == 0) stop("No peak files found in: ", peaks_dir)
peak_names <- gsub(".gz", "", basename(peak_files))

peak_gr <- lapply(peak_files, function(f) {
  df <- fread(f)[, 1:3]
  setnames(df, c("seqnames", "chromStart", "chromEnd"))
  makeGRangesFromDataFrame(df)
})
names(peak_gr) <- peak_names

# Histone bin files: columns 1-3 = chr, start, end; column 5 = score (signal in the 2000-bp bin).
hist_paths <- list(
  H3K27me3 = file.path(bins_dir, "H3K27me3_Mast_wt_bin_2000_mod.bed.gz"),
  H3K4me1  = file.path(bins_dir, "H3K4me1_Mast_wt_bin_2000_mod.bed.gz")
)

hist_gr <- lapply(hist_paths, function(p) {
  dt_raw  <- fread(p)
  file_dt <- dt_raw[, c(1, 2, 3, 5)]
  setnames(file_dt, c("seqnames", "chromStart", "chromEnd", "score"))
  file_dt$bed <- paste(file_dt$seqnames, file_dt$chromStart, file_dt$chromEnd, sep = "-")
  makeGRangesFromDataFrame(file_dt, keep.extra.columns = TRUE)
})

# ---- FEATURE MATRIX ---------------------------------------------------------------------------
# regions_gr: GRanges of the regions. Returns a matrix with one row per region and one column per feature.
build_X <- function(regions_gr) {
  # peak columns: number of overlaps with each peak file
  X <- sapply(peak_gr, function(g) countOverlaps(regions_gr, g))
  X <- matrix(X, nrow = length(regions_gr), dimnames = list(NULL, peak_names))

  # histone columns: sum of the scores of all bins overlapping each region (0 if none)
  for (h_name in names(hist_gr)) {
    f_gr <- hist_gr[[h_name]]
    hits <- findOverlaps(regions_gr, f_gr, type = "any")
    score_col <- rep(0, length(regions_gr))
    if (length(hits) > 0) {
      agg <- tapply(f_gr$score[subjectHits(hits)], queryHits(hits), sum)   # sum per region index
      score_col[as.integer(names(agg))] <- as.numeric(agg)
    }
    X <- cbind(X, score_col)
    colnames(X)[ncol(X)] <- h_name
  }
  as.matrix(X)
}

# ---- ONE ITERATION (one Y vector set) ------------------------------------------------------------
run_iteration <- function(y_vector_path) {
  Y_vector <- read.csv(y_vector_path)
  if ("X" %in% colnames(Y_vector)) Y_vector <- Y_vector[, -1]            # drop row-index column, if any
  if (!"chr" %in% colnames(Y_vector) && "seqnames" %in% colnames(Y_vector)) {
    colnames(Y_vector)[colnames(Y_vector) == "seqnames"] <- "chr"
  }
  if (!"bed" %in% colnames(Y_vector)) {
    Y_vector$bed <- paste(Y_vector$chr, Y_vector$start, Y_vector$end, sep = "-")
  }
  y_gr <- makeGRangesFromDataFrame(Y_vector, keep.extra.columns = TRUE)

  X <- build_X(y_gr)
  Y <- Y_vector$value                                                      # 1 = enhancer, 0 = not enhancer

  # 80/20 split (fixed seed for reproducibility)
  set.seed(1)
  in_train <- createDataPartition(Y, p = 0.80, list = FALSE)
  X_train  <- X[in_train, ]
  Y_train  <- Y[in_train]

  # LASSO logistic regression: lambda chosen by cross-validation, then the model is refitted at lambda.min
  cv_fit    <- cv.glmnet(X_train, Y_train, alpha = 1, family = "binomial")
  lasso_fit <- glmnet(X_train, Y_train, lambda = cv_fit$lambda.min, alpha = 1, family = "binomial")

  # features of the Tbx1 regions (same columns, same order as the training matrix)
  X_tbx1 <- build_X(tbx1_gr)
  stopifnot(identical(colnames(X_tbx1), colnames(X)))

  # If the single-lambda refit does not converge, glmnet returns an empty model (AUC = 0.5).
  # In that case the cv.glmnet fit at lambda.min is used instead; otherwise this is the plain refit.
  n_nonzero <- sum(as.numeric(coef(lasso_fit))[-1] != 0)
  if (n_nonzero == 0) {
    warning("Empty refit model in ", basename(y_vector_path), ": using cv.glmnet at lambda.min")
    auc_test  <- assess.glmnet(cv_fit, newx = X[-in_train, ], newy = Y[-in_train], s = "lambda.min")$auc
    pred_tbx1 <- predict(cv_fit, newx = X_tbx1, s = "lambda.min", type = "response")
  } else {
    auc_test  <- assess.glmnet(lasso_fit, newx = X[-in_train, ], newy = Y[-in_train])$auc
    pred_tbx1 <- predict(lasso_fit, newx = X_tbx1, type = "response")
  }
  message("Test AUC: ", round(auc_test, 4))

  pred_df <- as.data.frame(pred_tbx1)
  rownames(pred_df) <- tbx1_regions$bed

  list(pred = pred_df, auc = auc_test)
}

# ---- RUN ALL Y VECTOR SETS ----------------------------------------------------------------------------
# Missing files (e.g. Y_vector0_BL.csv in the original set) are skipped.
y_files <- file.path(y_dir, paste0("Y_vector", 0:9, "_BL.csv"))
predictions_list <- list()

for (i in seq_along(y_files)) {
  if (file.exists(y_files[i])) {
    message("Processing: ", basename(y_files[i]))
    res <- run_iteration(y_files[i])
    message("Predicted regions: ", nrow(res$pred))
    predictions_list[[i]] <- res$pred
  } else {
    message("File not found: ", y_files[i])
  }
}

predictions_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, predictions_list)

# ---- FINAL SCORE: mean over iterations ---------------------------------------------------------------------
if (length(predictions_list) > 0) {
  predictions_all  <- do.call(cbind, predictions_list)    # one column per iteration
  predictions_mean <- rowMeans(predictions_all)           # final enhancer probability of each region
  saveRDS(predictions_mean, file.path(out_dir, "predictions_mean.rds"))
  saveRDS(predictions_all,  file.path(out_dir, "predictions_all_iterations.rds"))
  print(predictions_mean)
  message("Done.")
}
