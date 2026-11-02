library(dplyr)
library(GenomicRanges)
library(data.table)
library(Biostrings)
library(glmnet)
library(caret)
library(ROCR)

ml_dir    <- "/media/olga/KINGSTON/MACHINE_LEARNING_ENHANCERS/github/data/Y_vector"
prova_dir <- "/media/olga/KINGSTON/MACHINE_LEARNING_ENHANCERS/bed_file/prova"
bins_dir  <- "/media/olga/KINGSTON/MACHINE_LEARNING_ENHANCERS/bed_file/bins_union"

# 14 regioni TBX1 definite direttamente senza file RData esterni
enhancers <- data.frame(
  seqnames = "chr16",
  start = c(18576567, 18577384, 18578371, 18580095, 18581868, 18586388, 18587866, 18590110, 18591397, 18601232, 18603959, 18607904, 18610483, 18617351),
  end   = c(18577247, 18578222, 18579005, 18580974, 18582636, 18587203, 18588576, 18590935, 18592314, 18602126, 18604876, 18608777, 18611385, 18618243)
)
enhancers$bed <- paste(enhancers$seqnames, enhancers$start, enhancers$end, sep = "-")
y.tbx1_base <- makeGRangesFromDataFrame(enhancers, keep.extra.columns = TRUE)

predizioni_list <- list()

for (i in 0:9) {
  y_path <- file.path(ml_dir, paste0("Y_vector", i, "_BL.csv"))
  if (!file.exists(y_path)) next
  
  message("Elaborazione: Y_vector", i, "_BL.csv")
  
  Y_vector <- read.csv(y_path)
  if ("X" %in% colnames(Y_vector) || colnames(Y_vector)[1] == "") Y_vector <- Y_vector[, -1]
  if (!"bed" %in% colnames(Y_vector)) {
    chr_col <- intersect(c("seqnames", "chr"), colnames(Y_vector))[1]
    Y_vector$bed <- paste(Y_vector[[chr_col]], Y_vector$start, Y_vector$end, sep = "-")
  }
  y <- makeGRangesFromDataFrame(Y_vector)
  
  old_wd <- getwd()
  setwd(prova_dir)
  filenames <- list.files(path = prova_dir)
  
  prova <- sapply(filenames, simplify = FALSE, function(file_name) {
    df <- fread(file_name)[, 1:3]
    setnames(df, c("seqnames", "chromStart", "chromEnd"))
    df1 <- makeGRangesFromDataFrame(df)
    res <- data.frame(countOverlaps(y, df1))
    colnames(res) <- gsub(".gz", "", file_name)
    res
  })
  setwd(old_wd)
  
  X <- do.call("cbind", prova)
  
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
  setwd(prova_dir)
  prova_tbx1 <- sapply(filenames, simplify = FALSE, function(file_name) {
    df <- fread(file_name)[, 1:3]
    setnames(df, c("seqnames", "chromStart", "chromEnd"))
    df1 <- makeGRangesFromDataFrame(df)
    res <- data.frame(countOverlaps(y.tbx1_base, df1))
    colnames(res) <- gsub(".gz", "", file_name)
    res
  })
  setwd(old_wd)
  
  X.tbx1 <- do.call("cbind", prova_tbx1)
  
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
  
  Y.lasso.pred.tbx1 <- predict(lasso.mod.train, newx = X.tbx1, type = "response")
  pred_df <- as.data.frame(Y.lasso.pred.tbx1)
  rownames(pred_df) <- enhancers$bed
  
  predizioni_list[[paste0("iter_", i)]] <- pred_df
}

predizioni_df <- do.call(cbind, predizioni_list)
predizioni_medie <- rowMeans(predizioni_df)
print(predizioni_medie)