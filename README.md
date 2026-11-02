# Enhancer prediction at the Tbx1 locus with LASSO logistic regression

Predicts the probability that 14 regions at the *Tbx1* locus (mouse, mm10, chr16:18,576,567-18,618,243)
are enhancers. A LASSO logistic regression (`glmnet`) is trained on known mouse enhancers
(VISTA Enhancer Browser) versus random non-enhancer regions, using as features the overlap with a
collection of epigenomic peak files plus two histone-mark signals. The final score of each region
is the mean over the models trained on 10 sets of negatives.

## Repository layout

```
scripts/
  01_generate_Y_vector.R        builds the labelled region sets (Y vectors)
  02_predict_tbx1_enhancers.R   builds features, trains LASSO, predicts the 14 Tbx1 regions
data/
  sources.xlsx                  origin of every input file (accession, URL, license)
  raw/                          inputs of script 01 (VISTA fasta, liftOver bed)
  peaks/                        peak files = columns of X (not versioned, see sources.xlsx)
  bins/                         the two 2000-bp histone bin files
  Y_vectors_original/           Y vectors used for the reported results
```

## How the Y vector is generated (script 01)

1. **Positives.** Mouse elements are downloaded from the VISTA Enhancer Browser
   (https://enhancer.lbl.gov/, "both" positive and negative elements, mouse) and saved as FASTA.
   Name fields give chromosome, start, end and the VISTA label (positive/negative).
   Coordinates are lifted over from mm9 to mm10 with UCSC liftOver. Only elements labelled
   *positive* are kept as enhancers (value = 1).
2. **Negatives.** Random regions with the same lengths as the positives, drawn with
   `bedtools shuffle` and **excluding all VISTA elements** (positive and negative) and using
   `mm10.chrom.sizes` (value = 0). Each Y vector set uses a different draw of negatives
   (seeds 1-9 in the regenerated sets; the original sets were not seeded, see below).
3. **Blacklist.** Regions overlapping the ENCODE mm10 blacklist (v2) are removed from the
   **final Y vector** (positives and negatives), equivalent to `bedtools intersect -v`.
4. Output: `Y_vector{i}_BL.csv` with columns `chr, start, end, bed, value`.

## Features and model (script 02)

- **X, peak columns:** for each file in `data/peaks/`, the number of peaks overlapping each region
  (`countOverlaps`), files in alphabetical order.
- **X, histone columns:** for `H3K27me3` and `H3K4me1`, the sum of the score (column 5) of the
  2000-bp bins overlapping each region.
- **Model:** `cv.glmnet` (alpha = 1, binomial), 80/20 split with `caret::createDataPartition`,
  `set.seed(1)`, lambda = `lambda.min`.
- **Prediction:** probability of enhancer for the 14 Tbx1 regions (hard-coded in the script);
  final score = mean over the available iterations.

## Reproducing the results

Requirements: R (>= 4.x) with `Biostrings`, `GenomicRanges`, `IRanges`, `rtracklayer`, `data.table`,
`dplyr`, `tidyr`, `glmnet`, `caret`, `ROCR`; `bedtools` (script 01 only).

**Exactly (original Y vectors):** copy the original `Y_vector{i}_BL.csv` into
`data/Y_vectors_original/`, put the peak and bin files in `data/peaks/` and `data/bins/`
(see `data/sources.xlsx`), set `Y_source <- "original"` in script 02 and run it from the repository root.

**From scratch:** put the VISTA fasta and its liftOver output in `data/raw/`, run script 01, set
`Y_source <- "generated"` in script 02 and run it. Negatives are random, so values will be close
but not identical to the original ones.

## Data availability

- **In this repository:** code, the original Y vectors (small CSV, genomic coordinates only) and
  `data/sources.xlsx`.
- **Not in this repository:** third-party peak files, VISTA downloads and the blacklist. They are
  public; accessions and URLs are in `data/sources.xlsx`. Check each license before redistributing.


## Citation and license

The original reference paper is Aurigemma, I., Lanzetta, O., Cirino, A. et al. Endothelial gene regulatory elements associated with cardiopharyngeal lineage differentiation. Commun Biol 7, 351 (2024). https://doi.org/10.1038/s42003-024-06017-8. This work was also used in Lanzetta, O., Bilio, M., Liebig, J., Jechow, K., Ten, F. W., Ferrentino, R., Aurigemma, I., Illingworth, E., Conrad, C., Lukassen, S., Angelini, C., Baldini, A. Tbx1 stabilizes differentiation of the cardiopharyngeal mesoderm and drives morphogenesis in the pharyngeal apparatus. Development (2025) 152 (12): dev204524. https://doi.org/10.1242/dev.204524.

GEO accession numbers for the peaks are listed in the sources.xlsx file.
