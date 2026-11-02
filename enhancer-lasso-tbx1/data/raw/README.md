Inputs for `scripts/01_generate_Y_vector.R`:
- `all_enhancer.fasta`  VISTA elements (mouse, positive and negative), mm9
- `ENHANCERS_mm10.bed`  the same elements after UCSC liftOver mm9 -> mm10
- `opposite_enhancers_senza_negativi.bed` (optional) original set of negative regions
The mm10 blacklist and `mm10.chrom.sizes` are downloaded automatically if missing.
