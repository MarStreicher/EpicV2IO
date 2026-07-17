#!/usr/bin/env Rscript
# Build a compact Parquet subset of the Peters et al. EPICv2 manifest.
#
# Upstream:
#   Peters TJ, et al. (2024) BMC Genomics 25:251.
#   https://doi.org/10.1186/s12864-024-10027-5  (CC BY 4.0; Additional file 4)
#   Bioconductor EPICv2manifest / AnnotationHub AH116484 (Artistic-2.0)
#
# The output is a derived column subset for EpicV2IO only. See NOTICE.

library(AnnotationHub)
library(arrow)

MANIFEST_COLUMNS <- c(
  "IlmnID",
  "CHR",
  "MismatchPos",
  "MissingPos",
  "CH_BLAT",
  "CH_WGBS_evidence"
)
OUTPUT_PATH <- "src/epicv2io/resources/peters_epicv2_manifest.parquet"

cat("Querying AnnotationHub...\n")
ah <- AnnotationHub()
cat("Downloading AH116484 (Peters EPICv2 manifest)...\n")
manifest <- ah[["AH116484"]]
cat("Converting to data.frame...\n")
df <- as.data.frame(manifest)

missing_cols <- setdiff(MANIFEST_COLUMNS, colnames(df))
if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing expected columns in manifest:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

df <- df[, MANIFEST_COLUMNS]

cat("Writing parquet...\n")
write_parquet(df, OUTPUT_PATH)
cat("Done!\n")
cat("Rows:", nrow(df), "\n")
cat("Cols:", ncol(df), "\n")
cat("Saved to:", OUTPUT_PATH, "\n")
