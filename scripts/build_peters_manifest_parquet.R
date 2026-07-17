#!/usr/bin/env Rscript

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
OUTPUT_PATH <- "src/packages/EpicV2IO/src/epicv2io/resources/peters_epicv2_manifest.parquet"

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