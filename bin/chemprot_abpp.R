#!/usr/bin/env Rscript
# Chemprot ABPP analysis (limma moderated t-test)

suppressPackageStartupMessages({
    library(tidyverse)
    library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: chemprot_abpp.R <msstats> <condition_col> <fdr_thresh> <out_prefix>")

msstats_file  <- args[1]
condition_col <- args[2]
fdr_thresh    <- as.numeric(args[3])
out_prefix    <- args[4]

quant <- read_tsv(msstats_file, show_col_types = FALSE)

wide <- quant %>%
    select(Protein, Sample, LogIntensity) %>%
    pivot_wider(names_from = Sample, values_from = LogIntensity) %>%
    column_to_rownames("Protein")

wide <- wide[apply(wide, 1, function(x) !all(is.na(x))), ]

sample_info <- quant %>% distinct(Sample, !!sym(condition_col))

# Auto-detect condition levels
design <- model.matrix(~ 0 + factor(sample_info[[condition_col]]))
colnames(design) <- levels(factor(sample_info[[condition_col]]))

# Build contrasts dynamically if exactly 2 levels
lvl <- levels(factor(sample_info[[condition_col]]))
if (length(lvl) == 2) {
    cont_matrix <- makeContrasts(contrast = design[,lvl[1]] - design[,lvl[2]], levels = design)
} else {
    cont_matrix <- makeContrasts(contrast = design[,1] - design[,2], levels = design)
}

fit  <- lmFit(wide, design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "none")
results <- results %>%
    rownames_to_column("Protein") %>%
    select(Protein, logFC, P.Value, adj.P.Val) %>%
    mutate(significant = adj.P.Val < fdr_thresh)

write_tsv(results, paste0(out_prefix, "_abpp_results.tsv"))
message("ABPP complete. Significant proteins: ", sum(results$significant))
