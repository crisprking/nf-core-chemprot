#!/usr/bin/env Rscript
# Chemprot Dose-Response analysis (4-parameter log-logistic)

suppressPackageStartupMessages({
    library(tidyverse)
    library(drc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: chemprot_dose.R <msstats> <conc> <min_r2> <min_pep> <out_prefix>")

msstats_file   <- args[1]
concentrations <- as.numeric(strsplit(args[2], ",")[[1]])
min_r2         <- as.numeric(args[3])
min_peptides   <- as.numeric(args[4])
out_prefix     <- args[5]

quant <- read_tsv(msstats_file, show_col_types = FALSE)
if ("NumPeptides" %in% colnames(quant)) quant <- quant %>% filter(NumPeptides >= min_peptides)

fit_drc <- function(df) {
    tryCatch({
        model <- drm(LogIntensity ~ Concentration, data = df,
                     fct = LL.4(names = c("Slope","Lower","Upper","EC50")),
                     control = drmc(errorm = FALSE, maxIt = 500))
        pred <- predict(model)
        r2 <- 1 - sum((df$LogIntensity - pred)^2) / sum((df$LogIntensity - mean(df$LogIntensity))^2)
        ec50 <- ED(model, 50, interval = "delta")
        tibble(EC50 = as.numeric(ec50[1]), EC50_lower = as.numeric(ec50[2]),
               EC50_upper = as.numeric(ec50[3]), Hill_slope = coef(model)[1], R2 = r2)
    }, error = function(e) tibble(EC50 = NA, EC50_lower = NA, EC50_upper = NA, Hill_slope = NA, R2 = NA))
}

results <- quant %>%
    group_by(Protein) %>%
    nest() %>%
    mutate(fit = map(data, fit_drc)) %>%
    unnest(fit) %>%
    filter(!is.na(EC50), R2 >= min_r2) %>%
    mutate(pEC50 = -log10(EC50))

write_tsv(results, paste0(out_prefix, "_dose_response.tsv"))
message("Dose-response complete: ", out_prefix, "_dose_response.tsv")
