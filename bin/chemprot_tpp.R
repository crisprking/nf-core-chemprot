#!/usr/bin/env Rscript
# Chemprot TPP analysis — robust sigmoid fitting with fallback

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: chemprot_tpp.R <msstats> <temps_csv> <min_r2> <min_pep> <out_prefix>")

msstats_file <- args[1]
temperatures <- as.numeric(strsplit(args[2], ",")[[1]])
min_r2       <- as.numeric(args[3])
min_peptides <- as.numeric(args[4])
out_prefix   <- args[5]

quant <- read_tsv(msstats_file, show_col_types = FALSE)
if ("NumPeptides" %in% colnames(quant)) {
    quant <- quant %>% filter(NumPeptides >= min_peptides)
}

fit_sigmoid <- function(df, temps) {
    bottom <- min(df$LogIntensity)
    top    <- max(df$LogIntensity)
    Tm     <- median(temps)
    slope  <- 1
    tryCatch({
        model <- nls(LogIntensity ~ bottom + (top - bottom) / (1 + exp((Tm - Temperature) / slope)),
                     data = df,
                     start = list(bottom = bottom, top = top, Tm = Tm, slope = slope),
                     control = nls.control(maxiter = 500, warnOnly = TRUE))
        pred <- predict(model)
        r2 <- 1 - sum((df$LogIntensity - pred)^2) / sum((df$LogIntensity - mean(df$LogIntensity))^2)
        list(Tm = coef(model)["Tm"], R2 = r2, success = TRUE)
    }, error = function(e) {
        list(Tm = NA, R2 = NA, success = FALSE)
    })
}

fits <- quant %>%
    group_by(Protein, Condition) %>%
    nest() %>%
    mutate(fit = map(data, ~ fit_sigmoid(.x, temperatures))) %>%
    mutate(Tm = map_dbl(fit, "Tm"),
           R2 = map_dbl(fit, "R2"),
           success = map_lgl(fit, "success"))

fits_good <- fits %>% filter(success, R2 >= min_r2)

if (nrow(fits_good) == 0) {
    write_tsv(tibble(Error = "No sigmoid fits met R2 threshold"),
              paste0(out_prefix, "_tpp_results.tsv"))
    message("No valid fits — check input data.")
    quit(status = 0)
}

tm_wide <- fits_good %>%
    select(Protein, Condition, Tm, R2) %>%
    pivot_wider(names_from = Condition, values_from = c(Tm, R2), names_sep = "_")

if (ncol(tm_wide) >= 3 && "Tm_Treated" %in% colnames(tm_wide) && "Tm_Control" %in% colnames(tm_wide)) {
    tm_wide <- tm_wide %>%
        mutate(delta_Tm = Tm_Treated - Tm_Control,
               significant = abs(delta_Tm) > 0.5 & !is.na(delta_Tm))
    tm_wide$p_value <- 2 * pnorm(-abs(tm_wide$delta_Tm) / sd(tm_wide$delta_Tm, na.rm = TRUE))
    tm_wide$padj <- p.adjust(tm_wide$p_value, method = "fdr")
} else {
    tm_wide$delta_Tm <- NA
    tm_wide$significant <- FALSE
    tm_wide$p_value <- NA
    tm_wide$padj <- NA
}

write_tsv(tm_wide, paste0(out_prefix, "_tpp_results.tsv"))
message("TPP complete: ", out_prefix, "_tpp_results.tsv")
