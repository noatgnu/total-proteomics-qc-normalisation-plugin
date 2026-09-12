library(SummarizedExperiment)
library(PRONE)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(svglite)
library(patchwork)

annotate  <- ggplot2::annotate
filter    <- dplyr::filter
select    <- dplyr::select
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
count     <- dplyr::count

detect_delimiter <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  if (ext == "csv") "," else "\t"
}

basename_of <- function(x) sub(".*[/\\\\]", "", x)

save_all_formats <- function(plot_obj, path_no_ext, width = 10, height = 6) {
  if (is.null(plot_obj)) return(invisible(NULL))
  if (inherits(plot_obj, "ggplot")) {
    ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = width, height = height, dpi = 300)
    ggsave(paste0(path_no_ext, ".png"), plot_obj, width = width, height = height, dpi = 300)
    ggsave(paste0(path_no_ext, ".svg"), plot_obj, width = width, height = height)
  } else {
    pdf(paste0(path_no_ext, ".pdf"), width = width, height = height); print(plot_obj); dev.off()
    png(paste0(path_no_ext, ".png"), width = width, height = height, units = "in", res = 300); print(plot_obj); dev.off()
    svglite::svglite(paste0(path_no_ext, ".svg"), width = width, height = height); print(plot_obj); dev.off()
  }
}

read_annotation <- function(annotation_file) {
  sep <- detect_delimiter(annotation_file)
  df <- read.table(annotation_file, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("Sample" %in% colnames(df)) || !("Condition" %in% colnames(df))) {
    stop("Annotation file must have Sample and Condition columns", call. = FALSE)
  }
  df$Sample <- trimws(df$Sample)
  df$Condition <- trimws(df$Condition)
  df
}

# Matches each annotation Sample against a pg_matrix column: exact string first, basename fallback.
match_sample_columns <- function(pg_colnames, annotation_df) {
  matched <- rep(NA_character_, nrow(annotation_df))
  for (i in seq_len(nrow(annotation_df))) {
    s <- annotation_df$Sample[i]
    if (s %in% pg_colnames) {
      matched[i] <- s
      next
    }
    hit <- pg_colnames[basename_of(pg_colnames) == basename_of(s)]
    if (length(hit) >= 1) matched[i] <- hit[1]
  }
  matched
}

run_qc_normalisation <- function(pg_matrix_file, stats_file, annotation_file, output_folder,
                                  min_unique_peptides = 2, contaminant_column = "Contaminant",
                                  cc_mapped_column = NULL) {

  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  for (sub in c("QC", "normalisation")) {
    dir.create(file.path(output_folder, sub), showWarnings = FALSE, recursive = TRUE)
  }

  message("Loading pg_matrix...")
  pg_sep <- detect_delimiter(pg_matrix_file)
  pg_raw <- read.table(pg_matrix_file, sep = pg_sep, header = TRUE,
                        na.strings = c("NA", "NaN", "N/A", "#VALUE!"),
                        check.names = FALSE, stringsAsFactors = FALSE)

  message("Loading annotation file...")
  annotation_df <- read_annotation(annotation_file)

  matched_cols <- match_sample_columns(colnames(pg_raw), annotation_df)
  if (any(is.na(matched_cols))) {
    stop(paste0("Could not match annotation samples to pg_matrix columns: ",
                paste(annotation_df$Sample[is.na(matched_cols)], collapse = ", ")), call. = FALSE)
  }
  sample_cols <- matched_cols
  group_levels <- unique(annotation_df$Condition)

  sample_annotation <- tibble::tibble(
    sample_name = sample_cols,
    group = factor(annotation_df$Condition, levels = group_levels),
    display = annotation_df$Condition
  )
  if ("BioReplicate" %in% colnames(annotation_df)) {
    sample_annotation$bioreplicate <- annotation_df$BioReplicate
  }
  message("Samples per group:")
  print(table(sample_annotation$group))
  readr::write_tsv(sample_annotation, file.path(output_folder, "sample_annotation.tsv"))

  # ---- optional contaminant filter (guarded) ----
  if (!is.null(contaminant_column) && contaminant_column != "" && contaminant_column %in% colnames(pg_raw)) {
    cflag <- stringr::str_trim(as.character(pg_raw[[contaminant_column]]))
    n_contam <- sum(cflag == "+", na.rm = TRUE)
    message(paste("Contaminants flagged '+':", n_contam))
    pg_nc <- pg_raw[is.na(cflag) | cflag != "+", ]
  } else {
    message("No contaminant column found - skipping contaminant filter")
    pg_nc <- pg_raw
  }

  # ---- minimum proteotypic peptide filter ----
  if (!("N.Proteotypic.Sequences" %in% colnames(pg_nc))) {
    stop("pg_matrix is missing the N.Proteotypic.Sequences column", call. = FALSE)
  }
  pg_clean <- pg_nc |> dplyr::filter(N.Proteotypic.Sequences >= min_unique_peptides)
  message(paste("Proteins before:", nrow(pg_nc), "after:", nrow(pg_clean)))

  # ---- optional curated-category annotation (guarded, never filtered on) ----
  gene_map <- data.frame(
    Protein.Group = pg_clean$Protein.Group,
    Genes = pg_clean$Genes,
    primary_gene = sub(";.*", "", pg_clean$Genes),
    stringsAsFactors = FALSE
  )
  cc_col <- if (!is.null(cc_mapped_column) && cc_mapped_column != "" && cc_mapped_column %in% colnames(pg_clean)) {
    cc_mapped_column
  } else {
    NA_character_
  }
  if (!is.na(cc_col)) {
    cc_val <- stringr::str_trim(as.character(pg_clean[[cc_col]]))
    gene_map$Category <- ifelse(!is.na(cc_val) & cc_val != "" & cc_val == gene_map$primary_gene, "CC.mapped", "None")
    message(paste("CC.mapped rows flagged:", sum(gene_map$Category == "CC.mapped")))
  } else {
    gene_map$Category <- "None"
  }

  # ---- log2 matrix ----
  int_mat <- pg_clean |> dplyr::select(dplyr::all_of(sample_cols)) |> as.matrix()
  mode(int_mat) <- "numeric"
  rownames(int_mat) <- pg_clean$Protein.Group
  int_mat[int_mat == 0] <- NA
  log2_mat <- log2(int_mat)

  annot <- as.data.frame(sample_annotation)
  rownames(annot) <- annot$sample_name

  readr::write_tsv(cbind(gene_map, as.data.frame(log2_mat)), file.path(output_folder, "log2_matrix.tsv"))

  # ---- run-level QC plots from the stats file (optional) ----
  if (!is.null(stats_file) && stats_file != "") {
    message("Loading stats file...")
    stats_sep <- detect_delimiter(stats_file)
    stats <- read.table(stats_file, sep = stats_sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

    stats$sample_match <- vapply(stats$File.Name, function(f) {
      hit <- sample_annotation$sample_name[basename_of(sample_annotation$sample_name) == basename_of(f)]
      if (length(hit) >= 1) hit[1] else NA_character_
    }, character(1))

    stats_q <- stats |> dplyr::filter(!is.na(sample_match)) |>
      dplyr::left_join(sample_annotation, by = c("sample_match" = "sample_name"))

    if (nrow(stats_q) == 0) {
      message("No stats rows matched annotation samples - skipping QC plots")
    } else {
      stats_q <- stats_q |>
        dplyr::arrange(match(sample_match, sample_annotation$sample_name)) |>
        dplyr::mutate(run_label = factor(sample_match, levels = sample_match))

      theme_prot <- function() {
        theme_bw(base_size = 11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 12),
                legend.position = "right")
      }

      p_proteins <- ggplot(stats_q, aes(run_label, Proteins.Identified, fill = group)) +
        geom_col(colour = "grey25", width = 0.75) +
        geom_text(aes(label = scales::comma(Proteins.Identified)), vjust = -0.3, size = 2.2, angle = 90, hjust = -0.05) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
        labs(title = "Protein Groups Identified per Run", x = NULL, y = "Protein Groups", fill = "Group") + theme_prot()

      p_precursors <- ggplot(stats_q, aes(run_label, Precursors.Identified, fill = group)) +
        geom_col(colour = "grey25", width = 0.75) +
        geom_text(aes(label = scales::comma(Precursors.Identified)), vjust = -0.3, size = 2.2, angle = 90, hjust = -0.05) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
        labs(title = "Precursors Identified per Run", x = NULL, y = "Precursors", fill = "Group") + theme_prot()

      p_missed <- ggplot(stats_q, aes(run_label, Average.Missed.Tryptic.Cleavages, colour = group, group = 1)) +
        geom_line(linewidth = 0.9, colour = "grey50") + geom_point(size = 2.6) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
        labs(title = "Average Missed Tryptic Cleavages per Run", x = NULL, y = "Avg Missed Cleavages", colour = "Group") + theme_prot()

      qc_panel <- (p_proteins / p_precursors / p_missed) +
        patchwork::plot_annotation(title = "DIA-NN Run-Level QC")
      save_all_formats(qc_panel, file.path(output_folder, "QC", "QC_run_summary"), width = 14, height = 16)
      save_all_formats(p_proteins,   file.path(output_folder, "QC", "QC_proteins_identified"),  width = 12, height = 5)
      save_all_formats(p_precursors, file.path(output_folder, "QC", "QC_precursors_identified"), width = 12, height = 5)
      save_all_formats(p_missed,     file.path(output_folder, "QC", "QC_missed_cleavages"),      width = 12, height = 5)
    }
  } else {
    message("No stats file provided - skipping run-level QC plots")
  }

  # ---- PRONE normalisation across all methods ----
  message("Running PRONE normalisation...")
  pg_for_prone <- pg_clean |> dplyr::select(Protein.Group, Genes, dplyr::all_of(sample_cols))
  meta_prone <- sample_annotation |>
    dplyr::select(sample_name, group) |>
    dplyr::rename(Sample = sample_name, Condition = group) |> as.data.frame()
  meta_prone$Column <- meta_prone$Sample

  se <- PRONE::load_data(
    data = pg_for_prone, md = meta_prone,
    protein_column = "Protein.Group", gene_column = "Genes",
    condition_column = "Condition"
  )

  all_methods <- c("Median", "Quantile", "VSN", "NormicsVSN")
  se_norm <- tryCatch(
    PRONE::normalize_se(se, methods = all_methods),
    error = function(e) {
      message(paste0("normalize_se failed with all methods (", e$message, "); retrying without NormicsVSN"))
      all_methods <<- setdiff(all_methods, "NormicsVSN")
      PRONE::normalize_se(se, methods = all_methods)
    }
  )

  pcv_plot <- tryCatch(PRONE::plot_intragroup_PCV(se_norm, ain = all_methods, condition = "Condition"),
                        error = function(e) { message(paste("PCV plot error:", e$message)); NULL })
  save_all_formats(pcv_plot, file.path(output_folder, "normalisation", "PRONE_intragroup_PCV"))

  corr_plot <- tryCatch(PRONE::plot_intragroup_correlation(se_norm, ain = all_methods, condition = "Condition"),
                         error = function(e) { message(paste("Correlation plot error:", e$message)); NULL })
  save_all_formats(corr_plot, file.path(output_folder, "normalisation", "PRONE_intragroup_correlation"))

  # ---- CV comparison across methods (the decision plot) ----
  cv_by_method <- purrr::map_dfr(all_methods, function(m) {
    mat <- as.matrix(SummarizedExperiment::assay(se_norm, m))[, sample_cols]
    rownames(mat) <- pg_clean$Protein.Group
    as.data.frame(mat) |>
      tibble::rownames_to_column("protein") |>
      tidyr::pivot_longer(-protein, names_to = "sample_name", values_to = "intensity") |>
      dplyr::left_join(annot |> dplyr::select(sample_name, group), by = "sample_name") |>
      dplyr::group_by(protein, group) |>
      dplyr::summarise(cv = 100 * sd(intensity, na.rm = TRUE) / mean(intensity, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(method = m)
  })

  method_colours <- setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(all_methods)),
    all_methods
  )

  p_cv <- ggplot(cv_by_method, aes(method, cv, fill = method)) +
    geom_violin(alpha = 0.7, colour = "grey30") +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
    geom_hline(yintercept = 20, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = method_colours) +
    scale_y_continuous(limits = c(0, NA)) +
    facet_wrap(~ group) +
    labs(title = "CV Comparison Across Normalisation Methods",
         subtitle = "Lower CV = tighter replicate agreement within group", x = NULL, y = "CV (%)") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  save_all_formats(p_cv, file.path(output_folder, "normalisation", "CV_comparison_all_methods"), width = 13, height = 7)

  message("Median CV per method (lower is better):")
  print(cv_by_method |> dplyr::group_by(method) |> dplyr::summarise(median_cv = median(cv, na.rm = TRUE)))

  for (m in all_methods) {
    mat <- as.matrix(SummarizedExperiment::assay(se_norm, m))[, sample_cols]
    rownames(mat) <- pg_clean$Protein.Group
    readr::write_tsv(cbind(gene_map, as.data.frame(mat)),
                      file.path(output_folder, paste0("normalized_", m, ".tsv")))
  }

  message("QC and normalisation complete.")
}

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (startsWith(arg, "--")) {
      key <- substring(arg, 3)
      if (i < length(args) && !startsWith(args[i + 1], "--")) {
        parsed[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        parsed[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  parsed
}

params <- parse_args(args)

pg_matrix_file <- params$pg_matrix_file
stats_file <- params$stats_file
annotation_file <- params$annotation_file
output_folder <- params$output_folder
min_unique_peptides <- ifelse(is.null(params$min_unique_peptides), 2, as.numeric(params$min_unique_peptides))
contaminant_column <- ifelse(is.null(params$contaminant_column), "Contaminant", params$contaminant_column)
cc_mapped_column <- params$cc_mapped_column

if (is.null(pg_matrix_file) || is.null(annotation_file) || is.null(output_folder)) {
  stop("Missing required arguments: pg_matrix_file, annotation_file, output_folder", call. = FALSE)
}

run_qc_normalisation(
  pg_matrix_file = pg_matrix_file,
  stats_file = stats_file,
  annotation_file = annotation_file,
  output_folder = output_folder,
  min_unique_peptides = min_unique_peptides,
  contaminant_column = contaminant_column,
  cc_mapped_column = cc_mapped_column
)
