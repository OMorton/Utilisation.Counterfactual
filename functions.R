
##################################
##    CONVENIENCE FUNCTIONS     ##
##################################

#  get_clean_formula_terms() ---------------------------------------------------

#   Prunes candidate BAM formula terms to prevent rank deficiency, matrix 
#   singularity, and mgcv model fitting crashes when running on subsets.

# Process:
#  1. Parametric Terms:
#     - Verifies terms exist in data and have non-zero variance (>1 unique value).
#     - Enforces cross-coverage for ordered factor threats against `used_f`
#       (drops threats nested exclusively within a single `used_f` level).
#     - Runs a QR matrix rank check (`model.matrix`) to detect and prune 
#       collinear parametric terms (e.g., alias conflicts with `-1 + used_f`).

#  2. Spatial Pre-checks:
#     - Calculates the total number of unique spatial coordinate pairs (`Latitude`, `Longitude`)
#       in the subset prior to smooth evaluation.

#  3. Smooth Terms (Safe Parsing & Marginality Enforcement):
#     - Safe Parsing: Safely extracts core covariates without breaking on nested function 
#       arguments (e.g., `d = c(2, 1)`, `bs = c('tp', 'cr')` inside spatial-temporal `ti()` terms).
#     - Random Effects (`bs = 're'`): Validates that random effect variables (e.g., `fyear`) 
#       exist as valid factors with >1 unique level.
#     - Non-'by' Smooths: Drops spatial splines (`s(Lat, Lon)`) and spatial-temporal 
#       tensors (`ti(Lat, Lon, year)`) if fewer than 5 unique coordinate pairs exist.
#     - 'by' Smooths & Marginality:
#         * Ordered factor smooths (`*_ord`) are dropped if their corresponding parametric term 
#           was pruned in Step 1.
#         * Numeric interaction smooths (`inter_*_num`) are pruned if they contain < 10 
#           non-zero observations OR if their parent ordered threat (`*_ord`) failed/was 
#           pruned in Step 1 (via `threat_lookup` mapping)

# Returns:
#   A vector of valid, clean formula term strings 

get_clean_formula_terms <- function(candidate_terms, data, threat_lookup = NULL) {
  is_smooth    <- grepl("^(s|ti|te)\\(", candidate_terms)
  para_terms   <- candidate_terms[!is_smooth]
  smooth_terms <- candidate_terms[is_smooth]
  
  # --- Step 1: Parametric Terms ---
  valid_para <- character(0)
  has_used_f <- "used_f" %in% names(data)
  
  if (has_used_f) {
    used_vec <- data$used_f
    used_levels <- unique(na.omit(used_vec))
    has_dual_used <- length(used_levels) > 1
  } else {
    has_dual_used <- FALSE
  }
  
  for (term in para_terms) {
    if (term %in% c("-1", "- 1")) {
      valid_para <- c(valid_para, "-1")
      next
    }
    
    if (!term %in% names(data)) next
    v <- na.omit(data[[term]])
    if (length(unique(v)) <= 1) next
    
    if (grepl("_ord$", term) && has_dual_used) {
      min_val <- min(v)
      has_threat_used0 <- any(v[used_vec == 0] > min_val, na.rm = TRUE)
      has_threat_used1 <- any(v[used_vec == 1] > min_val, na.rm = TRUE)
      
      if (has_threat_used0 && has_threat_used1) {
        valid_para <- c(valid_para, term)
      }
    } else {
      valid_para <- c(valid_para, term)
    }
  }
  
  # QR rank check on parametric terms
  if (length(valid_para) > 1) {
    test_form <- as.formula(paste("popchange ~", paste(valid_para, collapse = " + ")))
    X_para <- tryCatch(model.matrix(test_form, data = data), error = function(e) NULL)
    
    if (!is.null(X_para)) {
      qr_X <- qr(X_para, tol = 1e-7)
      if (qr_X$rank < ncol(X_para)) {
        bad_cols <- colnames(X_para)[qr_X$pivot[(qr_X$rank + 1):ncol(X_para)]]
        bad_vars <- unique(gsub("\\.L|\\.Q|\\.C|\\d+$", "", bad_cols))
        valid_para <- setdiff(valid_para, bad_vars)
      }
    }
  }
  
  # Spatial points check
  has_spatial_coords <- all(c("Latitude", "Longitude") %in% names(data))
  n_spatial_pts <- if (has_spatial_coords) nrow(unique(na.omit(data[, c("Latitude", "Longitude")]))) else 0
  
  # Map interaction terms back to their parent ord terms if lookup provided
  inter_to_ord <- list()
  if (!is.null(threat_lookup)) {
    for (t_info in threat_lookup) {
      inter_to_ord[[t_info$inter]] <- t_info$ord
    }
  }
  
  # --- Step 2: Smooth Terms ---
  valid_smooths <- character(0)
  
  for (sm in smooth_terms) {
    inner_str <- sub("^[^\\(]+\\((.*)\\)$", "\\1", sm)
    
    clean_str <- gsub("(bs|k|d|m|xt|sp|fx)\\s*=\\s*c\\([^\\)]+\\)", "", inner_str)
    clean_str <- gsub("(bs|k|d|m|xt|sp|fx)\\s*=\\s*[^,)]+", "", clean_str)
    clean_vars <- trimws(unlist(strsplit(clean_str, ",")))
    clean_vars <- clean_vars[clean_vars != "" & !grepl("by\\s*=", clean_vars)]
    
    # Random effects (bs = 're')
    if (grepl("bs\\s*=\\s*['\"]re['\"]", sm)) {
      if (all(clean_vars %in% names(data))) {
        fyear_var <- clean_vars[1]
        if (length(unique(na.omit(data[[fyear_var]]))) > 1) {
          valid_smooths <- c(valid_smooths, sm)
        }
      }
      next
    }
    
    # Non-'by' Smooths
    if (!grepl("by\\s*=", sm)) {
      if (!all(clean_vars %in% names(data))) next
      
      is_spatial_smooth <- any(c("Latitude", "Longitude") %in% clean_vars)
      if (is_spatial_smooth && n_spatial_pts < 5) next
      
      vars_ok <- all(sapply(clean_vars, function(v) length(unique(na.omit(data[[v]]))) > 1))
      if (vars_ok) valid_smooths <- c(valid_smooths, sm)
      next
    }
    
    # 'by' Interaction Smooths
    by_var <- trimws(sub(".*by\\s*=\\s*([^,)]+).*", "\\1", sm))
    
    if (by_var %in% names(data)) {
      val_vec <- na.omit(data[[by_var]])
      
      if (grepl("_ord$", by_var)) {
        # Main ordered factor smooths MUST have their parametric term survive Step 1
        if (by_var %in% valid_para) valid_smooths <- c(valid_smooths, sm)
      } else if (by_var == "used_f") {
        if (has_dual_used) valid_smooths <- c(valid_smooths, sm)
      } else {
        # Numeric interaction variables (inter_*_num)
        # MARGINALITY: Drop interaction smooth if parent _ord variable failed/was pruned!
        parent_ord <- inter_to_ord[[by_var]]
        
        if (!is.null(parent_ord) && !parent_ord %in% valid_para) {
          next # Skip this interaction smooth!
        }
        
        if (sum(val_vec != 0, na.rm = TRUE) >= 10) {
          valid_smooths <- c(valid_smooths, sm)
        }
      }
    }
  }
  
  c(valid_para, valid_smooths)
}

# class_att_plt()---------------------------------------------------------------

class_att_plt <- function(draws = NULL, summary = NULL, common.y = FALSE, taxa.cols = FALSE) {
  
  taxa <- unique(summary$taxa)
  
  taxa.col <- "darkred"
  blank.col <- "grey60"
  
  if (taxa.cols == T) {
    blank.col <- "grey95"
    if (taxa == "Actinopteri") {taxa.col <- "#2B5C70"}
    if (taxa == "Aves") {taxa.col <- "#D4A373"}
    if (taxa == "Elasmobranchii") {taxa.col <- "#4A8B99"}
    if (taxa == "Mammalia") {taxa.col <- "#7570B3"}
    if (taxa == "Reptilia") {taxa.col <- "#A08B76"}
  }
    
  p_att <- ggplot(draws, aes(x = year, y = att_ratio)) +
     geom_line(aes(group = .draw), alpha = 0.12, color = taxa.col, linewidth = 0.3) +
     geom_ribbon(data = summary, 
                 aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
                 fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
     geom_line(data = summary, aes(year, att_ratio_mn), colour = "black", linewidth = 0.7) +
     geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
     coord_cartesian(ylim = c(0, max(draws$att_ratio)),
                     xlim = c(1970, 2024), expand = FALSE) +
     labs(x = "Year", y = "Index") +
     theme_bw() +
    theme(plot.margin = margin(t = 0, unit = "pt"))
  
  if (common.y == TRUE) {
    p_att <- p_att + coord_cartesian(ylim = c(0, 4),xlim = c(1970, 2024), expand = F)
  }
  
  dat_counts <- mod_dat_full %>%
    filter(as.numeric(as.character(fyear)) >= 1969,
           Class == unique(summary$taxa)) %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    filter(year >= 1970) %>%
    group_by(year, used_f) %>%
    summarise(n_series = n_distinct(series), .groups = "drop") %>%
    mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
                               levels = c("Used", "Not used")))
  
  max_series_count <- dat_counts %>%
    group_by(year) %>%
    summarise(total_n = sum(n_series), .groups = "drop") %>%
    pull(total_n) %>%
    max()
  
  p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("Used" = taxa.col, "Not used" = blank.col)) +
    coord_cartesian(xlim = c(1970, 2024), expand = FALSE) +
    scale_y_continuous(breaks = c(max_series_count), limits = c(0, max_series_count)) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.y = element_line(),
          legend.position = "none",panel.grid = element_blank(),
          plot.margin = margin(b = 0, t = 5, unit = "pt"))
  
  return(list(bar = p_bar, att = p_att))
}

# region_att_plt()---------------------------------------------------------------

region_att_plt <- function(draws = NULL, summary = NULL, common.y = FALSE, reg.cols = T) {
  
  reg <- unique(summary$country)
  
  reg.col <- "darkred"
  blank.col <- "grey60"
  
  if (reg.cols == T) {
    blank.col <- "grey95"
    if (reg == "Africa") {reg.col <- "#C0392B"}
    if (reg == "Antarctic") {reg.col <- "#F39C12"}
    if (reg == "Asia") {reg.col <- "#00B4D8"}
    if (reg == "Europe") {reg.col <- "#388E3C"}
    if (reg == "North America") {reg.col <- "#8E44AD"}
    if (reg == "Oceania") {reg.col <- "#D81B60"}
    if (reg == "Latin America and Caribbean") {reg.col <- "#4E6E58"}
  }
  
  p_att <- ggplot(draws, aes(x = year, y = att_ratio)) +
    geom_line(aes(group = .draw), alpha = 0.12, color = reg.col, linewidth = 0.3) +
    geom_ribbon(data = summary, 
                aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
                fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
    geom_line(data = summary, aes(year, att_ratio_mn), colour = "black", linewidth = 0.7) +
    geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
    coord_cartesian(ylim = c(0, max(draws$att_ratio)),xlim = c(1970, 2024),
                    expand = FALSE) +
    scale_x_continuous(breaks = c(1980, 2000, 2020)) +
    labs(x = "Year", y = "Index") +
    theme_bw() +
    theme(plot.margin = margin(t = 0, r = 2, unit = "pt"))
  
  if (common.y == TRUE) {
    p_att <- p_att + coord_cartesian(ylim = c(0, 3.5),xlim = c(1970, 2024),  expand = F)
  }
  
  dat_counts <- mod_dat_full %>%
    filter(as.numeric(as.character(fyear)) >= 1969,
           Region == reg) %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    filter(year >= 1970) %>%
    group_by(year, used_f) %>%
    summarise(n_series = n_distinct(series), .groups = "drop") %>%
    mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
                               levels = c("Used", "Not used")))
  
  max_series_count <- dat_counts %>%
    group_by(year) %>%
    summarise(total_n = sum(n_series), .groups = "drop") %>%
    pull(total_n) %>%
    max()
  
  p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("Used" = reg.col, "Not used" = blank.col)) +
    coord_cartesian(xlim = c(1970, 2024), expand = FALSE) +
    scale_y_continuous(breaks = c(max_series_count), limits = c(0, max_series_count)) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.y = element_line(),
          legend.position = "none",panel.grid = element_blank(),
          plot.margin = margin(b = 0, t = 5, r = 2, unit = "pt"))
  
  return(list(bar = p_bar, att = p_att))
}

# country_att_plt()---------------------------------------------------------------

country_att_plt <- function(draws = NULL, summary = NULL, common.y = FALSE, reg.cols = T) {
  
  reg <- unique(summary$region)
  country <- unique(summary$country)
  target_iso <- unique(summary$iso)
  reg.col <- "darkred"
  blank.col <- "grey60"
  
  if (reg.cols == T) {
    blank.col <- "grey95"
    if (reg == "Africa") {reg.col <- "#C0392B"}
    if (reg == "Antarctic") {reg.col <- "#F39C12"}
    if (reg == "Asia") {reg.col <- "#00B4D8"}
    if (reg == "Europe") {reg.col <- "#388E3C"}
    if (reg == "North America") {reg.col <- "#8E44AD"}
    if (reg == "Oceania") {reg.col <- "#D81B60"}
    if (reg == "Latin America and Caribbean") {reg.col <- "#4E6E58"}
  }
  
  p_att <- ggplot(draws, aes(x = year, y = att_ratio)) +
    geom_line(aes(group = .draw), alpha = 0.12, color = reg.col, linewidth = 0.3) +
    geom_ribbon(data = summary, 
                aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
                fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
    geom_line(data = summary, aes(year, att_ratio_mn), colour = "black", linewidth = 0.7) +
    geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
    coord_cartesian(ylim = c(0, max(draws$att_ratio)),xlim = c(1970, 2024),
                    expand = FALSE) +
    scale_x_continuous(breaks = c(1980, 2000, 2020)) +
    labs(x = "Year", y = "Index") +
    theme_bw() +
    theme(plot.margin = margin(t = 0, r = 2, unit = "pt"))
  
  if (common.y == TRUE) {
    p_att <- p_att + coord_cartesian(ylim = c(0, 3.5),xlim = c(1970, 2024),  expand = F)
  }
  
  dat_counts <- mod_dat_full %>%
    filter(as.numeric(as.character(fyear)) >= 1969,
           Country == country) %>%
    mutate(year = as.numeric(as.character(fyear))) %>%
    filter(year >= 1970) %>%
    group_by(year, used_f) %>%
    summarise(n_series = n_distinct(series), .groups = "drop") %>%
    mutate(used_label = factor(ifelse(used_f == "1", "Used", "Not used"), 
                               levels = c("Used", "Not used")))
  
  max_series_count <- dat_counts %>%
    group_by(year) %>%
    summarise(total_n = sum(n_series), .groups = "drop") %>%
    pull(total_n) %>%
    max()
  
  p_bar <- ggplot(dat_counts, aes(x = year, y = n_series, fill = used_label)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("Used" = reg.col, "Not used" = blank.col)) +
    coord_cartesian(xlim = c(1970, 2024), expand = FALSE) +
    scale_y_continuous(breaks = c(max_series_count), limits = c(0, max_series_count)) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.y = element_line(),
          legend.position = "none",panel.grid = element_blank(),
          plot.margin = margin(b = 0, t = 5, r = 2, unit = "pt"))
  
  return(list(bar = p_bar, att = p_att))
}

# tidy formula tables() ----------------------------------------------------------
extract_model_info <- function(results_list) {
  map_dfr(results_list, function(x) {
    raw_formula <- paste(deparse(x$formula), collapse = " ")
    clean_formula <- gsub("\\s+", " ", raw_formula)
    tibble(formula = clean_formula, rho = round(x$rho_est, 3) )
  }, .id = "taxa")
}
