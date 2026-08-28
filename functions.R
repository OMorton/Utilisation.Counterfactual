
#  get_clean_formula_terms() ---------------------------------------------------

#   Prunes candidate BAM formula terms to prevent rank deficiency, matrix 
#   singularity, and mgcv model fitting crashes when running on subsets.

# Process:
#   1. Parametric Terms:
#      - Verifies terms have non-zero variance.
#      - Enforces cross-coverage for ordered factor threats against `used_f`
#        (drops threats nested exclusively within a single `used_f` level).
#      - Runs a QR matrix rank check (`model.matrix`) to detect and prune 
#        collinear parametric terms (e.g., `-1 + used_f` alias conflicts).

#   2. Smooth Terms:
#      - Non-'by' Smooths: Drops spatial terms (e.g., `s(Lat, Lon)`) or spatial 
#        tensors (`ti(Lat, Lon, year)`) if fewer than 5 unique coordinate pairs 
#        exist, ensuring sufficient spatial points for thin plate splines.
#      - 'by' Smooths: Enforces strict dependencies—ordered factor smooths are 
#        dropped if their corresponding parametric term was pruned; numeric 
#        interaction smooths are dropped if they contain no non-zero variance.

# Returns:
#   A vector of valid, clean formula term strings ready to be passed to 
#   as.formula() and mgcv::bam().

get_clean_formula_terms <- function(candidate_terms, data) {
  is_smooth    <- grepl("^(s|ti|te)\\(", candidate_terms)
  para_terms   <- candidate_terms[!is_smooth]
  smooth_terms <- candidate_terms[is_smooth]
  
  # cehck parametric coefs
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
    
    # Check ordered factor threats against used_f
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
  
  # QR rank check on remaining parametric terms via model.matrix
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
  
  # Pre-calculate spatial coordinate count
  has_spatial_coords <- all(c("Latitude", "Longitude") %in% names(data))
  if (has_spatial_coords) {
    n_spatial_pts <- nrow(unique(na.omit(data[, c("Latitude", "Longitude")])))
  } else {
    n_spatial_pts <- 0
  }
  
  # Check smooths
  valid_smooths <- character(0)
  
  for (sm in smooth_terms) {
    # Extract inner parameters
    inner_str  <- sub("^[^\\(]+\\((.*)\\)$", "\\1", sm)
    raw_args   <- unlist(strsplit(inner_str, ","))
    clean_args <- trimws(raw_args[!grepl("=", raw_args)]) # drop named params (bs='tp', k=20)
    
    
    if (!grepl("by\\s*=", sm)) {
      if (!all(clean_args %in% names(data))) next
      
      # Spatial check: >= 5 unique spatial points
      is_spatial_smooth <- any(c("Latitude", "Longitude") %in% clean_args)
      if (is_spatial_smooth && n_spatial_pts < 5) {
        next
      }
      
      # Verify all underlying variables in the smooth have variation (>1 unique value)
      vars_ok <- all(sapply(clean_args, function(v) length(unique(na.omit(data[[v]]))) > 1))
      if (vars_ok) {
        valid_smooths <- c(valid_smooths, sm)
      }
      next
    }
    
    # do the by smooths
    by_var <- trimws(sub(".*by\\s*=\\s*([^,)]+).*", "\\1", sm))
    
    if (by_var %in% names(data)) {
      val_vec <- na.omit(data[[by_var]])
      
      if (grepl("_ord$", by_var)) {
        # Ordered factor smooths MUST have their parametric counterpart survive Step 1
        if (by_var %in% valid_para) {
          valid_smooths <- c(valid_smooths, sm)
        }
      } else if (by_var == "used_f") {
        if (has_dual_used) valid_smooths <- c(valid_smooths, sm)
      } else {
        # Numeric interaction variables (inter_*_num)
        if (length(unique(val_vec)) > 1 && any(val_vec != 0)) {
          valid_smooths <- c(valid_smooths, sm)
        }
      }
    }
  }
  
  c(valid_para, valid_smooths)
}

# class_att_plt()---------------------------------------------------------------

class_att_plt <- function(draws = NULL, summary = NULL, common.y = FALSE) {
  
  p_att <- ggplot(draws, aes(x = year, y = att_ratio)) +
     geom_line(aes(group = .draw), alpha = 0.12, color = "darkred", linewidth = 0.3) +
     geom_ribbon(data = summary, 
                 aes(y = att_ratio_mn, ymin = att_ratio_lo, ymax = att_ratio_hi),
                 fill = NA, colour = "black", linetype = "dashed", linewidth = 0.5) +
     geom_line(data = summary, aes(year, att_ratio_mn), colour = "black", linewidth = 0.7) +
     geom_hline(yintercept = 1.0, linetype = "solid", color = "black", linewidth = 0.5) +
     coord_cartesian(ylim = c(0, max(draws$att_ratio)), expand = FALSE) +
     labs(x = "Year", y = "Index") +
     theme_bw() +
    theme(plot.margin = margin(t = 0, unit = "pt"))
  
  if (common.y == TRUE) {
    p_att <- p_att + coord_cartesian(ylim = c(0, 6), expand = F)
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
    scale_fill_manual(values = c("Used" = "darkred", "Not used" = "grey60")) +
    coord_cartesian(xlim = c(1970, max(aves_att_chained$year)), expand = FALSE) +
    scale_y_continuous(breaks = c(max_series_count), limits = c(0, max_series_count)) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.y = element_line(),
          legend.position = "none",panel.grid = element_blank(),
          plot.margin = margin(b = 0, t = 5, unit = "pt"))
  
  return(list(bar = p_bar, att = p_att))
}
