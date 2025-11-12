library(stringr)
library(dplyr)
library(tidyr)

# 1️⃣ Generate Fragments (|ACA, A|CA, |CA)
generate_modified_fragments <- function(seq_str, out_path = "Result/NIST_modified_MazF_cleavage_result.xlsx") {
  # Basic checks
  if (!is.character(seq_str) || length(seq_str) != 1) stop("seq_str must be a single character string.")
  if (grepl("[^ACGU]", seq_str)) stop("seq_str must contain only A/C/G/U (uppercase).")
  
  n <- nchar(seq_str)
  
  # Helper to get all (1-based) start indices of fixed pattern matches; returns integer(0) if none
  fixed_positions <- function(pattern, x) {
    pos <- gregexpr(pattern, x, fixed = TRUE)[[1]]
    if (identical(pos, -1L)) integer(0) else as.integer(pos)
  }
  
  # Boundary sets:
  # - Before "ACA"  → indices where "ACA" starts
  # - Before "CA"   → indices where "CA" starts
  pos_aca <- fixed_positions("ACA", seq_str)
  pos_ca  <- fixed_positions("CA",  seq_str)
  
  # Boundaries are 1 (5'-end), any rule-driven internal boundary, and n+1 (3'-end)
  boundaries <- sort(unique(c(1L, pos_aca, pos_ca, n + 1L)))
  
  # Label a boundary index for reporting (for start/end labels)
  boundary_label <- function(idx) {
    if (idx == 1L) return("5'-end")
    if (idx == n + 1L) return("3'-end")
    
    # Before ACA
    if (idx %in% pos_aca) return("|ACA")
    
    # Before CA (may or may not be A|CA depending on the previous base)
    if (idx %in% pos_ca) {
      if (idx > 1L && substr(seq_str, idx - 1L, idx - 1L) == "A") {
        return("A|CA")
      } else {
        return("|CA")
      }
    }
    # Fallback (should not happen with our current rule set)
    return("unknown")
  }
  
  # U counter
  count_U <- function(s) {
    if (nchar(s) == 0) return(0L)
    sum(strsplit(s, "", fixed = TRUE)[[1]] == "U")
  }
  
  # Enumerate all contiguous fragments between any two boundaries (i < j)
  res_list <- vector("list", length = (length(boundaries) * (length(boundaries) - 1)) %/% 2)
  k <- 0L
  for (i in seq_len(length(boundaries) - 1L)) {
    for (j in (i + 1L):length(boundaries)) {
      start_idx <- boundaries[i]
      end_idx   <- boundaries[j]
      frag_seq  <- if (end_idx - 1L >= start_idx) substr(seq_str, start_idx, end_idx - 1L) else ""
      
      k <- k + 1L
      res_list[[k]] <- list(
        Start_position   = start_idx,
        End_position = end_idx - 1L,
        Start_label   = boundary_label(start_idx),
        End_label     = boundary_label(end_idx),
        Fragment_size  = nchar(frag_seq),
        Fragment_sequence  = frag_seq,
        U_count       = count_U(frag_seq)
      )
    }
  }
  
  fragments_df <- do.call(rbind.data.frame, res_list)
  rownames(fragments_df) <- NULL
  
  # Write to xlsx
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required. Please install it with install.packages('writexl').")
  }
  writexl::write_xlsx(list(fragments = fragments_df), path = out_path)
  
  invisible(fragments_df)
}

# 2️⃣ Function to calculate theo mass for methylated U
add_fragment_mass_Me <- function(df) {
  
  # Base monoisotopic masses (RNA)
  base_mass <- c(
    A = 329.05252,
    C = 305.04129,
    G = 345.04743,
    U = 320.04095  # Methylated pseudo-U
  )
  
  # Get the fragment sequence column name
  frag_col <- grep("fragment.*sequence", names(df), ignore.case = TRUE, value = TRUE)
  if (length(frag_col) == 0) stop("No column matching 'Fragment sequence' found.")
  
  seqs <- as.character(df[[frag_col]])
  seqs <- trimws(seqs)
  
  compute_mass <- function(seq) {
    if (is.na(seq) || nchar(seq) == 0) return(NA_real_)
    seq <- toupper(seq)
    bases <- unlist(strsplit(seq, ""))
    
    if (!all(bases %in% names(base_mass))) {
      cat("Non-ACGU base found in:", seq, "\n")
      return(NA_real_)
    }
    
    sum(base_mass[bases])
  }
  
  # Compute masses
  df$theo_mass_modified <- sapply(seqs, compute_mass)
  
  # # Add chemical adjustments
  # df$mass <- df$mass + 97.9769 
  # if (nrow(df) > 1) {
  #   df$mass[-nrow(df)] <- df$mass[-nrow(df)] - 18.01
  # }
  # 
  return(df)
}


# 3️⃣ Match function (customizable column names)
match_by_ppm <- function(exp, df,
                         exp_mass_col = "monoisotopic_mass",
                         theo_mass_col = "theo_mass_modified",
                         intensity_col = "sum_intensity",
                         tol = 10) {
  fuzzy_inner_join(
    exp, df,
    by = setNames(theo_mass_col, exp_mass_col),
    match_fun = function(x, y) ppm_match(x, y, tol = tol)
  ) %>%
    group_by(.data[[exp_mass_col]]) %>%
    slice_max(order_by = .data[[intensity_col]], n = 1, with_ties = FALSE) %>%
    ungroup()
}

# 4️⃣ Unmatched finder (customizable column names)
find_unmatched <- function(df, df_matched,
                           theo_mass_col = "theo_mass_modified",
                           exp_mass_col = "monoisotopic_mass",
                           tol = 10) {
  fuzzy_join(
    df, df_matched,
    by = setNames(exp_mass_col, theo_mass_col),
    match_fun = function(x, y) ppm_match(x, y, tol = tol),
    mode = "anti"
  )
}
count_bases <- function(seq) {
  # Ensure input is a character string
  if (!is.character(seq)) stop("Input must be a character string.")
  
  # Convert to uppercase for consistency
  seq_upper <- toupper(seq)
  
  # Split into individual bases
  bases <- strsplit(seq_upper, "")[[1]]
  
  # Count occurrences of A, C, G, U
  counts <- table(factor(bases, levels = c("A", "C", "G", "U")))
  
  # Return as a named vector
  return(counts)
}

# 5️⃣ Mass reversion for modified
find_compositions_ppm_modified <- function(target_mass, tolerance_ppm = 10,
                                           max_length = 8,
                                           direction = "") {
  # Base monoisotopic masses (RNA) with U fully methylated
  base_mass <- c(
    A = 329.0525,
    C = 305.04129,
    G = 345.04743,
    U = 320.04095  # U is fully methylated
  )
  
  direction_mass <- switch(
    direction,
    "5"        =  97.9769,
    "3"        = -61.95579,
    "internal" =  18.01,
    0
  )
  
  tolerance <- target_mass * tolerance_ppm / 1e6
  results_list <- list(); counter <- 1
  
  for (len in 1:max_length) {
    for (A in 0:len) {
      for (C in 0:(len - A)) {
        for (G in 0:(len - A - C)) {
          U <- len - A - C - G
          total_mass <- A*base_mass["A"] + C*base_mass["C"] +
            G*base_mass["G"] + U*base_mass["U"] +
            direction_mass
          if (abs(total_mass - target_mass) <= tolerance) {
            results_list[[counter]] <- data.frame(
              length = len, A = A, C = C, G = G, U = U,
              direction = direction,
              total_mass = total_mass,
              mass_error_Da = total_mass - target_mass,
              mass_error_ppm = (total_mass - target_mass) / target_mass * 1e6,
              stringsAsFactors = FALSE
            )
            counter <- counter + 1
          }
        }
      }
    }
  }
  
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    rownames(results) <- NULL
    results <- results[order(abs(results$mass_error_ppm)), , drop = FALSE]
  } else {
    results <- data.frame(
      length = integer(),
      A = integer(), C = integer(), G = integer(), U = integer(),
      direction = character(),
      total_mass = numeric(), mass_error_Da = numeric(),
      mass_error_ppm = numeric(),
      stringsAsFactors = FALSE
    )
  }
  results
}


# 6️⃣ Find compositions in sequence
find_matching_compositions_modified <- function(seq, df_comp) {
  seq <- toupper(seq)
  seq_chars <- strsplit(seq, "")[[1]]
  n <- length(seq_chars)
  
  matches <- list()
  
  for (i in seq_len(nrow(df_comp))) {
    A_needed <- df_comp$A[i]
    C_needed <- df_comp$C[i]
    G_needed <- df_comp$G[i]
    U_needed <- df_comp$U[i]
    len <- df_comp$length[i]
    
    if (len > n) next
    
    # find all possible windows of the same length
    for (j in 1:(n - len + 1)) {
      frag <- seq_chars[j:(j + len - 1)]
      counts <- table(factor(frag, levels = c("A", "C", "G", "U")))
      
      if (counts["A"] == A_needed &&
          counts["C"] == C_needed &&
          counts["G"] == G_needed &&
          counts["U"] == U_needed) {
        matches[[length(matches) + 1]] <- cbind(
          df_comp[i, ],
          matched_seq = paste0(frag, collapse = ""),
          start = j,
          end = j + len - 1
        )
      }
    }
  }
  
  if (length(matches) > 0) {
    result <- do.call(rbind, matches)
    rownames(result) <- NULL
  } else {
    result <- data.frame()
  }
  
  return(result)
}

# 7️⃣ Match experimental masses to theoretical fragments by ±14 Da incrementss
ppm_diff <- function(observed, theoretical) {
  (observed - theoretical) / theoretical * 1e6
}

# ---- Core matcher ----
# exp_unmatched: data.frame with column `monoisotopic_mass`
# df: data.frame with at least `theo_mass_modified`; optional `n_U`, `sequence`, positions...
# tol_ppm: numeric ppm tolerance (default 10)
# k_global_range: used only if `n_U` not available in df; e.g., -6:6
# write_xlsx_path: optional path to write results
match_by_14_multiples <- function(exp_unmatched,
                                  df,
                                  tol_ppm = 10,
                                  k_global_range = -6:6,
                                  write_xlsx_path = "Result/HTC/Results_match_by_14.xlsx") {
  stopifnot("monoisotopic_mass" %in% names(exp_unmatched))
  stopifnot("theo_mass_modified" %in% names(df))
  
  df2 <- df %>%
    mutate(row_id = dplyr::row_number())
  
  # Build candidate k per fragment
  if ("n_U" %in% names(df2)) {
    # Per-fragment range informed by U count
    df_k <- df2 %>%
      mutate(k_min = -as.integer(n_U), k_max = as.integer(n_U)) %>%
      rowwise() %>%
      mutate(k_list = list(seq.int(k_min, k_max))) %>%
      ungroup() %>%
      select(-k_min, -k_max) %>%
      tidyr::unnest_longer(k_list, values_to = "k") %>%
      filter(abs(k) <= n_U)
  } else {
    # Global range if n_U is not available
    df_k <- df2 %>%
      mutate(k = list(k_global_range)) %>%
      tidyr::unnest_longer(k)
  }
  
  # Adjusted theoretical mass for each k (14 per methyl difference)
  df_k <- df_k %>%
    mutate(adjusted_theo_mass = theo_mass_modified + 14 * k)
  
  # Cross-join with experimental masses, then filter by ppm tolerance
  # (efficient enough for mid-size; for very large tables switch to data.table rolling join)
  exp2 <- exp_unmatched %>%
    mutate(exp_id = dplyr::row_number())
  
  # Cartesian join then filter (use dplyr >=1.1.0 syntax)
  candidates <- merge(exp2, df_k) %>%
    mutate(delta_da = monoisotopic_mass - adjusted_theo_mass,
           ppm = ppm_diff(monoisotopic_mass, adjusted_theo_mass),
           abs_ppm = abs(ppm)) %>%
    filter(abs_ppm <= tol_ppm)
  
  # If multiple hits per experimental mass, keep the best (smallest |ppm|), but
  # also return full list in case you want to inspect all
  best_hits <- candidates %>%
    group_by(exp_id) %>%
    slice_min(order_by = abs_ppm, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Still unmatched after allowing 14*k
  still_unmatched <- exp2 %>%
    anti_join(best_hits %>% select(exp_id), by = "exp_id")
  
  # Arrange columns nicely
  best_hits <- best_hits %>%
    select(
      exp_id, monoisotopic_mass,
      row_id, k,
      theo_mass_modified, adjusted_theo_mass,
      delta_da, ppm,
      dplyr::any_of(c("U_count","Fragment_sequence","Start_position","End_position","Fragment_size"))
    ) %>%
    arrange(abs(ppm))
  
  # Optionally write to xlsx
  if (!is.null(write_xlsx_path)) {
    writexl::write_xlsx(
      list(
        Matches = best_hits,
        All_candidates = res$all_candidates,
        Still_unmatched = still_unmatched %>% select(exp_id, monoisotopic_mass)
      ),
      write_xlsx_path
    )
  }
  
  list(
    matches = best_hits,
    all_candidates = candidates,          # full set of hits within tolerance (can be large)
    still_unmatched = still_unmatched
  )
}

