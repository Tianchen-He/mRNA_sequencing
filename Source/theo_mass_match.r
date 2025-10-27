library(fuzzyjoin)
library(tidyr)

# 1️⃣ Function to calculate theo mass
add_fragment_mass <- function(df) {
  # Base monoisotopic masses (RNA)
  base_mass <- c(
    A = 329.05252,
    C = 305.04129,
    G = 345.04743,
    U = 306.02530
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
  df$theo_mass <- sapply(seqs, compute_mass)
  
  # # Add chemical adjustments
  # df$mass <- df$mass + 97.9769 
  # if (nrow(df) > 1) {
  #   df$mass[-nrow(df)] <- df$mass[-nrow(df)] - 18.01
  # }
  # 
  return(df)
}

# 2️⃣ Match function (customizable column names)
match_by_ppm <- function(exp, df,
                         exp_mass_col = "monoisotopic_mass",
                         theo_mass_col = "theo_mass",
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

# 3️⃣ Unmatched finder (customizable column names)
find_unmatched <- function(df, df_matched,
                           theo_mass_col = "theo_mass",
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