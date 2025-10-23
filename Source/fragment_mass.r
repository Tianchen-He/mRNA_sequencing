find_compositions <- function(target_mass, tolerance = 0.01,
                              max_length = 8,  # search up to this many bases
                              max_methyl = 0,
                              direction = "") {
  # Base monoisotopic masses (RNA)
  base_mass <- c(
    A = 329.05252,
    C = 305.04129,
    G = 345.04743,
    U = 306.02530
  )
  methyl_mass <- 14.01565
  
  # Directional correction
  direction_mass <- switch(
    direction,
    "5" =  97.9769,
    "3" = -61.95579,
    "internal" = 18.01,
    0
  )
  
  results_list <- list()
  counter <- 1
  
  for (len in 1:max_length) {
    for (A in 0:len) {
      for (C in 0:(len - A)) {
        for (G in 0:(len - A - C)) {
          U <- len - A - C - G
          
          for (Me in 0:min(max_methyl, len)) {
            total_mass <- A * base_mass["A"] +
              C * base_mass["C"] +
              G * base_mass["G"] +
              U * base_mass["U"] +
              Me * methyl_mass +
              direction_mass
            
            if (abs(total_mass - target_mass) <= tolerance) {
              results_list[[counter]] <- data.frame(
                length = len,
                A = A, C = C, G = G, U = U, Me = Me,
                direction = direction,
                total_mass = total_mass,
                mass_error = total_mass - target_mass,
                stringsAsFactors = FALSE
              )
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
  
  # Combine results
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    rownames(results) <- NULL
    results <- results[order(abs(results$mass_error)), , drop = FALSE]
  } else {
    results <- data.frame(
      length = integer(),
      A = integer(), C = integer(), G = integer(), U = integer(),
      Me = integer(), direction = character(),
      total_mass = numeric(), mass_error = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}

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
  df$mass <- sapply(seqs, compute_mass)
  
  # # Add chemical adjustments
  # df$mass <- df$mass + 97.9769 
  # if (nrow(df) > 1) {
  #   df$mass[-nrow(df)] <- df$mass[-nrow(df)] - 18.01
  # }
  # 
  return(df)
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