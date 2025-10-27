# Find composition from mass with ppm tolerance
find_compositions_ppm <- function(target_mass, tolerance_ppm = 10,
                                  max_length = 8,    # search up to this many bases
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
  
  # Convert ppm tolerance to absolute mass tolerance (Da)
  tolerance <- target_mass * tolerance_ppm / 1e6
  
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
  }
  
  # Combine results
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    rownames(results) <- NULL
    results <- results[order(abs(results$mass_error_ppm)), , drop = FALSE]
  } else {
    results <- data.frame(
      length = integer(),
      A = integer(), C = integer(), G = integer(), U = integer(),
      Me = integer(), direction = character(),
      total_mass = numeric(), mass_error_Da = numeric(),
      mass_error_ppm = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}

# Find compositions in sequence
find_matching_compositions <- function(seq, df_comp) {
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