find_compositions <- function(length, target_mass, tolerance = 0.01,
                              max_methyl = 0, direction = "") {
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
    0
  )
  
  results_list <- list()
  counter <- 1
  
  for (A in 0:length) {
    for (C in 0:(length - A)) {
      for (G in 0:(length - A - C)) {
        U <- length - A - C - G
        
        for (Me in 0:min(max_methyl, length)) {
          total_mass <- A * base_mass["A"] + 
            C * base_mass["C"] +
            G * base_mass["G"] + 
            U * base_mass["U"] +
            Me * methyl_mass +
            direction_mass
          
          if (abs(total_mass - target_mass) <= tolerance) {
            results_list[[counter]] <- data.frame(
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
  
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    # 🧹 remove any row names completely
    rownames(results) <- NULL
    # sort by absolute mass error
    results <- results[order(abs(results$mass_error)), , drop = FALSE]
  } else {
    results <- data.frame(
      A = integer(), C = integer(), G = integer(), U = integer(),
      Me = integer(), direction = character(),
      total_mass = numeric(), mass_error = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}
