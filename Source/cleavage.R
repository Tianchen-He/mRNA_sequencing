# Simulate RNA cleavage for different enzymes
# Supported enzymes (extendable): "MazF", "RNaseT1", "RNaseA"

cleavage <- function(seq_ref, enzyme = "MazF") {
  stopifnot(is.character(seq_ref), length(seq_ref) == 1)
  n <- nchar(seq_ref)
  
  if (n == 0) {
    return(data.frame(
      Cleavage_site = integer(0),
      Cleave_position = integer(0),
      Fragment_size = integer(0),
      Fragment_sequence = character(0)
    ))
  }
  
  # Define recognition motif and where to cut relative to motif
  enzyme_rules <- list(
    MazF = list(motif = "ACA", cut_before = TRUE, offset = 0),
    RNaseT1 = list(motif = "G", cut_before = FALSE, offset = 0),
    RNaseA = list(motif = c("U", "C"), cut_before = FALSE, offset = 0)
  )
  
  if (!enzyme %in% names(enzyme_rules)) {
    stop(paste("Unsupported enzyme:", enzyme))
  }
  
  rule <- enzyme_rules[[enzyme]]
  
  # Find cleavage motif positions
  if (length(rule$motif) == 1) {
    pos <- gregexpr(rule$motif, seq_ref, fixed = TRUE)[[1]]
  } else {
    # multiple motifs
    pos_list <- lapply(rule$motif, function(m) gregexpr(m, seq_ref, fixed = TRUE)[[1]])
    pos <- sort(unlist(pos_list))
  }
  
  # No cleavage motif found → single full-length fragment
  if (identical(pos, -1L)) {
    return(data.frame(
      Cleavage_site = 0L,
      Cleave_position = 0L,
      Fragment_size = n,
      Fragment_sequence = substr(seq_ref, 1, n),
      stringsAsFactors = FALSE
    ))
  }
  
  # Determine cleavage boundary
  cleave_positions <- if (rule$cut_before) {
    pos + rule$offset
  } else {
    pos + nchar(if (length(rule$motif)==1) rule$motif else substr(seq_ref, pos, pos)) + rule$offset
  }
  
  # Build fragments
  starts <- c(1L, cleave_positions)
  ends <- c(cleave_positions - 1L, n)
  keep <- (ends - starts + 1L) > 0L
  starts <- starts[keep]
  ends   <- ends[keep]
  cleave_pos_for_frag <- c(0L, cleave_positions)[keep]
  
  data.frame(
    Cleavage_site = seq_along(starts) - 1L,
    Cleave_position = replace(cleave_pos_for_frag, 1, cleave_pos_for_frag[2]),
    Fragment_size = ends - starts + 1L,
    Fragment_sequence = substring(seq_ref, starts, ends),
    stringsAsFactors = FALSE
  )
}