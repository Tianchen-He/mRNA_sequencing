library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(readxl)
library(writexl)
library(formattable)
library(plotly)
library(ggrepel) #avoid label stacking when using ggplot
library(viridisLite) #turbo 256 color scale
library(patchwork) #for table presentation

# build the theoretical ladder based on ladder type and sample type
build_theo = function(known_sequence, reference, ladder_5, sample_type){
  return_df = data.frame(matrix(ncol = 2))
  colnames(return_df) = c("base_name", "theoretical_mass")
  if(sample_type == "natural_RNA"){
    if(ladder_5 == TRUE){
      return_df[1,1] = "5'"
      return_df[1,2] = 97.9769
    } else {
      return_df[1,1] = "3'"
      return_df[1,2] = -61.95579
    }
    i = 1
    while(i <= nchar(known_sequence)) {
      next_base = substr(known_sequence, i, i)
      for(j in 1:nrow(reference)){
        if(next_base == reference[j,1]){
          return_df[nrow(return_df) + 1,2] = return_df[nrow(return_df),2] + reference[j,2]
          return_df[nrow(return_df),1] = next_base
          break
        }
      }
      i = i+1
    }
    return_df = return_df %>% 
      mutate(n_position = row_number() - 1)
    if(ladder_5){
      return_df[nrow(return_df), 2] = return_df[nrow(return_df), 2] - 79.97107
    } else {
      return_df[nrow(return_df), 2] = return_df[nrow(return_df), 2] + 79.97107
    }
    return(return_df)
  } else if(sample_type == "synthetic_RNA"){
    if(ladder_5 == TRUE){
      return_df[1,1] = "5'"
      return_df[1,2] = 18.015
    } else {
      return_df[1,1] = "3'"
      return_df[1,2] = -61.95579
    }
    i = 1
    while(i <= nchar(known_sequence)) {
      next_base = substr(known_sequence, i, i)
      for(j in 1:nrow(reference)){
        if(next_base == reference[j,1]){
          return_df[nrow(return_df) + 1,2] = return_df[nrow(return_df),2] + reference[j,2]
          return_df[nrow(return_df),1] = next_base
          break
        }
      }
      i = i+1
    }
    return_df = return_df %>% 
      mutate(n_position = row_number() - 1)
    if(ladder_5){
      return_df[nrow(return_df), 2] = return_df[nrow(return_df), 2] - 79.97107
    }
    return(return_df)
  } else if(sample_type == "peptide"){
    return_df[1,1] = "Water"
    return_df[1,2] = 18.015
    i = 1
    while(i <= nchar(known_sequence)) {
      next_base = substr(known_sequence, i, i)
      for(j in 1:nrow(reference)){
        if(next_base == reference[j,1]){
          return_df[nrow(return_df) + 1,2] = return_df[nrow(return_df),2] + reference[j,2]
          return_df[nrow(return_df),1] = next_base
          break
        }
      }
      i = i+1
    }
    return_df = return_df %>% 
      mutate(n_position = row_number() - 1)
    return(return_df)
  }
}

internal_builder = function(input_sequence, dictionary){ 
  
  full_fragments = NULL
  
  for(i in 2:nchar(input_sequence)){
    temp_theo = build_theo(substr(input_sequence, i, nchar(input_sequence)), dictionary, TRUE, "synthetic_RNA")
    temp_theo = temp_theo[-1,]
    temp_theo = temp_theo %>% 
      mutate(starting_position = i, ending_position = i+row_number()-1)
    full_fragments = rbind(full_fragments, temp_theo)
  }
  
  full_fragments = full_fragments %>% 
    filter(ending_position != nchar(input_sequence)) %>% 
    rename(internal_fragment_length = n_position,  base_at_ending_position = base_name)
  
  return(full_fragments)
  
}

