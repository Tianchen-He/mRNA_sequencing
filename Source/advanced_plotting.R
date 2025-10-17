################################################################################
# This function plots figures showing partial modification branching in JACS 2024 style. 
  # file_directory: the location of your file, has to be in adduct finder output format
  # ladder_type_chosen: your choice of ladder type to plot, exactly the same as the "ladder_type" column in your file
  # scale_shape_chosen: select the shape for your point for each ladder type, i.e., c("i6A 3'" = 21, "ms2i6A 3'" = 22, "native 3'" = 23)
  # x_breaks: the break on your x axis, i.e., c(2500, 5000, 7500, 10000, 12500)
  # x_limits: the limit for your x axis, i.e., c(11000, 14000)
  # y_breaks: the break on your y axis, i.e., c(12, 14, 16, 18)
  # y_limits: the limit for your y axis, i.e., c(12, 18)

plot_2024JACS_style = function(file_directory, ladder_type_chosen, scale_shape_chosen, x_breaks, x_limits, y_breaks, y_limits){
  
  # load your file in adduct finder ourput format
  temp = read_xlsx(file_directory) %>% 
    select(monoisotopic_mass, apex_rt, ladder_type, sum_intensity) %>% 
    drop_na() %>% 
    filter(ladder_type %in% ladder_type_chosen)
  
  ggplot(temp, aes(x = monoisotopic_mass, y = apex_rt)) +
    geom_line(aes(group = ladder_type), color = "black", size = 0.7) +
    geom_point(
      aes(fill = sum_intensity, shape = ladder_type),
      color = "black",
      stroke = 0.5,
      size = 4,
      alpha = 0.9
    ) +
    theme_minimal() +
    scale_fill_gradientn(colors = c("#add8e6", "#bfaee8", "#c9a3eb", "#d49aef", "#ba55d3")) +
    scale_shape_manual(values = scale_shape_chosen) + 
    labs(
      x = "Monoisotopic Mass (Da)",
      y = "Retention Time (min)"
    ) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = y_limits
    )

}

################################################################################
# This function plots chromatogram figures based on cvs file directly extracted from freestyle software. 
  # file_location: the location of your file, has to be in the format of freestyle
  # output cvs, note: that there will be three useless rows in the begining.

plot_chromatogram = function(file_location){
  options(scipen = 0)
  
  chromatogram_df = read_csv(file_location, skip = 3) %>% 
    janitor::clean_names() %>% 
    mutate(percentage = relative_abundance/max(relative_abundance)*100)
  
  ggplot(chromatogram_df, aes(x = time_min, y = percentage)) +
    geom_line(color = "steelblue", size = 0.75) +  
    labs(
      x = "Time (min)",
      y = "Relative abundance"
    ) +
    theme_classic(base_size = 14) +  
    theme(
      panel.grid = element_blank(),     
      axis.ticks = element_line(),        
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
}

################################################################################
#This function plots the depth of the sequencing result, the x-axis will be the position of nucleotide, the y-axis will be
#the number of depth sequenced at that particular position.
  #first variable: the df you want to draw based on, this should be in the syntax of output by function theo_creater.
  #second variable: the positions that you don't want to be included in this depth plot.(Ex. the intact position)
advanced_depth = function(file_location, positions_to_eliminate){
  temp = read_xlsx(file_location)  %>%
    add_count(position_5, name = "count_appearances") %>% 
    distinct(position_5, .keep_all = TRUE) %>% 
    select(position_5, count_appearances) %>% 
    arrange(position_5) %>% 
    filter(!position_5 %in% positions_to_eliminate)
  
  temp = temp %>%
    rename(depth = count_appearances)
  
  ggplot(temp, aes(x = position_5, y = depth)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.5) + 
    labs(
      x = "Position",
      y = "Depth"
    ) +
    scale_x_continuous(
      breaks = c(20 ,40, 60, 80, 100, 120)
    ) +
    scale_y_continuous(
      breaks = c(0, 2, 4, 6, 8, 10, 12)
    ) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80"),
      axis.text.x = element_text(size = 28, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 28, face = "bold", colour = "black"),
      axis.title.x = element_text(size = 28, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 28, face = "bold", colour = "black")
    )
}

plot_rt_clean = function(df, time_limit) { 
  ggplot(df, aes(x=monoisotopic_mass, y=apex_rt)) + 
    geom_point(
      color="#69b3a2",
      alpha=0.9,
      size=2
    ) +
    theme_classic(base_size = 14) +  
    theme(
      panel.grid = element_blank(),     
      axis.ticks = element_line(),        
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    ) +
    scale_y_continuous(limits = c(0, time_limit), breaks = c(0, 4, 8, 12, 16, 20, 24)) +
    scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000)) + 
    theme(plot.title = element_text(
      size = 12,        # Font size
      face = "bold",    # Font style: "bold", "italic", "bold.italic"
      hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
    )
    )
}