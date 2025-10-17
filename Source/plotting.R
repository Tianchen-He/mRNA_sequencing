plot_mass_rt = function(df, label_it, time_limit){
  if(label_it){
    ggplot(df, aes(x=monoisotopic_mass, y=apex_rt)) + 
      geom_point(
        color="#69b3a2",
        alpha=0.9,
        size=2
      ) +
      theme_bw() +
      labs(
        title = paste("2D Mass Vs. RT plot for", file_location),
        x = "Monoisotopic Mass(Da)",
        y = "Retention Time(min)"
      ) +
      scale_y_continuous(limits = c(0, time_limit), breaks = seq(0, time_limit, by = 2)) +
      geom_text_repel(aes(label = as.character(round(monoisotopic_mass,2))),             
                      vjust = -1,                         
                      size = 3,                           
                      color = "black") +
      theme(plot.title = element_text(
        size = 12,        # Font size
        face = "bold",    # Font style: "bold", "italic", "bold.italic"
        hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
      )
      )
  } else {
    ggplot(df, aes(x=monoisotopic_mass, y=apex_rt)) + 
      geom_point(
        color="#69b3a2",
        alpha=0.9,
        size=2
      ) +
      theme_bw() +
      labs(
        title = paste("2D Mass Vs. RT plot for", file_location),
        x = "Monoisotopic Mass(Da)",
        y = "Retention Time(min)"
      ) +
      scale_y_continuous(limits = c(0, time_limit), breaks = seq(0, time_limit, by = 2)) +
      theme(plot.title = element_text(
        size = 12,        # Font size
        face = "bold",    # Font style: "bold", "italic", "bold.italic"
        hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
      )
      )
  }
}

plot_mass_intensity = function(df, label_it){
  if(label_it){
    ggplot(df, aes(x=monoisotopic_mass, y=sum_intensity)) + 
      geom_point(
        color="orange",
        alpha=0.9,
        size=2
      ) +
      theme_bw() +
      labs(
        title = paste("2D Mass Vs. Intensity plot for", file_location),
        x = "Monoisotopic Mass(Da)",
        y = "Sum Intensity"
      ) +
      geom_text_repel(aes(label = as.character(round(monoisotopic_mass,2))),             
                      vjust = -1,                         
                      size = 3,                           
                      color = "black") +
      theme(plot.title = element_text(
        size = 12,        
        face = "bold",    
        hjust = 0.5      
      )
      )
  } else {
    ggplot(df, aes(x=monoisotopic_mass, y=sum_intensity)) + 
      geom_point(
        color="orange",
        alpha=0.9,
        size=2
      ) +
      theme_bw() +
      labs(
        title = paste("2D Mass Vs. Intensity plot for", file_location),
        x = "Monoisotopic Mass(Da)",
        y = "Sum Intensity"
      ) +
      theme(plot.title = element_text(
        size = 12,        
        face = "bold",    
        hjust = 0.5      
      )
      )
  }
}

mass_match_table = function(prophet_df, reference_ladder_type){
  print(paste("The table below is the result of mass matching between the dataset under analysis and reference ladder", reference_ladder_type))
  formattable(prophet_df %>% 
                rename(observed_mass = monoisotopic_mass) %>% 
                mutate(ppm = abs((theoretical_mass - observed_mass)/theoretical_mass * 1000000)) %>% 
                select(n_position, base_name, observed_mass, theoretical_mass, ppm, apex_rt, sum_intensity), 
              align =c("l","c","c","c","c", "r"), 
              list(sum_intensity= color_tile("#DeF7E9", "#71CA97"),
                   ppm= color_tile("#ff7f7f", "red")))
}

plot_prophet = function(prophet_df, reference_ladder_type){
  ggplot(prophet_df, aes(x = monoisotopic_mass, y = apex_rt, fill = sum_intensity)) +
    geom_line(color = "gray", size = 1) +
    geom_point(color="gray",
               shape=21,
               alpha=0.9,
               size=3) +
    theme_bw() +
    labs(
      title = paste("2D Mass Matching Result for Dataset Under Analysis with reference ladder", reference_ladder_type),
      x = "Monoisotopic Mass(Da)",
      y = "Rentention Time(min)"
    ) +
    theme(plot.title = element_text(
      size = 12,        
      face = "bold",    
      hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
    )
    ) +
    geom_text(
      aes(label = base_name),             
      vjust = -1,                        
      size = 3,                          
      color = "black"                     
    ) + 
    scale_fill_gradientn(colors = turbo(256))
}

plot_position_abundance = function(prophet_df_1, prophet_df_2, ladder_name_1, ladder_name_2){
  print("the graph below shows the relative abundance distribution of matched nucleotides:")
  if(nrow(prophet_df_1) != 0 & nrow(prophet_df_2) != 0){
    mass_match_5 = ggplot(prophet_df_1, aes(x = n_position, y = relative_abundance, fill = sum_intensity)) +
      geom_bar(stat = "identity") +
      labs(title = paste0("Relative Intensity for Data Points Corresponding to Each Position on ", ladder_name_1, " Found During Mass Matching"),
           x = "Position",
           y = "Relative Intensity") +
      theme_bw() + 
      theme(plot.title = element_text(
        size = 8,        
        hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
      )
      ) +
      scale_fill_gradient(low = "darkblue", high = "red") +
      scale_x_continuous(breaks = seq(0, max(prophet_df_1$n_position, na.rm = TRUE), by = 2))
    mass_match_3 = ggplot(prophet_df_2, aes(x = n_position, y = relative_abundance, fill = sum_intensity)) +
      geom_bar(stat = "identity") +
      labs(title = paste0("Relative Intensity for Data Points Corresponding to Each Position on ", ladder_name_2, " Found During Mass Matching"),
           x = "Position",
           y = "Relative Intensity") +
      theme_bw() + 
      theme(plot.title = element_text(
        size = 8,        
        hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
      )
      ) +
      scale_fill_gradient(low = "darkblue", high = "red") +
      scale_x_continuous(breaks = seq(0, max(prophet_df_2$n_position, na.rm = TRUE), by = 2))
    mass_match_5/mass_match_3
  }
}

plot_branching = function(df_name){
  ggplot(df_name, aes(x = monoisotopic_mass, y = apex_rt, color = Log_Intensity, shape = ladder_type)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_gradientn(colors = turbo(256)) +
    labs(x = "Monoisotopic Mass(Da)", y = "Retention Time(min)") +
    geom_text_repel(aes(label = base_name), color = "black", vjust = -1, size = 5) +
    theme(
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 20),   
      plot.title = element_text(size = 16)   
    ) +
    theme_bw() +
    ggtitle("")
}

#This function is used to plot the native 3' and 5' ladders, the x axis is monoisotopic mass, the y-axis is retention time(apex_rt)
#the color represents the ladder type, which should only be native_5' and native_3'. The breaks for x and y can be adjusted as
#you wish.
#first variable: the df you want to draw based on, this should be in the syntax of output by function theo_creater.
#second variable: the breaks for x-axis(syntax: c(1000, 2000, 3000))
#third variable: the breaks for y-axis(syntax: c(2, 4, 6))
plot_natives = function(data_frame, x_breaks, y_breaks){
  temp = data_frame %>% 
    filter(ladder_type %in% c("native 5'", "native 3'"))
  ggplot(temp, aes(x = monoisotopic_mass, y = apex_rt, color = log10(sum_intensity), shape = ladder_type)) +
    geom_line(aes(group = ladder_type),color = "gray", size = 1) +
    geom_point(
      alpha=0.9,
      size=3) +
    theme_bw() +
    scale_color_gradientn(colors = turbo(256)) +
    labs(
      x = "Monoisotopic Mass (Da)",
      y = "Rentention Time (min)"
    ) +
    theme(plot.title = element_text(
      size = 14,        
      face = "bold",    
      hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
    ),
    panel.grid = element_blank()
    ) +
    geom_text_repel(
      aes(label = base_name),  
      size = 3,               
      color = "black",        
      max.overlaps = 70       
    ) +
    scale_x_continuous(
      breaks = x_breaks
    ) +
    scale_y_continuous(
      breaks = y_breaks
    )
}


#This function plots the depth of the sequencing result, the x-axis will be the position of nucleotides, the y-axis will be
#the number of depth sequenced at that particular position.
#first variable: the df you want to draw based on, this should be in the syntax of output by function theo_creater.
#second variable: the positions that you don't want to be included in this depth plot.(Ex. the intact position)
plot_depth = function(data_frame, positions_to_eliminate){
  temp = data_frame %>%
    add_count(position_5, name = "count_appearances") %>% 
    distinct(position_5, .keep_all = TRUE) %>% 
    select(position_5, count_appearances) %>% 
    arrange(position_5) %>% 
    filter(!position_5 %in% positions_to_eliminate)
  
  temp = temp %>%
    rename(depth = count_appearances)
  
  # Draw the bar plot
  ggplot(temp, aes(x = position_5, y = depth)) +
    geom_bar(stat = "identity", fill = "steelblue") + 
    labs(
      x = "Position",
      y = "Depth",
      title = "Sequencing Depth at Each Position"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      text = element_text(family = "Helvetica", face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )+
    scale_y_continuous(
      breaks = c(2, 4, 6, 8, 10, 12, 14)
    )+
    scale_x_continuous(
      breaks = c(20, 40, 60, 80)
    )
}


#This function is used to draw the chosen layers, user may want to do this if they want to see branching. The x axis is 
#monoisotopic mass, the y-axis is retention time(apex_rt), the color represents the layer type. The breaks and limits 
#for x and y can be adjusted as user wish.
#first variable: the df you want to draw based on, this should be in the syntax of output by function theo_creater.
#second variable: the layers user choose(syntax: c("layer1", "layer2"))
#third variable: the breaks for x-axis(syntax: c(1000, 2000, 3000))
#fourth variable: the breaks for y-axis(syntax: c(2, 4, 6))\
#fifth variable: the limits for x-axis(syntax: c(1000, 3000))
#sixth variable: the limits for y-axis(syntax: c(2, 6))
plot_chosen_layers = function(data_frame, chosen_layers, x_breaks, y_breaks, x_limits, y_limits){
  temp = data_frame %>% 
    select(monoisotopic_mass, apex_rt, ladder_type) %>% 
    drop_na() %>% 
    filter(ladder_type %in% chosen_layers) 
  
  ggplot(temp, aes(x = monoisotopic_mass, y = apex_rt, color = ladder_type)) +
    geom_line(aes(group = ladder_type),color = "gray", size = 1) +
    geom_point(
      alpha=0.9,
      size=3) +
    theme_bw() +
    labs(
      x = "Monoisotopic Mass (Da)",
      y = "Rentention Time (min)"
    ) +
    theme(plot.title = element_text(
      size = 14,        
      face = "bold",    
      hjust = 0.5       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
    )
    ) +
    theme_classic(base_size = 12) +
    scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = y_limits
    )
}

#This function is for plotting the homology search plot, user can later use lines to connect the dots in the same group to
#show their relationship.
#First variable: the data frame that this function is going to draw upon, this data frame needs to include columns named 
#monoisotopic_mass, sum_intensity, group(which mass group, Ex.M1), and label(Ex.M1IF1).
#Second variable: the colors that you want the groups to be.(Ex.c("M1" = "red", "M2" = "blue", "M3" = "green"))
plot_homology = function(data_frame, group_color){
  ggplot(data_frame, aes(x = as.numeric(monoisotopic_mass), y = as.numeric(sum_intensity), color = group)) +
    geom_point(size = 6) +  # Plot points with size 4
    geom_text_repel(aes(label = label), vjust = -1, size = 4, color = "black") +  # Add labels above points
    labs(
      x = "Monoisotopic mass (Da)",
      y = "Sum intensity",
      color = "Group"
    ) +
    scale_color_manual(values = group_color) +  # Define colors for groups
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    )
}