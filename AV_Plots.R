col <- c(
  # CP - discrete
  "CP1"  = "#E41A1C",  # Red
  "CP2"  = "#377EB8",  # Blue
  "CP3"  = "#4DAF4A",  # Green
  "CP4"  = "#984EA3",  # Purple
  "CP5"  = "#FF7F00",  # Orange
  "CP6"  = "#A65628"  # Brown
)

heatmap_theme<-theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.text.x = element_text(angle = 0, size = 8),
    strip.background = element_blank(),
    axis.ticks.x = element_line(),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_text(size=6)
  ) 

General_Theme<-theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size=8),
        axis.title = element_text(size = 10,face = "bold"),
        legend.text = element_text(size=8),
        plot.margin = margin(20,0,0,20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"))

score_theme<-General_Theme

CV_Theme<-General_Theme

# Define the title theme
title_theme <- theme(
  plot.title = element_text(face = "bold",size=12),
  legend.position = "none",
  plot.margin = margin(0, 0, 0, 0, "cm"),
  panel.spacing = unit(0, "cm")
)

include<-c("ENSG00000010292",
           "ENSG00000099985",
           "ENSG00000101350",
           "ENSG00000102572",
           "ENSG00000102804",
           "ENSG00000111728",
           "ENSG00000116133",
           "ENSG00000139970",
           "ENSG00000140694",
           "ENSG00000142089",
           "ENSG00000151208",
           "ENSG00000164111",
           "ENSG00000180881",
           "ENSG00000184557",
           "ENSG00000197451",
           "ENSG00000213585",
           "ENSG00000231389")


CVPlot<-function(Filt){
  ggp<-ggplot(Filt, 
         aes(y = Gene, x = value, color = Control)) +
    geom_point() +
    scale_color_manual(values = col) + 
    xlab("CV (%)") + 
    xlim(0, 21) + 
    geom_vline(xintercept = 15, linetype = "dashed") +
    geom_text(data = Filt %>% filter(value > 15 & value < 20), 
              aes(label = round(value, 1)), 
              hjust = -0.2, 
              size = 3, 
              show.legend = FALSE)+
    geom_text(data = Filt %>% filter(value > 20), 
              aes(label = round(value, 1)), 
              hjust = 1.3, 
              size = 3, 
              show.legend = FALSE)+
    CV_Theme
  return(ggp)
}

ScorePlot <- function(Full, error = F, ideal_range = NULL) {
  # Initialize the ggplot object
  ggp <- ggplot()
  
  # Add geom_rect first if ideal_range is specified
  if (!is.null(ideal_range)) {
    factor_levels <- levels(Full$Parameter)
    xmin_pos <- which(factor_levels %in% ideal_range[1])
    
    # Error check to ensure levels exist
    if (length(xmin_pos) == 0) {
      stop("One or both specified ideal_range factor levels do not exist in the data.")
    }
    
    ggp <- ggp + 
      geom_rect(aes(xmin = xmin_pos - 0.5, 
                    xmax = xmin_pos + 0.5, 
                    ymin = -Inf, ymax = Inf, fill = "Protocol"), 
                alpha = 0.2, inherit.aes = FALSE)
  }
  
  # Add the rest of the plot based on error flag
  if (error) {
    ggp <- ggp +
      geom_boxplot(data = Full, 
                   aes(x = Control, y = Score, shape = Parameter, color = Control, 
                       group = interaction(Control, Parameter)), 
                   position = position_dodge(width = 0.75), 
                   outlier.shape = NA, size = 0.1, width = 0.5) +
      geom_jitter(data = Full, 
                  aes(x = Control, y = Score, shape = Parameter, color = Control), 
                  position = position_dodge(width = 0.75), size = 0.5)
  } else {
    ggp <- ggp +
      geom_point(data = Full, aes(y = Score, x = Parameter, color = Control, group = Control)) +
      geom_line(data = Full, aes(y = Score, x = Parameter, color = Control, group = Control))
  }
  
  ggp <- ggp +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red", size = 0.8) +
    scale_color_manual(values = col) + 
    scale_fill_manual(name = "Protocol", values = "lightgreen") +
    ylim(0, 100) +
    annotate("text", x = Inf, y = 100, label = "High Risk", hjust = 1.1, vjust = 1, 
             size = 3, fontface = "bold", color = "black") +
    annotate("text", x = Inf, y = 0, label = "Low Risk", hjust = 1.1, vjust = 0, 
             size = 3, fontface = "bold", color = "black") +
    ylab("Tutivia\u2122 Risk Score") +
    score_theme +
    guides(color = guide_legend(order = 1, nrow = 1),
           fill = guide_legend(order = 2, override.aes = list(alpha = 0.2)))
  
  return(ggp)
}




CVPlot_New<-function(Data){
  ggp <- ggplot(Data, aes(y = Gene, x = Control, fill = value)) +
    geom_tile() +
    facet_grid(. ~ name, scales = "free", space = "free")+
    heatmap_theme + 
    scale_fill_gradientn(
      colours = c("navy", "white", "firebrick"),
      values = c(0, 15 / 20, 15 / 20, 1),  # Dynamic gradient with fixed cutoff at 15
      limits = c(0, 20),  # Set upper limit to 20 to ensure red cutoff
      oob = scales::squish,  # Ensure values above 15 stay red
      name = "CV (%)"
    ) +
    theme(axis.title.y = element_blank())
  return(ggp)
}



recode_genes<-function(list){
  ensembl_to_gene <- c(
    "ENSG00000099985" = "OSM",
    "ENSG00000101350" = "KIF3B",
    "ENSG00000102572" = "STK24",
    "ENSG00000102804" = "TSC22D1",
    "ENSG00000116133" = "DHCR24",
    "ENSG00000139970" = "RTN1",
    "ENSG00000151208" = "DLG5",
    "ENSG00000164111" = "ANXA5",
    "ENSG00000180881" = "CAPS2",
    "ENSG00000184557" = "SOCS3",
    "ENSG00000197451" = "HNRNPAB",
    "ENSG00000213585" = "VDAC1",
    "ENSG00000231389" = "HLA-DPA1",
    "ENSG00000111728" = "ST8SIA1",
    "ENSG00000010292" = "NCAPD2",
    "ENSG00000140694" = "PARN",
    "ENSG00000142089" = "IFITM3"
  )
  
  new <- ensembl_to_gene[list]
  
  return(new)
}