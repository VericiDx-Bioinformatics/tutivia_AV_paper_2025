col <- c(
  # Time - numeric
  "T0"   = "#440154",  # Dark purple
  "5D"   = "#3B528B",  # Deep blue
  "17D"  = "#21908C",  # Teal
  "1M"   = "#5DC863",  # Green
  "3M"   = "#AADC32",  # Lime
  "6M"   = "#FDE725",  # Yellow-green
  "12M"  = "#FEE08B",  # Soft yellow
  
  # W/WO - discrete
  "W"    = "#1F78B4",  # Blue
  "WO"   = "#33A02C",  # Green
  
  # LR - discrete
  "L1R1" = "#A6CEE3",  # Light blue
  "L1R2" = "#B2DF8A",  # Light green
  "L2R1" = "#FB9A99",  # Salmon pink
  "L2R2" = "#CAB2D6",  # Lavender
  
  # CP - discrete
  "CP1"  = "#E41A1C",  # Red
  "CP2"  = "#377EB8",  # Blue
  "CP3"  = "#4DAF4A",  # Green
  "CP4"  = "#984EA3",  # Purple
  "CP5"  = "#FF7F00",  # Orange
  "CP6"  = "#A65628",  # Brown
  
  # N - numeric
  "N25"  = "#1B9E77",  # Teal
  "N32"  = "#66A61E",  # Olive
  "N39"  = "#FFD92F",  # Yellow
  "N16"  = "#E6AB02",  # Amber
  "N64"  = "#D95F02",  # Burnt orange
  
  # RI - numeric (purple → green → yellow)
  "RI320" = "#482878",  # Deep purple
  "RI400" = "#414487",  # Indigo
  "RI500" = "#31688E",  # Blue-green
  "RI600" = "#35B779",  # Bright green
  "RI720" = "#FDE725",  # Yellow
  
  # LI - numeric (blue → light cyan → orange)
  "LI675" = "#2C7BB6",  # Blue
  "LI750" = "#99D8C9",  # Aqua
  "LI825" = "#FDAE61",  # Orange
  
  # FT - numeric
  "1FT" = "#9E3D9C",  # Aqua
  "2FT" = "#D73027",
  
  # Summary colors
  "RNA Stability - Storage Time" = "#5DC863",
  "Library Preparation - Interference" = "#1F78B4",
  "Reproducibility"   = "#FB9A99",
  "Repeatability" = "#CAB2D6",
  "Sequencing Characterization - Batch Size"    = "#E6AB02",
  "Library Preparation - RNA Input"   = "#35B779",
  "Sequencing Characterization - Library Input"   = "#8073AC",
  "RNA Stability - Freeze/Thaw"   = "#D73027"

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

General_Theme<-theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 10,face = "bold"),
        legend.text = element_text(size=8),
        plot.margin = margin(20,0,0,20))

Paired_Theme<-General_Theme+
  theme(panel.grid = element_blank())

CV_Theme<-General_Theme +
  theme(panel.grid.major.x = element_line(linewidth = 0.1,color="gray"),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1,color="gray"))


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

ScorePlot<-function(Full,error=F){
  if(error){
    ggp<-ggplot(Full, aes(x = Control, y = Score, shape = Parameter, color = Control)) +
      geom_boxplot(aes(group = interaction(Control, Parameter)), 
                   position = position_dodge(width = 0.75), 
                   outlier.shape = NA, size = 0.1,width = 0.5) +
      geom_jitter(position = position_dodge(width = 0.75), size = 0.5)
  }else{
    ggp<-ggplot(Full,
                aes(y = Score, x = Parameter, color = Control, group = Control))+
      geom_point() +
      geom_line()
  }
  ggp<-ggp +
    geom_hline(yintercept = 40, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = 60, linetype = "dashed", color = "red", size = 0.8) +
    scale_color_manual(values = col) + 
    ylim(0, 100) +
    annotate("text", x = Inf, y = 100, label = "High Risk", hjust = 1.1, vjust = 1, 
             size = 3, fontface = "bold", color = "black") +
    annotate("text", x = Inf, y = 0, label = "Low Risk", hjust = 1.1, vjust = 0, 
             size = 3, fontface = "bold", color = "black") +
    annotate("text", x = Inf, y = 50, label = "Intermediate Risk", hjust = 1.1, vjust = 0, 
               size = 3, fontface = "bold", color = "black")+ 
    ylab("Tutivia\u2122 Risk Score")+
    Paired_Theme+
    guides(color = guide_legend(nrow = 1))
  
  return(ggp)
  
}

CVPlot_New<-function(Data){
  
  ggp <- ggplot(Data, aes(y = Gene, x = Control, fill = value)) +
    geom_tile() +
    facet_grid(. ~ name, scales = "free", space = "free") +
    theme_minimal() +
    theme(
      legend.position = "right",
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      strip.text.x = element_text(angle = 0, size = 10),
      strip.background = element_blank(),
      axis.ticks.x = element_line(),
      panel.grid = element_blank()
    ) +
    scale_fill_gradientn(
      colours = c("navy", "white", "firebrick"),
      values = c(0, 15 / 20, 15 / 20, 1),  # Dynamic gradient with fixed cutoff at 15
      limits = c(0, 20),  # Set upper limit to 20 to ensure red cutoff
      oob = scales::squish,  # Ensure values above 15 stay red
      name = "CV (%)"
    ) +
    theme(axis.title = element_blank())
  return(ggp)
}