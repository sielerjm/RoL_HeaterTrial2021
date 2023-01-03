# Plot Settings



# Color Palette -----------------------------------------------------------
# Source: RcolorBrewer: https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/


## Sequential (Ordinal)--------------------------------------------------------------
#   Sequential palettes (first list of colors), which are suited to ordered data 
#     that progress from low to high (gradient).

pal.BuPu <- brewer.pal(9, "BuPu") 
pal.Greys <- brewer.pal(9, "Greys")


## Qualitative (categorical,  non-ordinal) ---------------------------------
#   Qualitative palettes (second list of colors), which are best suited to represent nominal or categorical data. 
#     They do not imply magnitude differences between groups. 

pal.Set1 <- brewer.pal(9, "Set1")  # Not colorblind friendly
pal.Set2 <- brewer.pal(8, "Set2")
pal.Dark2 <- brewer.pal(8, "Dark2")
pal.Paired <- brewer.pal(12, "Paired")


## Diverging (Numeric, spectrum) ------------------------------------------------------
#   Diverging palettes (third list of colors), which put equal emphasis on mid-range 
#     critical values and extremes at both ends of the data range.

pal.RdYlGn <- brewer.pal(11, "RdYlGn")
pal.Spectral <- brewer.pal(9, "Spectral")
pal.BrBg <- brewer.pal(11, "BrBG")



## Plot Specific Palettes --------------------------------------------------

col.Temp <- pal.Spectral[c(9,8,1)]
col.DPE <- RColorBrewer::brewer.pal(9, "YlOrRd")
col.Treat <- c(pal.Dark2[c(3,5)])





pal.TempDPE <- c(RColorBrewer::brewer.pal(9, "Blues")[c(1,3,5,7,9)],
                 RColorBrewer::brewer.pal(9, "Greens")[c(1,3,5,7,9)],
                 RColorBrewer::brewer.pal(9, "Reds")[c(1,3,5,7,9)]
                 )

TempDPE.breaks <- c("28°C_0", "32°C_0", "35°C_0", 
                    "28°C_14", "32°C_14", "35°C_14", 
                    "28°C_21", "32°C_21", "35°C_21", 
                    "28°C_28", "32°C_28", "35°C_28", 
                    "28°C_42", "32°C_42", "35°C_42"
                    )

col.TempDPE <- c(pal.TempDPE[1], pal.TempDPE[6], pal.TempDPE[11],
                 pal.TempDPE[2], pal.TempDPE[7], pal.TempDPE[12],
                 pal.TempDPE[3], pal.TempDPE[8], pal.TempDPE[13], 
                 pal.TempDPE[4], pal.TempDPE[9], pal.TempDPE[14], 
                 pal.TempDPE[5], pal.TempDPE[10], pal.TempDPE[15]
                 )



# List of color palettes --------------------------------------------------
## Use this in a loop to access color palletes for your treatment variables

col.list <- list(Temperature = col.Temp,
                 DPE = col.DPE,
                 Treatment = col.Treat)


# GGplot2 Theme -----------------------------------------------------------
#     Standardizing ggplot settings
#     (adapted from Keaton Stagaman)

my_theme <- theme_update(
  legend.position = "top",
  legend.box = "vertical",
  legend.box.just = "center",
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 12),
  legend.key = element_rect(fill = "white"),
  legend.key.size = unit(1, "line"), # legend symbol size
  legend.spacing.y = unit(0, 'cm'),
  
  strip.text = element_text(size = 12),
  
  plot.caption = element_text(hjust = 0, size = 12),
  
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 16),
  
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  panel.background = element_rect(fill = pal.Greys[1]),
  panel.grid.major = element_line(color = pal.Greys[3]),
  panel.grid.minor = element_line(color = pal.Greys[3]),
) 


## Alpha Plot Settings -----------------------------------------------------

alpha_plot_settings <- function(
    tmp.sig.labs, 
    tmp.y.pos = c(1, 1.075, 1.15)
    ){
  
  list(
    ggpubr::stat_pvalue_manual(tmp.sig.labs, 
                               label = "p.adj.signif", 
                               size = 6,
                               bracket.size = 1,
                               tip.length = 0,
                               y.position = tmp.y.pos,
                               hide.ns = T
                               ),
    scale_y_continuous(breaks = c(0, .25, .5, .75, 1), 
                         limits = c(0, 1.16)) 
  )
}



################################################################################
################################################################################