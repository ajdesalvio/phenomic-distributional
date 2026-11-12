library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(ggpubr)
library(ggrepel)
library(hms)
library(ggpmisc)
library(scales)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Flight dates and times
fdt <- read.csv(paste0(data_path, 'Flight_Date_DAT_Conversion_Time_of_Day.csv'))

# Remove low-quality DATs
DATs_exclude <- c(46, 71:72, 99:100)
fdt <- fdt %>% filter(!DAT %in% DATs_exclude)

# Variance components
varcomp <- fread(paste0(data_path, 'VarComp_Within_DAT_Standard_VIs.csv')) %>% as.data.frame()

# Left join the flight dates and times
varcomp <- varcomp %>% left_join(fdt, by = 'DAT')

# Save only rows associated with Pedigree
varcomp.ped <- varcomp %>% filter(grp == 'Pedigree')

# Choose a specific vegetation index to plot
varcomp.ped.plot <- varcomp.ped %>% filter(Vegetation.Index == 'PSRI_median')
varcomp.ped.plot <- varcomp.ped %>% filter(Vegetation.Index %in% c('PSRI_median', 'CVI_median', 'ExG2_median', 'ExG_TP_median'))

# Create a time-of-day object (class hms)
varcomp.ped.plot <- varcomp.ped.plot %>% mutate(FlightTime = as_hms(paste0(Flight_Time, ':00')))

# Remove flight date with heritability of 0
varcomp.ped.plot <- varcomp.ped.plot %>% filter(!DAT %in% c(35, 42))


p2 <- ggplot(varcomp.ped.plot, aes(x = FlightTime, y = Heritability)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "bottom",
           color = "red") +
  scale_x_time(
    breaks = scales::breaks_width("2 hours"),
    labels = scales::label_time("%H:%M")
  ) +
  facet_wrap(~Vegetation.Index) +
  theme_minimal(base_size = 12) +
  theme(
    # make the actual panel background white
    panel.background = element_rect(fill = "white", color = NA),
    # draw a thin black border around each panel
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 0.5),
    # give the facet label “tab” a light gray fill with a black outline
    strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 0.5),
    # make the strip text stand out
    strip.text       = element_text(color = "black"),
    # rotate the x labels
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
    # remove any panel grid lines if you like
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Flight Time (HH:MM)",
    y = "Heritability"
  )

p2

jpeg('FigureS2.jpg', height = 5, width = 8, units = 'in', res = 600)
p2
dev.off()
