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
fdt <- read.csv(paste0(data_path, 'Flight_Date_DAP_Conversion_G2F_Time_of_Day.csv'))
fdt$Env.DAP <- paste(fdt$Env, fdt$DAP, sep = '.')
fdt <- fdt %>% select(Env.DAP, Flight_Time)

# Variance components
varcomp <- fread(paste0(data_path, 'VarComp_per_Quantile_Within_DAP_RGB_VIs.csv')) %>% as.data.frame()
varcomp$Env.DAP <- paste(varcomp$Env, varcomp$DAP, sep = '.')

# Keep only the median
varcomp <- varcomp %>% filter(Probability == 0.50)

# Left join the flight dates and times
varcomp <- varcomp %>% left_join(fdt, by = 'Env.DAP')

# Save only rows associated with Pedigree
varcomp.ped <- varcomp %>% filter(grp == 'Pedigree')

# Choose a specific vegetation index to plot
varcomp.ped.plot <- varcomp.ped %>% filter(Vegetation.Index %in% c('ExGR', 'COM1', 'GminusR', 'CIVE'))

# Choose the environment(s) and remove NA or '' times
varcomp.ped.plot <- varcomp.ped.plot %>% 
                        filter(Env %in% c('TXH1.2020')) %>%
                        #rename(Flight_Time = Flight_Time.x) %>%
                        filter(!is.na(Flight_Time), Flight_Time != '')

# Create a time-of-day object (class hms)
varcomp.ped.plot <- varcomp.ped.plot %>% mutate(FlightTime = as_hms(paste0(Flight_Time, ':00')))

# Remove flight date with heritability of 0
varcomp.ped.plot <- varcomp.ped.plot %>% filter(!DAP %in% c(35, 42))


p1 <- ggplot(varcomp.ped.plot, aes(x = FlightTime, y = Heritability)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "bottom",
           color = "red") +
  scale_x_time(
    breaks = scales::breaks_width("1 hours"),
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

p1

jpeg('FigureS3.jpg', height = 5, width = 8, units = 'in', res = 600)
p1
dev.off()
