library(lubridate)
library(fda.usc)
library(lme4)
library(lmerTest)
library(ggplot2)
library(partR2)
library(dplyr)
library(tidyr)

data_path <-"C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/" # Include the forward slash at the end

# Choose the Vegetation Index name
vi_load <- 'ExG2'

# Read data from CSV files
data <- fread(paste0(data_path, 'CS24-STEL-MULTI-RTK-Distributions-', vi_load, '-FieldCoords.csv')) %>% as.data.frame() %>% arrange(DAT)

df_wide <- data %>%
  pivot_wider(
    id_cols = c(
      Flight_Date, Test, Tier, Row, Range, Row.Range, Row.Range0,
      Pedigree, Pltg_ID, Pltg.ID, Transplant_Date, DAT, Replicate
    ),
    names_from = probability,
    values_from = quantile,
    names_prefix = paste0(vi_load, "_Quant_")
  )

# Remove low-quality DATs
DATs_exclude <- c(46, 71:72, 99:100)
df_wide <- df_wide %>% filter(!DAT %in% DATs_exclude)

# Remove plant that died
df_wide <- df_wide %>% filter(!Pltg_ID == '5011_018')

df_wide = as.data.frame(df_wide)
df_wide$Y = df_wide[,14:113]
df_wide$id = df_wide$Pltg_ID
data2 = df_wide

plot(fdata(data2$Y))

# Extract fda objects for specific IDs
fda1 = fdata(data2$Y[data2$id == '5013_028', ])[25:28]
fda2 = fdata(data2$Y[data2$id == '5010_001', ])[25:28]
fda3 = fdata(data2$Y[data2$id == '5013_025', ])[25:28]
fda4 = fdata(data2$Y[data2$id == '5008_022', ])[25:28]

# Set colors and titles
colors = c("darkred", "red", "plum", "steelblue1", "blue", "darkblue")

# Plot individual fda objects
par(mfrow = c(2, 2))
plot(fda1, col = colors,  ylab = "Q(p)", xlab = "Prob, p" )
plot(fda2, col = colors, ylab = "Q(p)", xlab = "Prob, p")
plot(fda3, col = colors,  ylab = "Q(p)", xlab = "Prob, p")
plot(fda4, col = colors, ylab = "Q(p)", xlab = "Prob, p")

dev.off()

#### Fit Multilevel Distributional functional model ####
#Pedigree: fixed effect (this is the genotype identifier)
#Range: random (Y coordinate in the field)
#Row: random (X coordinate in the field)
#Replicate: random (the rep # of the pedigree)

data2$Pedigree= as.factor(data2$Pedigree) 
data2$Replicate= as.factor(data2$Replicate)
data2$Pedigree= relevel(data2$Pedigree, ref='5016')

data2$Y= as.matrix(data2$Y)

library(refund)
library(fastFMM)
fit_dti <- fui(
  Y ~ Pedigree+(Range+Row|Replicate),
  data = data2
)

source(paste0(data_path, 'Plot_FUI_Custom_V2.R'))

beta_vec <- c(5016, 5001:5015)
jpeg('fastFMM_MULTI_ExG2_V2.jpg', height = 7, width = 10, units = 'in', res = 600)
plot_fui_custom(fit_dti, num_row = 4, title_names = c(paste(5016, '(TM-1)'), 5001:5015), xlab = 'Quantile Level',
                beta_vec = beta_vec)
dev.off()

