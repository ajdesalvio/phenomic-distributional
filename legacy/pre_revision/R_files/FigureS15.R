library(ggplot2)
library(partR2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggh4x)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

varcomp <- fread(paste0(data_path, 'VarComp_Within_DAT_Across_Quantile_Revised_Herit.csv')) %>% as.data.frame()

# Rename "probability" to "Quantile Level" for figure
varcomp$grp <- gsub('probability', 'Quantile Level', varcomp$grp)
varcomp$grp <- gsub('Pedigree', 'Genotype', varcomp$grp)

varcomp$grp<- factor(varcomp$grp, levels = c('Genotype',
                                             'Genotype:Quantile Level',
                                             'Quantile Level',
                                             'Range:Quantile Level',
                                             'Row:Quantile Level',
                                             'Replicate:Quantile Level',
                                             'Residual'))

# Identify four flight dates from key points in the season
VIs_for_fig <- c('PSRI', 'CVI', 'ExG2', 'ExG_TP')
df.fig <- varcomp %>% filter(Vegetation.Index %in% VIs_for_fig)
df.fig$Vegetation.Index <- factor(df.fig$Vegetation.Index, levels = VIs_for_fig)

p <- df.fig %>% 
  ggplot(aes(x = DAT, y = Percent, fill = grp, text = grp)) +
  geom_area() +
  geom_line(aes(y = Heritability * 100, color = "Heritability"), linetype = 'solid', linewidth = 0.6) +
  geom_line(aes(y = R_squared * 100, color = "RSquared"), linetype = 'dashed', linewidth = 0.6) +
  scale_y_continuous("Variance explained by component (%)") +
  scale_x_continuous("DAT", guide = guide_axis(n.dodge = 1)) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      linewidth = 0.1,
      linetype = "solid"
    )
  ) +
  scale_color_manual(
    name   = NULL,
    values = c(Heritability = "black", RSquared = "black"),
    guide  = guide_legend(
      override.aes = list(
        linetype = c("solid", "dashed"),
        shape    = NA
      ),
      order = 1
    )
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c(Heritability = "solid", RSquared = "dashed"),
    guide  = "none"
  ) +
  guides(
    fill = guide_legend(title = "Variance \ncomponent", order = 2),
    color = guide_legend(order = 1)
  ) +
  facet_wrap(~ Vegetation.Index, nrow = 2)

p

jpeg("FigureS15.jpg", width = 7, height = 5, units = "in", res=600)
p
dev.off()
