library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggh4x)

#### ANOVA per quantile within DAT ####
# Set working directory and source external functions
data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Import variance components
varcomp <- fread(paste0(data_path, 'VarComp_per_Quantile_Within_DAT.csv')) %>% as.data.frame()
# Save only rows associated with Pedigree
varcomp_ped <- varcomp %>% filter(grp == 'Pedigree')
# Set factors for figures
varcomp$grp <- gsub('Pedigree', 'Genotype', varcomp$grp)
varcomp$grp<- factor(varcomp$grp, levels = c("Genotype",
                                             "Range",
                                             "Row",
                                             "Replicate",
                                             "Residual"))
varcomp$Probability <- as.numeric(varcomp$Probability)
# Identify four flight dates from key points in the season:
DATs_for_fig <- c(57, 86, 96, 113, 130)
VIs_for_fig <- c('PSRI', 'CVI', 'ExG2', 'ExG_TP')
df.fig <- varcomp %>% filter(DAT %in% DATs_for_fig, Vegetation.Index %in% VIs_for_fig)
df.fig$Vegetation.Index <- factor(df.fig$Vegetation.Index, levels = VIs_for_fig)
df.fig$DAT <- as.character(df.fig$DAT)
df.fig$DAT <- factor(df.fig$DAT, levels = DATs_for_fig)
df.fig$VI.Type <- ifelse(df.fig$Vegetation.Index %in% c('PSRI', 'CVI'), 'Multispectral VIs', 'RGB VIs')
df.fig$VI.Type <- factor(df.fig$VI.Type, levels = c('Multispectral VIs', 'RGB VIs'))

p1 <- df.fig %>%
  ggplot(aes(x = Probability, y = Percent, fill = grp)) +
  geom_area() +
  geom_line(aes(y = Heritability * 100, color = "Heritability"), linetype = 'solid', linewidth = 0.6) +
  geom_line(aes(y = R_squared * 100, color = "RSquared"), linetype = 'dashed', linewidth = 0.6) +
  scale_y_continuous("Variance explained by component (%)") +
  scale_x_continuous("Quantile Level within DAT",
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     guide = guide_axis(n.dodge = 1)) +
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
    color = guide_legend(order = 1) # This puts the Heritability line on top
  ) +
  facet_grid(Vegetation.Index ~ DAT) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))

p1

#### ANOVA within quantile across DAT ####
varcomp2 <- fread(paste0(data_path, 'VarComp_Within_Quantile_Across_DAT_Revised_Herit.csv')) %>% as.data.frame()
varcomp2$grp <- gsub('Pedigree', 'Genotype', varcomp2$grp)
varcomp2$grp<- factor(varcomp2$grp, levels = c('Genotype',
                                             'Genotype:DAT',
                                             'DAT',
                                             'Range:DAT',
                                             'Row:DAT',
                                             'Replicate:DAT',
                                             'Residual'))

VIs_for_fig <- c('PSRI', 'CVI', 'ExG2', 'ExG_TP')
df.fig2 <- varcomp2 %>% filter(Vegetation.Index %in% VIs_for_fig)
df.fig2$Vegetation.Index <- factor(df.fig2$Vegetation.Index, levels = VIs_for_fig)

p2 <- df.fig2 %>% 
  ggplot(aes(x=Probability, y=Percent, fill=grp, text=grp)) +
  geom_area() +
  geom_line(aes(x=Probability, y=Heritability*100), linetype='solid') +
  geom_line(aes(x=Probability, y=R_squared*100), linetype='dashed') +
  scale_y_continuous("V.E.C. (%)") +
  scale_x_continuous('Quantile Level', guide = guide_axis(n.dodge=1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() +
  theme(legend.position='right',
        legend.background =  element_rect(fill = alpha("white", 0.7) , linewidth=0.1, linetype="solid"))+
  guides(fill=guide_legend(title="Variance \ncomponent", ncol = 2)) +
  facet_wrap(~Vegetation.Index, nrow=1) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))

p2

#### TXH1.2020 G2F results for top 4 most heritable VIs ####
# Import variance components
varcomp3 <- fread(paste0(data_path, 'VarComp_per_Quantile_Within_DAP_RGB_VIs.csv')) %>% as.data.frame()
# Set factors for figures
varcomp3$grp <- gsub('Pedigree', 'Genotype', varcomp3$grp)
varcomp3$grp<- factor(varcomp3$grp, levels = c("Genotype", "Range", "Pass", "Replicate", "Residual"))
varcomp3$Probability <- as.numeric(varcomp3$Probability)
# Subset based on one Environment and VI
vi_for_fig <- c('ExGR','GminusR', 'COM1', 'EXR')
DAPs_for_fig <- c(47, 68, 88, 96, 111) # Based on the heritability heatmap figure
varcomp3.env <- varcomp3 %>% filter(Vegetation.Index %in% vi_for_fig,
                                  Env == 'TXH1.2020',
                                  DAP %in% DAPs_for_fig)
varcomp3.env$Vegetation.Index <- factor(varcomp3.env$Vegetation.Index, levels = vi_for_fig)

p3 <- varcomp3.env %>%
  ggplot(aes(x = Probability, y = Percent, fill = grp)) +
  geom_area() +
  # Map both color & linetype to the same labels
  geom_line(aes(y = Heritability * 100,
                color    = "Heritability",
                linetype = "Heritability"),
            linewidth = 0.6) +
  geom_line(aes(y = R_squared * 100,
                color    = "RSquared",
                linetype = "RSquared"),
            linewidth = 0.6) +
  # ONE manual color scale for the two lines
  scale_color_manual(
    name   = NULL,
    values = c(
      Heritability = "black",
      RSquared     = "black"
    )
  ) +
  # ONE manual linetype scale for the two lines
  scale_linetype_manual(
    name   = NULL,
    values = c(
      Heritability = "solid",
      RSquared     = "dashed"
    )
  ) +
  scale_fill_discrete(name = "Variance\ncomponent") +
  # Control legend order: lines first, then fills
  guides(
    color    = guide_legend(order = 1),
    linetype = guide_legend(order = 1),
    fill     = guide_legend(order = 2)) +
  scale_y_continuous("Variance explained by component (%)") +
  scale_x_continuous("Quantile Level within DAP",
                     breaks = c(0,0.2,0.4,0.6,0.8,1.0),
                     guide  = guide_axis(n.dodge = 1)) +
  theme_bw() +
  theme(
    legend.position   = "right",
    legend.box        = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.7),
                                     linewidth = 0.1,
                                     linetype  = "solid"),
    axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  facet_grid(Vegetation.Index ~ DAP)

p3

jpeg('Figure4.jpg', width = 9, height = 8.5, units = 'in', res = 600)
cowplot::plot_grid(p1, p3, p2, labels = c('A', 'B', 'C'),
                   ncol = 1, rel_heights = c(0.42, 0.42, 0.16))
dev.off()
