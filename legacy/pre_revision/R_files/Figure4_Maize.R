library(lubridate)
library(fda.usc)
library(lme4)
library(lmerTest)
library(ggplot2)
library(partR2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggh4x)

# Set working directory and source external functions
data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Import variance components
varcomp <- fread(paste0(data_path, 'VarComp_per_Quantile_Within_DAP_RGB_VIs.csv')) %>% as.data.frame()

# Save only rows associated with Pedigree
varcomp_ped <- varcomp %>% filter(grp == 'Pedigree')

# Useful summaries (grouping by VI, DAT, and Probability)
quantile_summary_vi <- varcomp %>% group_by(Vegetation.Index) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_dat <- varcomp %>% group_by(DAP) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_prob <- varcomp %>% group_by(Probability) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)

# Specific VI for paper conclusions
quantile_summary_dat_specific <- varcomp %>% filter(Vegetation.Index == 'ExGR') %>% group_by(DAP) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_prob_specific <- varcomp %>% filter(Vegetation.Index == 'ExGR') %>% group_by(Probability) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)

# Specific VI for paper conclusions, TXH1.2020 only
quantile_summary_dat_specific <- varcomp %>% filter(Vegetation.Index == 'ExGR', Env == 'TXH1.2020') %>% group_by(DAP) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_prob_specific <- varcomp %>% filter(Vegetation.Index == 'ExGR', Env == 'TXH1.2020') %>% group_by(Probability) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)

# Pedigree variance
ped_var_vi <- varcomp_ped %>% group_by(Vegetation.Index) %>% summarize(MedianPed = median(Percent)) %>% arrange(-MedianPed)

# Set factors for figures
varcomp$grp <- gsub('Pedigree', 'Genotype', varcomp$grp)
varcomp$grp<- factor(varcomp$grp, levels = c("Genotype",
                                             "Range",
                                             "Pass",
                                             "Replicate",
                                             "Residual"))

varcomp$Probability <- as.numeric(varcomp$Probability)

# Subset based on one Environment and VI
envs <- sort(unique(varcomp$Env))
env <- envs[2]
vi_for_fig <- 'ExGR'

for (env in envs) {
varcomp.env <- varcomp %>% filter(Vegetation.Index == vi_for_fig, Env == env)

p <- varcomp.env %>%
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
  
  # Your fill legend stays as before
  scale_fill_discrete(name = "Variance\ncomponent") +
  
  # Control legend order: lines first, then fills
  guides(
    color    = guide_legend(order = 1),
    linetype = guide_legend(order = 1),
    fill     = guide_legend(order = 2)
  ) +
  
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
                                     linetype  = "solid")
  ) +
  
  facet_wrap(~ DAP) +
  ggtitle(env)

p

jpeg(paste0('ANOVA_per_Quantile_Within_DAT_RGB_', vi_for_fig, '_', env, '.jpg'), width = 10, height = 6, units = "in", res=600)
print(p)
dev.off()
}

#### Subsetting for the paper ####
paper.envs <- c('IAH4.2021', 'MIH1.2020', 'MNH1.2020', 'MNH1.2021', 'TXH1.2020', 'TXH2.2021', 'WIH1.2021', 'WIH2.2020',
                'WIH2.2021', 'WIH3.2020', 'WIH3.2021')
varcomp.paper.envs <- varcomp %>% filter(Env %in% paper.envs)

quantile_summary_paper_vi <- varcomp.paper.envs %>% group_by(Vegetation.Index) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_paper_dat <- varcomp.paper.envs %>% group_by(DAP) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_paper_prob <- varcomp.paper.envs %>% group_by(Probability) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)

# Pedigree variance for locations in the paper
varcomp_ped_paper <- varcomp.paper.envs %>% filter(grp == 'Genotype')
ped_var_vi_paper <- varcomp_ped_paper %>% group_by(Vegetation.Index) %>% summarize(MedianPed = median(Percent)) %>% arrange(-MedianPed)

# Investigating TXH1.2020 results only
quantile_summary_txh.2020_vi <- varcomp %>% filter(Env == 'TXH1.2020') %>% group_by(Vegetation.Index) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_txh.2020_dat <- varcomp %>% filter(Env == 'TXH1.2020') %>% group_by(DAP) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_txh.2020_prob <- varcomp %>% filter(Env == 'TXH1.2020') %>% group_by(Probability) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)

# Pedigree variance for TXH1.2020
varcomp_ped_txh12020 <- varcomp %>% filter(Env == 'TXH1.2020') %>% filter(grp == 'Genotype')
ped_var_vi_txh12020 <- varcomp_ped_txh12020 %>% filter(Env == 'TXH1.2020') %>% group_by(Vegetation.Index) %>% summarize(MedianPed = median(Percent)) %>% arrange(-MedianPed)
