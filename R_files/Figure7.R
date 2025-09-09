library(tidyverse)
library(data.table)
library(tidyverse)
library(ggrepel)
library(cowplot)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Import variable importance scores from multiple filtering scenarios
varimp <- fread(paste0(data_path, 'Variable_Importance_Filtering_Scenarios_All_G2F_ExGR_V2.csv')) %>% as.data.frame()

# Create factors
varimp$Filtering_Params <- factor(varimp$Filtering_Params, levels = c('all', 'drop_01_1', 'drop_01_02_99_1',
                                                                      'drop_01_03_98_1', 'drop_01_05_96_1',
                                                                      'drop_01_10_91_1', 'drop_01_15_86_1',
                                                                      'drop_01_20_81_1', 'drop_01_25_76_1'))

# Visualize a new environment (something other than TXH1.2020)
p1 <- ggplot(varimp[varimp$Env=='IAH4.2021',], aes(x = Probability, y = Overall)) +
  geom_point(aes(color = Overall), size=2) +
  facet_wrap(~ DAP) +
  scale_color_viridis_c(option = "H", # Continuous viridis
                        name   = "Importance\nscore",
                        guide  = guide_colorbar(
                          barheight = unit(4, "cm"),
                          barwidth  = unit(0.4, "cm")
                        )
  ) +
  scale_y_continuous("Variable importance scores") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text       = element_text(color = "black", face = "plain", size = 10),
    axis.text.x      = element_text(angle = 90, vjust = 0.5)
  ) +
  labs(x = "Quantile Level") +
  facet_grid(Filtering_Params ~ DAP)

p1

jpeg('Elastic_Net_Sensitivity_Analysis_V1.jpg', width = 11, height = 10, units = 'in', res = 600)
p1
dev.off()

# Remove "drop" from the Filtering_Params column
varimp$Filtering_Params_V2 <- gsub('drop_', '', varimp$Filtering_Params)

# Create factors
varimp$Filtering_Params_V2 <- factor(varimp$Filtering_Params_V2, levels = c('all', '01_1', '01_02_99_1',
                                                                            '01_03_98_1', '01_05_96_1',
                                                                            '01_10_91_1', '01_15_86_1',
                                                                            '01_20_81_1', '01_25_76_1'))

p2 <- ggplot(varimp[varimp$Env=='IAH4.2021',], aes(x = Probability, y = Overall)) +
  geom_point(aes(color = Overall), size=2) +
  facet_wrap(~ DAP) +
  scale_color_viridis_c(option = "H",
                        name   = "Importance\nscore",
                        guide  = guide_colorbar(
                          barheight = unit(4, "cm"),
                          barwidth  = unit(0.4, "cm")
                        )
  ) +
  scale_y_continuous("Variable importance scores") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text       = element_text(color = "black", face = "plain", size = 10),
    axis.text.x      = element_text(angle = 90, vjust = 0.5)
  ) +
  labs(x = "Quantile Level") +
  facet_grid(Filtering_Params_V2 ~ DAP)

p2

jpeg('Elastic_Net_Sensitivity_Analysis_V2.jpg', width = 11, height = 9, units = 'in', res = 600)
p2
dev.off()

# Create a loop to save results from each environment
envs <- sort(unique(varimp$Env))
for (env in envs) {
  plot_loop <- 
    varimp %>% filter(!Filtering_Params == 'Median', Env == env) %>%
    ggplot(aes(x = Probability, y = Overall)) +
    geom_point(aes(color = Overall), size=2) +
    facet_wrap(~ DAP) +
    scale_color_viridis_c(option = "H",
                          name   = "Importance\nscore",
                          guide  = guide_colorbar(
                            barheight = unit(4, "cm"),
                            barwidth  = unit(0.4, "cm")
                          )
    ) +
    scale_y_continuous("Variable importance scores") +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      strip.text       = element_text(color = "black", face = "plain", size = 10),
      axis.text.x      = element_text(angle = 90, vjust = 0.5)
    ) +
    labs(x = "Quantile Level") +
    facet_grid(Filtering_Params_V2 ~ DAP) +
    ggtitle(env)
  
  jpeg(paste0('Elastic_Net_Sensitivity_Analysis', env, '.jpg'), width = 11, height = 9, units = 'in', res = 600)
  print(plot_loop)
  dev.off()
}


#### For the paper - filtering the varimp data ####
# Plot variable importance results for DAP 68
plot_dap <- varimp %>%
  filter(Env == 'TXH1.2020', DAP == 68, !Filtering_Params == 'Median') %>%
  ggplot(aes(x = Probability, y = Overall)) +
  geom_point(aes(color = Overall), size=2) +
  facet_wrap(~ DAP) +
  scale_color_viridis_c(option = "H",
                        name   = "Importance\nscore",
                        guide  = guide_colorbar(
                          barheight = unit(4, "cm"),
                          barwidth  = unit(0.4, "cm")
                        )
  ) +
  scale_y_continuous("Variable importance scores") +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text       = element_text(color = "black", face = "plain", size = 10),
    axis.text.x      = element_text(angle = 90, vjust = 0.5)
  ) +
  labs(x = "Quantile Level") +
  facet_grid(Filtering_Params_V2 ~ DAP) + 
  geom_vline(xintercept = 0.39, linetype = 'solid', linewidth = 1, color = '#dc267f') +
  ggtitle('TXH1.2020 DAP 68')
plot_dap

jpeg('Elastic_Net_Sensitivity_Analysis_TXH1.2020.DAP68.jpg', width = 11, height = 9, units = 'in', res = 600)
plot_dap
dev.off()

varimp.inspect3 <- varimp %>% filter(Env == 'IAH4.2021', DAP == 75) %>% arrange(Filtering_Params, -Overall)
