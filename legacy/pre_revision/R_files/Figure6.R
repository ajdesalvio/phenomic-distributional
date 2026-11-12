library(tidyverse)
library(ggplot2)
library(data.table)
library(cowplot)
library(RColorBrewer)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# RSquared and MAE statistics from within-environment elastic net yield prediction (the script that generated varimp data)
metrics <- read.csv(paste0(data_path, 'Elastic_Net_Yield_Pred_Metrics.csv'))
metrics$Scenario <- gsub('drop_', '', metrics$Scenario)

# Create an "environment" for all data so we can facet_wrap with one data frame
metrics_allenv <- metrics %>% group_by(Scenario) %>%
  summarize(RSquared = mean(RSquared),
            RMSE = mean(RMSE),
            MAE = mean(MAE),
            Bias_Mean = mean(Bias_Mean)) %>%
  mutate(Env = 'All Environments')

# rbind the data frames
metrics_combined <- bind_rows(metrics, metrics_allenv)

# Environments used in the distributional data paper plus All Environments
envs <- c('All Environments', 'IAH4.2021', 'MIH1.2020', 'MNH1.2020', 'MNH1.2021', 'TXH1.2020', 'TXH2.2021', 'WIH1.2021',
          'WIH2.2020', 'WIH2.2021', 'WIH3.2020', 'WIH3.2021')

# Filter for environments in the paper
metrics_combined <- metrics_combined %>% filter(Env %in% envs)

# Factor levels for filtering parameters and median-only
metrics_combined$Scenario <- factor(metrics_combined$Scenario, levels = c('Median',
                                                        'all', '01_1', '01_02_99_1',
                                                        '01_03_98_1', '01_05_96_1',
                                                        '01_10_91_1', '01_15_86_1',
                                                        '01_20_81_1', '01_25_76_1'))
metrics_combined$Env <- factor(metrics_combined$Env, levels = envs)

p1 <- metrics_combined %>%
  filter(Env %in% envs) %>%
  ggplot(aes(x = Scenario, y = RSquared, fill = Scenario)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Paired') +
  facet_wrap(~ Env) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label=round(RSquared, 3)), hjust=0.0, color="black",
            position = position_dodge(0.9), angle = 90,size=3.5) +
  ylim(c(0, 0.80))
p1

jpeg('Elastic.Net.Yield.Pred.V2.jpg', width = 8.5, height = 6, units = 'in', res = 900)
p1
dev.off()

metrics_combined %>% group_by(Scenario) %>% summarize(RSquared = mean(RSquared)) %>% arrange(RSquared)
