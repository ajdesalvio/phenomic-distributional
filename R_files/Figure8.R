library(tidyverse)
library(ggplot2)
library(data.table)
library(purrr)
library(abind)
library(cowplot)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Load BLUEs of maize quantile data
mz <- fread(paste0(data_path, 'Quantile_BLUEs_G2F_Distributional_ExGR.csv')) %>% as.data.frame()

# Average matrices with Fisher z-transformation
# Helper to keep only the quantile columns
get_quantile_mat <- function(tbl){
  tbl %>% # Still has Pedigree, DAP, etc.
    select(Pedigree, Probability, Quantile.BLUE) %>% # Keep essentials
    pivot_wider(names_from  = Probability,
                values_from = Quantile.BLUE,
                names_prefix = "q_") %>%   # q_0.01, q_0.02, - colnames
    select(starts_with("q_")) %>%          # Drop Pedigree, DAP, Estimate
    as.matrix()                            # Hand a pure numeric matrix to cor()
}

# Unique environments
envs <- c('IAH4.2021', 'MIH1.2020', 'MNH1.2020', 'MNH1.2021', 'TXH1.2020', 'TXH2.2021', 'WIH1.2021', 'WIH2.2020',
          'WIH2.2021', 'WIH3.2020', 'WIH3.2021')
# Remove erroneous DAPs
mz <- mz %>% filter(!DAP %in% c(-33, -9, -3, 0, 1, 5))

#### Fisher-adjusted mean correlations for each environment separately ####
# List to compile figures
fig_list <- list()

for (env in envs) {

  # One correlation matrix for each DAP
  cors_by_day <- mz %>% 
    filter(Vegetation.Index == 'ExGR', Env == env) %>% 
    group_split(DAP) %>% # One tibble per date
    map(~ cor(get_quantile_mat(.x),
              use = "pairwise.complete.obs"))

  # Fisher z-transform, average, back-transform
  # atanh() is Fisher z, tanh() is inverse
  z_array <- abind::abind(map(cors_by_day, atanh), along = 3)
  mean_z  <- apply(z_array, 1:2, mean, na.rm = TRUE)
  mean_z[is.infinite(mean_z)] <- NA
  mean_r  <- tanh(mean_z)

  # Tidy & plot
  if (is.null(colnames(mean_r))) {
    colnames(mean_r) <- rownames(mean_r) <- paste0("q_", sprintf("%05.2f", seq(0.01, 1, 0.01)))
  }

  # Put 1s on the diagonal
  diag(mean_r) <- 1
  print(env); print(range(mean_r, na.rm = T)) # For the results section in the paper

  # Wide to long
  heat_df <- mean_r %>% 
    as_tibble(rownames = "q1") %>% 
    pivot_longer(-q1, names_to = "q2", values_to = "r") %>% 
    # drop the "q_" prefix and convert to numeric so ggplot treats them as numbers
    mutate(across(c(q1, q2), ~ as.numeric(sub("^q_", "", .))))

  # Plot
  p_fisher <- ggplot(heat_df, aes(q1, q2, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(name  = 'Fisher-adjusted\nmean correlation\nacross DAPs',
                        low   = "#2166ac",
                        mid   = "white",
                        high  = "#b2182b",
                        limits = c(0, 1),
                        midpoint = 0.5) +
    coord_equal() +
    labs(x = "Quantile level", y = "Quantile level") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggtitle(env)
  p_fisher
  
  fig_list[[env]] <- p_fisher
  
  # Uncomment this to save the figures individually
  #jpeg(paste0('Fisher_Adj_Mean_Corr_', env, '.jpg'), width = 8, height = 6.5, units = 'in', res = 600)
  #print(p_fisher)
  #dev.off()

}

#### Fisher-adjusted mean correlations for all environments together ####
mz.paper <- mz %>% filter(Env %in% envs)
# Remove DAPs that don't contain real data
mz.paper <- mz.paper %>% filter(!DAP %in% c(-33, -9, -3, 0, 1, 5))
# Create an Env.DAP column since we need to create one correlation matrix per DAP but need to prevent creating an
# accidental combined matrix in case two locations share a flight DAP
mz.paper$Env.DAP <- paste(mz.paper$Env, mz.paper$DAP, sep = '.')
length(unique(mz.paper$Env.DAP)) # Number of correlation matrices that should be calculated in total

# One correlation matrix for each DAP
cors_by_day <- mz.paper %>% 
  filter(Vegetation.Index == 'ExGR') %>% 
  group_split(Env.DAP) %>% # One tibble per Env.DAP combination
  map(~ cor(get_quantile_mat(.x),
            use = "pairwise.complete.obs"))

# Fisher z-transform, average, back-transform
# atanh() is Fisher z, tanh() is inverse
z_array <- abind::abind(map(cors_by_day, atanh), along = 3)
mean_z  <- apply(z_array, 1:2, mean, na.rm = TRUE)
mean_z[is.infinite(mean_z)] <- NA
mean_r  <- tanh(mean_z)
dim(z_array) # 3rd dimension should match length(unique(mz.paper$Env.DAP))

# Tidy & plot
if (is.null(colnames(mean_r))) {
  colnames(mean_r) <- rownames(mean_r) <- paste0("q_", sprintf("%05.2f", seq(0.01, 1, 0.01)))
}

# Put 1s on the diagonal
diag(mean_r) <- 1
print(range(mean_r, na.rm = T)) # Again for the results section in the paper

# Wide to long
heat_df <- mean_r %>% 
  as_tibble(rownames = "q1") %>% 
  pivot_longer(-q1, names_to = "q2", values_to = "r") %>% 
  # drop the "q_" prefix and convert to numeric so ggplot treats them as numbers
  mutate(across(c(q1, q2), ~ as.numeric(sub("^q_", "", .))))

# Plot
p_fisher <- ggplot(heat_df, aes(q1, q2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(name  = 'Fisher-adjusted\nmean correlation\nacross DAPs',
                       low   = "#2166ac",
                       mid   = "white",
                       high  = "#b2182b",
                       limits = c(0, 1),
                       midpoint = 0.5) +
  coord_equal() +
  labs(x = "Quantile level", y = "Quantile level") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle('Combined')
p_fisher


# Add the all DAP figure to fig_list
fig_list[['All_Envs']] <- p_fisher
# Save one of the legends
legend <- get_legend(fig_list[[1]])
# Strip the legends from the plots in the list
fig_list_noleg <- lapply(fig_list, function(p) p + theme(legend.position = "none"))
# Append the legend as the 12th plot
fig_list_plus_leg <- c(fig_list_noleg, list(legend))
# Assemble the plot as a 3 x 4 matrix
combined <- plot_grid(plotlist = fig_list_plus_leg,
                      nrow = 3)
jpeg('Fisher_Adj_Mean_Corr_ByEnv_AllEnv_V3.jpg', width = 8, height = 6.5, units = 'in', res = 600)
combined
dev.off()
