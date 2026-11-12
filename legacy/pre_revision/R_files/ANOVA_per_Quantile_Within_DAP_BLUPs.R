library(lubridate)
library(fda.usc)
library(lme4)
library(lmerTest)
library(ggplot2)
library(partR2)
library(dplyr)
library(tidyr)
library(data.table)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Set this path to the inner folder "Individual_VIs_Distributional" of All_Distributional_Files
data_path_distrib <- paste0(data_path, 'Individual_VIs_Distributional/')
# Path where RData will be saved
RData_path <- paste0(data_path, 'RData/')
csvs <- list.files(path = data_path_distrib, pattern = '\\.csv$')
all_vis <- lapply(strsplit(as.character(csvs), '_'), '[', 11) %>% unlist()
all_vis <- gsub('.csv', '', all_vis)
CSV_VI <- data.frame(CSV = csvs, VI = all_vis)

# Make sure VIs and CSVs match
CSV_VI$Check <- lapply(strsplit(as.character(CSV_VI$CSV), '_'), '[', 11) %>% unlist()
CSV_VI$Check <- gsub('.csv', '', CSV_VI$Check)
identical(CSV_VI$Check, CSV_VI$VI)

# Initialize empty lists to store results
G_list <- list()
VC_list <- list()

# Remove MOH1.2020.13 due to it causing a failure with ANOVA

for (i in 1:length(csvs)) {
  
  current_vi <- all_vis[i]
  print(paste('Current VI:', current_vi))
  
  # Read data from CSV files
  print(paste('Loading VI data frame for', current_vi))
  df <- fread(paste0(data_path_distrib, '2020_2021_G2F_TIFs_Renamed_RGB_VIs_Distributions_Hue_FieldInfo_', current_vi, '.csv')) %>% as.data.frame() %>% arrange(DAP)
  print(paste0('VI data frame loaded for ', current_vi, '!'))
  
  # Env.DAP column
  df$Env.DAP <- paste(df$Env, df$DAP, sep = '.')
  
  # Basic plot for inspection of data
  sort(unique(df$DAP))
  df.plot <- df %>% filter(Env == 'WIH1.2020', DAP == 58, Pedigree %in% c('W10004_0741/PHP02', 'W10004_0691/PHP02'))
  plot(df.plot$probability, df.plot$quantile)
  
  # Create the Env.DAP.Prob column
  df$Env.DAP.Prob <- paste(df$Env.DAP, df$probability, sep = '.')
  # Vector containing Envs, DAPs, and quantiles to loop through
  env.dap.prob <- unique(df$Env.DAP.Prob)
  
  # Loop over each element of Env.DAP.Prob
  for (edp in env.dap.prob) {
    
    # Subset data for the current Env/DAP/Probability combination
    df.subset <- df[df$Env.DAP.Prob == edp, ]
    
    # Convert variables to appropriate data types
    df.subset$Range <- as.factor(df.subset$Range)
    df.subset$Pass <- as.factor(df.subset$Pass)
    df.subset$Replicate <- as.factor(df.subset$Replicate)
    df.subset$Pedigree <- as.factor(df.subset$Pedigree)

    # Progress report
    print(paste('ANOVA Model:', edp, '(Env.DAP.Probability)'))
        
    # Fit the linear mixed-effects model after setting "Inf" and NaN as NA
    response <- df.subset$quantile
    response[is.na(response) | response == 'Inf' | is.nan(response) | response == '-Inf'] <- NA
        
    # Attempt to fit the model; if it fails, catch the error and move to the next iteration
    anova <- tryCatch({
      lmer(response ~
             (1 | Pedigree) +
             (1 | Range) +
             (1 | Pass) +
             (1 | Replicate),
           data = df.subset)
    }, error = function(e) {
      message(paste("Error in lmer for", edp, ":", e$message))
      return(NULL)
    })
    
    # If an error occurred (anova is NULL), skip to the next iteration
    if (is.null(anova)) next
        
    # Calculate RMSE
    rmse_value <- sqrt(mean(residuals(anova)^2))
        
    # Calculate R-squared
    R <- MuMIn::r.squaredGLMM(anova)[, "R2c"]  # Conditional R-squared
        
    # Extract variance components
    VC <- as.data.frame(VarCorr(anova))
    VC$Percent <- round(VC$vcov / sum(VC$vcov) * 100, 2)
        
    # Calculate heritability
    nRep <- length(unique(df.subset$Replicate))
    Vg <- VC[VC$grp == "Pedigree", "vcov"]
    Ve <- VC[VC$grp == "Residual", "vcov"]
    heritability <- Vg / (Vg + (Ve / nRep))
    heritability <- round(heritability, 3)
        
    # Useful variables
    current_DAP <- df.subset$DAP[1]
    current_prob <- df.subset$probability[1]
    current_env <- df.subset$Env[1]
        
    # Add additional information to VC data frame
    VC$Rmse <- rmse_value
    VC$Heritability <- heritability
    VC$R_squared <- R
    VC$DAP <- current_DAP
    VC$Vegetation.Index <- current_vi
    VC$Probability <- current_prob
    VC$Env <- current_env
    print(VC)
        
    # Store VC in the list
    VC_list[[paste(current_env, current_DAP, current_prob, current_vi, sep = '_')]] <- VC
        
    # Extract random effects for 'Pedigree'
    G <- coef(anova)$Pedigree
    G$Pedigree <- rownames(G)
    G$DAP <- current_DAP
    G$Vegetation.Index <- current_vi
    G$Probability <- current_prob
    G$Env <- current_env
    colnames(G)[1] <- 'Quantile.BLUP'
        
    # Store G in the list
    G_list[[paste(current_env, current_DAP, current_prob, current_vi, sep = '_')]] <- G
  }
  # Save RData for current loop iteration
  save.image(paste0(RData_path, 'ANOVA_per_Quantile_Within_DAP_', current_vi, '.RData'))
}

# Combine the lists into data frames
G_df <- do.call(rbind, G_list)
VC_df <- do.call(rbind, VC_list)

# View the first few rows of the results
head(G_df)
head(VC_df)

# Data wrangling before exporting
VC_df_cols <- colnames(VC_df)
VC_df_cols <- VC_df_cols[!VC_df_cols %in% c('var1', 'var2')]
VC_df <- VC_df[,VC_df_cols]
VC_df$DAP_Prob_VI <- rownames(VC_df)
VC_df$Probability <- lapply(strsplit(as.character(VC_df$DAP_Prob_VI), '_'), '[', 2) %>% unlist() %>% as.numeric()

colnames(G_df)[1] <- 'Quantile.BLUP'
G_df$DAP_Prob_VI.Pedigree <- rownames(G_df)
G_df$VI.Pedigree <- lapply(strsplit(as.character(G_df$DAP_Prob_VI.Pedigree), '_'), '[', 3) %>% unlist()
G_df$Probability <- lapply(strsplit(as.character(G_df$DAP_Prob_VI.Pedigree), '_'), '[', 2) %>% unlist() %>% as.numeric()
G_df$Vegetation.Index <- lapply(strsplit(as.character(G_df$VI.Pedigree), '\\.'), '[', 1)

# Save the data frames
#fwrite(VC_df, paste0(data_path, 'VarComp_per_Quantile_Within_DAP_RGB_VIs.csv'))
#fwrite(G_df, paste0(data_path, 'Quantile_BLUPs_per_Quantile_Within_DAP_RGB_VIs.csv'))

# Save only rows associated with Pedigree
VC_df_ped <- VC_df %>% filter(grp == 'Pedigree')

# Useful summaries (grouping by VI, DAP, and Probability)
quantile_summary_vi <- VC_df %>% group_by(Vegetation.Index) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_dat <- VC_df %>% group_by(DAP) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)
quantile_summary_prob <- VC_df %>% group_by(Probability) %>% summarize(MedHerit = median(Heritability)) %>% arrange(-MedHerit)

# Pedigree variance
ped_var_vi <- VC_df_ped %>% group_by(Vegetation.Index) %>% summarize(MedianPed = median(Percent)) %>% arrange(-MedianPed)
