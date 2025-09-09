library(caret)
library(glmnet)
library(tidyverse)
library(ggplot2)
library(data.table)
library(doParallel)

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

# Filter BLUE data frame for most-heritable VI (ExGR for G2F)
blue.filt <- fread(paste0(data_path, 'Quantile_BLUEs_G2F_Distributional_ExGR.csv')) %>% as.data.frame()

# Unique environments
envs <- sort(unique(blue.filt$Env))

# Example visualization of the BLUE quantile data
blue.filt.example <- blue.filt %>% filter(DAP == 50, Pedigree == 'W10004_1018/PHZ51')
plot(blue.filt.example$Probability, blue.filt.example$Quantile.BLUE)

# Import yield BLUEs from 2020/2021 G2F
yield <- read.csv(paste0(data_path, 'Yield_DTA_DTS_BLUEs_G2F_2020_2021.csv'))

# Empty list to store results
varimp_results <- list()
elastic_net_models <- list()
elastic_results_list <- list()

# Running with out DEH1.2020 since all filtering scenarios failed
envs <- c('IAH4.2021', 'MIH1.2020', 'MNH1.2020', 'MNH1.2021',
          'TXH1.2020', 'TXH2.2021', 'WIH1.2021', 'WIH2.2020',
          'WIH2.2021', 'WIH3.2020', 'WIH3.2021')

for (env in envs) {

  # Filter yield for TXH1.2020
  yield.filt <- yield %>% filter(Env == env)

  # Pedigree overlap between yield and quantile data
  overlap <- Reduce(intersect, list(unique(blue.filt$Pedigree), unique(yield.filt$Pedigree)))

  # Further filter the yield and quantile data
  blue.filt.loop <- blue.filt %>% filter(Pedigree %in% overlap, Env == env)
  yield.filt.loop <- yield.filt %>% filter(Pedigree %in% overlap) %>% arrange(Pedigree)

  # Filtering parameters
  filter.list <- list(
    Median = c(seq(0.01, 0.49, 0.01), seq(0.51, 1.00, 0.01)),
    all = numeric(0),
    drop_01_1 = c(0.01, 1.00),
    drop_01_02_99_1 = c(0.01, 0.02, 0.99, 1.00),
    drop_01_03_98_1 = c(0.01, 0.02, 0.03, 0.98, 0.99, 1.00),
    drop_01_05_96_1 = c(seq(0.01, 0.05, 0.01), seq(0.96, 1.00, 0.01)),
    drop_01_10_91_1 = c(seq(0.01, 0.1, 0.01), seq(0.91, 1.00, 0.01)),
    drop_01_15_86_1 = c(seq(0.01, 0.15, 0.01), seq(0.86, 1.00, 0.01)),
    drop_01_20_81_1 = c(seq(0.01, 0.20, 0.01), seq(0.81, 1.00, 0.01)),
    drop_01_25_76_1 = c(seq(0.01, 0.25, 0.01), seq(0.76, 1.00, 0.01))
  )

  # Loop through the data filtering parameters
  for (scenario in names(filter.list)) {
    message('Running filtering scenario: ', env, ' ', scenario)
  
    # Filtering settings
    drops <- filter.list[[scenario]]
  
    # Filter the BLUE data - convert to integers to deal with floating-point errors
    drops_int <- round(filter.list[[scenario]] * 100)
  
    blue.scenario <- blue.filt.loop %>%
      mutate(ProbPct = round(Probability * 100)) %>%
      filter(!ProbPct %in% drops_int) %>%
      select(-ProbPct)
  
    # Make the quantile data into wide format
    blue.scenario.wide <- blue.scenario %>%
      pivot_wider(
        id_cols      = Pedigree,
        names_from   = c(DAP, Probability, Vegetation.Index),
        names_glue   = "{DAP}_{Probability}_{Vegetation.Index}",
        values_from  = Quantile.BLUE,
        values_fill  = NA    # in case some combinations are missing
      ) %>%
      arrange(Pedigree)
  
    # Join the data frames
    df <- yield.filt.loop %>% inner_join(blue.scenario.wide, by = 'Pedigree')
    
    # Replace any Inf, -Inf, NaN values with NA
    X <- as.matrix(df[, 24:ncol(df)])
    X[!is.finite(X)] <- NA
    df[, 24:ncol(df)] <- X
  
    # Remove rows where yield has an NA
    df <- df[complete.cases(df$Yield.t.ha.BLUE),]
  
    # Prepare to run elastic net
    cl <- makePSOCKcluster(parallel::detectCores() - 4)
    registerDoParallel(cl)
  
    # Set up k-fold cross-validation
    ctrl <- trainControl(
      method   = "repeatedcv",   # k‑fold CV repeated B times
      number   = 5,              # k = 5
      repeats  = 25,             # B = 25   → 125 folds total
      verboseIter = TRUE,        # print progress
      allowParallel = TRUE,      # let caret use the foreach backend
      returnResamp = "final",    # keep all resampled metrics
      savePredictions = "final"  # keep fold‑wise predictions if you need them
    )
  
    # Define the tuning grid
    # α is low‑dimensional and often exhibits a smooth performance profile
    # λ needs finer resolution to land near the sweet spot of shrinkage
    # λ is searched across a log-spaced grid
    tune_grid <- expand.grid(
      alpha  = seq(0, 1, by = 0.1),                       # 11 values
      lambda = exp(seq(log(1e-4), log(10), length = 100)) # 100 λ per α
    )
  
    # Fit the model
    set.seed(123)  # Reproducible CV splits
    message('Performing elastic net scenario: ', env, ' ', scenario)
    elastic_fit <- tryCatch(
      {
        train(
          x          = df[, 24:ncol(df)],
          y          = df$Yield.t.ha.BLUE,
          method     = "glmnet",
          trControl  = ctrl,
          tuneGrid   = tune_grid,
          preProcess = c("center", "scale"),
          na.action  = na.omit        # drop rows with any NA/Inf
        )
      },
      error = function(e) {
        message("Elastic net failed for ", env, ' ', scenario,
                ".  Reason: ", conditionMessage(e))
        return(NULL)                  # <-- key: propagate NULL to caller
      }
    )
    if (is.null(elastic_fit)) next
    
    # Check for all-NA metrics
    if (all(is.na(elastic_fit$results$RMSE))) {
      message("All RMSE values are NA for ", env, ' ', scenario, ". Skipping.")
      next
    }
  
    # Inspect results
    print(elastic_fit)
    plot(elastic_fit)          # profile of RMSE across (α, λ)
    best <- elastic_fit$bestTune
    cat("Optimal α  =", best$alpha,  "\n")
    cat("Optimal λ  =", signif(best$lambda, 3), "\n")
  
    # Extract the final model for predictions on new data
    final_model <- elastic_fit$finalModel
  
    # Variable importance
    varimp <- varImp(elastic_fit)$importance %>% arrange(-Overall)
    
    # Metrics
    pred_df <- elastic_fit$pred %>%
      dplyr::filter(alpha == best$alpha, lambda == best$lambda)
    obs  <- pred_df$obs
    pred <- pred_df$pred
    
    rmse  <- sqrt(mean((obs - pred)^2))
    mae   <- mean(abs(obs - pred))
    
    elastic_results <- tibble(
      Env = env,
      Scenario = scenario,
      RSquared = max(elastic_fit$results$Rsquared, na.rm = TRUE),
      RMSE = rmse,
      MAE = mae,
      Bias_Mean = mean(pred - obs),
      best_alpha = elastic_fit$bestTune$alpha,
      best_lambda = elastic_fit$bestTune$lambda
    )
    print(elastic_results)
  
    # Shut down the cluster
    stopCluster(cl)
    registerDoSEQ()
  
    # Data wrangling to view variable importance
    varimp$DAP_Probability_Vegetation.Index <- rownames(varimp)
    varimp$DAP <- lapply(strsplit(as.character(varimp$DAP_Probability_Vegetation.Index), '_'), "[", 1) %>% as.numeric %>% unlist()
    varimp$Probability <- lapply(strsplit(as.character(varimp$DAP_Probability_Vegetation.Index), '_'), "[", 2) %>% as.numeric %>% unlist()
    varimp$Vegetation.Index <- lapply(strsplit(as.character(varimp$DAP_Probability_Vegetation.Index), '_'), "[", 3) %>% unlist()
    varimp$Filtering_Params <- scenario
    varimp$Env <- env
  
    # Add results to lists
    varimp_results[[paste(env, scenario, sep = '_')]] <- varimp
    elastic_net_models[[paste(env, scenario, sep = '_')]] <- final_model
    elastic_results_list[[paste(env, scenario, sep = '_')]] <- elastic_results
    
  }
}

# Save all results (temporary in case of crash)
save.image(file = 'C:/Users/ajdes/Downloads/Quantile_BLUEs_G2F_Distributional_ExGR/20250821_Elastic_Net_V9_Home.RData')

# combine results from all filtering scenarios
varimp.all <- rbindlist(varimp_results)

# Combine metrics
metrics <- rbindlist(elastic_results_list, fill = TRUE)

# Variable importance plot (remove anything with an importance score of 0)
varimp.fig <- varimp.all %>% filter(Overall > 0) %>% arrange(Env, Filtering_Params, DAP, Probability)

# Save variable importance data
#write.csv(varimp.fig, paste0(data_path, 'Variable_Importance_Filtering_Scenarios_All_G2F_ExGR_V2.csv'), row.names = F)

# Save metrics
#write.csv(metrics, paste0(data_path, 'Elastic_Net_Yield_Pred_Metrics.csv'), row.names = F)

