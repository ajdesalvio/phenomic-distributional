library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

# Argument parsing from param_list_ExGR.txt (each index gets its own .txt file)
args <- commandArgs(trailingOnly = TRUE)

# Expected usage:
#   Rscript  my_script.R  <csv_file>  [<edp_start> <edp_end>]
#
# * <csv_file>     – required, file name only (or path) of the CSV to process
# * <edp_start>    – optional, first Env.DAP.Prob index to process
# * <edp_end>      – optional, last  Env.DAP.Prob index to process

if (length(args) < 1 || length(args) == 2) {
  stop(
    "\nUsage: Rscript my_script.R <csv_file> [<edp_start> <edp_end>]\n",
    "Provide either 1 argument (the CSV) or 3 arguments (CSV + range)."
  )
}

input_csv  <- args[1]
start_edp  <- if (length(args) >= 3) as.integer(args[2]) else NA_integer_
end_edp    <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

# Directories and extracting VI name
data_path <- '/scratch/user/aarondesalvio/G2F_Distributional/0.0_Data/'
data_path_server <- '/scratch/user/aarondesalvio/G2F_Distributional/0.0_Data/Individual_VIs_Distributional/'
RData_path <- paste0(data_path, 'RData/')
vi_name <- sub("^.*_([^_]*)\\.csv$", "\\1", basename(input_csv))

if (!file.exists(paste0(data_path_server, input_csv)))
  stop("Cannot find CSV: ", input_csv)

if (!is.na(start_edp) && (is.na(end_edp) || end_edp < start_edp))
  stop("EDP range is malformed: start=", start_edp, " end=", end_edp)

# Initialize empty lists to store results
G_list <- list()
VC_list <- list()

message(paste('Current VI:', vi_name))
  
# Read data from CSV files
message(paste('Loading VI data frame for', vi_name))
df <- fread(paste0(data_path_server, '2020_2021_G2F_TIFs_Renamed_RGB_VIs_Distributions_Hue_FieldInfo_Trimmed_', vi_name, '.csv')) %>% as.data.frame() %>% arrange(DAP)
message(paste0('VI data frame loaded for ', vi_name, '!'))
  
# Env.DAP column
df$Env.DAP <- paste(df$Env, df$DAP, sep = '.')
  
# Create the Env.DAP.Prob column
df$Env.DAP.Prob <- paste(df$Env.DAP, df$probability, sep = '.')

# Remove any rows associated with .NA.0.01, .NA.0.02, etc.
# All of these rows are either Michigan plots or Plot_IDs 1015203, 1015204, or 36301
# There should be 30900 rows removed
before <- nrow(df)
df <- df[!grepl('NA', df$Env.DAP.Prob), ]
before - nrow(df)

# Vector containing Envs, DAPs, and quantiles to loop through
env.dap.prob <- unique(df$Env.DAP.Prob)

# If an EDP range was provided, subset to that range (1‑based indices)
if (!is.na(start_edp)) {
  if (end_edp > length(env.dap.prob))
    stop("EDP end (", end_edp, ") exceeds available combinations (", length(env.dap.prob), ")")
    
  env.dap.prob <- env.dap.prob[start_edp:end_edp]
  message("Processing EDP indices ", start_edp, ":", end_edp,
          " (", length(env.dap.prob), " combinations)")
}
  
# Loop over each element of Env.DAP.Prob
for (edp in env.dap.prob) {
    
  # Subset data for the current Env/DAP/Probability combination
  df.subset <- df[df$Env.DAP.Prob == edp, ]
  
  # Convert variables to the appropriate data types
  df.subset$Range <- as.factor(df.subset$Range)
  df.subset$Pass <- as.factor(df.subset$Pass)
  df.subset$Replicate <- as.factor(df.subset$Replicate)
  df.subset$Pedigree <- as.factor(df.subset$Pedigree)

  # Progress report
  message(paste('ANOVA Model:', edp, '(Env.DAP.Probability)'))
  message(paste('Completing ANOVA', match(edp, env.dap.prob), 'of', length(env.dap.prob)))
    
  # Fit the linear mixed-effects model after setting "Inf" and NaN as NA
  response <- df.subset$quantile
  response[is.na(response) | response == 'Inf' | is.nan(response) | response == '-Inf'] <- NA
    
  # Attempt to fit the model; if it fails, catch the error and move to the next iteration
  anova <- tryCatch({
    lmer(response ~
            Pedigree +
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
    
  # Extract variance components of random effects
  VC <- as.data.frame(VarCorr(anova))
  VC$Percent <- round(VC$vcov / sum(VC$vcov) * 100, 2)
    
  # Useful variables
  current_DAP <- df.subset$DAP[1]
  current_prob <- df.subset$probability[1]
  current_env <- df.subset$Env[1]
    
  # Add additional information to VC data frame
  VC$Rmse <- rmse_value
  VC$R_squared <- R
  VC$DAP <- current_DAP
  VC$Vegetation.Index <- vi_name
  VC$Probability <- current_prob
  VC$Env <- current_env
  print(VC)
    
  # Store VC in the list
  VC_list[[paste(current_env, current_DAP, current_prob, vi_name, sep = '_')]] <- VC
    
  # Extract fixed effects for 'Pedigree'
  G <- coef(summary(anova)) %>% as.data.frame()
    
  # The baseline (reference) level for a factor in R is the first element
  # returned by the levels() function for that factor. Use this to retrieve the
  # "missing" Pedigree - the one that was set as the intercept
  G$Pedigree <- gsub('Pedigree', '', rownames(G))
  # Retrieve the Intercept pedigree
  G$Pedigree[1] <- levels(df.subset$Pedigree)[1]
    
  # To retrieve BLUEs, we need to add the reference level (intercept) to all
  # other estimate values that were returned
  G$Quantile.BLUE <- NA
  G$Quantile.BLUE[1] <- G$Estimate[1]
  G$Quantile.BLUE[2:nrow(G)] <- G$Quantile.BLUE[1] + G$Estimate[2:nrow(G)]
    
  # Useful variables for data wrangling
  G$DAP <- current_DAP
  G$Vegetation.Index <- vi_name
  G$Probability <- current_prob
  G$Env <- current_env
    
  # Store G in the list
  G_list[[paste(current_env, current_DAP, current_prob, vi_name, sep = '_')]] <- G
}


# Build unique tags: VI + start‑end indices (used for saving RData)
tag <- sprintf("%s_%05d_%05d",
               vi_name,
               start_edp, end_edp)

# Combine the lists into data frames
G_df <- rbindlist(G_list, fill = TRUE) %>% as.data.frame()
VC_df <- rbindlist(VC_list, fill = TRUE) %>% as.data.frame()

# Data wrangling before exporting
VC_df_cols <- colnames(VC_df)
VC_df_cols <- VC_df_cols[!VC_df_cols %in% c('var1', 'var2')]
VC_df <- VC_df[,VC_df_cols]

# Save RData
saveRDS(G_df,
        file = file.path(RData_path, paste0("G_",  tag, ".rds")),
        compress = "xz")

saveRDS(VC_df,
        file = file.path(RData_path, paste0("VC_", tag, ".rds")),
        compress = "xz")
