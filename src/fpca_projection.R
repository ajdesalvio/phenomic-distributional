truncate_fpca_inputs <- function(fpca_data, t_min = NULL, t_max = NULL, min_points = 2) {
    stopifnot(all(c("Ly", "Lt", "Lid") %in% names(fpca_data)))
    if (is.null(t_min) && is.null(t_max)) {
        stop("At least one of t_min or t_max must be provided.")
    }
    if (!is.null(t_min) && !is.null(t_max) && t_min > t_max) {
        stop("t_min must be less than or equal to t_max.")
    }
    Ly2 <- fpca_data$Ly
    Lt2 <- fpca_data$Lt
    Lid2 <- fpca_data$Lid
    keep_subject <- logical(length(Ly2))
    for (i in seq_along(Ly2)) {
        tt <- Lt2[[i]]
        yy <- Ly2[[i]]
        keep <- rep(TRUE, length(tt))
        if (!is.null(t_min)) {
            keep <- keep & tt >= t_min
        }
        if (!is.null(t_max)) {
            keep <- keep & tt <= t_max
        }
        Lt2[[i]] <- tt[keep]
        Ly2[[i]] <- yy[keep]
        keep_subject[i] <- sum(keep) >= min_points
    }
    fpca_data$Ly <- Ly2[keep_subject]
    fpca_data$Lt <- Lt2[keep_subject]
    fpca_data$Lid <- Lid2[keep_subject]
    fpca_data
}

project_fpca <- function(VI, envs, heldouts=envs, quantiles=sort(unique(VI$Probability))) {
VI <- VI %>% filter(Env %in% envs, !DAP %in% c(-33,-9,-3,0,1,5))
veg.indices <- sort(unique(VI$Vegetation.Index)); fpc_out <- list()
for (vi.name in veg.indices) {
    for (env in heldouts) {
        for (quant.level in quantiles) {
            wide <- VI %>% filter(Probability == quant.level) %>% arrange(Env, DAP, Pedigree) %>% filter(Vegetation.Index == 
                vi.name) %>% pivot_wider(id_cols = c(Pedigree, Env), names_from = DAP, values_from = "Quantile.BLUE", 
                names_glue = "Quantile.BLUE.{DAP}", values_fill = NA) %>% as.data.frame()
            wide.DAPs <- colnames(wide)[3:ncol(wide)]
            wide.DAPs <- gsub("Quantile.BLUE.", "", wide.DAPs) %>% as.numeric() %>% sort()
            df.dap <- paste0("Quantile.BLUE.", wide.DAPs)
            wide <- wide %>% select(Pedigree, Env, any_of(df.dap))
            rownames(wide) <- paste(wide$Pedigree, wide$Env, sep = ".")
            DAPs <- as.numeric(gsub("Quantile.BLUE.", "", colnames(wide)[3:ncol(wide)]))
            fpca_data <- MakeFPCAInputs(IDs = rep(rownames(wide), each = length(sort(DAPs))), tVec = rep(sort(DAPs), 
                length(rownames(wide))), t(wide[, -c(1:2)]))
            env.t <- setdiff(envs, env)
            env.h <- env
            train_ids <- rownames(wide)[wide$Env %in% env.t]
            test_ids <- rownames(wide)[wide$Env %in% env.h]
            fpca_train <- list(Ly = fpca_data$Ly[train_ids], Lt = fpca_data$Lt[train_ids], Lid = fpca_data$Lid[match(train_ids, 
                unlist(fpca_data$Lid))])
            fpca_obj <- tryCatch(FPCA(fpca_train$Ly, fpca_train$Lt, list(dataType = "Sparse", methodMuCovEst = "smooth", 
                methodBwCov = "GCV", methodBwMu = "GCV", plot = TRUE)), error = function(e) {
                message("FPCA failed for VI=", vi.name, " q=", quant.level, " heldout=", env, " : ", 
                  e$message)
                NULL
            })
            if (is.null(fpca_obj)) stop("FPCA failed")
            scores.t <- as.data.frame(fpca_obj$xiEst)
            nExtract <- ncol(scores.t)
            colnames(scores.t) <- paste0("FPC", seq_len(ncol(scores.t)))
            rownames(scores.t) <- names(fpca_train$Ly)
            scores.t$Pedigree.Env <- rownames(scores.t)
            scores.t <- scores.t %>% relocate(Pedigree.Env)
            cumFVE <- fpca_obj$cumFVE
            nFVE <- min(length(cumFVE), nExtract)
            fpca_variation <- numeric(nFVE)
            if (nFVE >= 1) 
                fpca_variation[1] <- cumFVE[1]
            if (nFVE > 1) {
                for (i in 2:nFVE) {
                  fpca_variation[i] <- cumFVE[i] - cumFVE[i - 1]
                }
            }
            for (i in 1:nFVE) {
                scores.t[[paste0("FPC", i, "_FVE")]] <- fpca_variation[i]
            }
            message(paste0("Current FPCA projection: ", vi.name, ", quantile level ", quant.level, ", held-out experiment ", 
                env))
            max.df <- data.frame(Lt = unlist(fpca_train$Lt), Ly = unlist(fpca_train$Ly))
            max.df <- max.df %>% drop_na(Ly)
            end.num <- max(max.df$Lt)
            start.num <- min(max.df$Lt)
            fpca_test <- list(Ly = fpca_data$Ly[test_ids], Lt = fpca_data$Lt[test_ids], Lid = fpca_data$Lid[match(test_ids, 
                unlist(fpca_data$Lid))])
            fpca_test_trunc <- truncate_fpca_inputs(fpca_test, t_max = end.num, t_min = start.num, min_points = 2)
            pred <- predict(fpca_obj, newLy = fpca_test_trunc$Ly, newLt = fpca_test_trunc$Lt)
            scores.h <- as.data.frame(pred$scores)
            colnames(scores.h) <- paste0("FPC", seq_len(ncol(scores.h)))
            heldout_ids_kept <- names(fpca_test_trunc$Ly)
            rownames(scores.h) <- heldout_ids_kept
            scores.h$Pedigree.Env <- rownames(scores.h)
            scores.h <- scores.h %>% relocate(Pedigree.Env)
            message("Cutoff ", end.num, ": n_train=", length(fpca_train$Ly), " n_heldout=", length(fpca_test_trunc$Ly), 
                " n_total=", nrow(wide))
            scores_all <- bind_rows(scores.t, scores.h)
            scores_all$Pedigree <- lapply(strsplit(as.character(scores_all$Pedigree.Env), "\\."), "[", 
                1) %>% unlist()
            scores_all$Env <- lapply(strsplit(as.character(scores_all$Pedigree.Env), "\\."), "[", 2) %>% 
                unlist()
            scores_all$Env.VI <- paste(scores_all$Env, vi.name, sep = ".")
            scores_all$Vegetation.Index <- vi.name
            scores_all$Quantile.Level <- quant.level
            scores_all$Heldout.Env <- env
            head(scores_all)
            scores_all <- scores_all %>% relocate(Pedigree.Env, Pedigree, Env, Env.VI, Vegetation.Index, 
                Quantile.Level, Heldout.Env) %>% arrange(Env, Pedigree)
            scores_all_filled <- scores_all %>% group_by(Vegetation.Index, Quantile.Level, Heldout.Env) %>% 
                fill(matches("FVE"), .direction = "downup") %>% ungroup()
            fpc_out[[paste0(vi.name, ".Quantile.", quant.level, ".Heldout.Env.", env)]] <- scores_all_filled
        }
    }
}
out <- rbindlist(fpc_out, fill=TRUE)
# Recover full environment IDs through the input key, not dotted-name splitting.
ids <- unique(VI[,c("Pedigree","Env")]); ids$Pedigree.Env <- paste(ids$Pedigree,ids$Env,sep=".")
stopifnot(!anyDuplicated(ids$Pedigree.Env))
ix <- match(out$Pedigree.Env,ids$Pedigree.Env); stopifnot(!anyNA(ix))
out$Pedigree <- ids$Pedigree[ix]; out$Env <- ids$Env[ix]
out$Env.VI <- paste(out$Env,out$Vegetation.Index,sep=".")
out
}
