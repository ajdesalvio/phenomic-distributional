plot_fui_custom <- function (fuiobj, num_row = NULL, xlab = "Functional Domain", 
          title_names = NULL, ylim = NULL, align_x = NULL, x_rescale = 1, 
          y_val_lim = 1.1, y_scal_orig = 0.05, beta_vec = NULL, return = FALSE) 
{
  num_var <- nrow(fuiobj$betaHat)
  plot_list <- res_list <- vector(length = num_var, "list")
  if (is.null(num_row)) 
    num_row <- ceiling(num_var/2)
  name = NULL
  align <- ifelse(is.null(align_x), 0, align_x * x_rescale)
  if (is.null(title_names)) 
    title_names <- rownames(fuiobj$betaHat)
  if (nrow(fuiobj$betaHat) != length(title_names)) 
    title_names <- rownames(fuiobj$betaHat)
  names(res_list) <- rownames(fuiobj$betaHat)
  for (r in 1:num_var) {
    beta_lab <- if (is.null(beta_vec)) r - 1 else beta_vec[r]
    if (is.null(fuiobj$betaHat.var)) {
      beta.hat.plt <- data.frame(s = fuiobj$argvals, beta = fuiobj$betaHat[r, 
      ])
      plot_list[[r]] <- ggplot() + theme_classic() + theme(plot.title = element_text(hjust = 0.5, 
                                                                                     face = "bold")) + geom_line(aes(x = s/x_rescale - 
                                                                                                                       align/x_rescale - 1/x_rescale, y = beta, color = "Estimate"), 
                                                                                                                 data = beta.hat.plt, alpha = 1, linewidth = 1) + 
        geom_hline(yintercept = 0, linetype = "dashed", 
                   color = "red") + scale_colour_manual(name = "", 
                                                        values = c(Estimate = "black")) + labs(x = xlab, 
                                                                                               y = bquote(beta[.(beta_lab)]*"(s)"), title = title_names[r]) + 
        theme(legend.position = "none")
    }
    else {
      beta.hat.plt <- data.frame(s = fuiobj$argvals, beta = fuiobj$betaHat[r, 
      ], lower = fuiobj$betaHat[r, ] - 2 * sqrt(diag(fuiobj$betaHat.var[, 
                                                                        , r])), upper = fuiobj$betaHat[r, ] + 2 * sqrt(diag(fuiobj$betaHat.var[, 
                                                                                                                                               , r])), lower.joint = fuiobj$betaHat[r, ] - fuiobj$qn[r] * 
        sqrt(diag(fuiobj$betaHat.var[, , r])), upper.joint = fuiobj$betaHat[r, 
        ] + fuiobj$qn[r] * sqrt(diag(fuiobj$betaHat.var[, 
                                                        , r])))
      plot_list[[r]] <- ggplot() + theme_classic() + theme(plot.title = element_text(hjust = 0.5, 
                                                                                     face = "bold")) + geom_ribbon(aes(x = s/x_rescale - 
                                                                                                                         align/x_rescale - 1/x_rescale, ymax = upper.joint, 
                                                                                                                       ymin = lower.joint), data = beta.hat.plt, fill = "gray20", 
                                                                                                                   alpha = 0.2) + geom_ribbon(aes(x = s/x_rescale - 
                                                                                                                                                    align/x_rescale - 1/x_rescale, ymax = upper, 
                                                                                                                                                  ymin = lower), data = beta.hat.plt, fill = "gray10", 
                                                                                                                                              alpha = 0.4) + geom_line(aes(x = s/x_rescale - 
                                                                                                                                                                             align/x_rescale - 1/x_rescale, y = beta, color = "Estimate"), 
                                                                                                                                                                       data = beta.hat.plt, alpha = 1, linewidth = 1) + 
        scale_colour_manual(name = "", values = c(Estimate = "black")) + 
        labs(x = xlab, y = bquote(beta[.(beta_lab)]*"(s)"), title = title_names[r]) + theme(legend.position = "none")
    }
    if (!is.null(ylim)) {
      plot_list[[r]] <- plot_list[[r]] + coord_cartesian(ylim = ylim)
      ylimit <- ylim
    }
    else {
      if (is.null(fuiobj$betaHat.var)) {
        ylimit <- c(min(beta.hat.plt$beta), max(beta.hat.plt$beta))
        y_adjust <- y_scal_orig * (max(beta.hat.plt$beta) - 
                                     min(beta.hat.plt$beta))
      }
      else {
        ylimit <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint))
        y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - 
                                     min(beta.hat.plt$lower.joint))
      }
      ylimit[1] <- ylimit[1] - y_adjust
    }
    xlim <- layer_scales(plot_list[[r]])$x$range$range
    x_range <- diff(xlim) * 0.1
    y_range <- diff(ylimit) * 0.1
    y_range_up <- diff(ylimit) * 0.02
    y_val_lim_vec <- c(1, y_val_lim)
    y_top <- (0.975) * diff(ylimit * y_val_lim_vec) + ylimit[1] * 
      y_val_lim_vec[1]
    plot_list[[r]] <- plot_list[[r]] + coord_cartesian(ylim = ylimit * 
                                                         y_val_lim_vec, xlim = xlim)
    if (!is.null(align_x)) {
      plot_list[[r]] <- plot_list[[r]] + geom_segment(aes_string(y = ylimit[1] - 
                                                                   y_range, yend = y_top, x = 0, xend = 0), inherit.aes = TRUE, 
                                                      color = "black", lwd = 0.5, alpha = 0.75, linetype = "dashed")
    }
    if (!is.null(fuiobj$betaHat.var)) {
      if (max(beta.hat.plt$upper.joint) > 0 & min(beta.hat.plt$lower.joint) < 
          0) {
        plot_list[[r]] <- plot_list[[r]] + geom_segment(aes_string(x = xlim[1] - 
                                                                     x_range, xend = xlim[2] + x_range, y = 0, yend = 0), 
                                                        inherit.aes = TRUE, color = "black", lwd = 0.5, 
                                                        alpha = 0.75, linetype = "dashed")
      }
      colnames(beta.hat.plt) <- c("s", "beta.hat", "CI.lower.pointwise", 
                                  "CI.upper.pointwise", "CI.lower.joint", "CI.upper.joint")
    }
    else {
      colnames(beta.hat.plt) <- c("s", "beta.hat")
    }
    res_list[[r]] <- beta.hat.plt
  }
  plot_return <- do.call(gridExtra::grid.arrange, c(plot_list, nrow = num_row))
  plot_return
  if (return == TRUE) {
    res_list$plot <- plot_return
    return(res_list)
  }
  else {
    return(plot_return)
  }
}