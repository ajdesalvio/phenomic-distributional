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

data_path <- 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files/'

#### ANOVA across all DATs standard VIs MULTI ####
varcomp <- fread(paste0(data_path, 'VarComp_Across_All_DATs_Standard_VIs_Revised_Herit.csv')) %>% as.data.frame()
varcomp$grp <- gsub('Pedigree', 'Genotype', varcomp$grp)
varcomp$grp<- factor(varcomp$grp, levels = c('Genotype',
                                             'Genotype:DAT',
                                             'DAT',
                                             'Range:DAT',
                                             'Row:DAT',
                                             'Replicate:DAT',
                                             'Residual'))

# Identify four flight dates from key points in the season
VIs_for_fig <- c('PSRI_median', 'CVI_median', 'ExG2_median', 'ExG_TP_median')
df.fig <- varcomp %>% filter(Vegetation.Index %in% VIs_for_fig)
df.fig$Vegetation.Index <- factor(df.fig$Vegetation.Index, levels = VIs_for_fig)

p1 <- ggplot(df.fig, aes(x = Vegetation.Index, y = Percent)) +
  geom_col(aes(fill = grp), width = 0.7) +
  geom_point(aes(y = Heritability * 100, shape = "Heritability"),
             fill  = "white", color = "black", size = 4, alpha = 0.5) +
  geom_point(aes(y = R_squared * 100,   shape = "RSquared"),
             size  = 4, alpha = 1) +
  scale_shape_manual(
    name   = NULL,
    values = c(Heritability = 1, RSquared = 2),
    guide  = guide_legend(
      nrow  = 2,
      byrow = TRUE,
      order = 1
    )
  ) +
  guides(
    fill = guide_legend(title = "Variance \ncomponent", order = 2, nrow = 4)
  ) +
  theme_bw() +
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.background = element_rect(
      fill     = alpha("white", 0.5),
      linewidth= 0.1,
      linetype = "solid"
    ),
    axis.ticks.x = element_blank()
  ) +
  xlab("Vegetation Index") +
  ylab("Variance explained by component (%)")

p1

# Get colors from the first plot
get_fill_colors <- function(p) {
  b  <- ggplot_build(p)
  sc <- b$plot$scales$get_scales("fill")
  stopifnot(!is.null(sc))
  lim <- sc$get_limits()
  setNames(sc$map(lim), lim)
}

cols <- get_fill_colors(p1)
cols

cols_for_p2 <- cols[c(1,4,5,6,7)]

#### ANOVA per DAT standard VIs MULTI ####
varcomp2 <- fread(paste0(data_path, 'VarComp_Within_DAT_Standard_VIs.csv')) %>% as.data.frame()
varcomp2$grp <- gsub('Pedigree', 'Genotype', varcomp2$grp)
varcomp2$grp<- factor(varcomp2$grp, levels = c("Genotype",
                                             "Range",
                                             "Row",
                                             "Replicate",
                                             "Residual"))
# VIs with highest overall heritability
VI.herit <- varcomp2 %>% group_by(Vegetation.Index) %>% summarize(Median.Herit = median(Heritability)) %>% arrange(-Median.Herit)
# Mean vs. median
VI.herit$Type <- lapply(strsplit(as.character(VI.herit$Vegetation.Index), '_'), '[', 2) %>% unlist()
# VIs with highest pedigree variance
VI.ped <- varcomp2 %>% filter(grp == 'Pedigree') %>% group_by(Vegetation.Index) %>% summarize(Median.Pedig.Var = median(Percent)) %>% arrange(-Median.Pedig.Var)
# DATs with low heritability
DAT.herit <- varcomp2 %>% group_by(DAT) %>% summarize(Median.Herit = median(Heritability)) %>% arrange(-Median.Herit)

VIs_for_fig <- c('PSRI_median', 'CVI_median', 'ExG2_median', 'ExG_TP_median')
df.fig2 <- varcomp2 %>% filter(Vegetation.Index %in% VIs_for_fig, !DAT %in% c(42, 45, 48))
df.fig2$Vegetation.Index <- factor(df.fig2$Vegetation.Index, levels = VIs_for_fig)

p2 <- df.fig2 %>% 
  ggplot(aes(x = DAT, y = Percent, fill = grp)) +
  geom_area() +
  geom_line(aes(y = Heritability * 100,
                colour    = "Heritability",
                linetype  = "Heritability"),
            linewidth = 0.6) +
  geom_line(aes(y = R_squared * 100,
                colour    = "RSquared",
                linetype  = "RSquared"),
            linewidth = 0.6) +
  scale_color_manual(
    name   = NULL,
    values = c(Heritability = "black", RSquared = "black"),
    guide  = guide_legend(
      nrow         = 2,
      byrow        = TRUE,
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
    fill = guide_legend(
      title = "Variance \ncomponent",
      nrow  = 3,
      byrow = TRUE,
      order = 2
    )
  ) +
  scale_y_continuous("Variance explained by component (%)") +
  scale_x_continuous(
    "DAT",
    guide = guide_axis(n.dodge = 1)
  ) +
  
  theme_bw() +
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.background = element_rect(
      fill     = alpha("white", 0.7),
      linewidth= 0.1,
      linetype = "solid"
    )
  ) +
  facet_wrap(~ Vegetation.Index, nrow = 2) +
  scale_fill_manual(values = as.vector(cols_for_p2))

p2

jpeg('FigureS1.jpg', width = 9.6, height = 6.5, units = 'in', res = 600)
cowplot::plot_grid(p1, p2, labels = c('A', 'B'))
dev.off()
