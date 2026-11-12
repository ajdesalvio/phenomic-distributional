# Rscript scripts/figures/render_correlations.R INPUT.csv NEW_OUTPUT
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cowplot)})
args<-commandArgs(TRUE);if(length(args)!=2L)stop('Usage: INPUT.csv NEW_OUTPUT')
out<-args[2];if(file.exists(out))stop('Use a new output directory')
data<-fread(args[1]);envs<-setdiff(unique(data$Panel),'Combined');panels<-c(envs,'Combined')
  stopifnot(nrow(data)==120000L,setequal(unique(data$Panel),panels),!anyDuplicated(data[,.(Panel,Q1,Q2)]),all(is.finite(data$r)),all(data$r>=-1 & data$r<=1))
  stopifnot(all(data[,.(N=.N,Rows=uniqueN(Q1),Columns=uniqueN(Q2)),by=Panel][,N==10000L & Rows==100L & Columns==100L]))
  stopifnot(setequal(data$Q1,1:100),setequal(data$Q2,1:100))
  dir.create(out,recursive=TRUE)
  warnings_seen<-character();render_checks<-list()
  limits<-if(min(data$r)<0)c(-1,1)else c(0,1)
  find_rasters<-function(g){
    answer<-list()
    if(inherits(g,'rastergrob'))answer<-list(g)
    if(!is.null(g$children))for(ch in g$children)answer<-c(answer,find_rasters(ch))
    if(!is.null(g$grobs))for(ch in g$grobs)answer<-c(answer,find_rasters(ch))
    answer
  }
  withCallingHandlers({
    plots<-lapply(panels,function(e){
      p<-ggplot2::ggplot(data[Panel==e],ggplot2::aes(Q1/100,Q2/100,fill=r))+
        ggplot2::geom_raster()+
        ggplot2::scale_fill_gradient2(name='Mean correlation',low='#2166ac',mid='white',high='#b2182b',limits=limits,midpoint=mean(limits),na.value='#eeeeee')+
        ggplot2::scale_x_continuous(breaks=c(.01,.25,.5,.75,1),expand=c(0,0))+
        ggplot2::scale_y_continuous(breaks=c(.01,.25,.5,.75,1),expand=c(0,0))+
        ggplot2::coord_equal(xlim=c(.005,1.005),ylim=c(.005,1.005),expand=FALSE)+
        ggplot2::labs(title=e,x='Quantile level',y='Quantile level')+
        ggplot2::theme_bw(base_size=10)+ggplot2::theme(plot.title=ggplot2::element_text(hjust=.5,size=11),legend.position='bottom',
          legend.key.width=grid::unit(.65,'cm'),legend.key.height=grid::unit(.3,'cm'),legend.margin=ggplot2::margin(0,0,0,0),
          axis.text=ggplot2::element_text(color='black'),panel.grid=ggplot2::element_blank())
      built<-ggplot2::ggplot_build(p)$data[[1]]
      stopifnot(nrow(built)==10000L,all(is.finite(built$x)),all(is.finite(built$y)),!anyNA(built$fill),sum(built$x==1)==100L,sum(built$y==1)==100L)
      grob<-ggplot2::ggplotGrob(p+ggplot2::theme(legend.position='none'))
      rasters<-find_rasters(grob)
      if(length(rasters)!=1L || !identical(dim(rasters[[1]]$raster),c(100L,100L)))stop('Incomplete raster panel: ',e)
      if(any(is.na(rasters[[1]]$raster)))stop('Unexpected missing raster colors: ',e)
      render_checks[[e]]<<-list(tiles=nrow(built),q100_row_tiles=sum(built$y==1),q100_column_tiles=sum(built$x==1),raster_dimensions=dim(rasters[[1]]$raster),raster_missing=sum(is.na(rasters[[1]]$raster)))
      p
    })
    legend<-cowplot::get_legend(plots[[1]])
    grid<-cowplot::plot_grid(plotlist=lapply(plots,function(p)p+ggplot2::theme(legend.position='none')),ncol=4)
    figure<-cowplot::plot_grid(grid,legend,ncol=1,rel_heights=c(1,.07))
    for(ext in c('png','jpg','pdf'))ggplot2::ggsave(file.path(out,paste0('Figure8.',ext)),figure,width=8,height=6.5,units='in',dpi=600,bg='white',limitsize=FALSE)
  },warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
if(length(warnings_seen))stop(paste(warnings_seen,collapse='; '))
jsonlite::write_json(render_checks,file.path(out,'render_checks.json'),pretty=TRUE,auto_unbox=TRUE)
