importance<-fread(inputs[['importance']]);stopifnot(all(importance$N_Outer_Fits==25L))
base_theme<-theme_gray(base_size=10,base_family="Arial")+
 theme(strip.background=element_rect(fill="#D9D9D9",colour=NA),strip.text=element_text(colour="black",size=9),
       panel.grid.major=element_line(colour="white",linewidth=.45),panel.grid.minor=element_line(colour="white",linewidth=.25),
       axis.title=element_text(colour="black"),axis.text=element_text(colour="#333333"),plot.margin=margin(5,7,5,5))
importance_plot<-function(en,day=NULL) {
 d<-copy(importance[Env==en]);if(!is.null(day))d<-d[DAP==day]
 d[,Scenario:=factor(Scenario,levels=c("Median","all"),labels=c("Median","All"))]
 title<-if(is.null(day))en else paste0(en,"  DAP ",day)
 ggplot(d,aes(Probability,Mean_Importance,colour=Mean_Importance))+
   geom_point(size=if(is.null(day)).85 else 1.5)+
   facet_grid(Scenario~DAP)+scale_colour_viridis_c(option="H",limits=c(0,100),breaks=c(0,25,50,75,100),
     name="Mean\nimportance",guide=guide_colourbar(barheight=grid::unit(2.4,"cm"),barwidth=grid::unit(.3,"cm")))+
   scale_x_continuous(limits=c(0,1),breaks=if(is.null(day))c(0,.5,1)else c(0,.25,.5,.75,1),expand=expansion(mult=c(.025,.025)))+
   scale_y_continuous(limits=c(0,103),breaks=c(0,25,50,75,100),expand=expansion(mult=c(.02,0)))+
   labs(x="Quantile level",y="Mean relative coefficient importance",title=title)+base_theme+
   theme(panel.border=element_rect(colour="black",fill=NA,linewidth=.35),
      strip.background=element_rect(fill="#D9D9D9",colour="black",linewidth=.35),
      strip.text=element_text(colour="black",size=8),
      strip.text.x=element_text(margin=margin(2,1,2,1)),strip.text.y=element_text(margin=margin(1,2,1,2)),
      axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=7.5),axis.text.y=element_text(size=8),
      panel.spacing=grid::unit(.055,"in"),plot.title=element_text(size=13),
      legend.title=element_text(size=9),legend.text=element_text(size=8))
}

dir.create(out,recursive=TRUE)
artifacts<-list()
export_figure<-function(plot,name,width,height,caption="") {
 files<-file.path(out,paste0(name,c(".png",".jpg",".pdf")))
 if(any(file.exists(files)))stop("Refusing figure overwrite")
 ragg::agg_png(files[1],width=width,height=height,units="in",res=600,background="white");print(plot);dev.off()
 ragg::agg_jpeg(files[2],width=width,height=height,units="in",res=600,quality=95,background="white");print(plot);dev.off()
 grDevices::cairo_pdf(files[3],width=width,height=height,family="Arial",onefile=TRUE);print(plot);dev.off()
 artifacts[[name]]<<-list(files=as.list(files),width_inches=width,height_inches=height,raster_dpi=600,caption=caption)
}
export_figure(importance_plot('TXH1.2020'),'Figure7',11,4.3)
supp_envs<-setdiff(envs,'TXH1.2020')
for(i in seq_along(supp_envs))export_figure(importance_plot(supp_envs[i]),paste0('FigureS',14L+i),11,4.3)
export_figure(importance_plot('TXH1.2020',68L),'FigureS25',6,5.5)
fwrite(importance,file.path(out,'Feature_Importance_Plotted.csv'))
