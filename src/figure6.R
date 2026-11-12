palette<-setNames(RColorBrewer::brewer.pal(10L,'Paired')[1:2],c('Median','all'))
loeo_raw<-fread(inputs[['loeo_metrics']])[Scaling=='WITHIN'&Scenario%in%c('Median','all')]
stopifnot(nrow(loeo_raw)==22L,all(loeo_raw$SourceStatus=='historical_not_refit'))
loeo<-loeo_raw[,.(Mean_Pearson=mean(Pearson),Heldouts=.N,CorrectedRefits=0L,HistoricalNotRefit=.N),by=.(Scaling,Scenario)]
within<-fread(inputs[['nested_summary']])
reps<-fread(inputs[['nested_repeats']])
stopifnot(nrow(within)==22L,nrow(reps)==110L,setequal(within$Env,envs),setequal(within$Scenario,c("Median","all")),
 !anyDuplicated(reps[,.(Env,Scenario,Outer_Repeat)]),all(within$N_Repeats==5L),all(is.finite(within$Mean_Pearson)))
recomputed<-reps[,.(Recomputed=mean(Pearson)),by=.(Env,Scenario)]
stopifnot(max(abs(merge(within,recomputed,by=c("Env","Scenario"))[,Mean_Pearson-Recomputed]))<1e-12)
within_all<-within[,.(Mean_Pearson=mean(Mean_Pearson)),by=Scenario][,Env:="All Environments"]
bars<-rbindlist(list(within[,.(Env,Scenario,Mean_Pearson)],within_all),use.names=TRUE)
repeat_all<-reps[,.(Pearson=mean(Pearson)),by=.(Scenario,Outer_Repeat)][,Env:="All Environments"]
points<-rbindlist(list(reps[,.(Env,Scenario,Outer_Repeat,Pearson)],repeat_all),use.names=TRUE)
bars[,`:=`(Env=factor(Env,levels=c("All Environments",envs)),Scenario=factor(Scenario,levels=c("Median","all")))]
points[,`:=`(Env=factor(Env,levels=c("All Environments",envs)),Scenario=factor(Scenario,levels=c("Median","all")))]
points[,X:=as.integer(Scenario)+(Outer_Repeat-3L)*.045]
labels<-merge(bars,points[,.(LabelY=max(Pearson)+.035),by=.(Env,Scenario)],by=c("Env","Scenario"))
base_theme<-theme_gray(base_size=10,base_family="Arial")+
 theme(strip.background=element_rect(fill="#D9D9D9",colour=NA),strip.text=element_text(colour="black",size=9),
       panel.grid.major=element_line(colour="white",linewidth=.45),panel.grid.minor=element_line(colour="white",linewidth=.25),
       axis.title=element_text(colour="black"),axis.text=element_text(colour="#333333"),plot.margin=margin(5,7,5,5))
pA<-ggplot(bars,aes(x=Scenario,y=Mean_Pearson,fill=Scenario))+
 geom_col(width=.7)+geom_point(data=points,aes(x=X,y=Pearson),inherit.aes=FALSE,shape=21,size=.8,stroke=.25,fill="white",colour="black")+
 geom_text(data=labels,aes(y=LabelY,label=sprintf("%.3f",Mean_Pearson)),size=2.7,family="Arial",vjust=0)+
 facet_wrap(~Env,ncol=4)+scale_fill_manual(values=palette[1:2],labels=c("Median","All"),drop=FALSE)+
 scale_x_discrete(labels=c(Median="Median",all="All"))+
 scale_y_continuous(limits=c(-.05,1.02),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(0,.02)))+
 labs(x="Scenario",y="Pearson correlation",fill="Scenario")+base_theme+
 theme(legend.position="none",axis.text.x=element_text(size=9),panel.spacing=grid::unit(.16,"in"))
# Keep the 11 held-out environment estimates for each of the two scenarios.
# These rows are historical validated LOEO fits, not newly fitted V15 models.
loeo_selected<-loeo[Scaling=="WITHIN" & Scenario %chin% c("Median","all")]
loeo_selected_raw<-loeo_raw[Scaling=="WITHIN" & Scenario %chin% c("Median","all")]
stopifnot(nrow(loeo_selected)==2L,nrow(loeo_selected_raw)==22L,
 all(loeo_selected$Heldouts==11L),all(loeo_selected$CorrectedRefits==0L),
 all(loeo_selected$HistoricalNotRefit==11L),setequal(loeo_selected_raw$Heldout,envs),
 all(loeo_selected_raw[,uniqueN(Heldout),by=Scenario]$V1==11L))
loeo_selected[,Scenario:=factor(Scenario,levels=c("Median","all"))]
loeo_selected[,Panel_Label:="Across environments"]
pB<-ggplot(loeo_selected,aes(x=Scenario,y=Mean_Pearson,fill=Scenario))+
 geom_col(width=.7)+geom_text(aes(label=sprintf("%.3f",Mean_Pearson)),vjust=-.3,size=2.7,family="Arial")+
 facet_wrap(~Panel_Label,ncol=1)+scale_fill_manual(values=palette[1:2],drop=FALSE)+
 scale_x_discrete(labels=c(Median="Median",all="All"))+
 scale_y_continuous(limits=c(0,.45),breaks=seq(0,.4,.1),expand=expansion(mult=c(0,.02)))+
 labs(x="Scenario",y="Pearson correlation")+base_theme+
 theme(legend.position="none",axis.text.x=element_text(size=9),strip.text=element_text(size=9))
figure6<-plot_grid(pA,pB,ncol=1,rel_heights=c(.66,.34),labels=c("A)","B)"),label_size=15,label_fontfamily="Arial",label_x=.012,label_y=.995,hjust=0,vjust=1)

dir.create(out,recursive=TRUE,showWarnings=FALSE)
files<-file.path(out,paste0("Figure6",c(".png",".jpg",".pdf")))
stopifnot(!any(file.exists(files)))
ragg::agg_png(files[1],width=7.25,height=7.25,units="in",res=600,background="white");print(figure6);dev.off()
ragg::agg_jpeg(files[2],width=7.25,height=7.25,units="in",res=600,quality=95,background="white");print(figure6);dev.off()
grDevices::cairo_pdf(files[3],width=7.25,height=7.25,family="Arial",onefile=TRUE);print(figure6);dev.off()
fwrite(bars,file.path(out,"Figure6A_Bar_Data_V5.csv"))
fwrite(points,file.path(out,"Figure6A_Repeat_Point_Data_V5.csv"))
fwrite(loeo_selected,file.path(out,"Figure6B_WITHIN_Median_All_Data_V5.csv"))
fwrite(loeo_selected_raw,file.path(out,"Figure6B_Heldout_Metrics_V5.csv"))
