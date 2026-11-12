dir.create(out,recursive=TRUE)
m<-fread(inputs[['cot_metrics']]);v<-fread(inputs[['cot_components']]);mz<-fread(inputs[['maize_components']])
vis<-c('PSRI','CVI','ExG2')
stopifnot(nrow(m)==30300L,nrow(v)==90900L)
palette<-c(Genotype="#F8766D",Block="#00B0F6",Range="#A3A500",Pass="#00BF7D",Replicate="#00B0F6",Residual="#E76BF3")
base_theme<-theme_bw(base_size=9)+theme(legend.position="right",legend.box="vertical",legend.title=element_text(size=9),
 legend.text=element_text(size=8),strip.background=element_rect(fill="grey85"),strip.text=element_text(size=9),
 axis.text=element_text(size=8),axis.title=element_text(size=9),panel.grid.minor=element_blank(),panel.spacing.x=grid::unit(4,"mm"),plot.margin=margin(3,4,3,3))
prepare<-function(d,day,groups,x="Probability") {
  d<-copy(d);d[grp=="Pedigree",grp:="Genotype"];d[,grp:=factor(grp,levels=groups)]
  d[,Vegetation.Index:=factor(Vegetation.Index,levels=if(all(unique(Vegetation.Index)%in%vis))vis else unique(Vegetation.Index))]
  keys<-c(intersect("Env",names(d)),day,"Vegetation.Index",x)
  stopifnot(!anyNA(d$grp),all(is.finite(d$H2Individual)),all(d$H2Individual>=0&d$H2Individual<=1))
  if("PercentUnrounded"%in%names(d))d[,PlotPercent:=PercentUnrounded]else d[,PlotPercent:=100*vcov/sum(vcov),by=keys]
  stopifnot(max(abs(d[,sum(PlotPercent),by=keys]$V1-100))<1e-9)
  d[,StackOrder:=match(as.character(grp),rev(groups))];setorderv(d,c(keys,"StackOrder"))
  d[,`:=`(Ymax=cumsum(PlotPercent),Ymin=shift(cumsum(PlotPercent),fill=0)),by=keys]
  curvekeys<-setdiff(keys,x)
  pts<-unique(d[,c(curvekeys,x),with=FALSE]);setorderv(pts,c(curvekeys,x))
  if(x=="Probability")pts[,Segment:=cumsum(c(TRUE,diff(get(x))>.010000001)),by=curvekeys] else pts[,Segment:=1L]
  d<-merge(d,pts,by=c(curvekeys,x),all.x=TRUE,sort=FALSE)
  d[,PlotGroup:=interaction(grp,Segment,drop=TRUE)];d
}
area_plot<-function(d,x,label) {
 ggplot(d,aes(x=.data[[x]]))+
  geom_ribbon(aes(ymin=Ymin,ymax=Ymax,fill=grp,group=PlotGroup),colour=NA)+
  geom_line(data=d[grp=="Genotype"],aes(y=H2Individual*100,group=Segment,colour="Heritability x 100"),linewidth=.45)+
  scale_fill_manual(values=palette,drop=TRUE,name="Variance\ncomponent")+
  scale_colour_manual(values=c("Heritability x 100"="black"),name=NULL)+
  guides(colour=guide_legend(order=1),fill=guide_legend(order=2))+
  scale_y_continuous("Variance component (%)",breaks=c(0,50,100),expand=expansion(mult=c(0,.01)))+
  labs(x=label)+coord_cartesian(ylim=c(0,100))+base_theme
}
cota<-prepare(v[DAT%in%c(57,86,96,113,130)],"DAT",c("Genotype","Block","Residual"))
pa<-area_plot(cota,"Probability","Quantile level within DAT")+facet_grid(Vegetation.Index~DAT)+scale_x_continuous(breaks=c(0,.5,1),limits=c(0,1))
mzb<-mz[Env=="TXH1.2020"&Vegetation.Index%in%c("ExGR","GminusR","COM1","EXR")&DAP%in%c(47,68,88,96,111)]
mzb[,Vegetation.Index:=factor(Vegetation.Index,levels=c("ExGR","GminusR","COM1","EXR"))]
mzb<-prepare(mzb,"DAP",c("Genotype","Range","Pass","Replicate","Residual"))
mzb[,Vegetation.Index:=factor(as.character(Vegetation.Index),levels=c("ExGR","GminusR","COM1","EXR"))]
pb<-area_plot(mzb,"Probability","Quantile level within DAP")+facet_grid(Vegetation.Index~DAP)+scale_x_continuous(breaks=c(0,.5,1),limits=c(0,1))
season<-m[,.(Median=median(H2Individual),Q25=quantile(H2Individual,.25),Q75=quantile(H2Individual,.75),Dates=.N),by=.(Vegetation.Index,Probability)]
stopifnot(all(season$Dates==101L));season[,Vegetation.Index:=factor(Vegetation.Index,levels=vis)]
pc<-ggplot(season,aes(Probability,Median))+geom_ribbon(aes(ymin=Q25,ymax=Q75),fill="grey80")+
 geom_line(linewidth=.45)+geom_vline(xintercept=.5,linetype="dotted",colour="grey40")+facet_wrap(~Vegetation.Index,nrow=1)+
 scale_x_continuous("Quantile level",breaks=c(0,.25,.5,.75,1),limits=c(0,1))+
 scale_y_continuous("Median H2",limits=c(0,1),breaks=c(0,.5,1))+base_theme+theme(legend.position="none")
fig4<-plot_grid(pa,pb,pc,ncol=1,labels=c("A","B","C"),label_size=15,rel_heights=c(.36,.43,.21),align="v",axis="lr")
save_plot_files<-function(p,name,width,height) {
  if(grepl("^Supporting_Figure",name))name<-sub("^Supporting_Figure([0-9]+)_.*$","FigureS\\1",name)
  for(ext in c("png","jpg","pdf"))ggsave(file.path(out,paste0(name,".",ext)),p,width=width,height=height,units="in",dpi=600,bg="white",limitsize=FALSE)
  # Smaller native-render preview is an additional output, not a replacement.
  ggsave(file.path(out,paste0(name,"_Preview.png")),p,width=width,height=height,units="in",dpi=160,bg="white")
}
save_plot_files(fig4,"Figure4",9,9)
area_plot<-function(d,x,label) {
 ggplot(d,aes(x=.data[[x]]))+
  geom_ribbon(aes(ymin=Ymin,ymax=Ymax,fill=grp,group=PlotGroup),colour=NA)+
  geom_line(data=d[grp=="Genotype"],aes(y=H2Individual*100,group=Segment,colour="Individual-unit H2 x 100"),linewidth=.45)+
  scale_fill_manual(values=palette,drop=TRUE,name="Variance\ncomponent")+
  scale_colour_manual(values=c("Individual-unit H2 x 100"="black"),name=NULL)+
  guides(colour=guide_legend(order=1),fill=guide_legend(order=2))+
  scale_y_continuous("Variance component (%)",breaks=c(0,50,100),expand=expansion(mult=c(0,.01)))+
  labs(x=label)+coord_cartesian(ylim=c(0,100))+base_theme
}

cotmed<-prepare(v[QuantileID==50L],"DAT",c("Genotype","Block","Residual"),x="DAT")
s1<-area_plot(cotmed,"DAT","Days after transplanting")+facet_wrap(~Vegetation.Index,ncol=1)+scale_x_continuous(breaks=c(35,60,90,120,140))
save_plot_files(s1,"Supporting_Figure1_Cotton_Median_Within_Date_V1",8,6.5)
times<-fread(inputs[["cotton_flight"]]);stopifnot(!anyDuplicated(times$DAT))
ft<-merge(m[QuantileID==50L],times[,.(DAT,Flight_Time)],by="DAT",all.x=TRUE)
ft<-ft[!is.na(Flight_Time)&nzchar(Flight_Time)]
stopifnot(all(grepl("^[0-9]{1,2}:[0-9]{2}$",ft$Flight_Time)))
ft[,TimeHours:=as.numeric(sub(":.*","",Flight_Time))+as.numeric(sub(".*:","",Flight_Time))/60]
ft[,Vegetation.Index:=factor(Vegetation.Index,levels=vis)]
ftime_plot<-function(d,ncol) {
 z<-d[,{test<-cor.test(TimeHours,H2Individual);.(N=.N,Pearson=unname(test$estimate),P=test$p.value)},by=Vegetation.Index]
 z[,Label:=sprintf("n = %d; r = %.2f; p = %.3f",N,Pearson,P)]
 p<-ggplot(d,aes(TimeHours,H2Individual))+geom_point(size=1.5)+geom_smooth(method="lm",se=FALSE,colour="#3465A4",linewidth=.6)+
 geom_text(data=z,aes(x=-Inf,y=Inf,label=Label),inherit.aes=FALSE,hjust=-.03,vjust=1.4,size=2.8)+
 facet_wrap(~Vegetation.Index,ncol=ncol)+scale_x_continuous("Flight time (hour of day)")+scale_y_continuous("Individual-unit H2",limits=c(0,1))+
 base_theme+theme(legend.position="none")
 list(plot=p,statistics=z)
}
ct<-ftime_plot(ft,3);save_plot_files(ct$plot,"Supporting_Figure2_Cotton_Flight_Time_V1",9,3.2)
mt<-ftime_plot(fread(inputs[["maize_flight"]]),2);save_plot_files(mt$plot,"Supporting_Figure3_Maize_Flight_Time_V1",8,5.5)
for(i in seq_along(envs)) {
 e<-envs[i];d<-prepare(mz[Env==e&Vegetation.Index=="ExGR"],"DAP",c("Genotype","Range","Pass","Replicate","Residual"))
 p<-area_plot(d,"Probability","Quantile level within DAP")+facet_wrap(~DAP,ncol=5)+scale_x_continuous(breaks=c(0,.5,1),limits=c(0,1))+ggtitle(e)
 save_plot_files(p,sprintf("Supporting_Figure%d_Maize_ExGR_%s_V1",i+3L,e),10,6.5)
}
fwrite(season,file.path(out,"Figure4C_Across_Date_Quantile_Summary_V1.csv"));fwrite(ct$statistics,file.path(out,"Cotton_Flight_Time_Statistics_V1.csv"));fwrite(ft,file.path(out,"Cotton_Flight_Time_Plot_Data_V1.csv"))
fwrite(mt$statistics,file.path(out,"Maize_Flight_Time_Statistics_V1.csv"))
fwrite(cota,file.path(out,"Figure4A_Plot_Data_V1.csv"));fwrite(mzb,file.path(out,"Figure4B_Plot_Data_V1.csv"))
