# Rscript scripts/maize/prepare_feature_figures.R V15_RUN_ROOT NEW_OUTPUT
library(data.table);setDTthreads(1)
a<-commandArgs(TRUE);stopifnot(length(a)==2L);prod<-a[1];out<-a[2]
if(file.exists(out))stop('Use a new output directory')
envs<-c('IAH4.2021','MIH1.2020','MNH1.2020','MNH1.2021','TXH1.2020','TXH2.2021','WIH1.2021','WIH2.2020','WIH2.2021','WIH3.2020','WIH3.2021')
feature_summaries<-list();validation<-list()
for(en in envs) {
 c<-fread(file.path(prod,en,"Outer_Coefficients_V15.csv"));q<-fread(file.path(prod,en,"Input_QA_V15.csv"))
 stopifnot(all(c$Env==en),setequal(c$Scenario,c("Median","all")),!anyDuplicated(c[,.(Scenario,Outer_Repeat,Outer_Fold,Feature)]))
 c<-c[Feature!="(Intercept)"]
 stopifnot(all(is.finite(c$Coefficient_Standardized)),all(is.finite(c$Importance)),all(c$Selected==(c$Coefficient_Standardized!=0)),
           all(grepl("^DAP-?[0-9]+_Q[0-9]{3}_ExGR$",c$Feature)))
 c[,Recomputed:=if(max(abs(Coefficient_Standardized))>0)100*abs(Coefficient_Standardized)/max(abs(Coefficient_Standardized))else 0,
   by=.(Scenario,Outer_Repeat,Outer_Fold)]
 err<-max(abs(c$Importance-c$Recomputed));stopifnot(err<1e-9)
 c[,`:=`(DAP=as.integer(sub("^DAP(-?[0-9]+)_.*$","\\1",Feature)),ProbPct=as.integer(sub(".*_Q([0-9]{3})_ExGR$","\\1",Feature)))]
 stopifnot(all(c$ProbPct%in%1:100),all(c[Scenario=="Median",ProbPct]==50L))
 fitcounts<-c[,.N,by=.(Scenario,Outer_Repeat,Outer_Fold)]
 stopifnot(nrow(fitcounts)==50L,all(fitcounts[Scenario=="all",N]==q$N_Features),all(fitcounts[Scenario=="Median",N]==q$N_Features/100))
 s<-c[,.(Mean_Importance=mean(Recomputed),SD_Importance=sd(Recomputed),Selection_Frequency_Pct=100*mean(Selected),
           Mean_Abs_Standardized_Coefficient=mean(abs(Coefficient_Standardized)),N_Outer_Fits=.N),
       by=.(Env,Scenario,Feature,DAP,ProbPct)]
 stopifnot(all(s$N_Outer_Fits==25L),all(s$Mean_Importance>=0 & s$Mean_Importance<=100),
           all(s[Scenario=="all",.N,by=DAP]$N==100L))
 s[,Probability:=ProbPct/100];feature_summaries[[en]]<-s
 validation[[en]]<-data.table(Env=en,Outer_Fits=nrow(fitcounts),Candidates_All=q$N_Features,
   Summary_Rows=nrow(s),Max_Importance_Reconstruction_Error=err,Each_Feature_Fits=25L)
}
importance<-rbindlist(feature_summaries)
dir.create(out,recursive=TRUE)
fwrite(importance,file.path(out,'Feature_Importance_Plotted.csv'))
fwrite(rbindlist(validation),file.path(out,'Input_Validation.csv'))
