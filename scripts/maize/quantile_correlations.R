# Rscript scripts/maize/quantile_correlations.R CONFIG.json NEW_OUTPUT
suppressPackageStartupMessages(library(data.table));setDTthreads(1L)
args<-commandArgs(TRUE);if(length(args)!=2L)stop('Usage: CONFIG.json NEW_OUTPUT')
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)),winslash='/')
config<-jsonlite::read_json(args[1],simplifyVector=TRUE)
resolve<-function(entry){x<-entry$path;if(!grepl('^([A-Za-z]:[/\\\\]|/)',x))x<-file.path(dirname(args[1]),x);stopifnot(unname(tools::md5sum(x))==entry$md5);x}
predictor<-resolve(config$inputs$blue);support<-resolve(config$inputs$support)
policy<-file.path(dirname(script),'../../src/nested_cv.R');policy_env<-new.env();sys.source(policy,policy_env)
out<-args[2];if(file.exists(out))stop('Use a new output directory');dir.create(out,recursive=TRUE)
  envs<-c('IAH4.2021','MIH1.2020','MNH1.2020','MNH1.2021','TXH1.2020','TXH2.2021','WIH1.2021','WIH2.2020','WIH2.2021','WIH3.2020','WIH3.2021')
  invalid_dates<-c(-33,-9,-3,0,1,5)
  message('Read corrected genotype BLUEs and observed raw-support means')
  x<-fread(predictor);x<-x[Env %in% envs & !DAP %in% invalid_dates & Vegetation.Index=='ExGR']
  s<-fread(support);s<-s[Env %in% envs & !DAP %in% invalid_dates & Vegetation.Index=='ExGR']
  key<-c('Env','Pedigree','DAP','Probability','Vegetation.Index')
  stopifnot(!anyDuplicated(x[,..key]),!anyDuplicated(s[,..key]))
  d<-merge(x,s,by=key,all.x=TRUE,sort=FALSE);rm(x,s);gc()
  if(anyNA(d$N.Raw.Finite))stop('Missing raw support for a genotype BLUE')
  d[,Q:=as.integer(round(Probability*100))]
  stopifnot(all(abs(d$Q/100-d$Probability)<1e-10),all(d$Q %in% 1:100),setequal(unique(d$Env),envs))
  groups<-unique(d[,.(Env,DAP)])[order(match(Env,envs),DAP)]
  matrices<-vector('list',nrow(groups));names(matrices)<-paste(groups$Env,groups$DAP,sep='|')
  diagnostics<-screening<-list();eps<-1e-12
  for(i in seq_len(nrow(groups))){
    g<-groups[i];a<-d[Env==g$Env & DAP==g$DAP]
    blue<-dcast(a,Pedigree~Q,value.var='Quantile.BLUE')
    raw<-dcast(a,Pedigree~Q,value.var='Quantile.Raw.Mean')
    stopifnot(identical(blue$Pedigree,raw$Pedigree),all(as.character(1:100)%in%names(blue)))
    b<-as.matrix(blue[,as.character(1:100),with=FALSE]);r<-as.matrix(raw[,as.character(1:100),with=FALSE])
    b[!is.finite(b)]<-NA_real_;r[!is.finite(r)]<-NA_real_
    bad_blue<-vapply(seq_len(100),function(j)policy_env$v15_numeric_constant(b[,j]),logical(1))
    ss<-t(vapply(seq_len(100),function(j)policy_env$v15_support_stats(r[,j]),numeric(7)))
    keep<-!bad_blue & !as.logical(ss[,'Drop']) & colSums(is.finite(b))>=3L
    b[,!keep]<-NA_real_
    n<-crossprod(is.finite(b));c<-suppressWarnings(stats::cor(b,use='pairwise.complete.obs'))
    c[n<3L]<-NA_real_;c[!is.finite(c)]<-NA_real_
    stopifnot(all(c[is.finite(c)]>=-1-1e-12),all(c[is.finite(c)]<=1+1e-12))
    endpoints<-sum(is.finite(c)&abs(c)>=1-eps & row(c)!=col(c))/2
    z<-c;z[]<-atanh(pmax(-1+eps,pmin(1-eps,as.numeric(c))));diag(z)<-NA_real_
    matrices[[i]]<-list(r=c,z=z,n=n,keep=keep)
    diagnostics[[i]]<-data.table(Env=g$Env,DAP=g$DAP,N_Genotypes=nrow(b),Retained_Quantiles=sum(keep),Excluded_Quantiles=sum(!keep),Endpoint_Pairs_Clipped=endpoints,Finite_OffDiagonal_Pairs=sum(is.finite(z))/2,Min_Pair_N=if(any(is.finite(z)))min(n[is.finite(z)]) else NA_integer_)
    screening[[i]]<-data.table(Env=g$Env,DAP=g$DAP,Q=1:100,Probability=(1:100)/100,Retained=keep,BLUE_Constant=bad_blue,Raw_N=ss[,'N'],Raw_UniquePct=ss[,'UniquePct'],Raw_FrequencyRatio=ss[,'FrequencyRatio'],Raw_Support_Drop=as.logical(ss[,'Drop']))
    message('Correlations ',i,'/',nrow(groups),' ',g$Env,' DAP',g$DAP,' retained ',sum(keep))
  }
  average<-function(indices){
    array<-simplify2array(lapply(matrices[indices],`[[`,'z'))
    if(length(indices)==1L)array<-array(array,dim=c(100,100,1))
    z<-apply(array,1:2,function(v)if(all(!is.finite(v)))NA_real_ else mean(v[is.finite(v)]))
    result<-tanh(z);validq<-Reduce(`|`,lapply(matrices[indices],`[[`,'keep'));diag(result)<-ifelse(validq,1,NA_real_)
    data.table(Q1=rep(1:100,100),Q2=rep(1:100,each=100),r=as.vector(result),N_Date_Matrices=as.vector(apply(is.finite(array),1:2,sum)))
  }
  data<-rbindlist(lapply(c(envs,'Combined'),function(e){v<-average(if(e=='Combined')seq_len(nrow(groups)) else which(groups$Env==e));v[,Panel:=e];v}))
  stopifnot(nrow(data)==120000L,!anyDuplicated(data[,.(Panel,Q1,Q2)]))
  ranges<-data[,.(Min=min(r,na.rm=TRUE),Max=max(r,na.rm=TRUE),Missing_Cells=sum(!is.finite(r)),Matrices=max(N_Date_Matrices)),by=Panel]
  limits<-if(min(data$r,na.rm=TRUE)<0)c(-1,1)else c(0,1)
  fwrite(data,file.path(out,'Figure8_Plotted_Correlations_V5.csv'))
  fwrite(ranges,file.path(out,'Figure8_Panel_Ranges_V5.csv'))
  fwrite(rbindlist(diagnostics),file.path(out,'Figure8_Date_Diagnostics_V5.csv'))
  fwrite(rbindlist(screening),file.path(out,'Figure8_Quantile_Support_V5.csv'))
  saveRDS(list(groups=groups,matrices=matrices),file.path(out,'Figure8_Per_Date_Matrices_V5.rds'))
jsonlite::write_json(list(status='PASS',dates=nrow(groups),minimum_pair_count=3L,endpoint_clip_epsilon=eps,environments=envs,source_md5=as.list(tools::md5sum(c(predictor,support,policy,script)))),file.path(out,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
