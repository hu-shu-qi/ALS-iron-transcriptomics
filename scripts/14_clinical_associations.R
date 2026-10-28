a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=2)stop('Usage: Rscript 14_clinical_associations.R PRIVATE_WORKBOOK NEW_OUTPUT_DIR')
if(file.exists(a[2]))stop('Output exists')
stopifnot(file.exists(a[1]))
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
input<-normalizePath(a[1]);dir.create(a[2],recursive=TRUE);out<-normalizePath(a[2])
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(readxl))
raw<-as.data.frame(read_excel(input,sheet='Data1',col_types='text'));names(raw)<-trimws(names(raw))
pick<-function(p){j<-grep(p,names(raw));stopifnot(length(j)==1);raw[[j]]}
num<-function(x)suppressWarnings(as.numeric(x))
flag<-function(x)ifelse(is.na(x)|trimws(x)=='',0L,ifelse(x=='No',0L,ifelse(grepl('^Yes',x),1L,NA_integer_)))
d<-data.frame(ferritin=num(pick('^Serum_ferritin \\(')),iron=num(pick('^serum_iron')),tsat=num(pick('^transferrin_saturation')),age=num(pick('^Blood_collection_age')),sex=factor(trimws(pick('^Sex ')),levels=c('F','M')),onset=factor(trimws(pick('^Onset_site')),levels=c('Spinal','Bulbar','Mixture')),alsfrs=num(pick('^ALSFRS-R_baseline$')),followup=num(pick('^ALSFRS-R_follow_up$')),days=num(pick('^follow-up_interval')),anemia=flag(pick('^anemia$')),inflammation=flag(pick('^Inflammation$')),supplements=flag(pick('^iron_supplements$')))
id<-trimws(pick('^Number$'));stopifnot(!anyNA(id),!anyDuplicated(id))
d$y<-ifelse(is.finite(d$days)&d$days>0,(d$alsfrs-d$followup)/d$days*30.4,NA_real_)
recorded<- -num(pick('ALSFRS-R/month'))
stopifnot(all(is.na(d$y)==is.na(recorded)),max(abs(d$y-recorded),na.rm=TRUE)<1e-8,all(d$ferritin>0,na.rm=TRUE),all(d$iron>=0,na.rm=TRUE),all(d$tsat>=0 & d$tsat<=100,na.rm=TRUE))
write<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na='')
strict<-!is.na(d$anemia)&!is.na(d$inflammation)&!is.na(d$supplements)&d$anemia==0&d$inflammation==0&d$supplements==0
pdf(file.path(out,'model_diagnostics.pdf'),width=9,height=7)
par(mfrow=c(2,2));for(v in c('ferritin','iron','tsat','y'))hist(d[[v]],main=v,xlab=v,col='grey80')
set.seed(20260916);res<-list()
for(g in c('main','strict_sensitivity')){
 s<-if(g=='main')d else d[strict,]
 for(v in c('ferritin','iron','tsat')){
  ok<-is.finite(s[[v]])&is.finite(s$y);x<-s[[v]][ok];y<-s$y[ok];n<-length(x)
  ct<-cor.test(x,y,method='spearman',exact=FALSE)
  b<-replicate(5000,{j<-sample.int(n,n,replace=TRUE);cor(x[j],y[j],method='spearman')});stopifnot(length(b)==5000,all(is.finite(b)));ci<-quantile(b,c(.025,.975))
  res[[length(res)+1]]<-data.frame(analysis=g,marker=v,n=n,rho=unname(ct$estimate),ci_low=ci[1],ci_high=ci[2],p=ct$p.value)
 }
}
r<-do.call(rbind,res);r$q_BH_three<-ave(r$p,r$analysis,FUN=function(p)p.adjust(p,'BH'));write(r,'correlations.csv')
adj<-list();diag<-list()
for(v in c('ferritin','iron','tsat')){
 s<-d;s$x<-switch(v,ferritin=log2(s$ferritin),iron=s$iron/5,tsat=s$tsat/10)
 s<-s[complete.cases(s[,c('y','x','age','sex','alsfrs','onset')]),]
 fit<-lm(y~x+age+sex+alsfrs+onset,data=s)
 X<-model.matrix(fit);stopifnot(qr(X)$rank==ncol(X));bread<-solve(crossprod(X));h<-hatvalues(fit);e<-residuals(fit)
 V<-bread%*%crossprod(X,X*(e/(1-h))^2)%*%bread
 se<-sqrt(V['x','x']);beta<-coef(fit)['x'];df<-df.residual(fit);crit<-qt(.975,df)
 adj[[v]]<-data.frame(marker=v,n=nrow(s),beta=unname(beta),se_HC3=se,ci_low=unname(beta-crit*se),ci_high=unname(beta+crit*se),p=2*pt(-abs(beta/se),df),df=df,unit=switch(v,ferritin='per doubling',iron='per 5 umol/L',tsat='per 10 percentage points'))
 aux<-lm(I(e^2)~X[,-1]);bp<-length(e)*summary(aux)$r.squared
 diag[[v]]<-data.frame(marker=v,n=nrow(s),r_squared=summary(fit)$r.squared,max_leverage=max(h),n_cook_above_4_n=sum(cooks.distance(fit)>4/nrow(s)),max_cook=max(cooks.distance(fit)),residual_shapiro_p=shapiro.test(e)$p.value,variance_auxiliary_p=pchisq(bp,ncol(X)-1,lower.tail=FALSE))
 par(mfrow=c(2,2));plot(fit,which=c(1,2,3,5),id.n=0,sub.caption=paste(v,'model diagnostics'))
}
dev.off()
a<-do.call(rbind,adj);a$q_BH_three<-p.adjust(a$p,'BH');write(a,'adjusted_HC3.csv');write(do.call(rbind,diag),'model_diagnostics.csv')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(r,row.names=FALSE);print(a,row.names=FALSE);print(do.call(rbind,diag),row.names=FALSE)
flag_counts<-data.frame(field=c('anemia','inflammation','supplements'),missing=vapply(d[c('anemia','inflammation','supplements')],function(z)sum(is.na(z)),0L))
write.csv(flag_counts,file.path(out,'flag_missingness.csv'),row.names=FALSE)
write.csv(data.frame(seed=20260916,B=5000,unit='patient paired resampling',CI='percentile 2.5 and 97.5 percent; R quantile default type 7',rng=paste(RNGkind(),collapse='/'),id_unique=!anyDuplicated(id),clinical_rows=nrow(d),endpoint_rows=sum(is.finite(d$y))),file.path(out,'bootstrap_spec.csv'),row.names=FALSE)
write.csv(data.frame(path=c(input,normalizePath(script)),md5=unname(tools::md5sum(c(input,script)))),file.path(out,'run_inputs.csv'),row.names=FALSE)
