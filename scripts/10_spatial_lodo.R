a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=6)stop('Usage: Rscript 10_spatial_lodo.R REGION_RDS META_CSV CANDIDATES_CSV SPATIAL_RESULTS_CSV SUPPORTED_PAIRS_CSV NEW_OUTPUT_DIR')
suppressPackageStartupMessages({library(data.table);library(edgeR)})
stopifnot(all(file.exists(a[1:5])))
if(file.exists(a[6]))stop('Output already exists')
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
out<-a[6];dir.create(out,recursive=TRUE)
for(z in c('lodo','cache'))dir.create(file.path(out,z))
inputs<-normalizePath(a[1:5]);wr<-function(x,n)fwrite(x,file.path(out,'lodo',paste0(n,'.csv')))
hash <- tools::md5sum(inputs)
lst <- readRDS(inputs[1]); md <- fread(inputs[2],na.strings=c('','NA')); md <- md[match(names(lst),file)]
md[,group:=factor(diagnosis,levels=c('Control','ALS'))]
cand <- fread(inputs[3]); reference <- fread(inputs[4]); supported <- fread(inputs[5])
stopifnot(!anyDuplicated(cand$gene),!anyDuplicated(md$donor_id),!anyNA(md$donor_id))
if(!'FDR_candidates464'%in%names(reference))setnames(reference,'FDR_candidates','FDR_candidates464')
regs <- c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
forms <- list(primary=~group,sex_complete=~sex+group,age_sex_complete=~age+sex+group,depth=~depth_z+group,exclude_comorbidity=~group)
jobs <- list()
for(mod in names(forms)) {
 ix <- switch(mod,sex_complete=which(!is.na(md$sex)),age_sex_complete=which(complete.cases(md[,.(age,sex)])),exclude_comorbidity=which(md$sensitivity_eligible),seq_len(nrow(md)))
 jobs[[length(jobs)+1]] <- list(job=paste0(mod,'__full'),model=mod,ix=ix,omit='none')
 for(i in ix) jobs[[length(jobs)+1]] <- list(job=paste0(mod,'__drop_',md$donor_id[i]),model=mod,ix=setdiff(ix,i),omit=md$donor_id[i])
}
fprint <- c(hash,tools::md5sum(script))
res <- list(); designs <- list(); filters <- list()
for(reg in regs) {
 x <- sapply(lst,function(z)z[,reg]); rownames(x)<-rownames(lst[[1]])
 stopifnot(all(vapply(lst,function(z)identical(rownames(z),rownames(x)),TRUE)),all(x>=0),all(x==round(x)))
 mm0 <- copy(md); mm0[,depth_z:=as.numeric(scale(log10(colSums(x))))]
 keep <- filterByExpr(DGEList(x),design=model.matrix(~group,mm0))
 filters[[reg]] <- data.table(region=reg,gene=rownames(x),included=keep)
 xx <- x[keep,,drop=FALSE]
 for(j in jobs) {
  dest <- file.path(out,'cache',paste0(reg,'_',j$job,'.rds'))
  if(file.exists(dest)) { z<-readRDS(dest);stopifnot(identical(z$fingerprint,fprint)) } else {
   mm<-mm0[j$ix];d<-model.matrix(forms[[j$model]],mm)
   ok<-nrow(d)==nrow(mm) && qr(d)$rank==ncol(d) && all(table(mm$group)>=3) && nrow(d)>ncol(d)
   dg<-data.table(region=reg,job=j$job,model=j$model,omitted_donor=j$omit,n_Control=sum(mm$group=='Control'),n_ALS=sum(mm$group=='ALS'),rank=qr(d)$rank,df=nrow(d)-ncol(d),testable=ok,genes=sum(keep))
   if(ok) {
    y<-calcNormFactors(DGEList(xx[,j$ix,drop=FALSE])); y<-estimateDisp(y,d,robust=TRUE); fit<-glmQLFit(y,d,robust=TRUE)
    tt<-as.data.table(topTags(glmQLFTest(fit,coef='groupALS'),n=Inf,sort.by='none')$table,keep.rownames='gene')[gene%in%cand$gene]
   } else tt<-data.table(gene=character(),logFC=numeric(),PValue=numeric(),FDR=numeric())
   z<-list(table=tt,design=dg,fingerprint=fprint);saveRDS(z,dest)
   message(reg,' / ',j$job)
  }
  tt<-copy(z$table);tt[,`:=`(region=reg,job=j$job)];res[[length(res)+1]]<-tt;designs[[length(designs)+1]]<-z$design
 }
}
info<-rbindlist(lapply(jobs,function(j)data.table(job=j$job,model=j$model,omitted_donor=j$omit)))
a<-merge(CJ(gene=cand$gene,region=regs,job=info$job),rbindlist(res,fill=TRUE),by=c('gene','region','job'),all.x=TRUE)
a<-merge(a,info,by='job');a<-merge(a,cand[,.(gene,bulk_direction)],by='gene')
a[,FDR464:=p.adjust(ifelse(is.finite(PValue),PValue,1),'BH'),by=job];a[!is.finite(PValue),FDR464:=NA_real_]
a[,same_as_bulk:=sign(logFC)==bulk_direction]
baseline<-a[omitted_donor=='none']
ck<-merge(baseline,reference,by=c('gene','region','model'),suffixes=c('_lodo','_full'))
stopifnot(nrow(ck)==nrow(cand)*length(regs)*length(forms),identical(is.na(ck$PValue_lodo),is.na(ck$PValue_full)),max(abs(ck$logFC_lodo-ck$logFC_full),na.rm=TRUE)<1e-8,max(abs(ck$PValue_lodo-ck$PValue_full),na.rm=TRUE)<1e-8,max(abs(ck$FDR464-ck$FDR_candidates464),na.rm=TRUE)<1e-8)
lo<-merge(a[omitted_donor!='none'],baseline[,.(gene,region,model,baseline_logFC=logFC,baseline_FDR=FDR464)],by=c('gene','region','model'))
lo[,delta:=logFC-baseline_logFC]
su<-lo[,{
 ii<-which(is.finite(logFC)&is.finite(FDR464));k<-length(ii)
 list(n_planned=.N,n_tested=k,baseline_logFC=baseline_logFC[1],baseline_FDR=baseline_FDR[1],min_logFC=if(k)min(logFC[ii])else NA_real_,max_logFC=if(k)max(logFC[ii])else NA_real_,direction_retained_n=if(k)sum(sign(logFC[ii])==sign(baseline_logFC[ii]))else NA_integer_,bulk_direction_retained_n=if(k)sum(same_as_bulk[ii])else NA_integer_,FDR_retained_n=if(k)sum(FDR464[ii]<.05)else NA_integer_,max_FDR=if(k)max(FDR464[ii])else NA_real_,max_abs_delta=if(k)max(abs(delta[ii]))else NA_real_,most_influential_donor=if(k)omitted_donor[ii[which.max(abs(delta[ii]))]]else NA_character_)
 },by=.(gene,region,model)]
focus<-merge(su,supported[,.(gene,region)],by=c('gene','region'))
summary<-focus[,.(supported_pairs=.N,all_runs_tested=sum(n_tested==n_planned),all_directions_retained=sum(direction_retained_n==n_planned),all_FDR_retained=sum(FDR_retained_n==n_planned),min_FDR_retention=min(FDR_retained_n/n_planned),max_abs_logFC_change=max(max_abs_delta)),by=model]
wr(a,'all_candidate_results');wr(rbindlist(designs),'designs');wr(rbindlist(filters),'fixed_expression_filter');wr(lo,'leave_one_donor_out');wr(su,'all_pair_summary');wr(focus,'supported_24_pair_summary');wr(summary,'supported_summary')
wr(focus[n_tested<n_planned|direction_retained_n<n_planned|FDR_retained_n<n_planned],'sensitive_supported_pairs')
stopifnot(nrow(a)==nrow(cand)*length(regs)*length(jobs),nrow(focus)==nrow(supported)*length(forms),all(unname(tools::md5sum(inputs))==unname(hash)))
wr(data.table(path=inputs,md5=unname(hash)),'input_manifest')
writeLines(capture.output(sessionInfo()),file.path(out,'lodo','sessionInfo.txt'))
print(summary);message('LODO COMPLETE')
