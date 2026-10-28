a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4)stop('Usage: Rscript 9_spatial_expression.R REGION_RDS SAMPLE_CSV CANDIDATES_CSV NEW_OUTPUT_DIR')
suppressPackageStartupMessages({library(data.table);library(edgeR)})
stopifnot(all(file.exists(a[1:3])))
if(file.exists(a[4]))stop('Output already exists')
lst<-readRDS(a[1]);meta<-fread(a[2],na.strings=c('','NA'));cand<-fread(a[3])
stopifnot(!anyDuplicated(meta$donor_id),!anyDuplicated(meta$file),!anyDuplicated(cand$gene),!anyNA(match(names(lst),meta$file)))
if(!'bulk_direction'%in%names(cand))cand[,bulk_direction:=sign(logFC)]
regs<-c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
forms<-list(primary=~group,sex_complete=~sex+group,age_sex_complete=~age+sex+group,depth=~depth_z+group,exclude_comorbidity=~group)
stopifnot(all(vapply(lst,function(x)all(regs%in%colnames(x)),TRUE)))
common<-Reduce(intersect,lapply(lst,rownames));dir.create(a[4],recursive=TRUE);out<-normalizePath(a[4])
wr<-function(x,n)fwrite(x,file.path(out,paste0(n,'.csv')))
de<-list();diagnostics<-list();profiles<-list()
for(reg in regs){
 md<-copy(meta[match(names(lst),file)]);md[,group:=factor(diagnosis,levels=c('Control','ALS'))]
 md[,profile_id:=paste(donor_id,reg,sep='___')]
 x<-sapply(lst,function(z)z[common,reg]);rownames(x)<-common;colnames(x)<-md$profile_id
 stopifnot(all(is.finite(x)),all(x>=0),all(x==round(x)))
 md[,`:=`(library_size=colSums(x),detected_genes=colSums(x>0))];stopifnot(all(md$library_size>0))
 md[,depth_z:=as.numeric(scale(log10(library_size)))];md[,region:=reg];profiles[[reg]]<-md
 y0<-DGEList(x);keep<-filterByExpr(y0,design=model.matrix(~group,md));y0<-calcNormFactors(y0[keep,,keep.lib.sizes=FALSE])
 for(mod in names(forms)){
  ix<-switch(mod,sex_complete=which(!is.na(md$sex)),age_sex_complete=which(complete.cases(md[,.(age,sex)])),exclude_comorbidity=which(md$sensitivity_eligible),seq_len(nrow(md)))
  mm<-md[ix];d<-model.matrix(forms[[mod]],mm)
  stopifnot(nrow(d)==nrow(mm),qr(d)$rank==ncol(d),all(table(mm$group)>=3),nrow(d)>ncol(d))
  y<-calcNormFactors(y0[,ix]);y<-estimateDisp(y,d,robust=TRUE)
  q<-glmQLFTest(glmQLFit(y,d,robust=TRUE),coef='groupALS')
  z<-as.data.table(topTags(q,n=Inf,sort.by='none')$table,keep.rownames='gene');z[,`:=`(region=reg,model=mod)];de[[length(de)+1]]<-z
  diagnostics[[length(diagnostics)+1]]<-data.table(region=reg,model=mod,n_Control=sum(mm$group=='Control'),n_ALS=sum(mm$group=='ALS'),df=nrow(d)-ncol(d),genes=nrow(y))
 }
 message('Completed ',reg)
}
all<-rbindlist(de);all[,FDR_all_regions:=p.adjust(PValue,'BH'),by=model]
fwrite(all,file.path(out,'DE_all_genes.csv.gz'))
cg<-merge(CJ(gene=cand$gene,region=regs,model=names(forms)),all,by=c('gene','region','model'),all.x=TRUE)
cg<-merge(cg,cand[,.(gene,bulk_direction)],by='gene')
adj<-function(p){q<-p.adjust(ifelse(is.finite(p),p,1),'BH');q[!is.finite(p)]<-NA_real_;q}
cg[,FDR_candidates:=adj(PValue),by=model];cg[,same_as_bulk:=sign(logFC)==bulk_direction]
st<-cg[,.(n_test=sum(is.finite(PValue)),same_direction=all(same_as_bulk%in%TRUE),max_FDR=max(FDR_candidates)),by=.(gene,region)]
st[,supported:=n_test==length(forms)&same_direction&is.finite(max_FDR)&max_FDR<.05]
wr(cg,'spatial_candidates_all_models');wr(st,'spatial_stability');wr(st[supported==TRUE],'spatial_supported_pairs')
wr(rbindlist(diagnostics),'model_design');wr(rbindlist(profiles),'donor_region_metadata')
wr(data.table(input=normalizePath(a[1:3]),md5=unname(tools::md5sum(a[1:3]))),'input_manifest')
wr(data.table(candidate_genes=nrow(cand),regions=length(regs),family_per_model=nrow(cand)*length(regs),supported_pairs=sum(st$supported),supported_genes=uniqueN(st[supported==TRUE,gene])),'summary')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
