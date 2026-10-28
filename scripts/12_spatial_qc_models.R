a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4)stop('Usage: Rscript 12_spatial_qc_models.R QC_DIR CONFIG_CSV CANDIDATES_CSV NEW_OUTPUT_DIR')
if(file.exists(a[4]))stop('Output exists')
suppressPackageStartupMessages({library(data.table);library(edgeR)})
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
cfg<-fread(a[2],na.strings=c('','NA'));cand<-fread(a[3]);qc_dir<-normalizePath(a[1]);out<-a[4]
stopifnot(!anyDuplicated(cfg$donor_id),!anyDuplicated(cand$gene),all(cand$bulk_direction%in%c(-1,1)))
paths<-file.path(qc_dir,'cache',paste0(cfg$donor_id,'.rds'))
res<-lapply(paths,readRDS);dr<-fread(file.path(qc_dir,'donor_region_retention.csv'))
regs<-c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
qc_names<-c('as_supplied','mt20','mt20_mad3','mt20_mad4')
stopifnot(all(vapply(res,function(z)identical(rownames(z$matrices$as_supplied),rownames(res[[1]]$matrices$as_supplied)),TRUE)))
dir.create(out,recursive=TRUE);wr<-function(x,n)fwrite(x,file.path(out,paste0(n,'.csv')))
forms<-list(primary=~group,sex_complete=~sex+group,age_sex_complete=~age+sex+group,depth=~depth_z+group,exclude_comorbidity=~group)
de<-list();designs<-list()
for(reg in regs){
 base<-sapply(res,function(z)z$matrices$as_supplied[,reg]);rownames(base)<-rownames(res[[1]]$matrices$as_supplied)
 md0<-copy(cfg);md0[,group:=factor(diagnosis,levels=c('Control','ALS'))]
 keepgenes<-filterByExpr(DGEList(base),design=model.matrix(~group,md0))
 for(qc in qc_names){xx<-sapply(res,function(z)z$matrices[[qc]][,reg]);rownames(xx)<-rownames(base)
 qname<-qc;rr<-dr[region==reg & qc==qname];md<-copy(md0);md[,n_spots:=rr$n_kept[match(donor_id,rr$donor_id)]]
 md[,depth_z:=as.numeric(scale(log10(colSums(xx))))]
 for(mod in names(forms)){ix<-which(md$n_spots>=20);if(mod=='sex_complete')ix<-intersect(ix,which(!is.na(md$sex)));if(mod=='age_sex_complete')ix<-intersect(ix,which(complete.cases(md[,.(age,sex)])));if(mod=='exclude_comorbidity')ix<-intersect(ix,which(md$sensitivity_eligible))
 mm<-md[ix];d<-model.matrix(forms[[mod]],mm);ok<-length(ix)>ncol(d) & qr(d)$rank==ncol(d) & all(table(mm$group)>=3)
 designs[[length(designs)+1]]<-data.table(qc=qc,region=reg,model=mod,n_Control=sum(mm$group=='Control'),n_ALS=sum(mm$group=='ALS'),testable=ok)
 if(!ok)next
 y<-calcNormFactors(DGEList(xx[keepgenes,ix,drop=FALSE]));y<-estimateDisp(y,d,robust=TRUE);f<-glmQLFit(y,d,robust=TRUE);zz<-as.data.table(topTags(glmQLFTest(f,coef='groupALS'),n=Inf,sort.by='none')$table,keep.rownames='gene');zz[,`:=`(qc=qc,region=reg,model=mod)];de[[length(de)+1]]<-zz
 }}}
a<-rbindlist(de);fwrite(a,file.path(out,'DE_all_genes.csv.gz'));wr(rbindlist(designs),'model_design')
grid<-CJ(gene=cand$gene,region=regs,model=names(forms),qc=qc_names)
cg<-merge(grid,a,by=c('gene','region','model','qc'),all.x=TRUE)
cg[,FDR_candidates:=p.adjust(ifelse(is.na(PValue),1,PValue),'BH'),by=.(qc,model)]
cg[is.na(PValue),FDR_candidates:=NA_real_]
cg<-merge(cg,cand[,.(gene,bulk_direction)],by='gene')
cg[,same_as_bulk:=sign(logFC)==bulk_direction]
st<-cg[,.(n_test=sum(is.finite(PValue)),same_direction=all(same_as_bulk%in%TRUE),max_FDR=max(FDR_candidates)),by=.(gene,region,qc)]
st[,supported:=n_test==length(forms)&same_direction&is.finite(max_FDR)&max_FDR<.05]
wr(cg,'candidate_all_models');wr(st,'candidate_stability')
wr(st[,.(supported_pairs=sum(supported),supported_genes=uniqueN(gene[supported==TRUE])),by=qc],'summary')
inputs<-normalizePath(c(paths,file.path(qc_dir,'donor_region_retention.csv'),commandArgs(trailingOnly=TRUE)[2:3],script))
wr(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),'input_manifest')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(fread(file.path(out,'summary.csv')))
