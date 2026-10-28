a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=5)stop('Usage: Rscript 8_gene_set_models.R BULK_RDS REGION_RDS META_CSV SN_DONOR_SCORES_CSV NEW_OUTPUT_DIR')
suppressPackageStartupMessages({library(data.table);library(Matrix);library(edgeR);library(UCell);library(AUCell);library(singscore);library(BiocParallel)})
if(file.exists(a[5]))stop('Output already exists')
stopifnot(all(file.exists(a[1:4])))
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
pkg<-dirname(dirname(normalizePath(script)))
source(file.path(pkg,'scripts/lib/scoring.R'))
parts<-strsplit(readLines(file.path(pkg,'resources/iron_gene_sets.gmt')),'\t');sets<-setNames(lapply(parts,function(x)unique(x[-c(1,2)])),vapply(parts,`[`,'',1));sets<-sets[sort(names(sets))]
regs<-c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
forms<-list(sex=~sex+group)
spforms<-list(primary=~group,sex_complete=~sex+group,age_sex_complete=~age+sex+group,depth=~depth_z+group,exclude_comorbidity=~group)
out<-a[5];dir.create(out,recursive=TRUE);wr<-function(x,n)fwrite(x,file.path(out,paste0(n,'.csv')))
allsc<-list();tests<-list();coverage<-list()
b<-readRDS(a[1]);collect(b$expression,as.data.table(b$metadata),'bulk','Bulk',forms,20260910)
lst<-readRDS(a[2]);meta<-fread(a[3],na.strings=c('','NA'));stopifnot(!anyDuplicated(meta$file),!anyNA(match(names(lst),meta$file)))
common<-Reduce(intersect,lapply(lst,rownames))
for(reg in regs){md<-copy(meta[match(names(lst),file)]);md[,group:=factor(diagnosis,levels=c('Control','ALS'))];x<-sapply(lst,function(z)z[common,reg]);rownames(x)<-common;colnames(x)<-paste(md$donor_id,reg,sep='___');md[,depth_z:=as.numeric(scale(log10(colSums(x))))];collect(cpm(DGEList(x),log=TRUE,prior.count=2),md,reg,'Spatial',spforms,20260911+match(reg,regs));message('Scored ',reg)}
sn<-fread(a[4]);cts<-unique(sn[eligible_profile==TRUE,.(cell_type,donor_id,group)])[,.(ok=sum(group=='Control')>=3&sum(group=='ALS')>=3),by=cell_type][ok==TRUE,as.character(cell_type)]
for(ct in cts)for(meth in c('UCell','AUCell','singscore'))for(pan in names(sets)){
 z<-sn[as.character(cell_type)==ct&method==meth&panel==pan&eligible_profile==TRUE];t<-fit_score(z,~sex+group);t[,`:=`(layer='snRNA',context=ct,method=meth,panel=pan,model='sex')];tests[[length(tests)+1]]<-t
}
t<-rbindlist(tests);t[,FDR:=adj(P),by=.(layer,method,model)];t[,FDR_HC3:=adj(P_HC3),by=.(layer,method,model)]
wr(t,'score_tests');wr(rbindlist(allsc,fill=TRUE),'bulk_spatial_donor_scores');wr(rbindlist(coverage),'coverage')
wr(data.table(path=normalizePath(a[1:4]),md5=unname(tools::md5sum(a[1:4]))),'input_manifest')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
