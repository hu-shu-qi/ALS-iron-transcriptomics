a<-commandArgs(trailingOnly=TRUE)
if(!length(a)%in%c(3,4))stop('Usage: Rscript 7_single_nucleus_scores.R ANNOTATED_RDS PSEUDOBULK_RDS NEW_OUTPUT_DIR [BPCELLS_ROOT]')
stopifnot(all(file.exists(a[1:2])))
if(file.exists(a[3]))stop('Output already exists')
suppressPackageStartupMessages({library(data.table);library(Matrix);library(edgeR);library(UCell);library(AUCell);library(singscore);library(BiocParallel);library(Seurat);library(BPCells)})
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
pkg<-dirname(dirname(normalizePath(script)))
source(file.path(pkg,'scripts/lib/scoring.R'))
gmt<-file.path(pkg,'resources/iron_gene_sets.gmt');parts<-strsplit(readLines(gmt),'\t');sets<-setNames(lapply(parts,function(x)unique(x[-c(1,2)])),vapply(parts,`[`,'',1));sets<-sets[sort(names(sets))]
forms<-list(sex=~sex+group)
root<-a[3];dir.create(root,recursive=TRUE)
for(z in c('results','chunks','logs'))dir.create(file.path(root,z))
file.copy(gmt,file.path(root,'iron_gene_sets.gmt'))
wr<-function(x,n)fwrite(x,file.path(root,'results',paste0(n,'.csv')))
annotated_input<-normalizePath(a[1]);pseudobulk_input<-normalizePath(a[2])
backing_root<-if(length(a)==4)normalizePath(a[4],mustWork=TRUE) else NULL
helper<-file.path(pkg,'scripts/lib/bpcells.R')
source(helper)
read_counts<-function(obj){
 x<-LayerData(obj,assay='RNA',layer='counts')
 if(!is.null(backing_root)){
  relocated<-relocate_bpcells(x,backing_root)
  wr(relocated$mapping,'BPCells_path_mapping');x<-relocated$counts
 }
 x
}
obj<-readRDS(annotated_input)
md<-as.data.table(obj@meta.data,keep.rownames='cell_id');md[,snRNA_GSM:=sub('_.*','',orig.ident)]
pb<-readRDS(pseudobulk_input);pm<-as.data.table(pb$metadata)
don<-unique(pm[,.(snRNA_GSM,donor_id)]);stopifnot(!anyDuplicated(don$snRNA_GSM))
md[,donor_id:=don$donor_id[match(snRNA_GSM,don$snRNA_GSM)]];md[,profile_id:=paste(donor_id,cell_type,sep='___')]
counts<-read_counts(obj);stopifnot(identical(colnames(counts),md$cell_id),nrow(counts)==32869,ncol(counts)==100265,!anyNA(md$donor_id))
wr(data.table(panel=names(sets),requested=lengths(sets),detected=lengths(lapply(sets,intersect,y=rownames(counts)))),'sn_coverage')
blocks<-split(seq_len(ncol(counts)),ceiling(seq_len(ncol(counts))/500))
for(k in seq_along(blocks)){
 dest<-file.path(root,'chunks',sprintf('chunk_%03d.rds',k))
 if(file.exists(dest))next
 ix<-blocks[[k]];x<-as.matrix(counts[,ix,drop=FALSE]);sc<-score3(x,20260910+k)
 saveRDS(list(cell_id=colnames(x),scores=sc,set_md5=unname(tools::md5sum(file.path(root,'iron_gene_sets.gmt')))),dest)
 rm(x,sc);gc(verbose=FALSE);message('SN chunk ',k,'/',length(blocks),' ',Sys.time())
}
chunks<-lapply(seq_along(blocks),function(k)readRDS(file.path(root,'chunks',sprintf('chunk_%03d.rds',k))))
stopifnot(identical(unlist(lapply(chunks,`[[`,'cell_id'),use.names=FALSE),md$cell_id),all(vapply(chunks,function(z)identical(z$set_md5,unname(tools::md5sum(file.path(root,'iron_gene_sets.gmt')))),TRUE)))
ags<-list()
for(meth in c('UCell','AUCell','singscore')){
 mat<-do.call(rbind,lapply(chunks,function(z)z$scores[[meth]]))
 w<-cbind(md[,.(cell_id,profile_id)],as.data.table(mat));long<-melt(w,id.vars=c('cell_id','profile_id'),variable.name='panel',value.name='score')
 a<-long[,.(score=mean(score),median_score=median(score),n_cells=.N),by=.(profile_id,panel)];a[,method:=meth];ags[[meth]]<-a
}
a<-merge(rbindlist(ags),pm[,!c('n_cells'),with=FALSE],by='profile_id',all.x=TRUE);wr(a,'sn_donor_scores')
cts<-pm[eligible_profile==TRUE,.(ok=sum(group=='Control')>=3 & sum(group=='ALS')>=3),by=cell_type][ok==TRUE,as.character(cell_type)]
t<-list()
for(ct in cts)for(meth in c('UCell','AUCell','singscore'))for(pan in names(sets))for(mod in names(forms)){
 z<-a[as.character(cell_type)==ct & method==meth & panel==pan & eligible_profile==TRUE]
 v<-fit_score(z,forms[[mod]]);v[,`:=`(layer='snRNA',context=ct,method=meth,panel=pan,model=mod)];t[[length(t)+1]]<-v
}
t<-rbindlist(t);t[,FDR:=adj(P),by=.(method,model)];t[,FDR_HC3:=adj(P_HC3),by=.(method,model)]
wr(t,'sn_score_tests')
writeLines(capture.output(sessionInfo()),file.path(root,'logs/sn_sessionInfo.txt'));message('COMPLETE snRNA three-method scoring')
inputs<-c(annotated_input,pseudobulk_input,gmt,script,helper)
wr(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),'input_manifest')
