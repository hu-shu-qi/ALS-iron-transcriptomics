a<-commandArgs(trailingOnly=TRUE)
if(!length(a)%in%c(3,4))stop('Usage: Rscript 4_aggregate_single_nucleus.R ANNOTATED_RDS DONOR_CSV NEW_OUTPUT_DIR [BPCELLS_ROOT]')
stopifnot(all(file.exists(a[1:2])))
if(file.exists(a[3]))stop('Output exists')
suppressPackageStartupMessages({library(Seurat);library(BPCells);library(Matrix);library(data.table)})
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
obj<-readRDS(a[1]);clinical<-fread(a[2])
stopifnot(!anyDuplicated(clinical$snRNA_GSM),!anyDuplicated(clinical$donor_id))
md<-as.data.table(obj@meta.data,keep.rownames='cell_id')
md[,snRNA_GSM:=sub('_.*','',orig.ident)]
ix<-match(md$snRNA_GSM,clinical$snRNA_GSM)
stopifnot(!anyNA(ix),all(as.character(md$group)==clinical$group[ix]),!anyNA(md$cell_type))
md[,donor_id:=as.character(clinical$donor_id[ix])]
md[,profile_id:=paste(donor_id,cell_type,sep='___')]
counts<-LayerData(obj,assay='RNA',layer='counts')
mapping<-NULL
helper<-file.path(dirname(normalizePath(script)),'lib/bpcells.R')
if(length(a)==4){source(helper);r<-relocate_bpcells(counts,a[4]);counts<-r$counts;mapping<-r$mapping}
stopifnot(identical(colnames(counts),md$cell_id),!anyDuplicated(rownames(counts)))
keys<-factor(md$profile_id,levels=unique(md$profile_id))
pb<-as.matrix(BPCells::pseudobulk_matrix(counts,keys,method='sum'))
stopifnot(all(is.finite(pb)),all(pb>=0),all(pb==round(pb)),setequal(colnames(pb),md$profile_id))
pm<-md[,.(n_cells=.N),by=.(profile_id,donor_id,snRNA_GSM,group,cell_type)]
pm<-merge(pm,clinical[,.(snRNA_GSM,age,sex)],by='snRNA_GSM',all.x=TRUE,sort=FALSE)
pm[,group:=factor(group,levels=c('Control','ALS'))]
pm[,eligible_profile:=n_cells>=20]
cov<-pm[eligible_profile==TRUE,.(n_Control=sum(group=='Control'),n_ALS=sum(group=='ALS')),by=cell_type]
cov<-merge(unique(pm[,.(cell_type)]),cov,by='cell_type',all.x=TRUE)
cov[,eligible_type:=!is.na(n_Control)&n_Control>=3&n_ALS>=3]
dir.create(a[3],recursive=TRUE)
saveRDS(list(counts=pb,metadata=pm),file.path(a[3],'sn_pseudobulk.rds'))
fwrite(pm,file.path(a[3],'sn_profile_metadata.csv'));fwrite(cov,file.path(a[3],'sn_celltype_eligibility.csv'))
if(!is.null(mapping))fwrite(mapping,file.path(a[3],'BPCells_path_mapping.csv'))
inputs<-normalizePath(c(a[1:2],script,helper))
fwrite(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(a[3],'input_manifest.csv'))
writeLines(capture.output(sessionInfo()),file.path(a[3],'sessionInfo.txt'))
cat('COMPLETE:',nrow(pb),'genes,',ncol(pb),'profiles,',nrow(md),'nuclei\n')
