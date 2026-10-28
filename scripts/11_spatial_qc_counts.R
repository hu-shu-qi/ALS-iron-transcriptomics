a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4)stop('Usage: Rscript 11_spatial_qc_counts.R INPUT_DIR CONFIG_CSV REGION_RDS NEW_OUTPUT_DIR')
if(file.exists(a[4]))stop('Output exists')
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table)})
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
cfg<-fread(a[2],na.strings=c('','NA'))
stopifnot(!anyDuplicated(cfg$file),!anyDuplicated(cfg$donor_id))
input_dir<-normalizePath(a[1]);baseline_counts<-readRDS(a[3]);out<-a[4]
stopifnot(all(file.exists(file.path(input_dir,cfg$file))),setequal(names(baseline_counts),cfg$file))
regs<-c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
qc_names<-c('as_supplied','mt20','mt20_mad3','mt20_mad4')
dir.create(file.path(out,'cache'),recursive=TRUE)
wr<-function(x,n)fwrite(x,file.path(out,paste0(n,'.csv')))
for(i in seq_len(nrow(cfg))){
 dest<-file.path(out,'cache',paste0(cfg$donor_id[i],'.rds'));input<-file.path(input_dir,cfg$file[i])
 stamp<-file.info(input)[,c('size','mtime')]
 if(file.exists(dest)){cached<-readRDS(dest);stopifnot(identical(cached$stamp,stamp));next}
 message('QC extract ',i,'/18 ',cfg$donor_id[i]);cc<-gzcon(gzfile(input,'rb'));o<-readRDS(cc);close(cc)
 md<-as.data.table(o@meta.data,keep.rownames='spot_id')
 layers<-grep('^counts($|\\.)',Layers(o[['Spatial']]),value=TRUE)
 x<-do.call(cbind,lapply(layers,function(l)LayerData(o,assay='Spatial',layer=l)))
 stopifnot(!anyDuplicated(colnames(x)),setequal(colnames(x),md$spot_id));x<-x[,md$spot_id]
 stopifnot(all(Matrix::colSums(x)==md$nCount_Spatial),all(Matrix::colSums(x>0)==md$nFeature_Spatial),all(is.finite(md$percent.mt)),all(md$percent.mt>=0 & md$percent.mt<=100))
 co<-rbindlist(lapply(names(o@images),function(im){z<-as.data.table(GetTissueCoordinates(o[[im]]),keep.rownames='spot_id');z[,.(spot_id,image_id=im)]}))
 stopifnot(!anyDuplicated(co$spot_id),setequal(co$spot_id,md$spot_id));md[,image_id:=co$image_id[match(spot_id,co$spot_id)]]
 md[,`:=`(donor_id=cfg$donor_id[i],group=cfg$diagnosis[i],region=as.character(manual_3))]
 masks<-list(as_supplied=rep(TRUE,nrow(md)),mt20=md$percent.mt<=20);th<-list()
 for(k in c(3,4)){ok<-md$percent.mt<=20;for(field in c('nCount_Spatial','nFeature_Spatial')){
 v<-log1p(md[[field]]);med<-median(v);dev<-mad(v);lo<-if(dev>0)med-k*dev else -Inf;hi<-if(dev>0)med+k*dev else Inf
 ok<-ok & v>=lo & v<=hi;th[[length(th)+1]]<-data.table(donor_id=cfg$donor_id[i],k=k,field=field,log1p_lower=lo,log1p_upper=hi,mad=dev)};masks[[paste0('mt20_mad',k)]]<-ok}
 mats<-list();ret<-list()
 for(q in qc_names){keep<-masks[[q]];mats[[q]]<-sapply(regs,function(r)Matrix::rowSums(x[,keep & !is.na(md$region) & md$region==r,drop=FALSE]));rownames(mats[[q]])<-rownames(x)
 md[,kept:=keep];ret[[q]]<-md[,.(n_before=.N,n_kept=sum(kept),retention=mean(kept)),by=.(donor_id,group,image_id,region)];ret[[q]][,qc:=q]}
 stopifnot(isTRUE(all.equal(mats$as_supplied,baseline_counts[[cfg$file[i]]][rownames(x),regs],check.attributes=FALSE,tolerance=0)))
 mapping<-md[,.(n_spots=.N,sample_values=paste(unique(sample),collapse=';'),condition_values=paste(unique(condition),collapse=';'),mt_min=min(percent.mt),mt_max=max(percent.mt),n_mt_gt20=sum(percent.mt>20)),by=.(donor_id,image_id)]
 z<-list(stamp=stamp,matrices=mats,retention=rbindlist(ret),thresholds=rbindlist(th),mapping=mapping)
 saveRDS(z,dest);rm(o,x,md,z);gc()
}
res<-lapply(cfg$donor_id,function(d)readRDS(file.path(out,'cache',paste0(d,'.rds'))));ret<-rbindlist(lapply(res,`[[`,'retention'))
wr(ret,'image_region_retention');wr(rbindlist(lapply(res,`[[`,'thresholds')),'qc_thresholds');wr(rbindlist(lapply(res,`[[`,'mapping')),'image_mapping')
dr<-ret[region%in%regs,.(n_before=sum(n_before),n_kept=sum(n_kept)),by=.(donor_id,group,region,qc)];dr[,retention:=n_kept/n_before];wr(dr,'donor_region_retention')
inputs<-normalizePath(c(file.path(input_dir,cfg$file),a[2:3],script))
wr(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),'input_manifest')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('Spatial QC counts prepared.')
