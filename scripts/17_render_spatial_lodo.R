a <- commandArgs(trailingOnly=TRUE)
if (length(a)!=2L) stop('Usage: Rscript 17_render_spatial_lodo.R LODO_DIRECTORY NEW_OUTPUT_DIR')
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
stopifnot(!file.exists(a[2]))
input <- file.path(a[1],'all_candidate_results.csv')
mods <- c('primary','sex_complete','age_sex_complete','depth','exclude_comorbidity')
regs <- c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
if(file.exists(input)) {
  x <- fread(input)
stopifnot(all(c('gene','region','model','job','omitted_donor','logFC','PValue','FDR464','bulk_direction') %in% names(x)), !anyDuplicated(x[,.(gene,region,job)]))
stopifnot(setequal(unique(x$model),mods),setequal(unique(x$region),regs))
stopifnot(all(x[,.(n=.N),by=job]$n==uniqueN(x$gene)*length(regs)))
x[,checked_FDR:=p.adjust(ifelse(is.finite(PValue),PValue,1),'BH'),by=job]
x[!is.finite(PValue),checked_FDR:=NA_real_]
stopifnot(isTRUE(all.equal(x$checked_FDR,x$FDR464,tolerance=1e-10)))
baseline <- x[omitted_donor=='none']
stopifnot(!anyDuplicated(baseline[,.(gene,region,model)]))
supported <- baseline[,.(supported=.N==length(mods) && all(is.finite(PValue)) && all(is.finite(logFC)) && all(sign(logFC)==bulk_direction) && all(FDR464<0.05)),by=.(gene,region)][supported==TRUE]
lo <- merge(x[omitted_donor!='none'],baseline[,.(gene,region,model,baseline_logFC=logFC,baseline_FDR=FDR464)],by=c('gene','region','model'))
su <- lo[,{
  k <- is.finite(logFC)&is.finite(FDR464)
  list(n_planned=.N,n_tested=sum(k),baseline_logFC=baseline_logFC[1],baseline_FDR=baseline_FDR[1],min_logFC=if(any(k))min(logFC[k])else NA_real_,max_logFC=if(any(k))max(logFC[k])else NA_real_,direction_retained_n=sum(sign(logFC[k])==sign(baseline_logFC[k])),FDR_retained_n=sum(FDR464[k]<0.05),max_FDR=if(any(k))max(FDR464[k])else NA_real_)
},by=.(gene,region,model)]
focus <- merge(su,supported[,.(gene,region)],by=c('gene','region'))
stopifnot(nrow(focus)>0L,all(focus$n_planned>0L),all(focus$n_tested==focus$n_planned))

} else {
  input <- file.path(a[1],'FigureS3_source_data.csv')
  focus <- fread(input)
  stopifnot(!anyDuplicated(focus[,.(gene,region,model)]),setequal(unique(focus$model),mods),all(focus$n_tested==focus$n_planned))
}
rshort <- setNames(c('AH','GM','WM','PWM'),regs)
pdata <- copy(focus); pdata[,pair:=paste(gene,rshort[region],sep=' | ')]
ord <- sort(unique(pdata$pair)); pdata[,pair:=factor(pair,levels=rev(ord))]
pdata[,model:=factor(model,levels=mods)]
pdata[,`:=`(retention=FDR_retained_n/n_planned,label=paste0(FDR_retained_n,'/',n_planned),ypos=as.numeric(pair)+ifelse(model=='primary',0.12,-0.12))]
theme_set(theme_classic(base_size=7.5,base_family='Arial')+theme(axis.text=element_text(size=7,colour='#272727'),axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),legend.text=element_text(size=7),legend.title=element_text(size=7),plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=7),strip.background=element_blank(),strip.text=element_text(size=7.5,face='bold'),plot.margin=margin(3,4,3,4),legend.key.height=unit(3,'mm'),legend.key.width=unit(4,'mm')))
p1 <- ggplot(pdata[model %in% c('primary','age_sex_complete')],aes(baseline_logFC,ypos,colour=model))+geom_vline(xintercept=0,linetype=2,linewidth=0.3,colour='#BBBBBB')+geom_segment(aes(x=min_logFC,xend=max_logFC,yend=ypos),linewidth=0.45)+geom_point(size=1.6)+scale_y_continuous(breaks=seq_along(ord),labels=rev(ord),limits=c(0.4,length(ord)+0.6),expand=c(0,0))+scale_colour_manual(values=c(primary='#286E86',age_sex_complete='#B57845'),labels=c('Group','Age + sex'))+labs(title='a  Effect direction and LODO range',x=expression('ALS - Control log'[2]*' fold change'),y=NULL,colour=NULL)+theme(legend.position='bottom')
ml <- c(primary='Group',sex_complete='Sex',age_sex_complete='Age + sex',depth='Depth',exclude_comorbidity='No comorbidity')
p2 <- ggplot(pdata,aes(model,pair))+geom_tile(aes(fill=retention),colour='white',linewidth=0.2)+geom_text(aes(label=label),size=2.3,family='Arial')+scale_fill_gradient(low='#F6EADB',high='#78B2B5',limits=c(0,1),guide='none')+scale_x_discrete(labels=ml)+labs(title='b  FDR retained after donor deletion',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(angle=40,hjust=1,size=7))
p <- (p1|p2)+plot_annotation(caption='a  Range: minimum–maximum LODO effects\nb  Runs with FDR < 0.05 / planned runs',theme=theme(plot.caption=element_text(size=7,hjust=0)))
dir.create(a[2],recursive=TRUE)
for(ext in c('pdf','svg','png','tiff')) {
  fn <- file.path(a[2],paste0('FigureS3_LODO.',ext))
  if(ext=='pdf') grDevices::cairo_pdf(fn,width=170/25.4,height=170/25.4,family='Arial')
  if(ext=='svg') svglite::svglite(fn,width=170/25.4,height=170/25.4)
  if(ext=='png') ragg::agg_png(fn,width=170,height=170,units='mm',res=300)
  if(ext=='tiff') ragg::agg_tiff(fn,width=170,height=170,units='mm',res=600,compression='lzw')
  print(p); dev.off()
}
fwrite(focus,file.path(a[2],'FigureS3_source_data.csv'))
fwrite(pdata,file.path(a[2],'FigureS3_plot_data.csv'))
script <- gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
inputs <- c(input,script)
fwrite(data.table(file=basename(inputs),md5=unname(tools::md5sum(inputs))),file.path(a[2],'input_manifest.csv'))
writeLines(capture.output(sessionInfo()),file.path(a[2],'sessionInfo.txt'))
message('Rendered Figure S3 from ',uniqueN(focus[,.(gene,region)]),' supported pairs. Effect ranges are not confidence intervals.')
