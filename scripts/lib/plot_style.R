options(stringsAsFactors=FALSE)
try(Sys.setlocale('LC_CTYPE','en_US.UTF-8'),silent=TRUE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrastr)})
rd<-function(n)fread(file.path(src,paste0(n,'.csv')),encoding='UTF-8')
theme_set(theme_classic(base_size=7.5,base_family='Arial')+theme(axis.text=element_text(size=7,colour='#272727'),axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),legend.text=element_text(size=7),legend.title=element_text(size=7),plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=7),strip.background=element_blank(),strip.text=element_text(size=7.5,face='bold'),plot.margin=margin(3,4,3,4),legend.key.height=unit(3,'mm'),legend.key.width=unit(4,'mm')))
savefig<-function(p,n,h,w=170) {
 for(ext in c('pdf','svg','png','tiff')) {
  fn<-file.path(plots,paste0(n,'.',ext))
  if(ext=='pdf')grDevices::cairo_pdf(fn,width=w/25.4,height=h/25.4,family='Arial')
  if(ext=='svg')svglite::svglite(fn,width=w/25.4,height=h/25.4)
  if(ext=='png')ragg::agg_png(fn,width=w,height=h,units='mm',res=300)
  if(ext=='tiff')ragg::agg_tiff(fn,width=w,height=h,units='mm',res=600,compression='lzw')
  print(p);dev.off()
 }
 message('EXPORTED ',n)
}
divfill<-function(title,lim=NULL)scale_fill_gradient2(low='#3876A6',mid='white',high='#B74E45',midpoint=0,limits=if(is.null(lim))function(x)c(-max(abs(x)),max(abs(x)))else lim,breaks=function(x)c(x[1],0,x[2]),labels=function(x)format(signif(x,2),trim=TRUE,scientific=FALSE),na.value='#D4D4D4',name=title,guide=guide_colorbar(direction='horizontal',barwidth=unit(26,'mm'),barheight=unit(2.5,'mm'),title.position='top',label.position='bottom'))
ct<-c('Oligodendrocytes','Astrocytes','Microglia','OPCs','VLMC','Lymphocytes','PVM','Endothelial','Neurons','Pericytes','Ependymal','Schwann')
ctshort<-setNames(c('Oligo','Astro','Micro','OPC','VLMC','Lymph','PVM','Endo','Neuron','Pericytes','Ependymal','Schwann'),ct)
pal<-setNames(c('#4178A4','#D2944A','#388E83','#9174A8','#AA715C','#C8708B','#8D9650','#68AFC0','#D56053','#636DA2','#6E8B64','#A8A8A8'),ct)
key<-rd('Gene_set_label_key');key[,label:=paste(set_id,short_label)];setlabs<-setNames(key$label,key$panel)
reg<-c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter');rshort<-setNames(c('AH','GM','WM','PWM'),reg)
