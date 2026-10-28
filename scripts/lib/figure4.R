options(stringsAsFactors=FALSE)
try(Sys.setlocale('LC_CTYPE','en_US.UTF-8'),silent=TRUE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggtext)})
a<-fread(file.path(src,'Figure4_gene_matrix.csv'))
b<-fread(file.path(src,'Figure4_UCell.csv'))
key<-fread(file.path(src,'Gene_set_label_key.csv'))
ct<-c('Oligodendrocytes','Astrocytes','Microglia','OPCs','VLMC','Lymphocytes','PVM','Endothelial','Neurons','Pericytes','Ependymal')
ctl<-setNames(c('Oligo','Astro','Micro','OPC','VLMC','Lymph','PVM','Endo','Neuron','Pericytes','Ependymal'),ct)
stopifnot(nrow(a)==275,nrow(b)==154,uniqueN(a$gene)==25,all(is.finite(a[support=='sex_model_support',FDR])),all(a[support=='sex_model_support',FDR]>0))
s<-copy(a[support=='sex_model_support'])
stopifnot(nrow(s)==34,all(s$FDR<.05),all(abs(s$logFC)>=.25))
s[,priority:=match(cell_type,ct)];setorder(s,gene,FDR,priority)
ordering<-s[,.(primary_type=cell_type[1],primary_FDR=FDR[1]),by=gene]
ordering[,priority:=match(primary_type,c('Astrocytes','Microglia','Oligodendrocytes','OPCs','Endothelial'))]
setorder(ordering,priority,primary_FDR,gene)
a[,gene:=factor(gene,levels=rev(ordering$gene))];a[,cell_type:=factor(cell_type,levels=ct)]
b[,context:=factor(context,levels=ct)];b[,panel:=factor(panel,levels=rev(key$panel))]
sig<-b[is.finite(FDR)&FDR<.05,.N,by=context]
counts<-setNames(rep(0L,length(ct)),ct);counts[as.character(sig$context)]<-sig$N
stopifnot(sum(counts)==9,counts['Microglia']==6,counts['Astrocytes']==2,counts['PVM']==1)
blabels<-setNames(paste0(ctl,' (',counts,')'),ct)
set_colours<-c(GOBP='#008577',GOMF='#3376B8',HALLMARK='#AE7900',MODULE='#9151A5',REACTOME='#D06030')
setlabels<-setNames(sprintf("<span style='color:%s'><b>%s</b></span> %s",set_colours[sub('_.*','',key$panel)],key$set_id,key$short_label),key$panel)
theme_set(theme_classic(base_size=7.5,base_family='Arial')+theme(
 axis.text=element_text(size=7,colour='#272727'),axis.line=element_blank(),axis.ticks=element_blank(),
 axis.text.x=element_text(angle=45,hjust=1,face='bold'),plot.title=element_text(size=8,face='bold'),
 plot.subtitle=element_text(size=7),legend.text=element_text(size=7),legend.title=element_text(size=7),
 plot.margin=margin(4,4,4,4),legend.position='right',legend.box='vertical',legend.spacing.y=unit(5,'mm')))
fc_limit<-max(abs(a$logFC),na.rm=TRUE)
pa<-ggplot(a,aes(cell_type,gene))+
 geom_point(data=a[support!='not_tested'&support!='sex_model_support'],colour='#D4D4D4',size=.65)+
 geom_point(data=a[support=='sex_model_support'],aes(fill=logFC,size=-log10(FDR)),shape=21,colour='transparent',stroke=0)+
 scale_x_discrete(labels=ctl,drop=FALSE)+
 scale_y_discrete(labels=function(x)paste0(ifelse(x%in%c('C3','CTSB','HAMP','SLC11A1','TF'),paste0('<b>',x,'</b>'),x),ifelse(x=='SLC11A1','\u2020','')),drop=FALSE)+
 scale_fill_gradient2(low='#17639B',mid='white',high='#AA3830',midpoint=0,limits=c(-fc_limit,fc_limit),transform=scales::pseudo_log_trans(sigma=.35),breaks=c(-5,-1,0,1,5),
   name='log2 fold change',guide=guide_colourbar(order=1,barwidth=unit(24,'mm'),barheight=unit(3,'mm'),direction='horizontal',title.position='top',label.position='bottom'))+
 scale_size_area(max_size=3.6,breaks=c(1.5,2,3),name=expression(-log[10]('Within-type FDR')),
   guide=guide_legend(order=2,override.aes=list(fill='#8D9AA8',colour='transparent')))+
 labs(title='a  Cell-type gene support',subtitle='25 genes | 34 supported gene–cell-type pairs',x=NULL,y=NULL,
   caption='Colour: supported pair   |   Grey dot: other tested pair   |   Blank: not tested')+
 theme(axis.text.y=element_markdown(size=7),panel.grid.major=element_line(colour='#F0F0F0',linewidth=.2),plot.caption=element_text(size=6.5,hjust=0))
score_limit<-max(abs(b$effect),na.rm=TRUE)
pb<-ggplot(b,aes(context,panel))+
 geom_tile(aes(fill=effect),colour='white',linewidth=.3)+
 geom_point(data=b[is.finite(FDR)&FDR<.05],shape=8,size=1.6)+
 scale_x_discrete(labels=blabels,drop=FALSE)+scale_y_discrete(labels=setlabels,drop=FALSE)+
 scale_fill_gradient2(low='#3876A6',mid='white',high='#B74E45',midpoint=0,limits=c(-score_limit,score_limit),
  breaks=c(-score_limit,0,score_limit),labels=function(x)sprintf('%.3f',x),name='UCell effect',
  guide=guide_colourbar(barwidth=unit(24,'mm'),barheight=unit(3,'mm'),direction='horizontal',title.position='top',label.position='bottom'))+
 labs(title='b  Iron-related UCell scores',subtitle='14 sets | * BH FDR < 0.05 across 154 tests',x=NULL,y=NULL,
  caption='Numbers in parentheses: significant sets per cell type')+
 theme(axis.text.y=element_markdown(size=7),plot.caption=element_text(size=6.5,hjust=0))
pa<-pa+labs(subtitle=NULL,caption=NULL)
pb<-pb+labs(subtitle=NULL,caption=NULL)
p<-(pa/pb+plot_layout(heights=c(1.35,1)))+
 plot_annotation(caption=paste0('a  Colour: supported pair | Grey dot: other tested pair | Blank: not tested\n',
 'Bold: cross-layer candidates | Colour scale: pseudo-log spacing; values are log2 fold changes\n',
 'b  * BH FDR < 0.05 | Parentheses: significant sets per cell type'),
 theme=theme(plot.caption=element_text(size=6.5,hjust=0,lineheight=1.2)))

savefig(p,'Figure4_cell_type_support',205)
