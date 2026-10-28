suppressPackageStartupMessages(library(ggtext))
loc<-rd('Figure5_region_locator')
loc[,group:=factor(group,levels=c('Control','ALS'))]
loc[,display_region:=ifelse(region%in%reg,rshort[region],'Other')]
loc[,`:=`(plot_x=x-min(x),plot_y=y-min(y)),by=image_id]
rp<-c(AH='#C28A53',GM='#8F789D',WM='#61989E',PWM='#8CA763',Other='#D8D8D8')
locator<-ggplot(loc,aes(plot_x,-plot_y,colour=display_region))+
 geom_point_rast(size=.28,raster.dpi=600)+facet_wrap(~group,nrow=1)+
 scale_colour_manual(values=rp,breaks=names(rp))+coord_equal()+theme_void(base_size=7.5,base_family='Arial')+
 labs(title='a  Original region annotations',colour=NULL)+
 theme(plot.title=element_text(size=8,face='bold'),strip.text=element_text(size=7.5,face='bold'),
 legend.position='right',legend.text=element_text(size=7),legend.key.height=unit(3,'mm'))+
 guides(colour=guide_legend(override.aes=list(size=2)))
a<-rd('Figure5_gene_matrix');b<-rd('Figure5_UCell')
stopifnot(nrow(a)==56,uniqueN(a$gene)==14,nrow(b)==56,sum(a$five_model_support,na.rm=TRUE)==24,
 all(a[five_model_support==TRUE,FDR_candidates464]>0),all(is.finite(a[five_model_support==TRUE,FDR_candidates464])),
 sum(b$FDR<.05,na.rm=TRUE)==0)
a[,region:=factor(region,levels=reg)];a[,gene:=factor(gene,levels=rev(sort(unique(gene))))]
gl<-setNames(paste0(ifelse(levels(a$gene)%in%c('C3','CTSB','HAMP','SLC11A1','TF'),paste0('<b>',levels(a$gene),'</b>'),levels(a$gene)),ifelse(levels(a$gene)=='SLC11A1','\u2020','')),levels(a$gene))
a[,stars:=fifelse(is.na(FDR_candidates464),'',fifelse(FDR_candidates464<.001,'***',fifelse(FDR_candidates464<.01,'**',fifelse(FDR_candidates464<.05,'*',''))))]
stopifnot(all(a[stars!='',is.finite(logFC)]))
plim<-max(abs(a$logFC),na.rm=TRUE)
p1<-ggplot(a,aes(region,gene))+
 geom_tile(aes(fill=logFC),colour='white',linewidth=.3)+
 geom_text(aes(label=stars),size=3.2,vjust=.7,colour=ifelse(is.finite(a$logFC)&abs(a$logFC)>1.7,'white','#202020'))+
 scale_fill_gradient2(low='#3876A6',mid='white',high='#B74E45',midpoint=0,limits=c(-plim,plim),na.value='#D4D4D4',
  breaks=c(-plim,0,plim),labels=function(x)sprintf('%.1f',x),name='log2 fold change',
  guide=guide_colourbar(order=1,direction='horizontal',barwidth=unit(25,'mm'),barheight=unit(3,'mm'),title.position='top',label.position='bottom'))+
 scale_x_discrete(labels=rshort,drop=FALSE)+scale_y_discrete(labels=gl,drop=FALSE)+
 labs(title='b  Regional gene expression\n    (primary model)',x=NULL,y=NULL)+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_markdown(size=7),
  legend.position='bottom',legend.justification='top',legend.box='vertical',legend.spacing.y=unit(1,'mm'),legend.box.spacing=unit(1,'mm'),axis.text.x=element_text(face='bold'))
b[,context:=factor(context,levels=reg)];b[,panel:=factor(panel,levels=rev(key$panel))]
set_colours<-c(GOBP='#008577',GOMF='#3376B8',HALLMARK='#AE7900',MODULE='#9151A5',REACTOME='#D06030')
labs_set<-setNames(sprintf("<span style='color:%s'><b>%s</b></span>",set_colours[sub('_.*','',key$panel)],key$set_id),key$panel)
set_key<-copy(b[1:5]);set_key[,collection:=factor(names(set_colours),levels=names(set_colours))]
p2<-ggplot(b,aes(context,panel))+geom_tile(aes(fill=effect),colour='white',linewidth=.3)+
 geom_point(data=set_key,aes(colour=collection),alpha=0,size=0,show.legend=TRUE)+
 scale_colour_manual(values=set_colours,labels=c('GO-BP','GO-MF','Hallmark','Module','Reactome'),name='Set',
  guide=guide_legend(order=2,nrow=2,byrow=TRUE,title.position='top',override.aes=list(alpha=1,size=3,shape=15)))+
 scale_x_discrete(labels=rshort,drop=FALSE)+scale_y_discrete(labels=labs_set,drop=FALSE)+divfill('UCell effect')+
 labs(title='c  Regional set scores',x=NULL,y=NULL)+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_markdown(size=7),
  axis.text.x=element_text(face='bold'),legend.position='bottom',legend.box='vertical',
  legend.spacing.y=unit(1,'mm'),legend.box.spacing=unit(1,'mm'))+
 guides(fill=guide_colourbar(order=1,direction='horizontal',barwidth=unit(26,'mm'),barheight=unit(2.5,'mm'),title.position='top',label.position='bottom'))
p<-locator/(p1|p2)+plot_layout(heights=c(.8,1.5))+
 plot_annotation(caption='b  Stars: primary-model FDR (* < 0.05; ** < 0.01; *** < 0.001)\n24 five-model supported pairs: Figure S3 | Grey: not tested | Bold: cross-layer candidates\nc  0/56 comparisons with BH FDR < 0.05 | Set definitions: Figure S4',
  theme=theme(plot.caption=element_text(size=6.5,hjust=0,lineheight=1.2)))

savefig(p,'Figure5_spatial_support',185)
