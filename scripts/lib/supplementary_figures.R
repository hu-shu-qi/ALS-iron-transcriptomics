library(ggtext)
colset<-c(GOBP='#008577',GOMF='#3376B8',HALLMARK='#AE7900',MODULE='#9151A5',REACTOME='#D06030')
setlabs<-setNames(sprintf("<span style='color:%s'><b>%s</b></span> %s",colset[sub('_.*','',key$panel)],key$set_id,key$short_label),key$panel)

gmt<-strsplit(readLines(file.path(pkg,'resources/iron_gene_sets.gmt')),'\t')
members<-rbindlist(lapply(gmt,function(x)data.table(panel=x[1],gene=unique(x[-c(1,2)]))))
stopifnot(nrow(members)==792,uniqueN(members$gene)==545)
definitions<-merge(key,members[,.(genes=.N),by=panel],by='panel',sort=FALSE)
setorder(definitions,set_id)
definitions[,collection:=sub('_.*','',panel)]
definitions[,y:=15-seq_len(.N)]
definitions[,full:=vapply(gsub('_','_ ',panel),function(x)paste(strwrap(x,width=52),collapse='\n'),character(1))]
definitions[,full:=gsub('_ ','_',full,fixed=TRUE)]
definitions[,label:=paste0('<b>',set_id,'</b>  ',short_label)]
kp<-ggplot(definitions)+
 geom_rect(aes(xmin=0,xmax=1,ymin=y-.47,ymax=y+.47,fill=collection),alpha=.11,colour=NA)+
 geom_richtext(aes(x=.015,y=y,label=label,colour=collection),hjust=0,size=2.45,label.colour=NA,fill=NA)+
 geom_text(aes(x=.46,y=y,label=full),hjust=0,size=2.3,lineheight=1.05,family='Arial')+
 geom_text(aes(x=.97,y=y,label=genes),hjust=1,size=2.5,family='Arial')+
 scale_fill_manual(values=colset,guide='none')+scale_colour_manual(values=colset,guide='none')+
 annotate('text',x=c(.015,.46,.97),y=14.8,label=c('Set ID and short name','Human MSigDB v2026.1.Hs name','Genes'),hjust=c(0,0,1),size=2.65,fontface='bold')+
 coord_cartesian(xlim=c(0,1),ylim=c(.35,15.1),expand=FALSE,clip='off')+
 labs(title='Iron-related gene-set definitions',caption='792 memberships across 545 unique genes | Memberships overlap across sets')+
 theme_void(base_family='Arial')+theme(plot.title=element_text(size=9,face='bold'),plot.caption=element_text(size=7,hjust=0),plot.margin=margin(6,6,6,6))
savefig(kp,'FigureS4_gene_set_key',210)

a<-rd('Gene_set_score_results')[layer=='Bulk']
stopifnot(nrow(a)==42)
a[,panel:=factor(panel,levels=rev(key$panel))]
a[,method:=factor(method,levels=c('UCell','AUCell','singscore'))]
a[,significant:=is.finite(FDR)&FDR<.05]
concord<-a[,.(concordant=all(status=='tested'&is.finite(FDR)&FDR<.05)&&uniqueN(sign(effect))==1),by=panel]
stopifnot(sum(concord$concordant)==4)
a<-merge(a,concord,by='panel',sort=FALSE)
a[,stars:=fifelse(significant,'*','')]
p<-ggplot(a[status=='tested'],aes(effect,panel))+
 geom_vline(xintercept=0,linetype=2,linewidth=.3,colour='#9A9A9A')+
 geom_segment(aes(x=CI_low,xend=CI_high,yend=panel),linewidth=.5,colour='#4178A4')+
 geom_point(aes(shape=concordant),size=1.8,colour='#4178A4')+
 geom_text(aes(x=CI_high,label=stars),hjust=-.4,size=3)+
 geom_text(data=a[status!='tested'],aes(x=0,label='Not testable'),size=2.4,colour='#777777')+
 facet_wrap(~method,nrow=1,scales='free_x')+
 scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),breaks=c('TRUE','FALSE'),labels=c('Significant in all three algorithms; same direction','Other tested set'),name=NULL)+
 scale_y_discrete(labels=setlabs,drop=FALSE)+scale_x_continuous(n.breaks=3,expand=expansion(mult=c(.10,.22)))+
 labs(x='ALS - Control score effect (95% CI)',y=NULL,
 caption='Sex-adjusted models | 8 ALS and 4 control donors\n* BH FDR < 0.05 across 14 sets per algorithm; interval scales differ between algorithms')+
 theme(axis.text.y=element_markdown(size=7),panel.spacing.x=unit(4,'mm'),legend.position='bottom',legend.box='vertical',plot.caption=element_text(size=7,hjust=0))+
 guides(shape=guide_legend(ncol=1))
savefig(p,'FigureS5_bulk_three_algorithms',150)

a<-rd('Gene_set_score_results')[layer=='snRNA'&method!='UCell']
a[,context:=factor(context,levels=ct[1:11])];a[,panel:=factor(panel,levels=rev(key$panel))]
pp<-lapply(c('AUCell','singscore'),function(mm)ggplot(a[method==mm],aes(context,panel))+
 geom_tile(aes(fill=effect),colour='white',linewidth=.15)+
 geom_point(data=a[method==mm&is.finite(FDR)&FDR<.05],shape=8,size=1.3)+
 scale_x_discrete(labels=ctshort)+scale_y_discrete(labels=setlabs)+divfill(paste(mm,'effect'))+
 labs(title=mm,x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right',axis.text.y=element_markdown(size=7),axis.text.x=element_text(angle=45,hjust=1,size=7)))
savefig(wrap_plots(pp,ncol=1)+plot_annotation(caption='* BH FDR < 0.05 across 154 comparisons per algorithm',theme=theme(plot.caption=element_text(size=7,hjust=0))),'FigureS2_all_type_scores',165)
