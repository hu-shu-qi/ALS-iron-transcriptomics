suppressPackageStartupMessages(library(ggrepel))
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(axis.text=element_text(size=7,colour='#272727'),plot.title=element_text(size=8,face='bold'),legend.text=element_text(size=7),plot.margin=margin(5,5,5,5)))
savefig<-function(p,n,h){for(ext in c('pdf','svg','png','tiff')){
 f<-file.path(plots,paste0(n,'.',ext))
 if(ext=='pdf')cairo_pdf(f,width=170/25.4,height=h/25.4,family='Arial')
 if(ext=='svg')svglite::svglite(f,width=170/25.4,height=h/25.4)
 if(ext=='png')ragg::agg_png(f,width=170,height=h,units='mm',res=300)
 if(ext=='tiff')ragg::agg_tiff(f,width=170,height=h,units='mm',res=600,compression='lzw')
 print(p);dev.off()
}}
b<-fread(file.path(source_root,'Figure2_all_genes.csv'))
b[,status:=fcase(FDR<.01&logFC>.585,'Up',FDR<.01&logFC<(-.585),'Down',default='Other')]
i<-b[iron==TRUE];i[,status:=fcase(FDR<.05&logFC>=.585,'Up',FDR<.05&logFC<=(-.585),'Down',default='Other')]
stopifnot(nrow(i)==510,sum(b$status!='Other')==3004,sum(i$status!='Other')==116)
vol<-function(z,title,q,labels=FALSE){
 stopifnot(all(is.finite(z$FDR)),all(z$FDR>0))
 p<-ggplot(z,aes(logFC,-log10(FDR)))+geom_point(aes(colour=status),size=.55,alpha=.7)+scale_colour_manual(values=c(Down='#3876A6',Other='#C5C5C5',Up='#B74E45'),breaks=c('Down','Other','Up'))+geom_hline(yintercept=-log10(q),linetype=2,linewidth=.3)+geom_vline(xintercept=c(-.585,.585),linetype=2,linewidth=.3)+labs(title=title,x='ALS - Control log2 fold change',y=expression(-log[10]('BH FDR')),colour=NULL)+theme(legend.position='bottom')
 p<-p+guides(colour=guide_legend(override.aes=list(size=2.2,alpha=1)))
 if(labels)p<-p+geom_text_repel(data=z[gene%in%c('C3','CTSB','HAMP','SLC11A1','TF')],aes(label=gene),seed=20260918,size=2.5,fontface='bold',min.segment.length=0,max.overlaps=Inf)
 p
}
pa<-vol(b,'a  All tested genes',.01)+labs(subtitle='17,960 genes | 1,362 up | 1,642 down')
pb<-vol(i,'b  Iron-related genes',.05,TRUE)+labs(subtitle=paste0('510 detected | ',sum(i$status=='Up'),' up | ',sum(i$status=='Down'),' down'))
en<-fread(file.path(source_root,'Conditional_enrichment_all_28.csv'))[stage=='bulk'];key<-rd('Gene_set_label_key');en<-merge(en,key,by='panel');en[,label:=paste(set_id,short_label)];en[,label:=factor(label,levels=rev(en$label[order(en$set_id)]))]
stopifnot(nrow(en)==14,all(is.finite(en$FDR28)),all(en$FDR28>0 & en$FDR28<=1))
set_colours<-c(GOBP='#008577',GOMF='#3376B8',HALLMARK='#AE7900',MODULE='#9151A5',REACTOME='#D06030')
set_backgrounds<-setNames(vapply(set_colours,function(z){v<-round(col2rgb(z)*.20+255*.80);rgb(v[1],v[2],v[3],maxColorValue=255)},character(1)),names(set_colours))
label_source<-sub('_.*','',en$panel)
stopifnot(all(label_source %in% names(set_colours)))
axis_labels<-setNames(sprintf("<span style='color:%s'><b>%s</b></span> <b>%s</b>",set_colours[label_source],en$set_id,en$short_label),as.character(en$label))
source_key<-paste(sprintf("<span style='color:%s;font-size:10pt'>■</span> %s",set_backgrounds,c('GO-BP','GO-MF','Hallmark','Module','Reactome')),collapse=' &nbsp; ')
pc<-ggplot(en,aes(fold_enrichment,label))+
 geom_rect(data=data.frame(y=seq_along(levels(en$label)),label=levels(en$label)),
   aes(xmin=-Inf,xmax=Inf,ymin=y-.5,ymax=y+.5),inherit.aes=FALSE,
   fill=unname(set_backgrounds[label_source[match(levels(en$label),as.character(en$label))]]),colour='white',linewidth=.3)+
 geom_vline(xintercept=1,linetype=2,linewidth=.3,colour='#999999')+
 geom_point(aes(size=overlap_n,fill=-log10(FDR28)),shape=21,colour='#526071',stroke=.3)+
 geom_text(data=en[FDR28<.05],label='*',nudge_x=.1,size=3.5,vjust=.75)+
 scale_fill_gradientn(colours=c('#C6DFEA','#9186B5','#A63660'),name=expression(-log[10]('BH FDR')),guide=guide_colourbar(order=1,barheight=unit(26,'mm'),barwidth=unit(3,'mm')))+
 scale_size_area(max_size=4,name='Candidate overlap',guide=guide_legend(order=2,override.aes=list(fill='#B8C6D5')))+
 scale_y_discrete(labels=axis_labels)+
 labs(title='c  Conditional enrichment of bulk candidates',subtitle='116 candidates / 510 detected iron-related genes | BH across 28 tests',caption=paste0('&#42; BH FDR &lt; 0.05<br><b>Set:</b> ',source_key),x='Fold enrichment',y=NULL)+
 theme(legend.position='right',legend.spacing.y=unit(5,'mm'),axis.text.y=ggtext::element_markdown(size=7,colour='#272727'),plot.caption=ggtext::element_markdown(size=7,hjust=0,lineheight=1.3))
pa<-pa+labs(subtitle=NULL,caption='FDR < 0.01; |log2 FC| > 0.585')+theme(plot.caption=element_text(size=6.5,hjust=0))
pb<-pb+labs(subtitle=NULL,caption='FDR < 0.05; |log2 FC| ≥ 0.585')+theme(plot.caption=element_text(size=6.5,hjust=0))
pc<-pc+labs(subtitle=NULL)
savefig(wrap_elements(full=(pa|pb))/pc+plot_layout(heights=c(1,1.1)),'Figure2_bulk_candidates',180)
