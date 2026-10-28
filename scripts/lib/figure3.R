 u<-rd('Figure3_umap');u[,cell_type:=factor(cell_type,levels=ct)];u[,group:=factor(group,levels=c('Control','ALS'))]
 set.seed(20260918);u<-u[order(runif(.N))]
 umbase<-function(z)ggplot(z,aes(UMAP1,UMAP2,colour=cell_type))+geom_point_rast(size=.12,alpha=.6,raster.dpi=600)+scale_colour_manual(values=pal,drop=FALSE)+coord_equal(xlim=range(u$UMAP1),ylim=range(u$UMAP2))+theme_void(base_size=7.5,base_family='Arial')+theme(plot.title=element_text(size=8,face='bold'),strip.text=element_text(size=7.5),legend.text=element_text(size=7),legend.key.height=unit(3,'mm'),plot.margin=margin(3,4,3,4))
 ua<-umbase(u)+labs(title='a  Annotated single-nucleus atlas',subtitle='100,265 nuclei | 12 donors',colour=NULL)+guides(colour=guide_legend(ncol=3,override.aes=list(size=1.8,alpha=1)))+theme(legend.position='bottom')
 ub<-umbase(u)+facet_wrap(~group,ncol=2)+labs(title='b  Shared embedding by group',subtitle='4 control donors | 8 ALS donors')+theme(legend.position='none')
 m<-rd('Figure3_marker_dotplot');mk<-rd('Figure3_marker_key');m[,cell_type:=factor(cell_type,levels=rev(ct))];m[,gene:=factor(gene,levels=mk$gene)]
 pc<-ggplot(m,aes(gene,cell_type))+geom_point(aes(size=pct_expressing,fill=z_expression),shape=21,stroke=0)+scale_size_area(max_size=2.8,limits=c(0,100),breaks=c(25,50,100),name='Expressing (%)')+guides(size=guide_legend(override.aes=list(fill='#777777')))+divfill('Mean expression\n(z per gene)')+scale_y_discrete(labels=ctshort)+labs(title='c  Canonical annotation markers',x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=60,hjust=1,size=7),axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right')
 n<-rd('Figure3_donor_composition');nd<-unique(n[,.(donor_id,group)]);setorder(nd,group,donor_id);nd[,display:=paste0(ifelse(group=='Control','C','A'),sprintf('%02d',seq_len(.N))),by=group];n<-merge(n,nd,by=c('donor_id','group'));n[,display:=factor(display,levels=nd$display)];n[,cell_type:=factor(cell_type,levels=ct)];n[,group:=factor(group,levels=c('Control','ALS'))]
 pd<-ggplot(n,aes(display,fraction,fill=cell_type))+geom_col(width=.8)+facet_grid(~group,scales='free_x',space='free_x')+scale_fill_manual(values=pal,drop=FALSE)+scale_y_continuous(labels=function(x)paste0(x*100,'%'),breaks=seq(0,1,.25),expand=c(0,0))+coord_cartesian(ylim=c(0,1))+labs(title='d  Recovered nuclei by donor',x=NULL,y='Nuclear fraction')+theme(legend.position='none',axis.text.x=element_text(size=7))
 f3<-((ua|ub)+plot_layout(widths=c(1.12,1)))/pc/pd+plot_layout(heights=c(1.4,1.35,1))

p<-f3
top<-p[[1]]
ua<-top[[1]];ub<-top[[2]]
ua<-ua+labs(subtitle=NULL);ub<-ub+labs(subtitle=NULL)
ua<-ua+labs(title='a  Single-nucleus atlas')+guides(colour=guide_legend(ncol=4,byrow=TRUE,title=NULL,override.aes=list(size=1.8,alpha=1)))+
 theme(legend.position='bottom',legend.text=element_text(size=7),legend.key.width=unit(4,'mm'),legend.key.height=unit(3,'mm'),legend.key.spacing.x=unit(4,'mm'))
cairo_pdf(tempfile(fileext='.pdf'),width=170/25.4,height=200/25.4,family='Arial')
gt<-ggplotGrob(ua)
dev.off()
legend<-gt$grobs[[which(gt$layout$name=='guide-box-bottom')]]
stopifnot(inherits(legend,'gtable'))
top<-(ua+theme(legend.position='none')|ub+theme(legend.position='none'))+plot_layout(widths=c(1,2))
c<-p[[2]]
cell_scale<-c$scales$get_scales('y')
cell_scale$labels[c('Pericytes','Ependymal')]<-c('Pericytes','Ependymal')
c$layers[[1]]$aes_params$colour<-'transparent'
c$layers[[1]]$aes_params$stroke<-0
c<-c+guides(size=guide_legend(title='Expressing (%)',order=1,override.aes=list(fill='#777777',colour='transparent',stroke=0)),
 fill=guide_colorbar(title='Mean expression\n(z per gene)',direction='horizontal',order=2,
 barwidth=unit(23,'mm'),barheight=unit(3,'mm'),title.position='top',label.position='bottom'))+
 theme(legend.position='right',legend.box='vertical',legend.spacing.y=unit(7,'mm'),
 legend.title=element_text(size=7,lineheight=1.15,margin=margin(b=4)),
 legend.text=element_text(size=7),legend.box.margin=margin(0,3,0,5))
p[[2]]<-c
d<-p[[3]]
d$layers[[1]]$position<-position_stack(reverse=TRUE)
d<-d+guides(fill=guide_legend(ncol=4,byrow=TRUE,title=NULL))+
  theme(legend.position='bottom',legend.text=element_text(size=7),
        legend.key.height=unit(3,'mm'),legend.key.width=unit(4,'mm'),
        legend.margin=margin(2,0,0,0))
p<-wrap_plots(top,wrap_elements(full=legend),c,d,ncol=1,heights=c(1.4,.32,1.35,1.3))
savefig(p,'Figure3_single_nucleus_atlas',200)
