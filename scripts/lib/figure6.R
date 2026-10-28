suppressPackageStartupMessages(library(ggtext))

blue<-'#3376B8'; coral<-'#D06030'
marker_fill<-c(Ferritin='#F3F9F7',`Serum iron`='#FDF5F1',TSAT='#F8F4FA')

flow<-rd('Figure6_cohort_flow')
cohort_n<-flow[metric=='analysis_cohort',value]
marker_n<-setNames(flow$value,flow$metric)
stopifnot(length(cohort_n)==1L,all(is.finite(marker_n[c('main_ferritin','main_iron','main_tsat')])))
boxes<-data.table(
 x=c(2,1,2,3),y=c(2.02,1,1,1),
 label=c(sprintf('<b>%d patients</b><br>with evaluable follow-up functional decline',cohort_n),
         sprintf('<b>Ferritin</b><br>n = %d',marker_n['main_ferritin']),sprintf('<b>Serum iron</b><br>n = %d',marker_n['main_iron']),sprintf('<b>TSAT</b><br>n = %d',marker_n['main_tsat'])),
 fill=c('#E8F0F8',rep('#F3F6F9',3)),edge=rep('#69849B',4))
pf<-ggplot()+
 geom_segment(aes(x=2,y=1.70,xend=2,yend=1.47),linewidth=.45,colour='#7A858F')+
 geom_segment(aes(x=1,y=1.47,xend=3,yend=1.47),linewidth=.45,colour='#7A858F')+
 geom_segment(data=data.table(x=1:3),aes(x=x,y=1.47,xend=x,yend=1.27),linewidth=.45,colour='#7A858F')+
 ggtext::geom_richtext(data=boxes,aes(x,y,label=label,fill=fill),colour='#272727',label.colour='#69849B',size=2.75,family='Arial',
            label.size=.3,label.r=unit(2,'mm'),label.padding=unit(2.2,'mm'),lineheight=1.1)+
 scale_fill_identity()+scale_colour_identity()+
 coord_cartesian(xlim=c(.35,3.65),ylim=c(.52,2.52),clip='off')+theme_void(base_family='Arial')+
 labs(title='a  Clinical analysis sets',caption='Overlapping marker-specific subsets')+
 theme(plot.caption=element_text(size=7,hjust=.5),plot.title=element_text(size=8,face='bold'),plot.margin=margin(3,4,1,4))

a<-rd('Figure6_correlations')
stopifnot(nrow(a)==6,all(is.finite(as.matrix(a[,.(rho,ci_low,ci_high,n)]))),all(a$ci_low<=a$rho),all(a$ci_high>=a$rho))
a[,marker_label:=factor(fcase(marker=='ferritin','Ferritin',marker=='iron','Serum iron',default='TSAT'),
                        levels=c('TSAT','Serum iron','Ferritin'))]
a[,y:=as.numeric(marker_label)+ifelse(analysis=='main',.13,-.13)]
a[,analysis:=factor(analysis,levels=c('main','strict_sensitivity'))]
a[,value_label:=sprintf('%+.3f [%+.3f, %+.3f]',rho,ci_low,ci_high)]
bands<-data.table(y=1:3,fill=unname(marker_fill[c('TSAT','Serum iron','Ferritin')]))
pr<-ggplot(a,aes(rho,y,colour=analysis,shape=analysis))+
 geom_rect(data=bands,aes(xmin=-.36,xmax=.76,ymin=y-.42,ymax=y+.42,fill=fill),
           inherit.aes=FALSE,colour=NA)+scale_fill_identity()+
 geom_vline(xintercept=0,linetype=2,linewidth=.35,colour='#8D969E')+
 geom_segment(aes(x=ci_low,xend=ci_high,yend=y),linewidth=.75,lineend='round')+
 geom_point(size=2.5,stroke=.7)+
 geom_text(aes(x=.285,label=paste0('n = ',n)),hjust=0,colour='#272727',size=2.5,family='Arial',show.legend=FALSE)+
 geom_text(aes(x=.415,label=value_label),hjust=0,colour='#272727',size=2.5,family='Arial',show.legend=FALSE)+
 scale_colour_manual(values=c(main=blue,strict_sensitivity=coral),labels=c('Main analysis','Strict sensitivity'))+
 scale_shape_manual(values=c(main=16,strict_sensitivity=18),labels=c('Main analysis','Strict sensitivity'))+
 scale_y_continuous(breaks=1:3,labels=c('TSAT','Serum iron','Ferritin'),limits=c(.5,3.62),expand=c(0,0))+
 scale_x_continuous(breaks=c(-.3,-.2,-.1,0,.1,.2),labels=function(x)sprintf('%.1f',x),limits=c(-.36,.76),expand=c(0,0))+
 labs(title='b  Association with follow-up functional decline',x='Spearman rho (95% bootstrap CI)',y=NULL,colour=NULL,shape=NULL)+
 annotate('text',x=.285,y=3.48,label='Patients',hjust=0,size=2.5,fontface='bold',family='Arial')+
 annotate('text',x=.415,y=3.48,label='rho [95% CI]',hjust=0,size=2.5,fontface='bold',family='Arial')+
 theme(legend.position='bottom',legend.direction='horizontal',legend.box.spacing=unit(1,'mm'),
       panel.background=element_rect(fill='white',colour=NA),axis.line.y=element_blank(),axis.ticks.y=element_blank(),
       axis.text.y=element_text(face='bold'),plot.margin=margin(2,4,2,4))

p<-pf/pr+plot_layout(heights=c(.78,1.42))

savefig(p,'Figure6_clinical_associations',115)
