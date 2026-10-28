a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=3)stop('Usage: 2_prepare_bulk.R NORMALIZED_CSV SAMPLE_CSV NEW_OUTPUT_DIR')
suppressPackageStartupMessages(library(data.table))
stopifnot(all(file.exists(a[1:2])),!file.exists(a[3]))
d<-fread(a[1]);m<-fread(a[2])
required<-c('donor_id','bulk_GSM','expression_sample_id','group','sex')
stopifnot(all(required%in%names(m)),all(c('gene_ID','gene_name')%in%names(d)))
stopifnot(nrow(m)==12,!anyNA(m[,..required]),all(nzchar(as.matrix(m[,..required]))))
for(k in c('donor_id','bulk_GSM','expression_sample_id'))stopifnot(!anyDuplicated(m[[k]]))
stopifnot(all(m$group%in%c('Control','ALS')),sum(m$group=='Control')==4,sum(m$group=='ALS')==8)
stopifnot(all(m$sex%in%c('Female','Male')),all(m$expression_sample_id%in%names(d)))
setorder(m,bulk_GSM);cols<-m$expression_sample_id
m[,group:=factor(group,levels=c('Control','ALS'))];m[,sex:=factor(sex,levels=c('Female','Male'))]
stopifnot(all(vapply(d[,..cols],is.numeric,TRUE)))
x<-as.matrix(d[,..cols]);stopifnot(all(is.finite(x)),all(x>=0))
d[,gene:=toupper(trimws(gene_name))];d[,median_expression:=apply(x,1,median)]
setorder(d,gene,-median_expression,gene_ID)
d[,retained:=seq_len(.N)==1,by=gene]
symbols<-d[,.(gene_ID,gene,median_expression,retained)]
d<-d[retained & !is.na(gene) & gene!='']
x<-as.matrix(d[,..cols]);rownames(x)<-d$gene
keep<-rowSums(x>1)>=3 & apply(x,1,var)>0
filters<-data.table(gene=rownames(x),retained=keep)
x<-log2(x[keep,,drop=FALSE]+1)
stopifnot(nrow(x)>0,!anyDuplicated(rownames(x)),identical(colnames(x),m$expression_sample_id))
design<-model.matrix(~sex+group,m);stopifnot(qr(design)$rank==ncol(design))
dir.create(a[3],recursive=TRUE);out<-normalizePath(a[3])
saveRDS(list(expression=x,metadata=m),file.path(out,'bulk_expression.rds'))
fwrite(m,file.path(out,'bulk_metadata.csv'))
fwrite(symbols,file.path(out,'bulk_symbol_selection.csv'))
fwrite(filters,file.path(out,'bulk_expression_filter.csv'))
fwrite(data.table(input_rows=nrow(symbols),unique_symbols=nrow(filters),retained_genes=nrow(x),donors=ncol(x),controls=sum(m$group=='Control'),ALS=sum(m$group=='ALS')),file.path(out,'summary.csv'))
fwrite(data.table(input_mode='DESeq2-normalized expression',transform='log2(x+1)',deduplication='highest median; gene_ID ascending tie break',filter='expression>1 in >=3 samples; variance>0'),file.path(out,'processing.csv'))
script<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
inputs<-normalizePath(c(a[1:2],script))
fwrite(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_manifest.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('COMPLETE: ',nrow(x),' genes, ',ncol(x),' donors')
