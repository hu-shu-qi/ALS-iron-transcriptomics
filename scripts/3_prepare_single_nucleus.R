a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2L) stop('Usage: Rscript 3_prepare_single_nucleus.R PUBLIC_H5_DIRECTORY NEW_OUTPUT_DIR')
suppressPackageStartupMessages({library(Seurat);library(BPCells);library(data.table)})
stopifnot(dir.exists(a[1]),!file.exists(a[2]))
script <- gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
config <- file.path(dirname(dirname(normalizePath(script))),'config')
annotation <- fread(file.path(config,'sn_frozen_annotation.csv.gz'))
genes <- fread(file.path(config,'sn_frozen_features.csv'))$gene
gene_map <- fread(file.path(config,'sn_input_feature_names.csv'))
stopifnot(!anyDuplicated(annotation$cell_id),!anyDuplicated(genes),nrow(annotation)==100265L)
samples <- sort(unique(annotation$orig.ident))
files <- file.path(a[1],paste0(samples,'_filtered_feature_bc_matrix.h5'))
stopifnot(all(file.exists(files)))
dir.create(a[2],recursive=TRUE);out<-normalizePath(a[2])
dir.create(file.path(out,'bpcells_raw'))
objects <- list()
for(i in seq_along(samples)) {
  sample <- samples[i]; message('Rebuilding ',sample)
  counts <- open_matrix_10x_hdf5(files[i],feature_type='Gene Expression')
  stopifnot(identical(rownames(counts),gene_map$feature_id))
  rownames(counts) <- make.unique(gene_map$symbol)
  stopifnot(!is.null(counts),!anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)))
  backing <- file.path(out,'bpcells_raw',sample)
  bp <- write_matrix_dir(convert_matrix_type(counts,type='uint32_t'),dir=backing)
  obj <- CreateSeuratObject(counts=bp,project=sample,min.cells=3,min.features=200)
  ann <- annotation[orig.ident==sample]
  barcodes <- substring(ann$cell_id,nchar(sample)+2L)
  stopifnot(all(barcodes %in% colnames(obj)))
  obj <- obj[,barcodes]
  objects[[sample]] <- obj
  rm(counts,bp,obj);gc()
}
o <- merge(objects[[1]],y=objects[-1],add.cell.ids=samples)
o <- JoinLayers(o)
stopifnot(setequal(colnames(o),annotation$cell_id),setequal(rownames(o),genes))
o <- o[genes,annotation$cell_id]
stopifnot(identical(colnames(o),annotation$cell_id),identical(rownames(o),genes))
stopifnot(max(abs(o$nCount_RNA-annotation$nCount_RNA))==0,max(abs(o$nFeature_RNA-annotation$nFeature_RNA))==0)
meta <- as.data.frame(annotation[,!c('cell_id','UMAP_1','UMAP_2')]);rownames(meta)<-annotation$cell_id
factor_levels <- fread(file.path(config,'sn_factor_levels.csv'))
for(column_name in unique(factor_levels$field)) {
  lev <- factor_levels[field==column_name][order(order)]
  meta[[column_name]] <- factor(meta[[column_name]],levels=lev$level,ordered=lev$ordered[1])
  stopifnot(!anyNA(meta[[column_name]]))
}
o@meta.data <- meta
embedding <- as.matrix(annotation[,.(UMAP_1,UMAP_2)]);rownames(embedding)<-annotation$cell_id
o[['umap']] <- CreateDimReducObject(embeddings=embedding,key='UMAP_',assay='RNA')
Idents(o) <- o$cell_type
saveRDS(o,file.path(out,'GSE287257_annotated.rds'))
inputs <- c(files,file.path(config,c('sn_frozen_annotation.csv.gz','sn_frozen_features.csv','sn_input_feature_names.csv','sn_factor_levels.csv')),script)
fwrite(data.table(file=basename(inputs),md5=unname(tools::md5sum(inputs))),file.path(out,'input_manifest.csv'))
fwrite(data.table(check=c('frozen_nuclei','frozen_features','per_nucleus_counts','per_nucleus_detected_genes'),value=c(ncol(o),nrow(o),max(abs(o$nCount_RNA-annotation$nCount_RNA)),max(abs(o$nFeature_RNA-annotation$nFeature_RNA)))),file.path(out,'checks.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('Single-nucleus input prepared with the study annotations and UMAP.')
