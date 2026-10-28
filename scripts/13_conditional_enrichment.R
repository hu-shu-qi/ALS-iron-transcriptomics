a <- commandArgs(trailingOnly=TRUE)
if (!length(a) %in% c(3L,4L)) stop('Usage: Rscript 13_conditional_enrichment.R BULK_DE_CSV SPATIAL_CANDIDATES_CSV NEW_OUTPUT_DIR [GMT]')
suppressPackageStartupMessages(library(data.table))
script <- gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
code <- dirname(normalizePath(script))
gmt <- if (length(a)==4L) a[4] else file.path(dirname(code),'resources/iron_gene_sets.gmt')
stopifnot(all(file.exists(c(a[1:2],gmt))), !file.exists(a[3]))
lines <- strsplit(readLines(gmt,warn=FALSE),'\t',fixed=TRUE)
stopifnot(length(lines)==14L, all(lengths(lines)>2L))
sets <- setNames(lapply(lines,function(z) unique(z[-(1:2)])),vapply(lines,`[`,character(1),1L))
stopifnot(!anyDuplicated(names(sets)), length(unique(unlist(sets)))==545L)
b <- fread(a[1]); sp <- fread(a[2])
stopifnot(all(c('gene','FDR','logFC') %in% names(b)), !anyDuplicated(b$gene))
required <- c('gene','region','model','PValue','logFC')
stopifnot(all(required %in% names(sp)), !anyDuplicated(sp[,.(gene,region,model)]))
iron_symbols <- unique(unlist(sets))
u_bulk <- sort(b[gene %in% iron_symbols & is.finite(FDR) & is.finite(logFC),gene])
cand <- b[gene %in% u_bulk & FDR<0.05 & abs(logFC)>=0.585,.(gene,bulk_direction=sign(logFC))]
regions <- c('Anterior_horns','Grey_matter','White_matter','Posterior_white_matter')
models <- c('primary','sex_complete','age_sex_complete','depth','exclude_comorbidity')
stopifnot(setequal(unique(sp$region),regions), setequal(unique(sp$model),models))
grid <- merge(CJ(gene=cand$gene,region=regions,model=models),sp[,..required],by=c('gene','region','model'),all.x=TRUE)
grid <- merge(grid,cand,by='gene')
stopifnot(all(!is.finite(grid$PValue) | (grid$PValue>=0 & grid$PValue<=1)))
grid[,FDR_candidates464:=p.adjust(ifelse(is.finite(PValue),PValue,1),'BH'),by=model]
grid[!is.finite(PValue),FDR_candidates464:=NA_real_]
support <- grid[,.(supported=.N==length(models) && all(is.finite(PValue)) && all(is.finite(logFC)) && all(sign(logFC)==bulk_direction) && all(FDR_candidates464<0.05)),by=.(gene,region)]
u_spatial <- sort(unique(grid[is.finite(PValue),gene]))
q_spatial <- sort(unique(support[supported==TRUE,gene]))
stopifnot(all(q_spatial %in% u_spatial))
results <- list(); memberships <- list()
for (stage in c('bulk','spatial')) {
  universe <- if (stage=='bulk') u_bulk else u_spatial
  selected <- if (stage=='bulk') sort(cand$gene) else q_spatial
  stopifnot(length(selected)>0L, all(selected %in% universe))
  memberships[[stage]] <- data.table(stage=stage,gene=universe,selected=universe %in% selected)
  for (panel in names(sets)) {
    gs <- intersect(sets[[panel]],universe); hit <- sort(intersect(selected,gs))
    counts <- c(length(hit),length(selected)-length(hit),length(gs)-length(hit),length(universe)-length(selected)-length(gs)+length(hit))
    stopifnot(all(counts>=0L),sum(counts)==length(universe))
    tab <- matrix(counts,nrow=2L,byrow=TRUE)
    ft <- fisher.test(tab,alternative='greater'); ci <- fisher.test(tab,alternative='two.sided')
    results[[length(results)+1L]] <- data.table(stage=stage,panel=panel,selected_n=length(selected),universe_n=length(universe),set_n=length(gs),overlap_n=length(hit),fold_enrichment=if(length(gs)) (length(hit)/length(selected))/(length(gs)/length(universe)) else NA_real_,OR=unname(ft$estimate),CI_low=ci$conf.int[1],CI_high=ci$conf.int[2],P=ft$p.value,genes=paste(hit,collapse=';'),selected_in_set=counts[1],selected_outside_set=counts[2],unselected_in_set=counts[3],unselected_outside_set=counts[4])
  }
}
en <- rbindlist(results); en[,FDR28:=p.adjust(P,'BH')]
stopifnot(nrow(en)==28L)
dir.create(a[3],recursive=TRUE)
fwrite(en,file.path(a[3],'Conditional_enrichment_all_28.csv'))
fwrite(en[stage=='bulk'],file.path(a[3],'Figure2_enrichment.csv'))
fwrite(rbindlist(memberships),file.path(a[3],'enrichment_background_membership.csv'))
fwrite(support[supported==TRUE],file.path(a[3],'spatial_supported_pairs.csv'))
inputs <- normalizePath(c(a[1:2],gmt,script))
fwrite(data.table(file=basename(inputs),md5=unname(tools::md5sum(inputs))),file.path(a[3],'input_manifest.csv'))
writeLines(capture.output(sessionInfo()),file.path(a[3],'sessionInfo.txt'))
print(unique(en[,.(stage,selected_n,universe_n)]))
