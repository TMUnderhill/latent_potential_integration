#!/usr/bin/env Rscript
library(argparse)
library(Matrix)
parser <- ArgumentParser(description='pipeline')
parser$add_argument("--rdata", dest = "rdata")
parser$add_argument("--celltype", dest = "celltype")
parser$add_argument("--tssranges", dest = "tssranges")
parser$add_argument("--fragments", dest = "fragments")
parser$add_argument("--overlapping", dest = "overlapping")
parser$add_argument("--nbatches", dest = "nbatches", type = "integer")
parser$add_argument("--cpm_cutoff", dest = "cpm_cutoff", type = "double")
parser$add_argument("--srand", dest = "srand", type = "integer")
# args <- parser$parse_args(strsplit("--rdata bat.Rdata --celltype FAP --tssranges wstssranges.Rdata --fragments fragments_bat.tsv.gz --overlapping overlapingTSS500.csv --nbatches 10 --cpm_cutoff 2.5 --srand 0", split = " ")[[1]])
args <- parser$parse_args()

message(paste("Loading RData file: ", args$rdata))
load(args$rdata)
message(paste("Loading TSSranges file: "), args$tssranges)
load(args$tssranges)

# tissue name
tissue_name <- strsplit(args$rdata, '\\.')[[1]][1]
celltype_name <- args$celltype
filename_prefix <- paste(tissue_name, celltype_name, sep = "_")

sample_rna <- get(paste(tissue_name, "rna", sep = '_'))
sample_atac <- get(tissue_name)

# some necessary redirects here to local fragment file 
sample_atac@assays$peaks@fragments[[1]]@path <- args$fragments 

message("Constructing ATAC FeatureMatrix")
counts_atac <- FeatureMatrix(fragments = sample_atac@assays$peaks@fragments, features = wstssranges, cells = colnames(sample_atac))

message("Constructing ATAC counts matrix")
library(Seurat)
# construct lookup table GRanges <--> corresponding TSS
# wstssranges contains TSS +/- 200bp
genes <- wstssranges$gene_name
names(genes) <- GRangesToString(grange = wstssranges)
rownames(counts_atac) <- make.unique(genes[rownames(counts_atac)])
counts_atac <- counts_atac[rownames(counts_atac)!="",]

sample_atac[['TSS']] <- CreateAssayObject(counts = counts_atac)
# pick TSS with most reads and assign that to be the 'ATAC counts'
tot_reads <- rowSums(counts_atac)
tss_names <- sapply(names(tot_reads), function(x) strsplit(x, split = " ")[[1]][[1]])
# now pick the TSS with the most reads
# if there are more than one such TSS, pick the first one (??)
tss_max_idx <- sapply(unique(tss_names), function(x) which(tot_reads == max(tot_reads[tss_names == x]) & tss_names == x)[[1]])
counts_atac <- counts_atac[tss_max_idx, ]
rownames(counts_atac) <- names(tss_max_idx)
# now remove overlapping 500bp
overlap_tss <- read.csv(args$overlapping)$x
counts_atac <- counts_atac[!(rownames(counts_atac) %in% overlap_tss), ]

message("Constructing RNA counts matrix")
cellidx_rna <- which(grepl(args$celltype, sample_rna@active.ident))
counts_rna <- sample_rna@assays$RNA@counts[, cellidx_rna]

message("Aggregating counts")
set.seed(args$srand)
batch_idx_rna <- sample.int(args$nbatches, ncol(counts_rna), replace = T)
batch_idx_atac <- sample.int(args$nbatches, ncol(counts_atac), replace = T)

# form pseudoreplicates
counts_rna_agg <- lapply(1:args$nbatches, function(x) rowSums(counts_rna[, batch_idx_rna == x]))
counts_atac_agg <- lapply(1:args$nbatches, function(x) rowSums(counts_atac[, batch_idx_atac == x]))
# find common genes
# should replace this with outer join and remove genes that got 0 reads for ATAC
library(plyr)
# genes_common <- intersect(names(counts_rna_agg[[1]]), names(counts_atac_agg[[1]]))
# drop blacklisted genes
# genes_common <- genes_common[!grepl("Rik", genes_common) & !grepl(paste0("Gm","[0-9]"), genes_common) & !grepl(paste0("Mir","[0-9]"), genes_common) & !grepl("mt-", genes_common) & !grepl("RP23-", genes_common) & !grepl("RP24-", genes_common) & !grepl("os$", genes_common) & !grepl(paste0("Olfr","[0-9]"), genes_common)& !grepl(paste0("AC","[0-9]"), genes_common)]
# rna_atac_summ <- lapply(1:args$nbatches, function(i) data.frame(rna = counts_rna_agg[[i]][genes_common], atac = counts_atac_agg[[i]][genes_common]))

counts_rna_agg_df <- data.frame(counts_rna_agg);
names(counts_rna_agg_df) <- paste0("rna_", c(1:args$nbatches))
counts_rna_agg_df$genename <- rownames(counts_rna_agg_df)
rownames(counts_rna_agg_df) <- NULL
counts_atac_agg_df <- data.frame(counts_atac_agg);
names(counts_atac_agg_df) <- paste0("atac_", c(1:args$nbatches))
counts_atac_agg_df$genename <- rownames(counts_atac_agg_df)
rownames(counts_atac_agg_df) <- NULL
# merge and remove blacklisted genes 
rna_atac_summ <- join(counts_rna_agg_df, counts_atac_agg_df, by = "genename", type = "full")
rownames(rna_atac_summ) <- rna_atac_summ$genename
rna_atac_summ$genename <- NULL
keep_genes <- !grepl("Rik", rownames(rna_atac_summ)) & !grepl(paste0("Gm","[0-9]"), rownames(rna_atac_summ)) & !grepl(paste0("Mir","[0-9]"), rownames(rna_atac_summ)) & !grepl("mt-", rownames(rna_atac_summ)) & !grepl("RP23-", rownames(rna_atac_summ)) & !grepl("RP24-", rownames(rna_atac_summ)) & !grepl("os$", rownames(rna_atac_summ)) & !grepl(paste0("Olfr","[0-9]"), rownames(rna_atac_summ))& !grepl(paste0("AC","[0-9]"), rownames(rna_atac_summ))
rna_atac_summ <- rna_atac_summ[keep_genes, ]
# remove genes that are NA in either RNA or ATAC.
# optionally retain genes that got NA in RNA (and replace with 0 later)
# rna_atac_summ <- rna_atac_summ[!is.na(rowSums(rna_atac_summ[, grepl("rna", names(rna_atac_summ))])), ] 
rna_atac_summ <- rna_atac_summ[!is.na(rowSums(rna_atac_summ)), ] 
# rna_atac_summ[is.na(rna_atac_summ)] <- 0 # if we want to keep RNA

library(limma)
library(edgeR)
library(dplyr)
library(svglite)

d0 <- DGEList(rna_atac_summ)
d0 <- calcNormFactors(d0)
drop <- which(apply(cpm(d0), 1, max) < args$cpm_cutoff)
d <- d0[-drop, ]
assay <- sapply(names(rna_atac_summ), function(x) strsplit(x, split = "_")[[1]][[1]])

mm <- model.matrix(~ 0 + assay)
# do this twice to save as png and svg 
png(filename=paste0(filename_prefix, "_voom_before.png"))
y <- voom(d, mm, plot = T, normalize.method = "quantile")
dev.off()
svglite(paste0(filename_prefix, "_voom_before.svg"))
y <- voom(d, mm, plot = T, normalize.method = "quantile")
dev.off()

fit <- lmFit(y, mm)
contr <- makeContrasts(assayatac - assayrna, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# do this twice to save as png and svg 
png(filename=paste0(filename_prefix, "_voom_after.png"))
plotSA(tmp, main = "Final model: Mean-variance trend")
dev.off()
svglite(paste0(filename_prefix, "_voom_after.svg"))
plotSA(tmp, main = "Final model: Mean-variance trend")
dev.off()
top <- topTable(tmp, sort.by = "logFC", n = Inf)
top <- top[order(top$logFC, decreasing = T), ]

write.table(top, file = paste0(filename_prefix, "_toptable.csv"))
