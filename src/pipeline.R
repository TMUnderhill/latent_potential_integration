#!/usr/bin/env Rscript

# Load dependencies
library(argparse)
library(Matrix)
library(Seurat)
library(plyr)
library(limma)
library(edgeR)
library(dplyr)
library(svglite)
library(tools)
library(Signac)

parser <- ArgumentParser(description = "Pipeline for generating latent potential RNA vs. ATAC plots")
parser$add_argument("--rdata", dest = "rdata", help = "Path to .RData file containing RNA and ATAC sample")
parser$add_argument("--celltype", dest = "celltype", help = "Cell-type key")
parser$add_argument("--tssranges", dest = "tssranges", help = "Path to TSS ranges .RData file")
parser$add_argument("--fragments", dest = "fragments", help = "Path to fragments .tsv.gz file corresponding to .RData file")
parser$add_argument("--overlapping", dest = "overlapping", help = "Path to .csv file listing overlapping TSS")
parser$add_argument("--nbatches", dest = "nbatches", type = "integer", help = "Number of pseudoreplicates")
parser$add_argument("--cpm_cutoff", dest = "cpm_cutoff", type = "double", help = "CPM cutoff")
parser$add_argument("--srand", dest = "srand", type = "integer", help = "Random seed")
args <- parser$parse_args()
# args <- parser$parse_args(strsplit('--rdata ../data/bat.Rdata --celltype FAP
# --tssranges ../data/wstssranges.Rdata --fragments
# ../data/fragments_bat.tsv.gz --overlapping ../data/overlapingTSS500.csv
# --nbatches 10 --cpm_cutoff 2.5 --srand 0', split = ' ')[[1]])

message("Loading RData file: ", args$rdata)
load(args$rdata)
message("Loading TSSranges file: ", args$tssranges)
load(args$tssranges)

# Fetch tissue and celltype name
tissue_name <- file_path_sans_ext(basename(args$rdata))
celltype_name <- args$celltype
# build output file name prefix as {tissue}_{celltype}
filename_prefix <- paste(tissue_name, celltype_name, sep = "_")
message("Output filename prefix: ", filename_prefix)

# Fetch RNA/ATAC samples now that we have loaded the .RData files
sample_rna <- get(paste(tissue_name, "rna", sep = "_"))
sample_atac <- get(tissue_name)

# A necessary redirect here to local fragment file instead of what was
# originally saved in the RData
sample_atac@assays$peaks@fragments[[1]]@path <- args$fragments

message("Constructing ATAC feature matrix")
counts_atac <- FeatureMatrix(fragments = sample_atac@assays$peaks@fragments, features = wstssranges,
    cells = colnames(sample_atac))

message("Constructing ATAC counts matrix")
# construct lookup table: GRanges range <-> corresponding TSS wstssranges
# contains TSS site +/- 200bp range
genes <- wstssranges$gene_name
names(genes) <- GRangesToString(grange = wstssranges)
rownames(counts_atac) <- make.unique(genes[rownames(counts_atac)])
counts_atac <- counts_atac[rownames(counts_atac) != "", ]

sample_atac[["TSS"]] <- CreateAssayObject(counts = counts_atac)
# pick TSS with most reads and assign that to be the 'ATAC counts'
tot_reads <- rowSums(counts_atac)
tss_names <- sapply(names(tot_reads), function(x) strsplit(x, split = " ")[[1]][[1]])
# now pick the TSS with the most reads if there are more than one such TSS,
# pick the first one
tss_max_idx <- sapply(unique(tss_names), function(x) which(tot_reads == max(tot_reads[tss_names ==
    x]) & tss_names == x)[[1]])
counts_atac <- counts_atac[tss_max_idx, ]
rownames(counts_atac) <- names(tss_max_idx)
# now remove overlapping TSSs (within 500bp)
overlap_tss <- read.csv(args$overlapping)$x
counts_atac <- counts_atac[!(rownames(counts_atac) %in% overlap_tss), ]

message("Constructing RNA counts matrix")
cellidx_rna <- which(grepl(args$celltype, sample_rna@active.ident))
counts_rna <- sample_rna@assays$RNA@counts[, cellidx_rna]

message("Aggregating counts into pseudobulk")
set.seed(args$srand)
batch_idx_rna <- sample.int(args$nbatches, ncol(counts_rna), replace = T)
batch_idx_atac <- sample.int(args$nbatches, ncol(counts_atac), replace = T)

# form pseudoreplicates of ATAC and RNA counts
counts_rna_agg <- lapply(1:args$nbatches, function(x) rowSums(counts_rna[, batch_idx_rna ==
    x]))
counts_atac_agg <- lapply(1:args$nbatches, function(x) rowSums(counts_atac[, batch_idx_atac ==
    x]))
counts_rna_agg_df <- data.frame(counts_rna_agg)
names(counts_rna_agg_df) <- paste0("rna_", c(1:args$nbatches))
counts_rna_agg_df$genename <- rownames(counts_rna_agg_df)
rownames(counts_rna_agg_df) <- NULL
counts_atac_agg_df <- data.frame(counts_atac_agg)
names(counts_atac_agg_df) <- paste0("atac_", c(1:args$nbatches))
counts_atac_agg_df$genename <- rownames(counts_atac_agg_df)
rownames(counts_atac_agg_df) <- NULL

# merge and remove blacklisted genes
rna_atac_summ <- join(counts_rna_agg_df, counts_atac_agg_df, by = "genename", type = "full")
rownames(rna_atac_summ) <- rna_atac_summ$genename
rna_atac_summ$genename <- NULL

# filter genes
keep_genes <- !grepl("Rik", rownames(rna_atac_summ)) & !grepl(paste0("Gm", "[0-9]"),
    rownames(rna_atac_summ)) & !grepl(paste0("Mir", "[0-9]"), rownames(rna_atac_summ)) &
    !grepl("mt-", rownames(rna_atac_summ)) & !grepl("RP23-", rownames(rna_atac_summ)) &
    !grepl("RP24-", rownames(rna_atac_summ)) & !grepl("os$", rownames(rna_atac_summ)) &
    !grepl(paste0("Olfr", "[0-9]"), rownames(rna_atac_summ)) & !grepl(paste0("AC",
    "[0-9]"), rownames(rna_atac_summ))
rna_atac_summ <- rna_atac_summ[keep_genes, ]

# remove genes that are NA in either RNA or ATAC.  alternatively, retain genes
# that got NA in RNA (and replace with 0 later)
rna_atac_summ <- rna_atac_summ[!is.na(rowSums(rna_atac_summ)), ]
# if we want to treat RNA NAs as zeros.  rna_atac_summ <-
# rna_atac_summ[!is.na(rowSums(rna_atac_summ[, grepl('atac',
# names(rna_atac_summ))])), ] rna_atac_summ[is.na(rna_atac_summ)] <- 0

# Pseudobulk DE using RNA vs. ATAC counts.
message("Performing differential expression testing")
d0 <- DGEList(rna_atac_summ)
d0 <- calcNormFactors(d0)
drop <- which(apply(cpm(d0), 1, max) < args$cpm_cutoff)
d <- d0[-drop, ]
assay <- sapply(names(rna_atac_summ), function(x) strsplit(x, split = "_")[[1]][[1]])
# Enforce no intercept in model
mm <- model.matrix(~0 + assay)

# Save figures prior to running voom: do this twice to save as png and svg
png(filename = paste0(filename_prefix, "_voom_before.png"))
y <- voom(d, mm, plot = T, normalize.method = "quantile")
dev.off()
svglite(paste0(filename_prefix, "_voom_before.svg"))
y <- voom(d, mm, plot = T, normalize.method = "quantile")
dev.off()

# voom fit
fit <- lmFit(y, mm)
contr <- makeContrasts(assayatac - assayrna, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# Save figures after running voom: do this twice to save as png and svg
png(filename = paste0(filename_prefix, "_voom_after.png"))
plotSA(tmp, main = "Final model: Mean-variance trend")
dev.off()
svglite(paste0(filename_prefix, "_voom_after.svg"))
plotSA(tmp, main = "Final model: Mean-variance trend")
dev.off()

# Write table of top DE genes
top <- topTable(tmp, sort.by = "logFC", n = Inf)
top <- top[order(top$logFC, decreasing = T), ]

write.table(top, file = paste0(filename_prefix, "_toptable.csv"))

