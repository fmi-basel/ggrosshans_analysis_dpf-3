#!/usr/bin/env Rscript

########################
### LOAD PACKAGES	 ###
########################
rm(list=ls())
library(optparse)
library("edgeR")
library("gplots")
library("RColorBrewer")
library("rtracklayer")

########################
### PREPROCESSING	 ###
########################

option_list <- list(
  make_option(c("--conditions"), action="store", type="character", help="condition names"),
  make_option(c("--counts"), action="store", type="character", help="expression table"),
  make_option(c("--outfolder"), action="store", type="character", help="output folder")
 )
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS]",
                           option_list = option_list, add_help_option=TRUE)
opt <- parse_args(opt_parser)

# Estimated number of reads for all replicates and conditions
x <- read.csv(opt$counts, header=TRUE, row.names=1)

# Condition names
conds <- read.table(opt$conditions, sep=',')
conds <- c(t(conds))

names <- unique(conds)
c1 <- length(conds[conds %in% names[1]])
c2 <- length(conds[conds %in% names[2]])
conds_comma <- paste(conds, collapse=",")

group <- factor(c(rep(names[1],c1),rep(names[2],c2)))
group <- relevel(group, ref=names[1])

#design <- model.matrix(~0 + group)
design <- model.matrix(~group)

y <- DGEList(counts=x, group=group)

# Normalization
y <- calcNormFactors(y)

# create MDS plot
pdf(paste(opt$outfolder, "MDS.pdf", sep='/'))
plotMDS(y)
dev.off()

png(paste(opt$outfolder, "MDS.png", sep='/'))
plotMDS(y)
dev.off()

y <- estimateDisp(y, design, robust=TRUE)

# create BCV plot
pdf(paste(opt$outfolder, "BCV.pdf", sep='/'))
plotBCV(y)
dev.off()

png(paste(opt$outfolder, "BCV.png", sep='/'))
plotBCV(y)
dev.off()


# Testing for DE genes
fit <- glmFit(y, design)

# Conduct likelihood ratio tests for KO vs WT and show the top genes:
lrt <- glmLRT(fit)
lrt_df <- lrt$table
#print(group)
#print(design)
#print(head(lrt_df))
#quit()


# Create a table with standard p-values
pst <- data.frame(lrt_df$PValue)

# Create a table with adjusted p-values (FDR)
pval <- data.frame(FDR=p.adjust(lrt_df$PValue, method="BH"))
rownames(pval) <- rownames(lrt_df)

# keep <- pval$FDR<=0.05
# pval_sig <- lrt_df[keep,]
# create a big table

lrt_df_with_FDR <- merge(lrt_df,pval,by="row.names",all.x=TRUE)

colnames(lrt_df_with_FDR) <- c("id", "logFC", "logCPM", "LR", "PValue", "FDR")

# ensembl_gtf <- import(opt$gtf)
# ensembl_gene_name_type <- as.data.frame(ensembl_gtf[,c("transcript_id","source","gene_id","gene_biotype","type")])
# ensembl_gene_name_type_small <- ensembl_gene_name_type[,c("gene_id","gene_biotype","type","seqnames","start","end")]
# ensembl_gene_name_type_small_only_genes <- ensembl_gene_name_type_small[ensembl_gene_name_type_small$type=="gene",]
# ensembl_gene_name_type_small_only_genes_unique <- unique(ensembl_gene_name_type_small_only_genes)
# colnames(ensembl_gene_name_type_small_only_genes_unique) <- c("id","gene_biotype","type","chromosome","start","end")
# final_table <- merge(lrt_df_with_FDR, ensembl_gene_name_type_small_only_genes_unique, by="id", all.x=TRUE)

write.table(lrt_df_with_FDR,
			file = paste(opt$outfolder,"final_table.tsv", sep='/'),
			quote = FALSE,
			sep = "\t",
			eol = "\n",
			row.names = FALSE,
			col.names = TRUE)


final_table_FDR_low <- lrt_df_with_FDR[lrt_df_with_FDR$FDR<0.05,]
write.table(final_table_FDR_low,
			file = paste(opt$outfolder,"final_table_FDR_low.tsv", sep='/'),
			quote = FALSE,
			sep = "\t",
			eol = "\n",
			row.names = FALSE,
			col.names = TRUE)

pdf(paste(opt$outfolder,"volcano_plot.pdf", sep='/'))
with(lrt_df_with_FDR, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", col='gray'))
with(subset(lrt_df_with_FDR, FDR<0.05), points(logFC, -log10(FDR), pch=20, main="Volcano plot", col='black'))
legend("topright", c("All genes", "Significant (FDR <0.05)"), fill=c("gray","black"))
dev.off()

png(paste(opt$outfolder,"volcano_plot.png", sep='/'))
with(lrt_df_with_FDR, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", col='gray'))
with(subset(lrt_df_with_FDR, FDR<0.05), points(logFC, -log10(FDR), pch=20, main="Volcano plot", col='black'))
legend("topright", c("All genes", "Significant (FDR <0.05)"), fill=c("gray","black"))
dev.off()

#
pdf(paste(opt$outfolder,"ma_plot.pdf", sep='/'))
with(lrt_df_with_FDR, plot(logCPM, logFC, pch=20, main="MA plot", col='gray'))
with(subset(lrt_df_with_FDR, FDR<0.05), points(logCPM, logFC, pch=20, col='black'))
legend("bottomright", c("All genes", "Significant (FDR <0.05)"), fill=c("gray","black"))
dev.off()

png(paste(opt$outfolder,"ma_plot.png", sep='/'))
with(lrt_df_with_FDR, plot(logCPM, logFC, pch=20, main="MA plot", col='gray'))
with(subset(lrt_df_with_FDR, FDR<0.05), points(logCPM, logFC, pch=20, col='black'))
legend("bottomright", c("All genes", "Significant (FDR <0.05)"), fill=c("gray","black"))
dev.off()
