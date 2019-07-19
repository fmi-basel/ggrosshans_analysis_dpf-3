args <- commandArgs(TRUE)
input_gtf <- args[1]
input_tr_id_gene_id <- args[2]
output_gtf <- args[3]

# for testing
# input_gtf <- "c_elegans.WS220.annotations.trs.exon.corrected.gtf"
# input_tr_id_gene_id <- "/work/gbioinfo/DB/WormBase/WS220/mart_export_WS220_WBtoTR.txt"
# output_gtf <- "c_elegans.WS220.annotations.trs.exon.corrected.with_gene_id.gtf"

# read gtf file
gtf <- read.delim(input_gtf, header=FALSE, as.is=TRUE)
names(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# read gene id to transcript id
gene_id_transcript_id <- read.delim(input_tr_id_gene_id, header=TRUE, as.is=TRUE)
names(gene_id_transcript_id) <- c("gene_id", "transcript_id")

# create a transcirpt id column in the main gtf file
gtf[,"transcript_id"] <- sub("transcript_id ","", unlist(lapply(strsplit(gtf[,"attribute"], ";"), '[[', 1)))
gtf[,"transcript_biotype"] <- sub(" transcript_biotype ","", unlist(lapply(strsplit(gtf[,"attribute"], ";"), '[[', 2)))
gtf[,"gene_id"] <- "None"


for (row in 1:nrow(gtf)) {
  transcript_id = gtf[row, "transcript_id"]
  gene_id = gene_id_transcript_id[gene_id_transcript_id[,"transcript_id"] == transcript_id,][,'gene_id']
  # print(transcript_id)
  # print(gene_id)
  if (length(gene_id) > 1){
    gene_id = unlist(strsplit(gene_id, " ")[1])
  }
  gtf[row, "gene_id"] = gene_id
}

# merge the two dataframes
# final_gtf <- merge(gtf, gene_id_transcript_id, by="transcript_id")


gtf[,"attribute"] <- paste("gene_id \"", gtf[,"gene_id"], "\";",
      " transcript_id \"", gtf[,"transcript_id"], "\";",
      " transcript_biotype \"", gtf[,"transcript_biotype"], "\";", sep="")

final_gtf <- gtf[ , !(names(gtf) %in% c("gene_id", "transcript_id", "transcript_biotype"))]

write.table(final_gtf, file=output_gtf, quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
