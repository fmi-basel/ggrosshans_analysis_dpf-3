args <- commandArgs(TRUE)
input_gff <- args[1]
output_gtf <- args[2]

exonsGFF <- read.delim(input_gff, header=FALSE, as.is=TRUE)
names(exonsGFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# rename chromosomes to match ce10 (UCSC)
exonsGFF[,1] <- sub("CHROMOSOME_","chr",exonsGFF[,1])
exonsGFF[,1] <- sub("chrMtDNA","chrM",exonsGFF[,1])
# keep only miRNA_mature_transcript and remove curated_miRNA
exonsGFF <- exonsGFF[exonsGFF[,2] != "curated_miRNA",]
# remove the suffix "_mature_transcript" from the annotation keywords
exonsGFF[,2] <- sub("_mature_transcript","",exonsGFF[,2])

# fix attributes column
exonsGFF[,9] <- sub("Transcript", "transcript_id", exonsGFF[,9])
exonsGFF[,9] <- sub(" ", " \"", exonsGFF[,9])
exonsGFF[,9] <- paste(exonsGFF[,9], "\";", sep = "")

# transcript biotype
exonsGFF[,2] <- paste("transcript_biotype \"", exonsGFF[,2], "\"", sep = "")

# fix attributes
exonsGFF[,9] <- paste(exonsGFF[,9], exonsGFF[,2], sep = " ")

# fix trabscript source
exonsGFF[,2] <- "c_elegans.WS220"

# write output file
write.table(exonsGFF, file=output_gtf, quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
