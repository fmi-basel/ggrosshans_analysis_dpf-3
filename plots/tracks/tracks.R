setwd("/tungstenfs/scratch/ggrossha/gypafoiv/projects/GROUP_PROJECTS/rajani/small_RNAs_Rajani/plots/tracks")

library(Gviz)
library(rtracklayer)
library(GenomicFeatures)

options(ucscChromosomeNames=FALSE)

# annotation
gtf_path <- "/tungstenfs/scratch/ggrossha/GROUP/annotation/ce11/wormbase/WS270/c_elegans.PRJNA13758.WS270.canonical_geneset.sorted.gtf"

# repeats
repeats_path <- "../../00_annotation/results/annotation/ce_11_repeats.filtered.bed"

# small RNAs

WT_A_plus_path <- "../../01_small_RNA_seq_15_C/results/wild_type_A/bowtie/alignment.sorted.filtered.unique_mappers.plus.bw"
WT_A_minus_path<- "../../01_small_RNA_seq_15_C/results/wild_type_A/bowtie/alignment.sorted.filtered.unique_mappers.minus.bw"
Dpf_3_null_A_plus_path <- "../../01_small_RNA_seq_15_C/results/Dpf_3_null_A/bowtie/alignment.sorted.filtered.unique_mappers.plus.bw"
Dpf_3_null_A_minus_path <- "../../01_small_RNA_seq_15_C/results/Dpf_3_null_A/bowtie/alignment.sorted.filtered.unique_mappers.minus.bw"
mut_2_A_plus_path <- "../../01_small_RNA_seq_15_C/results/mut_2_A/bowtie/alignment.sorted.filtered.unique_mappers.plus.bw"
mut_2_A_minus_path <- "../../01_small_RNA_seq_15_C/results/mut_2_A/bowtie/alignment.sorted.filtered.unique_mappers.minus.bw"
mut_7_A_plus_path <- "../../01_small_RNA_seq_15_C/results/mut_7_A/bowtie/alignment.sorted.filtered.unique_mappers.plus.bw"
mut_7_A_minus_path <- "../../01_small_RNA_seq_15_C/results/mut_7_A/bowtie/alignment.sorted.filtered.unique_mappers.minus.bw"
# total RNAs
total_rna_WT_A_str1_path <- "../../02_total_RNAseq_15_C/results/wild_type_A/STAR/wild_type_A_Signal.Unique.str1.out.bw"
total_rna_WT_A_str2_path <- "../../02_total_RNAseq_15_C/results/wild_type_A/STAR/wild_type_A_Signal.Unique.str2.out.bw"
total_rna_Dpf_3_null_A_str1_path <- "../../02_total_RNAseq_15_C/results/Dpf_3_null_A/STAR/Dpf_3_null_A_Signal.Unique.str1.out.bw"
total_rna_Dpf_3_null_A_str2_path <- "../../02_total_RNAseq_15_C/results/Dpf_3_null_A/STAR/Dpf_3_null_A_Signal.Unique.str2.out.bw"
total_rna_mut_2_A_str1_path <- "../../02_total_RNAseq_15_C/results/mut_2_A/STAR/mut_2_A_Signal.Unique.str1.out.bw"
total_rna_mut_2_A_str2_path <- "../../02_total_RNAseq_15_C/results/mut_2_A/STAR/mut_2_A_Signal.Unique.str2.out.bw"
total_rna_mut_7_A_str1_path <- "../../02_total_RNAseq_15_C/results/mut_7_A/STAR/mut_7_A_Signal.Unique.str1.out.bw"
total_rna_mut_7_A_str2_path <- "../../02_total_RNAseq_15_C/results/mut_7_A/STAR/mut_7_A_Signal.Unique.str2.out.bw"

# read objects
repeats <- AnnotationTrack(repeats_path, name="repeats")

txdb <- makeTxDbFromGFF(gtf_path)
geneTrack <- GeneRegionTrack(txdb, name="genes")

WT_A_plus <- DataTrack(WT_A_plus_path, name="WT A +")
WT_A_minus <- DataTrack(WT_A_minus_path, name="WT A -", ylim=c(0,40))
Dpf_3_null_A_plus <- DataTrack(Dpf_3_null_A_plus_path, name="Dpf-3 null A +")
Dpf_3_null_A_minus <- DataTrack(Dpf_3_null_A_minus_path, name="Dpf-3 null A -", ylim=c(0,40))
mut_2_A_plus <- DataTrack(mut_2_A_plus_path, name="mut-2 A +")
mut_2_A_minus <- DataTrack(mut_2_A_minus_path, name="mut-2 A -", ylim=c(0,40))
mut_7_A_plus <- DataTrack(mut_7_A_plus_path, name="mut-7 A +")
mut_7_A_minus <- DataTrack(mut_7_A_minus_path, name="mut-7 A -", ylim=c(0,40))

total_rna_WT_A_str1 <- DataTrack(total_rna_WT_A_str1_path, name = "WT A", ylim=c(0,1.2))
total_rna_WT_A_str2 <- DataTrack(total_rna_WT_A_str2_path, name = "WT A", ylim=c(0,1.2))
total_rna_Dpf_3_null_A_str1 <- DataTrack(total_rna_Dpf_3_null_A_str1_path, name = "Dpf-3 null A", ylim=c(0,1.2))
total_rna_Dpf_3_null_A_str2 <- DataTrack(total_rna_Dpf_3_null_A_str2_path, name = "Dpf-3 null A", ylim=c(0,1.2))
total_rna_mut_2_A_str1 <- DataTrack(total_rna_mut_2_A_str1_path, name = "mut-2 A", ylim=c(0,1.2))
total_rna_mut_2_A_str2 <- DataTrack(total_rna_mut_2_A_str2_path, name = "mut-2 A", ylim=c(0,1.2))
total_rna_mut_7_A_str1 <- DataTrack(total_rna_mut_7_A_str1_path, name = "mut-7 A", ylim=c(0,1.2))
total_rna_mut_7_A_str2 <- DataTrack(total_rna_mut_7_A_str2_path,name = "mut-7 A", ylim=c(0,1.2))


# ZC15.3
chromosome <- "V"
start <- 20303588
end <- 20304908

pdf("ZC15.3.pdf")
plotTracks(
  list(repeats, 
       geneTrack,
       WT_A_minus,
       Dpf_3_null_A_minus,
       mut_2_A_minus,
       mut_7_A_minus,
       total_rna_WT_A_str2,
       total_rna_Dpf_3_null_A_str2,
       total_rna_mut_2_A_str2,
       total_rna_mut_7_A_str2),
  from=start,
  to=end,
  chromosome=chromosome,
  type="histogram",
  col.histogram = "darkblue",
  fill.histogram = "darkblue",
  fontcolor = "black",
  shape = "box",
  extend.left=100,
  extend.right=100,
  col.axis="black",
  col=NULL,
  fontsize = 10,
  col.frame="white",
  background.title="white",
  col.title="black"
)
dev.off()


WT_A_plus <- DataTrack(WT_A_plus_path, name="WT A +")
WT_A_minus <- DataTrack(WT_A_minus_path, name="WT A -", ylim=c(0,150))
Dpf_3_null_A_plus <- DataTrack(Dpf_3_null_A_plus_path, name="Dpf-3 null A +")
Dpf_3_null_A_minus <- DataTrack(Dpf_3_null_A_minus_path, name="Dpf-3 null A -", ylim=c(0,150))
mut_2_A_plus <- DataTrack(mut_2_A_plus_path, name="mut-2 A +")
mut_2_A_minus <- DataTrack(mut_2_A_minus_path, name="mut-2 A -", ylim=c(0,150))
mut_7_A_plus <- DataTrack(mut_7_A_plus_path, name="mut-7 A +")
mut_7_A_minus <- DataTrack(mut_7_A_minus_path, name="mut-7 A -", ylim=c(0,150))

total_rna_WT_A_str1 <- DataTrack(total_rna_WT_A_str1_path, name = "WT A", ylim=c(0,2))
total_rna_WT_A_str2 <- DataTrack(total_rna_WT_A_str2_path, name = "WT A", ylim=c(0,2))
total_rna_Dpf_3_null_A_str1 <- DataTrack(total_rna_Dpf_3_null_A_str1_path, name = "Dpf-3 null A", ylim=c(0,2))
total_rna_Dpf_3_null_A_str2 <- DataTrack(total_rna_Dpf_3_null_A_str2_path, name = "Dpf-3 null A", ylim=c(0,2))
total_rna_mut_2_A_str1 <- DataTrack(total_rna_mut_2_A_str1_path, name = "mut-2 A", ylim=c(0,2))
total_rna_mut_2_A_str2 <- DataTrack(total_rna_mut_2_A_str2_path, name = "mut-2 A", ylim=c(0,2))
total_rna_mut_7_A_str1 <- DataTrack(total_rna_mut_7_A_str1_path, name = "mut-7 A", ylim=c(0,2))
total_rna_mut_7_A_str2 <- DataTrack(total_rna_mut_7_A_str2_path,name = "mut-7 A", ylim=c(0,2))



# c38d9.2: 
chromosome <- "V"
start <- 17567513
end <- 17571742

pdf("C38D9.2.pdf")
plotTracks(
  list(repeats, 
       geneTrack,
       WT_A_minus,
       Dpf_3_null_A_minus,
       mut_2_A_minus,
       mut_7_A_minus,
       total_rna_WT_A_str2,
       total_rna_Dpf_3_null_A_str2,
       total_rna_mut_2_A_str2,
       total_rna_mut_7_A_str2),
  from=start,
  to=end,
  chromosome=chromosome,
  type="histogram",
  col.histogram = "darkblue",
  fill.histogram = "darkblue",
  fontcolor = "black",
  shape = "box",
  extend.left=200,
  extend.right=200,
  col.axis="black",
  col=NULL,
  fontsize = 10,
  col.frame="white",
  background.title="white",
  col.title="black"
)
dev.off()

