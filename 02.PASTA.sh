# polyApipe: get polyA
BAM="bam/"
OUTPUT="polyApipe_result/"
Sarray=("$@")
num=${#Sarray[@]}
for ((i=0;i<num;i++)); do
    time polyApipe.py \
             -i ${BAM}/${Sarray[i]}_merged_sorted.bam \
             -o ${OUTPUT}/${Sarray[i]} \
             --depth_threshold 5 \
             --region_size 300 \
             --minpolyA 5 \
             -t 6 \
             --no_peaks
done


# polyApipe: peak gff
OUTPUT="polyApipe_result/"
time polyApipe.py \
         -i ${OUTPUT} \
         -o ${OUTPUT}/GSM \
         --polyA_bams \
         --depth_threshold 5 \
         --region_size 300 \
         --minpolyA 5 \
         -t 7 \
         --no_count


# polyApipe: count
BAM="bam/"
OUTPUT="polyApipe_result/"
Sarray=("$@")
num=${#Sarray[@]}
for ((i=0;i<num;i++)); do
        # count_in_polyA: {OUT_ROOT}_counts/{sample}_counts.tab.gz
        time polyApipe.py \
             -i ${BAM}/${Sarray[i]}_merged_sorted.bam \
             -o ${OUTPUT}/${Sarray[i]} \
             --peaks_gff ${OUTPUT}/GSM_polyA_peaks_chr.gff \
             --depth_threshold 5 \
             --region_size 300 \
             --minpolyA 5 \
             -t 6
done


# PASTA: filter polyA site
library(PASTA)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(dplyr)

dir = "meta/"
setwd(dir)

counts.file <- "polyApipe_result/impute_polyA_site/1006_1007_counts_modified.tab.gz"
peak.file <- "/polyApipe_result/GSM_polyA_peaks_chr.gff"
polyAdb.file <- "human_PAS_hg38.txt"
peak_short_name_file = read.table("polyApipe_result/GSM_polyA_peaks_chr.tsv",sep='\t')
colnames(peak_short_name_file) <- c('peak','peak_short_name','strand')

options(Seurat.object.assay.calcn = TRUE)
polyA.counts = ReadPolyApipe(counts.file = counts.file, peaks.file = peak.file, filter.chromosomes = TRUE)

polyA.assay = CreatePolyAAssay(counts = polyA.counts, genome = "hg38", validate.fragments = FALSE)
# add annotations to polyA.assay
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(polyA.assay) <- annotations
save(annotations, file="polyA_assay_annotations.Rdata")

seurat_obj <- CreateSeuratObject(counts = polyA.assay,
                                 assay = "polyA")
seurat_obj <- GetPolyADbAnnotation(seurat_obj, assay = "polyA", polyAdb.file = polyAdb.file, max.dist = 50)

meta <- seurat_obj[["polyA"]][[]] %>%
  filter(
    Intron.exon_location == "3'_most_exon",
    !is.na(Gene_Symbol)
  ) %>%
  group_by(Gene_Symbol) %>%
  filter(n() > 1) %>%
  ungroup()

meta=as.data.frame(meta)
meta <- merge(meta, peak_short_name_file[, c("peak", "peak_short_name", "strand")],
                      by = c("peak", "strand"),
                      all.x = TRUE)
write.table(meta,file="metadata_tandem.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

meta <- seurat_obj[["polyA"]][[]] %>%
  filter(
    Intron.exon_location %in% c("3'_most_exon",'Intron'),
    !is.na(Gene_Symbol)
  ) %>%
  group_by(Gene_Symbol) %>%
  filter(n() > 1) %>%
  ungroup()

meta=as.data.frame(meta)
meta <- merge(meta, peak_short_name_file[, c("peak", "peak_short_name", "strand")],
                      by = c("peak", "strand"),
                      all.x = TRUE)
write.table(meta,file="metadata_intronic_and_tandem.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

meta=seurat_obj[['polyA']][[]]
meta <- merge(meta, peak_short_name_file[, c("peak", "peak_short_name", "strand")],
                      by = c("peak", "strand"),
                      all.x = TRUE)
write.table(meta,file="metadata_all_polyA_site.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

timeend <- Sys.time()
runningtime <- timeend-timestart
print(runningtime)


# PASTA: background (tandem)
library(PASTA)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(data.table)
library(MGLM)
library(gplm)
library(stats)

counts.file = "polyApipe_result/tandem_merged_filtered_counts_modified.background.tab"
output_dir = "PASTA_result/background/"
setwd(output_dir)

peak.file <- "polyApipe_result/GSM_polyA_peaks_chr.gff"
polyAdb.file <- "human_PAS_hg38.txt"
outputfile <- "tandem_background.Rdata"

metadata <- read.table("barcode_celltype.tsv",header=T,row.names="barcode",sep="\t")

polyA_meta <- fread(metadata_tandem.txt")
polyA_meta <- as.data.frame(polyA_meta)
rownames(polyA_meta) <- polyA_meta$peak
load("polyA_assay_annotations.Rdata")

options(Seurat.object.assay.calcn = TRUE)
polyA.counts = ReadPolyApipe(counts.file = counts.file, peaks.file = peak.file, filter.chromosomes = TRUE)
polyA.assay = CreatePolyAAssay(counts = polyA.counts, genome = "hg38", validate.fragments = FALSE)

# add annotations to polyA.assay
Annotation(polyA.assay) <- annotations

seurat_obj <- CreateSeuratObject(counts = polyA.assay,
                                 assay = "polyA",
                                 meta.data = metadata)

seurat_obj[['polyA']] <- AddMetaData(seurat_obj[['polyA']],polyA_meta)
features <- rownames(seurat_obj[['polyA']])
Idents(seurat_obj) <- "cellTypePaper"

message("Calculating PolyA Residuals...")
source("fit.background.R")
seurat_obj <- CalcPolyAResiduals(seurat_obj,
                                assay = "polyA",
                                features = features,
                                gene.names = "Gene_Symbol",
                                sample.n = 10000,
                                verbose = TRUE)


# PASTA: calculate polyA residuals (tandem)
library(PASTA)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(data.table)
library(MGLM)
library(gplm)
library(stats)

args<-commandArgs(trailingOnly = TRUE)
counts.file = args[1]
output_dir = args[2]

setwd(output_dir)
peak.file <- "polyApipe_result/GSM_polyA_peaks_chr.gff"
polyAdb.file <- "human_PAS_hg38.txt"
gene <- gsub("^[^.]+\\.([^.]+)\\..*$", "\\1", basename(counts.file))
outputfile <- paste0("seurat_obj_",gene,".Rdata")

metadata <- read.table("barcode_celltype_subset.tsv",header=T,row.names="barcode",sep="\t")
polyA_meta <- fread("metadata_tandem.txt")
polyA_meta <- as.data.frame(polyA_meta)
rownames(polyA_meta) <- polyA_meta$peak
load("polyA_assay_annotations.Rdata")

options(Seurat.object.assay.calcn = TRUE)
polyA.counts = ReadPolyApipe(counts.file = counts.file, peaks.file = peak.file, filter.chromosomes = TRUE)
polyA.assay = CreatePolyAAssay(counts = polyA.counts, genome = "hg38", validate.fragments = FALSE)

# add annotations to polyA.assay
Annotation(polyA.assay) <- annotations

seurat_obj <- CreateSeuratObject(counts = polyA.assay,
                                 assay = "polyA",
                                 meta.data = metadata)

seurat_obj[['polyA']] <- AddMetaData(seurat_obj[['polyA']],polyA_meta)
features <- rownames(seurat_obj[['polyA']])
Idents(seurat_obj) <- "cellTypePaper"

message("Calculating PolyA Residuals...")
source("fit_revise.R")
seurat_obj <- CalcPolyAResiduals(seurat_obj,
                                assay = "polyA",
                                features = features,
                                gene.names = "Gene_Symbol",
                                sample.n = 10000,
                                verbose = TRUE,
                                grid.var = "PASTA_result/background/tandem_background.Rdata")

save(seurat_obj, file = outputfile)

exp = as.matrix(seurat_obj[['polyA']]$scale.data)
outputfile = paste0("polyA_residuals_",gene,".csv")
write.csv(exp, file = outputfile, quote=FALSE)
