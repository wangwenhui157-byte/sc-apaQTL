######################### input prepare ######################### 
#!/usr/bin/Rscript

suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("yaml"))
suppressMessages(library("stringr"))
suppressMessages(library("optparse"))
suppressMessages(library("R.utils"))

VERSION = "V1.1.1 2024-12-23"

# Function
if(TRUE){
  checkPath = function(path, mode="file", exists=FALSE){
    path = normalizePath(path, winslash="/", mustWork=exists)
    if(mode=="file"){
      dirPath = dirname(path)
    }else if(mode=="dir"){
      dirPath = path
    }
    if(!dir.exists(dirPath) & !exists){dir.create(dirPath, recursive=TRUE)}
    
    return(path)
  }
}

# Paraments
if(TRUE){
  parser = list()
  parser[["filter"]] = make_option("--filter", dest="filter", type="integer", default=0, help="")
  parser[["type"]] = make_option("--type", dest="type", type="character", default="sigSNP", help="")
  parser[["del"]] = make_option("--del", dest="del", type="character", default="", help="")
  parser[["config"]] = make_option("--config", dest="config", type="character", default="config.yml", help="")
  parser[["qtl"]] = make_option("--qtl", dest="qtl", type="character", default=".", help="")
  parser[["sigSNP"]] = make_option("--sigSNP", dest="sigSNP", type="character", default=".", help="")
  parser[["check"]] = make_option("--check", dest="check", type="character", default=".", help="")
  parser[["out"]] = make_option("--out", dest="out", type="character", default=".", help="")
  
  args = parse_args(OptionParser(option_list=parser))
  FILTERMARKER = args$filter
  TYPE = args$type
  DELCELLTYPE = args$del
  CONFIG = args$config
  QTLPATH = args$qtl
  SIGSNP = args$sigSNP
  CHECKPATH = args$check
  OUTPATH = args$out
}

# Process paraments
if(TRUE){
  FILTERMARKER = ifelse(FILTERMARKER==1, TRUE, FALSE)
  DELCELLTYPE = unlist(strsplit(DELCELLTYPE, ','))
  setDTthreads(30)
  CONFIG = checkPath(CONFIG, mode="file", exists=TRUE)
  QTLPATH = checkPath(QTLPATH, mode="dir", exists=TRUE)
  SIGSNP = checkPath(SIGSNP, mode="file", exists=TRUE)
  CHECKPATH = checkPath(CHECKPATH, mode="dir", exists=TRUE)
  OUTPATH = checkPath(OUTPATH, mode="dir", exists=FALSE)
}

# Print paraments
if(TRUE){
  cat('\n', rep('=', 47), "CONFIG", rep('=', 47), '\n', sep='')
  cat("[Version]", VERSION, '\n', sep='')
  cat("[Date]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n', sep='')
  cat("[Parament]filter: ", FILTERMARKER, '\n', sep='')
  cat("[Parament]type: ", TYPE, '\n', sep='')
  cat("[Parament]del: ", DELCELLTYPE, '\n', sep='')
  cat("[Parament]config: ", CONFIG, '\n', sep='')
  cat("[Parament]qtl: ", QTLPATH, '\n', sep='')
  cat("[Parament]sigSNP: ", SIGSNP, '\n', sep='')
  cat("[Parament]check: ", CHECKPATH, '\n', sep='')
  cat("[Parament]out: ", OUTPATH, '\n', sep='')
  cat(rep('-', 48), "LOG", rep('-', 49), '\n', sep='')
}

# main
if(TRUE){
  # Loading
  if(TRUE){
    cat("[Progress]Loading file...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    cat("\t[Info]Loading file:", CONFIG, '\n')
    config = yaml::read_yaml(CONFIG)
    cellTypeV = config[["CELLTYPEV"]]
    # Read sigSNP data and filter pairs that are significant in at least one cell type
    cat("\t[Info]Loading file:", SIGSNP, '\n')
    dtSig = fread(SIGSNP, sep='\t', quote=FALSE) %>%
      .[!cellType %in% DELCELLTYPE]
    #sigTag = dtSig[, APATag] %>%
    #  unique()
    #sigSNP = dtSig[, snps] %>%
    #  unique()
    sigPairs = dtSig[, pairs] %>%
      unique()
    cat(sprintf("\t[Info]There are %s pairs significant in at least one cell type.", length(sigPairs)), '\n')
    # For QTL data: Read Matrix eQTL results and filter pairs significant in ≥1 cell type
    if(FALSE){
      qtl = list.files(QTLPATH, pattern="_cis.tsv") %>%
        setNames(., stringr::str_remove(., pattern="_cis.tsv")) %>%
        .[names(.) %in% cellTypeV] %>%
        .[!names(.) %in% DELCELLTYPE] %>%
        lapply(function(filename){
          filepath = file.path(QTLPATH, filename)
          cat("\t[Info]Loading file:", filepath, '\n')
          fread(filepath, sep='\t', quote=FALSE, select=c(1,2,3,5)) %>%
            .[, .(gene, SNP, beta, `p-value`)] %>%
            .[gene %in% sigTag] %>%
            .[SNP %in% sigSNP] %>%
            .[, pairs:=paste0(gene, '@', SNP)] %>%
            .[pairs %in% sigPairs]
        })
    }
    if(FILTERMARKER){
      qtl = list.files(QTLPATH, pattern="_cis.tsv") %>%
        setNames(., stringr::str_remove(., pattern="_cis.tsv")) %>%
        .[names(.) %in% cellTypeV] %>%
        .[!names(.) %in% DELCELLTYPE] %>%
        lapply(function(filename){
          filepath =file.path(QTLPATH, filename)
          cat("\t[Info]Loading file:", filepath, '\n')
          fread(filepath, sep='\t', quote=FALSE, select=c("gene", "SNP", "pairs", "beta", "se", "p-value")) %>%
            .[pairs %in% sigPairs] %>%
            setkey(pairs)
        })
    }
    if(!FILTERMARKER){
      qtl = list.files(QTLPATH, pattern="_cis.tsv") %>%
        setNames(., stringr::str_remove(., pattern="_cis.tsv")) %>%
        .[names(.) %in% cellTypeV] %>%
        .[!names(.) %in% DELCELLTYPE] %>%
        lapply(function(filename){
          filepath = file.path(QTLPATH, filename)
          cat("\t[Info]Loading file:", filepath, '\n')
          fread(filepath, sep='\t', quote=FALSE, select=c("gene", "SNP", "pairs", "beta", "se", "p-value")) %>%
            setkey(pairs)
        })
    }
    
    filepath = file.path(CHECKPATH, sprintf("mash_qtl_1_%s_%s.rds", ifelse(FILTERMARKER, "filter", "nofilter"), TYPE))
    cat(sprintf("\t[Info]Saving check file: %s", filepath), '\n')
    saveRDS(qtl, filepath, compress=FALSE)
  }
  # Process
  if(TRUE){
    cat("[Progress]Formatting data...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    # Select intersection of pairs present in every cell type  
    cat("\t[Info]Filtering APATag...\n")
    intersectPairs = qtl %>%
      lapply(function(dtCellType){
        cat('+')
        dtCellType[, pairs]
      }) %>%
      Reduce(intersect, .) %>%
      sort()
    cat('\n')
    cat(sprintf("\t[Info]There are %s pairs shared in %s cellTypes.", length(intersectPairs), length(qtl)), '\n')
    # Process Matrix eQTL results, retain only pairs common to all cell type datasets
    cat("\t[Info]Filtering summary table...\n")
    qtl = qtl %>% 
      lapply(function(dtCellType){
        cat('+')
        dtCellType[pairs %in% intersectPairs] %>%
          setorder(pairs)
      })
    cat('\n')
    # Saving file
    filepath = file.path(CHECKPATH, sprintf("mash_qtl_2_%s_%s.rds", ifelse(FILTERMARKER, "filter", "nofilter"), TYPE))
    cat(sprintf("\t[Info]Saving check file: %s", filepath), '\n')
    saveRDS(qtl, filepath, compress=FALSE)
    # del pairs
    qtl = purrr::map(qtl, ~ .x[, pairs:=NULL])
    # Change APATag and SNP IDs to fix bugs
    if(TRUE){
      idMap = list("APATag"=qtl[[1]][["gene"]] %>% unique() %>% data.table("APATag"=.),
                   "snps"=qtl[[1]][["SNP"]] %>% unique() %>% data.table("snps"=.))
      idMap[["APATag"]][["APATagMapped"]] = paste0("tag", seq_len(nrow(idMap[["APATag"]])))
      idMap[["snps"]][["snpsMapped"]] = paste0("snp", seq_len(nrow(idMap[["snps"]])))
      
      qtl = qtl %>%
        lapply(function(dtCellType){
          cat('+')
          dtCellType = dtCellType %>%
            merge.data.table(idMap[["snps"]], by.x="SNP", by.y="snps") %>%
            setkey(NULL) %>%
            merge.data.table(idMap[["APATag"]], by.x="gene", by.y="APATag") %>%
            setkey(NULL) %>%
            .[, gene:=NULL] %>%
            .[, SNP:=NULL] %>%
            setnames(c("APATagMapped", "snpsMapped"), c("gene", "SNP")) %>%
            setorder(gene, SNP)
        })
      cat('\n')
      # Saving file
      filepath = file.path(OUTPATH, sprintf("mash_idMap_%s_%s.rds", ifelse(FILTERMARKER, "filter", "nofilter"), TYPE))
      cat("\t[Info]Saving file:", filepath, '\n')
      saveRDS(idMap, filepath)
      filepath = file.path(CHECKPATH, sprintf("mash_qtl_3_%s_%s.rds", ifelse(FILTERMARKER, "filter", "nofilter"), TYPE))
      cat(sprintf("\t[Info]Saving check file: %s", filepath), '\n')
      saveRDS(qtl, filepath)
    }
    #qtl = qtl %>%
    #  purrr::map(~ .x[, gene:=stringr::str_replace_all(gene, pattern="\\.", replacement="_")]) %>%  # APATag中的.替换为_
    #  purrr::map(~ .x[, SNP:=stringr::str_replace_all(SNP, pattern=":", replacement="_")])
    #qtl = qtl %>%
    #  purrr::map(~ .x[, se:=sqrt(((beta)^2)/qchisq(`p-value`, 1,lower.tail=F))])
    qtl = qtl %>%
      purrr::map(~ setnames(.x,
                            old=c("gene", "SNP", "beta", "p-value", "se"),
                            new=c("gene_id", "variant_id", "slope", "pval_nominal", "slope_se"))) %>%
      purrr::map(~ .x[, .(gene_id, variant_id, slope, pval_nominal, slope_se)])
    
    allAPATag = qtl[[1]][, gene_id] %>%
      unique() %>%
      as.data.table()
  }
  # Saving
  if(TRUE){
    cat("[Progress]Saving file...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')

    outAPATag = file.path(OUTPATH, sprintf("apatag.tsv"))
    cat("\t[Info]Saving file:", outAPATag, '\n')
    fwrite(allAPATag, outAPATag, sep='\t', quote=FALSE, col.names=FALSE)

    sampleFile = c()
    for(i in names(qtl)){
      outfile = file.path(OUTPATH, sprintf("%s.tsv", i))
      cat("\t[Info]Saving file:", outfile, '\n')
      fwrite(qtl[[i]], outfile, sep='\t', quote=FALSE, col.names=TRUE)
      R.utils::gzip(outfile)
      sampleFile = c(sampleFile, sprintf("%s.tsv.gz", i))
    }

    sampleFile = as.data.table(sampleFile)
    outfile = file.path(OUTPATH, sprintf("sample.tsv"))
    cat("\t[Info]Saving file:", outfile, '\n')
    fwrite(sampleFile, outfile, sep='\t', quote=FALSE, col.names=FALSE)
  }
  cat(rep('-', 30), paste0('[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ']', "All parts finished"), rep('-', 30), '\n', sep='')
}


######################### run ######################### 
#Linux
#transform summary data to mash format
#!/bin/bash
sos run fastqtl_to_mash --data-list sample.tsv --gene-list apatag.tsv -j 20 --cols 3 5 4

#/bin/bash
sos run mashr_flashr_workflow.ipynb mash -j 10 --data fastqtl_to_mash_output/sample.mash.rds -v3


######################### plot ######################### 
#magnitude
out <- readRDS("fastqtl_to_mash_output/sample.mash.rds")
maxb <- out$strong.b
maxz <- out$strong.z
out <- readRDS("mashr_flashr_workflow_output/sample.mash.EZ.posterior.rds")
pm.mash <- out$PosteriorMean
lfsr.all <- out$lfsr
standard.error <- maxb/maxz
pm.mash.beta <- pm.mash*standard.error

thresh <- 0.05
pm.mash.beta <- pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
lfsr.mash <- lfsr.all[rowSums(lfsr.all<0.05)>0,]
shared.fold.size <- matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size) <- rownames(shared.fold.size) <- colnames(maxz)
for (i in 1:ncol(lfsr.mash)){
  for (j in 1:ncol(lfsr.mash)) {
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    shared.fold.size[i,j] = mean(quotient > 0.5 & quotient < 2, na.rm = TRUE)
  }
}

rownames(shared.fold.size) <- sub("\\.tsv$", "", rownames(shared.fold.size))
colnames(shared.fold.size) <- sub("\\.tsv$", "", colnames(shared.fold.size))
lat <- shared.fold.size

celltypes_order <- c("Mono_NC", "Mono_C", "DC", "B_Mem", "B_IN", "NK_R", "CD4_CTL", "MAIT",
                     "gdT", "T_Reg", "CD8_TCM", "CD4_TEM", "CD8_Naive", "NK", "CD8_TEM",
                     "CD4_Naive", "CD4_TCM")
celltypes_order <- rev(celltypes_order)
lat <- lat[celltypes_order, celltypes_order]
lat[lower.tri(lat)] <- NA

clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(100)

celltype_colors <- celltype_colors[celltypes_order]

annotation_col <- data.frame(cellType = factor(colnames(lat), levels = celltypes_order))
rownames(annotation_col) <- colnames(lat)

annotation_colors_col <- list()
annotation_colors_col[[1]] <- celltype_colors[celltypes_order]
names(annotation_colors_col) <- NULL

p <- pheatmap(lat,
              color = clrs,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              display_numbers = FALSE,
              fontsize_row = 12,
              fontsize_col = 10,
              angle_col = 45,
              na_col = "white",
              border_color = "white",
              cellwidth = 16,
              cellheight = 16,
              labels_col = rep("", ncol(lat)),
              annotation_row = annotation_col,
              annotation_col = annotation_col,
              annotation_colors = list(cellType = celltype_colors),
              annotation_names_row = FALSE,
              annotation_names_col = FALSE,
              annotation_legend = FALSE
             )

#sign
out <- readRDS("fastqtl_to_mash_output/sample.mash.rds")
maxb <- out$strong.b
maxz <- out$strong.z
out <- readRDS("mashr_flashr_workflow_output/sample.mash.EZ.posterior.rds")
pm.mash <- out$PosteriorMean
lfsr.all <- out$lfsr
standard.error <- maxb/maxz
pm.mash.beta <- pm.mash*standard.error

thresh=0.05
pm.mash.beta=pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]
shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(maxz)
for(i in 1:ncol(lfsr.mash)){
  for(j in 1:ncol(lfsr.mash)){
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    shared.fold.size[i,j] = mean(quotient > 0, na.rm = TRUE)
  }
}

rownames(shared.fold.size) <- sub("\\.tsv$", "", rownames(shared.fold.size))
colnames(shared.fold.size) <- sub("\\.tsv$", "", colnames(shared.fold.size))
lat <- shared.fold.size

celltypes_order <- c("Mono_NC", "Mono_C", "DC", "B_Mem", "B_IN", "NK_R", "CD4_CTL", "MAIT",
                      "gdT", "T_Reg", "CD8_TCM", "CD4_TEM", "CD8_Naive", "NK", "CD8_TEM",
                      "CD4_Naive", "CD4_TCM")
celltypes_order <- rev(celltypes_order)
lat <- lat[celltypes_order, celltypes_order]
lat[lower.tri(lat)] <- NA

clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(100)

celltype_colors <- celltype_colors[celltypes_order]

annotation_col <- data.frame(cellType = factor(colnames(lat), levels = celltypes_order))
rownames(annotation_col) <- colnames(lat)

annotation_colors_col <- list()
annotation_colors_col[[1]] <- celltype_colors[celltypes_order]
names(annotation_colors_col) <- NULL


p <- pheatmap(lat,
              color = clrs,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              display_numbers = FALSE,
              fontsize_row = 12,
              fontsize_col = 10,
              angle_col = 45,
              na_col = "white",
              border_color = "white",
              cellwidth = 16,
              cellheight = 16,
              labels_col = rep("", ncol(lat)),
              annotation_row = annotation_col,
              annotation_col = annotation_col,
              annotation_colors = list(cellType = celltype_colors),
              annotation_names_row = FALSE,
              annotation_names_col = FALSE,
              annotation_legend = FALSE
              )
