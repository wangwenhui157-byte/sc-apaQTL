library(data.table)
library(dplyr)
library(moloc)
library(purrr)

diseaseV <- c("MS.GCST003566", "PSO.GCST90243956", "SLE.GCST90011866", "T1DM.GCST90014023",
              "CD.ibdgc", "IBD.ibdgc", "RA.jenger", "Celiac.phecode557")
cellTypeV <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD8_Naive", "CD8_TEM",
               "T_Reg", "NK", "B_Mem", "B_IN", "Mono_C", "Mono_NC")

gwas_pt <- paste0(dis, ".ma")
gwas <- fread(gwas_pt)

individualQTL_pt <- "individuals2.rds"
individualSummary = readRDS(individualQTL_pt) %>%
  purrr::map(~ length(.x)) %>%
  as.data.table() %>%
  data.table::transpose(keep.names="cellType") %>%
  setnames("V1", "num") %>%
  .[cellType %in% cellTypeV] %>%
  .[, .(cellType, num)]

pqtl_info_pt <- "00.pQTL.info.sub.new.list"
pqtl_info <- fread(pqtl_info_pt)
pqtl_genes <- unique(pqtl_info$gene_name)

sig_pt <- "COLOC.SMR.utr.txt"
sig <- fread(sig_pt)
sig <- sig[disease %in% diseaseV & cellType %in% cellTypeV]
sig <- subset(sig, COLOC == "sig.coloc")
sig <- sig[Gene_Symbol %in% pqtl_genes]
head(sig)

sig <- subset(sig, disease == dis)
celltypes <- unique(sig$cellType)

out_file <- paste0(dis, ".moloc.res.txt")
if (file.exists(out_file)) file.remove(out_file)

for (ct in celltypes){
  message("Processing Disease: ", dis, "; Celltype: ", ct)
  N.apaqtl <- individualSummary[cellType==ct, num]

  sig_sub <- sig[cellType==ct, ]
  apas <- unique(sig_sub$APATag)

  apaqtl_pt <- paste0("QTL/", ct, "_cis.tsv")
  apaqtl <- fread(apaqtl_pt)
  esi_pt <- paste0("fromMatrix/", ct, ".esi")
  esi <- fread(esi_pt)
  colnames(esi) <- c("chr", "SNP", "probe", "pos", "A1", "A2", "maf")
  esi <- esi[, c("SNP", "A1", "A2", "maf")]
  apaqtl <- merge(apaqtl, esi, by = "SNP")

  apaqtl <- apaqtl[gene %in% apas, ]
  apasnps <- unique(apaqtl$SNP)
  gwas_ct <- gwas[SNP %in% apasnps, ]

  common_genes <- unique(sig_sub$Gene_Symbol)
  for (gene in common_genes){
    message("Processing Disease: ", dis, "; Celltype: ", ct, "; Gene: ", gene)

    pqtl_file <- unique(pqtl_info[pqtl_info$gene_name == gene, ]$file)
    pqtl_pt <- paste0("raw_data/", pqtl_file)
    pqtl <- fread(pqtl_pt)
    pqtl <- pqtl[rsids %in% apasnps, ]

    apas <- unique(sig_sub[sig_sub$Gene_Symbol == gene, ]$APATag)
    for (apa in apas){
      message("Processing Disease: ", dis, "; Celltype: ", ct, "; Gene: ", gene, "APATag: ", apa)

      apaqtl_moloc <- apaqtl[gene==apa, ]

      common_snp <- Reduce(intersect, list(apaqtl_moloc$SNP, pqtl$rsids, gwas_ct$SNP))
      if (length(common_snp) == 0) {
        message("No common SNPs. For Disease: ", dis, "; Celltype: ", ct, "; Gene: ", gene, "; APATag: ", apa)
        next
      }

      apaqtl_moloc <- apaqtl_moloc[match(common_snp, apaqtl_moloc$SNP), ]
      pqtl_moloc <- pqtl[match(common_snp, pqtl$rsids), ]
      gwas_moloc <- gwas_ct[match(common_snp, gwas_ct$SNP), ]

      gwas_moloc <- gwas_moloc[, c("SNP", "A1", "A2", "b", "se", "p", "N", "freq")]
      colnames(gwas_moloc) <- c("SNP", "A1", "A2", "BETA", "SE", "PVAL", "N", "MAF")

      apaqtl_moloc[, N := N.apaqtl]
      apaqtl_moloc[, se := beta / `t-stat`]
      apaqtl_moloc <- apaqtl_moloc[, c("SNP", "A1", "A2", "beta", "se", "p-value", "N", "maf")]
      colnames(apaqtl_moloc) <- c("SNP", "A1", "A2", "BETA", "SE", "PVAL", "N", "MAF")

      pqtl_moloc <- pqtl_moloc[, c("rsids", "effectAllele", "otherAllele", "Beta", "SE", "Pval", "N", "ImpMAF")]
      colnames(pqtl_moloc) <- c("SNP", "A1", "A2", "BETA", "SE", "PVAL", "N", "MAF")

      # allel harmonise
      merged <- merge(gwas_moloc, apaqtl_moloc, by = "SNP", suffixes = c(".gwas", ".apaqtl"))
      cols_pqtl <- setdiff(colnames(pqtl_moloc), "SNP")
      setnames(pqtl_moloc, cols_pqtl, paste0(cols_pqtl, ".pqtl"))
      merged <- merge(merged, pqtl_moloc, by = "SNP")

      ok_idx <- merged$A1.apaqtl == merged$A1.gwas & merged$A2.apaqtl == merged$A2.gwas
      flip_idx <- merged$A1.apaqtl == merged$A2.gwas & merged$A2.apaqtl == merged$A1.gwas
      merged$BETA.apaqtl <- ifelse(flip_idx, -merged$BETA.apaqtl, merged$BETA.apaqtl)
      merged$A1.apaqtl <- ifelse(flip_idx, merged$A1.gwas, merged$A1.apaqtl)
      merged$A2.apaqtl <- ifelse(flip_idx, merged$A2.gwas, merged$A2.apaqtl)
      merged$BETA.apaqtl <- ifelse(!(ok_idx | flip_idx), NA, merged$BETA.apaqtl)

      ok_idx <- merged$A1.pqtl == merged$A1.gwas & merged$A2.pqtl == merged$A2.gwas
      flip_idx <- merged$A1.pqtl == merged$A2.gwas & merged$A2.pqtl == merged$A1.gwas
      merged$BETA.pqtl <- ifelse(flip_idx, -merged$BETA.pqtl, merged$BETA.pqtl)
      merged$A1.pqtl <- ifelse(flip_idx, merged$A1.gwas, merged$A1.pqtl)
      merged$A2.pqtl <- ifelse(flip_idx, merged$A2.gwas, merged$A2.pqtl)
      merged$BETA.pqtl <- ifelse(!(ok_idx | flip_idx), NA, merged$BETA.pqtl)

      setDT(merged)
      merged <- merged[!is.na(BETA.apaqtl) & !is.na(BETA.pqtl)]
      all(merged$A1.apaqtl == merged$A1.gwas & merged$A2.apaqtl == merged$A2.gwas)
      all(merged$A1.pqtl   == merged$A1.gwas & merged$A2.pqtl   == merged$A2.gwas)
      head(merged)

      print(paste(
        "Disease: ", dis,
        "; Celltype: ", ct,
        "; Gene:", gene,
        "; APATag: ", apa,
        "; Have common snps:", nrow(merged), "\n"
      ))

      if (nrow(merged) == 0) {
        message("No SNPs after harmonization. For Disease: ", dis, "; Celltype: ", ct, "; Gene: ", gene, "; APATag: ", apa)
        next
      }

      moloc_input <- list(
        apaqtl = data.frame(
          SNP  = merged$SNP,
          BETA = merged$BETA.apaqtl,
          varbeta = merged$SE.apaqtl^2,
          pvalues = merged$PVAL.apaqtl,
          SE   = merged$SE.apaqtl,
          MAF  = merged$MAF.apaqtl,
          N    = merged$N.apaqtl,
          type    = "apaqtl"
        ),
        pqtl = data.frame(
          SNP  = merged$SNP,
          BETA = merged$BETA.pqtl,
          varbeta = merged$SE.pqtl^2,
          pvalues = merged$PVAL.pqtl,
          SE   = merged$SE.pqtl,
          MAF  = merged$MAF.pqtl,
          N    = merged$N.pqtl,
          type    = "pqtl"
        ),
        gwas = data.frame(
          SNP  = merged$SNP,
          BETA = merged$BETA.gwas,
          varbeta = merged$SE.gwas^2,
          pvalues = merged$PVAL.gwas,
          SE   = merged$SE.gwas,
          MAF  = merged$MAF.gwas,
          N    = merged$N.gwas,
          type    = "gwas"
        )
      )

      tryCatch({
        moloc_res <- moloc_test(moloc_input, prior_var = "default", save.SNP.info = FALSE)
        moloc_res

        if (!is.null(moloc_res$best_snp)) {
          c_res0 <- c(t(moloc_res$best_snp))
          names(c_res0) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b",
                             "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c",
                             "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc",
                             "pp_abc", "best_SNP_abc")
          moloc_combine_res0 <- data.frame(celltype=ct,
                                           disease=dis,
                                           APATag=apa,
                                           nsnps=moloc_res$nsnps,
                                           t(c_res0))
          fwrite(moloc_combine_res0, out_file, sep="\t", append=TRUE)
        } else {
          message(paste("Moloc analysis did not return best_snp for ---",
                        "Disease:", dis,
                        "; Celltype:", ct,
                        "; APATag:", apa))
        }

      }, error = function(e) {
        message(paste("Moloc analysis failed for ---",
                      "Disease:", dis,
                      "; Celltype:", ct,
                      "; APATag:", apa))
        message("Error message:", e$message)
      })

      rm(merged, apaqtl_moloc, pqtl_moloc, gwas_moloc, moloc_input, res)
      gc()
    } # end APATag loop

    rm(pqtl)
    gc()
  } # end Gene loop

  rm(apaqtl)
  gc()
} # end celltype loop
