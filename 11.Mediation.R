library(survey)
library(gsmr2)
library(TwoSampleMR)
library(plinkbinr)
library(ieugwasr)
library(data.table)

diseaseV <- c("MS.GCST003566", "PSO.GCST90243956", "SLE.GCST90011866", "T1DM.GCST90014023",
              "CD.ibdgc", "IBD.ibdgc", "RA.jenger", "Celiac.phecode557")
cellTypeV <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD8_Naive", "CD8_TEM",
               "T_Reg", "NK", "B_Mem", "B_IN", "Mono_C", "Mono_NC")

pqtl_info_pt <- "00.pQTL.info.sub.new.list"
pqtl_info <- fread(pqtl_info_pt)
pqtl_genes <- unique(pqtl_info$gene_name)


######################### apaQTL-pQTL-MR ##########################
sig_pt <- "COLOC.SMR.utr.txt"
sig <- fread(sig_pt)
sig <- sig[disease %in% diseaseV & cellType %in% cellTypeV]
sig <- subset(sig, COLOC == "sig.coloc")
sig <- sig[Gene_Symbol %in% pqtl_genes]
head(sig)

celltypes <- unique(sig$cellType)

or_file <- "apaQTL.pQTL.MR/02.apaQTL.pQTL.MR.res"
if (file.exists(or_file)) {
  file.remove(or_file)
}

dat_file <- "apaQTL.pQTL.MR/02.apaQTL.pQTL.MR.dat"
if (file.exists(dat_file)) {
  file.remove(dat_file)
}

for (ct in celltypes){
  message("Processing Celltype: ", ct)

  sig_sub <- sig[cellType==ct, ]
  apas <- unique(sig_sub$APATag)

  apaqtl_pt <- paste0("QTL/", ct, "_cis.tsv")
  apaqtl <- fread(apaqtl_pt)
  esi_pt <- paste0("fromMatrix/", ct, ".esi")
  esi <- fread(esi_pt)
  colnames(esi) <- c("chr", "SNP", "probe", "pos", "A1", "A2", "maf")
  esi <- esi[, c("SNP", "chr", "A1", "A2", "maf")]
  apaqtl <- merge(apaqtl, esi, by = "SNP")
  apaqtl[, se := abs(beta / `t-stat`)]
  apaqtl[, F_value := (beta^2) / (se^2)]
  apaqtl <- apaqtl[apaqtl$FDR < 0.05, ]
  apaqtl <- apaqtl[apaqtl$F_value > 10, ]

  common_genes <- unique(sig_sub$Gene_Symbol)
  for (gene in common_genes){
    message("Processing Celltype: ", ct, "; Gene: ", gene)

    pqtl_file <- pqtl_info[pqtl_info$gene_name == gene, ]
    pqtl_file <- pqtl_file[order(-pqtl_file$nonnorm_smpnorm_corr), "file"][1]
    pqtl_pt <- paste0("raw_data/", pqtl_file)
    pqtl <- fread(pqtl_pt)

    apas <- unique(sig_sub[sig_sub$Gene_Symbol == gene, ]$APATag)
    for (apa in apas){
      message("Processing Celltype: ", ct, "; Gene: ", gene, "; APATag: ", apa)

      apaqtl_exp <- apaqtl[gene == apa, ]

      if (nrow(apaqtl_exp) == 0) {
        message("  No apaQTL SNP for ", apa, " — skip")
        next
      }

      apaqtl_exp$rsid <- apaqtl_exp$SNP
      apaqtl_exp$pval <- apaqtl_exp$`p-value`
      chr_number <- unique(apaqtl_exp$chr)

      apaqtl_exp <- tryCatch({
        message("  Running ld_clump for ", apa, " ...")
        ld_clump(apaqtl_exp,
                 plink_bin = get_plink_exe(),
                 bfile =paste0("1000G.EUR.QC.",chr_number),
                 clump_kb=100,
                 clump_r2 =0.01)
      }, error = function(e) {
        message("  [ERROR] ld_clump failed for ", apa, " in ", ct)
        message("  Error message: ", e$message)
        return(NULL)
      })

      if (is.null(apaqtl_exp)) {
        message("  Skip ", apa, " due to ld_clump failure")
        next
      }

      if (nrow(apaqtl_exp) == 0) {
        message("  No SNP left after LD clumping for ", apa, " — skip")
        next
      }

      apaqtl_exp$rsid <- NULL
      apaqtl_exp$pval <- NULL
      apaqtl_exp <- data.frame(apaqtl_exp)
      apaqtl_exp <- format_data(
        apaqtl_exp,
        type = "exposure",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        pval_col = "p-value",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col = "maf"
      )

      snp_overlap <- sum(pqtl$rsids %in% apaqtl_exp$SNP)
      message("  Overlap SNPs: ", snp_overlap)
      if(snp_overlap == 0){
        message("  No overlapping SNPs for ", apa, " — skip")
        next
      }

      pqtl_out <- pqtl[pqtl$rsids %in% apaqtl_exp$SNP, ]
      pqtl_out <- data.frame(pqtl_out)
      pqtl_out <- format_data(
        pqtl_out,
        type = "outcome",
        snp_col = "rsids",
        beta_col = "Beta",
        se_col = "SE",
        pval_col = "Pval",
        effect_allele_col = "effectAllele",
        other_allele_col = "otherAllele",
        eaf_col = "ImpMAF"
      )

      dat1 <- harmonise_data(apaqtl_exp, pqtl_out)
      dat1$id.exposure <- apa
      dat1$id.outcome <- gene
      dat1$exposure <- "apaQTL"
      dat1$outcome <- "pQTL"

      mr1 <- mr(dat1)
      print(mr1)

      if (!"b" %in% colnames(mr1)) {
        message("  No column 'b' in MR result for ", apa, " — skip")
        next
      }

      or <- generate_odds_ratios(mr1)

      heterogeneity_dat <- mr_heterogeneity(dat1)
      ivw_het <- heterogeneity_dat[heterogeneity_dat$method == "Inverse variance weighted", ]
      or$Q <- if (nrow(ivw_het) == 0) NA else ivw_het$Q
      or$heterogeneity_pval <- if (nrow(ivw_het) == 0) NA else ivw_het$Q_pval
      or$Q_df <- if (nrow(ivw_het) == 0) NA else ivw_het$Q_df

      pleiotropy_dat <- mr_pleiotropy_test(dat1)
      or$pleiotropy_pval <- if (is.null(pleiotropy_dat$pval)) NA else pleiotropy_dat$pval

      not_ivw <- or$method != "Inverse variance weighted"
      or$Q[not_ivw] <- NA
      or$heterogeneity_pval[not_ivw] <- NA
      or$Q_df[not_ivw] <- NA
      or$pleiotropy_pval[not_ivw] <- NA
      or$celltype <- ct

      or_dt <- data.table(or)[, .(celltype, id.exposure, id.outcome, outcome, exposure, method, nsnp, b, se, pval, or,
                                  or_lci95, or_uci95, Q, heterogeneity_pval, Q_df, pleiotropy_pval)]
      fwrite(or_dt, file = or_file, sep = "\t", append = TRUE, col.names = !file.exists(or_file))

      dat1 <- data.table(dat1)
      fwrite(dat1, file = dat_file, sep = "\t", append = TRUE, col.names = !file.exists(dat_file))
    }
  }
}


######################### pQTL-GWAS-MR ######################### 
######################### apaQTL-GWAS-MR ######################### 


######################### result merge ######################### 
apaqtl_pqtl <- fread("apaQTL.pQTL.MR/02.apaQTL.pQTL.MR.res")
apaqtl_pqtl <- apaqtl_pqtl[
  (nsnp == 1 & method == "Wald ratio") |
    (nsnp > 1 & method == "Inverse variance weighted")
]
apaqtl_pqtl[, apaQTL.to.pQTL := ifelse(pval < 0.05, "sig.apaQTL.to.pQTL", NA_character_)]
apaqtl_pqtl <- apaqtl_pqtl[, c("celltype", "id.exposure", "id.outcome", "b", "se", "pval", "apaQTL.to.pQTL")]
colnames(apaqtl_pqtl) <- c("cellType", "APATag", "Gene_Symbol", "b.apaQTL.to.pQTL", "se.apaQTL.to.pQTL", "pval.apaQTL.to.pQTL", "apaQTL.to.pQTL")

pqtl_gwas <- fread("pQTL.GWAS.MR/03.pQTL.GWAS.MR.res")
pqtl_gwas <- pqtl_gwas[
  (nsnp == 1 & method == "Wald ratio") |
    (nsnp > 1 & method == "Inverse variance weighted")
]
pqtl_gwas[, pqtl.to.gwas := ifelse(pval < 0.05, "sig.pqtl.to.gwas", NA_character_)]
pqtl_gwas <- pqtl_gwas[, c("id.exposure", "id.outcome", "b", "se", "pval", "pqtl.to.gwas")]
colnames(pqtl_gwas) <- c("Gene_Symbol", "disease", "b.pqtl.to.gwas", "se.pqtl.to.gwas", "pval.pqtl.to.gwas", "pqtl.to.gwas")

apaqtl_gwas <- fread("apaQTL.GWAS.MR/04.apaQTL.GWAS.MR.res")
apaqtl_gwas <- apaqtl_gwas[
  (nsnp == 1 & method == "Wald ratio") |
    (nsnp > 1 & method == "Inverse variance weighted")
]
apaqtl_gwas[, apaqtl.to.gwas := ifelse(pval < 0.05, "sig.apaqtl.to.gwas", NA_character_)]
apaqtl_gwas <- apaqtl_gwas[, c("celltype", "id.exposure", "id.outcome", "b", "se", "pval", "apaqtl.to.gwas")]
colnames(apaqtl_gwas) <- c("cellType", "APATag", "disease", "b.apaqtl.to.gwas", "se.apaqtl.to.gwas", "pval.apaqtl.to.gwas", "apaqtl.to.gwas")

sig <- merge(sig, apaqtl_pqtl, by = c("cellType", "APATag", "Gene_Symbol"), all.x = TRUE, sort = FALSE)
sig <- merge(sig, pqtl_gwas, by = c("Gene_Symbol", "disease"), all.x = TRUE, sort = FALSE)
sig <- merge(sig, apaqtl_gwas, by = c("cellType", "APATag", "disease"), all.x = TRUE, sort = FALSE)

sig[, mediation := ifelse(
  pval.apaQTL.to.pQTL < 0.05 &
    pval.pqtl.to.gwas < 0.05 &
    pval.apaqtl.to.gwas < 0.05,
  "Yes", "No"
)]
sig[mediation == "Yes", indirect := b.apaQTL.to.pQTL * b.pqtl.to.gwas]
sig[mediation == "Yes", total := b.apaqtl.to.gwas]
sig[mediation == "Yes", PM := indirect / total]

sig$indirect_effect <- sig$b.apaQTL.to.pQTL * sig$b.pqtl.to.gwas
sig$direct_effect <- sig$b.apaqtl.to.gwas - sig$indirect_effect
sig$mediation_ratio <- sig$indirect_effect / sig$b.apaqtl.to.gwas
