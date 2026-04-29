library(parallel)
library(data.table)
library(coloc)
library(dplyr)

cellTypeV <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_CTL", "CD8_Naive", "CD8_TCM", "CD8_TEM",
               "T_Reg", "gdT", "MAIT", "NK", "NK_R", "B_Mem", "B_IN", "Mono_C", "Mono_NC", "DC")
diseaseV <- c("MS.GCST003566", "PSO.GCST90243956", "SLE.GCST90011866", "T1DM.GCST90014023",
              "CD.ibdgc", "IBD.ibdgc", "RA.jenger", "Celiac.phecode557")

individualQTL_pt <- "individuals2.rds"
individualSummary = readRDS(individualQTL_pt) %>%
  purrr::map(~ length(.x)) %>%
  as.data.table() %>%
  data.table::transpose(keep.names="cellType") %>%
  setnames("V1", "num") %>%
  .[cellType %in% cellTypeV] %>%
  .[, .(cellType, num)]

window <- 1e6
for (ct in cellTypeV){
  message("Celltype: ", ct)

  N.qtl <- individualSummary[cellType==ct, num]

  for (dis in diseaseV){
    message("Celltype: ", ct, "; Disease: ", dis)

    gwas_sig_pt <- paste0(dis, ".GWAS.summary", ".sig")
    message("Celltype: ", ct, "; Disease: ", dis, "; GWAS_sig_pt: ", gwas_sig_pt)
    gwas_sig <- fread(gwas_sig_pt)

    clump_pt <- paste0(dis, ".leadSNP.clump.txt")
    message("Celltype: ", ct, "; Disease: ", dis, "; Clump_pt: ", clump_pt)
    clump_snp <- fread(clump_pt)
    clump_snp <- clump_snp[, c("rsid", "bp", "id")]
    setnames(clump_snp, "rsid", "SNP")
    clump_snp <- merge(clump_snp, gwas_sig, by = "SNP")
    rm(gwas_sig)
    gc(FALSE)

    chrs <- unique(clump_snp$chr.gwas)
    for (chr in chrs){
      message("Celltype: ", ct, "; Disease: ", dis, "; Chr: ", chr)

      qtl_pt <- paste0(ct, ".QTL.summary", ".", chr)
      message("Celltype: ", ct, "; Disease: ", dis, "; QTL_pt: ", qtl_pt)
      qtl <- fread(qtl_pt, select=c("SNP","chr.qtl","pos.qtl","gene.qtl",
                                    "beta.qtl","t-stat.qtl","A1.qtl","A2.qtl","p-value.qtl"))
      qtl[, se.qtl := abs(beta.qtl / `t-stat.qtl`)]
      qtl[, chr.qtl := paste0("chr", chr.qtl)]

      gwas_pt <- paste0(dis, ".GWAS.summary", ".", chr)
      message("Celltype: ", ct, "; Disease: ", dis, "; Gwas_pt: ", gwas_pt)
      gwas <- fread(gwas_pt, select=c("SNP","chr.gwas","pos.gwas",
                                      "b.gwas","se.gwas","p.gwas","freq.gwas","A1.gwas","A2.gwas","N.gwas"))

      clump_snp_chr <- clump_snp[chr.gwas == chr]

      for (i in seq_len(nrow(clump_snp_chr))) {
        lead <- clump_snp_chr[i]
        message("Celltype: ", ct, "; Disease: ", dis, "; GWAS loci: ", lead$SNP)
        chr0 <- lead$chr.gwas
        pos0 <- as.numeric(lead$pos.gwas)

        gwas_win <- gwas[chr.gwas == chr0 &
                           pos.gwas %between% c(pos0 - window, pos0 + window)]
        qtl_win <- qtl[chr.qtl == chr0 &
                         pos.qtl %between% c(pos0 - window, pos0 + window)]

        if (nrow(gwas_win) == 0 || nrow(qtl_win) == 0) next
        gwas_win <- gwas_win[order(p.gwas), .SD[1], by = SNP]

        genes <- unique(qtl_win$gene.qtl)
        for (gene0 in genes) {
          message("Celltype: ", ct, "; Disease: ", dis, "; GWAS loci: ", lead$SNP, "; Gene: ", gene0)

          qtl_gene <- qtl_win[gene.qtl == gene0]

          common_snp <- intersect(gwas_win$SNP, qtl_gene$SNP)
          if (length(common_snp) == 0) next

          gwas2 <- gwas_win[SNP %in% common_snp]
          qtl2  <- qtl_gene[SNP %in% common_snp]

          setkey(gwas2,SNP)
          setkey(qtl2,SNP)

          m <- qtl2[gwas2, nomatch=0]
          m[, allele_match := fifelse(
            A1.gwas == A1.qtl & A2.gwas == A2.qtl, "same",
            fifelse(
              A1.gwas == A2.qtl & A2.gwas == A1.qtl, "flip",
              "drop"
            )
          )]
          m <- m[allele_match != "drop"]
          m <- m[
            complete.cases(
              SNP,
              b.gwas,
              se.gwas,
              p.gwas,
              freq.gwas,
              beta.qtl,
              se.qtl,
              `p-value.qtl`
            )
          ]

          if (nrow(m) == 0) next
          m[allele_match=="flip", beta.qtl := -beta.qtl]
          m[, allele_match := NULL]

          m <- m[m$se.gwas != 0 & m$se.qtl != 0, ]
          coloc_res <- coloc.abf(
            dataset1 = list(
              snp = m$SNP, beta = m$b.gwas, varbeta = m$se.gwas^2,
              pvalues = m$p.gwas, N = unique(m$N.gwas), type = "cc"
            ),
            dataset2 = list(
              snp = m$SNP, beta = m$beta.qtl, varbeta = m$se.qtl^2,
              pvalues = m$`p-value.qtl`, N = N.qtl, type="quant"
            ),
            MAF = m$freq.gwas
          )

          coloc_all <- as.data.table(coloc_res$results)
          setorder(coloc_all, -SNP.PP.H4)
          coloc_all[, `:=`(
            celltype = ct,
            disease  = dis,
            chr      = chr0,
            lead_snp = lead$SNP,
            lead_pos = pos0,
            gene     = gene0,
            nsnp     = nrow(m)
          )]
          setcolorder(
            coloc_all,
            c("celltype", "disease", "chr", "lead_snp", "lead_pos", "gene", "nsnp")
          )

          coloc_all_pt <- paste0(ct, "_", dis, ".coloc.all.txt")
          fwrite(coloc_all, file = coloc_all_pt, sep = "\t",
                 append = file.exists(coloc_all_pt),
                 col.names = !file.exists(coloc_all_pt))

          coloc_sum <- data.table(
            celltype = ct,
            disease  = dis,
            chr      = chr0,
            lead_snp = lead$SNP,
            lead_pos = pos0,
            gene     = gene0,
            nsnp     = nrow(m),
            PP.H0 = coloc_res$summary["PP.H0.abf"],
            PP.H1 = coloc_res$summary["PP.H1.abf"],
            PP.H2 = coloc_res$summary["PP.H2.abf"],
            PP.H3 = coloc_res$summary["PP.H3.abf"],
            PP.H4 = coloc_res$summary["PP.H4.abf"]
          )

          coloc_sum_pt <- paste0(ct, "_", dis, ".coloc.sum.txt")
          fwrite(coloc_sum, file = coloc_sum_pt, sep = "\t",
                 append = file.exists(coloc_sum_pt),
                 col.names = !file.exists(coloc_sum_pt))

          rm(coloc_res, coloc_all, coloc_sum, m, gwas2, qtl2, qtl_gene, common_snp)
        }
      }
    }
  }
  message("==== Finished: ", ct, " ====")
}
