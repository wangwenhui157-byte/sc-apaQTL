library("reticulate")
library("yaml")
library("tidyverse")
library("Seurat")
library("future")
library("future.apply")
library("phateR")
library("slingshot")
library("dittoSeq")
library("ggpubr")
library("Matrix")
library("data.table")
library("stringr")
library("preprocessCore")
library("impute")
library("utils")
library("lme4")
library("qvalue")
library("ComplexHeatmap")
library("rsvg")
library("qpdf")

# Paraments-debug
if(TRUE){
  FDRCUTOFF = 0.05
  stageMax = 6
  leadSNPMarker = TRUE
  CLASSIFY = "utr"
  
  CONFIG = "config.yml"
  FILEPATH = "onek1kPasta"
  SEURATFILE = file.path(FILEPATH, "pseudoB/seuratB.rds")
  SIGSNPFILE = file.path(FILEPATH, ifelse(CLASSIFY=="utr","QTL","QTLIntron"), "resultCombined_sigSNP_cis.tsv")
  LINKSNPFILE = file.path(FILEPATH, ifelse(CLASSIFY=="utr","QTL","QTLIntron"), "resultCombined_leadSNP_cis.tsv")
  APATAGMAP = file.path(FILEPATH, ifelse(CLASSIFY=="utr","pheno","phenoIntron"), "tagRef.tsv")
  MATRIXEQTLPATH = file.path(FILEPATH, "pseudoB", CLASSIFY, "QTL", sprintf("Q%s", stageMax))
  COVPATH = file.path(FILEPATH, "pseudoB", CLASSIFY, "pheno", sprintf("Q%s", stageMax), "peer")
  GWASFILE = "GWAS/All_associations_v1.0.2.tsv"
  
  PHENOPATH = file.path(FILEPATH, "pseudoB", CLASSIFY, "pheno", sprintf("Q%s", stageMax))
  GENOPATH = "hg19hg38/split"
  SNPREFFILE = "hg19hg38/snpTable_uniq.tsv"
  
  CHECKPATH = file.path(FILEPATH, "pseudoB", CLASSIFY, "check")
  
  PLOTPATH = file.path(FILEPATH, "plot")
}

# Function
if(TRUE){
  mergeTagSNP = function(tag, snpId, pheno, genotype, cov, regression=FALSE){
    tag = tag
    snpId = snpId
    
    pheno = pheno
    genotype = genotype
    cov = cov
    regression = regression
    

    if(!tag %in% colnames(pheno)){cat("\t[Info]", tag, "not existed in pheno.\n");return(NA)}

    phenoTemp = pheno[, .SD, .SDcols=c("newId", "sampleId", tag)] %>%
      na.omit() %>%
      setnames(tag, "tag") %>%
      .[, sampleId:=as.character(sampleId)] %>%
      .[, tag:=as.numeric(tag)]
    

    genotypeTemp = genotype[, snpId, drop=FALSE] %>%
      as.matrix() %>%
      as.data.frame() %>%
      setnames("snps") %>%
      rownames_to_column("sampleId") %>%
      mutate(sampleId=as.character(sampleId),
             snps=as.integer(snps))
    

    df = merge.data.table(phenoTemp, genotypeTemp, by="sampleId", all=FALSE)
    df[, sampleId:=NULL]

    df = merge.data.table(df, cov, by="newId", all=FALSE) %>%
      na.omit() %>%
      .[, stage:=stringr::str_extract(stage, pattern="[0-9]")] %>%
      setcolorder(c("newId", setdiff(names(.), "newId")))
    

    if(regression){

      df = df %>%
        mutate(across(c("tag", "snps", "stage", "Gender", "Age"), as.numeric)) %>%
        mutate(across(starts_with("pf") | starts_with("pc"), as.numeric))
      formula = "tag ~ Gender + Age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2" 
      model = lm(formula , data=df)
      residuals = resid(model)
      df[["tag_regression"]] = residuals
    }
    

    colCharacter = c("newId", "sampleId")
    colNumeric = c("tag", "snps", "Age", "stage", "Gender",
                   names(df)[stringr::str_detect(names(df), pattern="pc|pf")])
    df[, (colCharacter):=lapply(.SD, as.character), .SDcols=colCharacter]
    df[, (colNumeric):=lapply(.SD, as.numeric), .SDcols=colNumeric]
    
    return(df)
  }
  fitLmm = function(tag, snpId, pheno, genotype, cov, progress=NULL){
    tag = tag
    snpId = snpId
    pheno = pheno
    genotype = genotype
    cov = cov
    
    if(!is.null(progress)){
      progress <<- progress +1
      cat(progress, '\r')
    }
    
    df = mergeTagSNP(tag=tag, snpId=snpId, pheno=pheno, genotype=genotype, cov=cov)
    if(is.na(df)){return(NA)}
    df[, snps_stage:=snps*stage]
    
    if(length(unique(df[["stage"]]))==1){
      return(NA)
    }
    # Fit null model
    fit0 = lme4::lmer(tag ~ snps + Age + Gender + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2 + (1 | sampleId) + stage,
                      data=df,
                      REML=FALSE)
    # Fit augmented model
    fit = lme4::lmer(tag ~ snps + Age + Gender + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2 + (1 | sampleId) + stage + snps_stage,
                     data=df,
                     REML=FALSE)
    # Likelihood Ratio Test
    anovaResult = anova(fit, fit0)
    

    singularMarker = isSingular(fit)

    fit@pp=merPredD()
    fit@frame = data.frame()
    fit@resp = lmerResp()
    
    funcResult = list("anova"=anovaResult,
                      "singular"=singularMarker)
    return(funcResult)
  }
  fitLmmSq = function(tag, snpId, pheno, genotype, cov, progress=NULL){
    tag = tag
    snpId = snpId
    pheno = pheno
    genotype = genotype
    cov = cov
    
    if(!is.null(progress)){
      progress <<- progress +1
      cat(progress, '\r')
    }
    
    df = mergeTagSNP(tag=tag, snpId=snpId, pheno=pheno, genotype=genotype, cov=cov)
    if(is.na(df)){return(NA)}
    
    if(length(unique(df[["stage"]]))==1){
      return(NA)
    }
    
    df[, stage2:=stage*stage]
    df[, snps_stage:=snps*stage]
    df[, snps_stage2:=snps*stage*stage]
    
    # 空模型
    fit0 = lme4::lmer(tag ~ snps + Age + Gender + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2 + (1 | sampleId) + stage + stage2 + snps_stage,
                      data=df,
                      REML=FALSE)
    # Fit augmented model
    fit = lme4::lmer(tag ~ snps + Age + Gender + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pf1 + pf2 + (1 | sampleId) + stage + stage2 + snps_stage + snps_stage2,
                     data=df,
                     REML=FALSE)
    # Perform LRT
    anovaResult = anova(fit, fit0)
    

    singularMarker = isSingular(fit)

    fit@frame = data.frame()
    fit@flist=list()
    fit@pp=merPredD()
    
    funcResult = list("anova"=anovaResult,
                      "singular"=singularMarker)
    return(funcResult)
  }
  makeStagePlot = function(tag, geneName, snpId, rsid, major, minor, 
                           pheno=pheno, genotype=genotype, cov=cov, regression=FALSE, reverse=FALSE){
    tag = tag
    geneName = geneName
    snpId = snpId
    rsid = rsid
    if(is.na(rsid)){rsid=snpId}
    
    major = major
    minor = minor
    
    pheno = pheno
    genotype = genotype
    cov = cov
    
    regression = regression
    reverse = reverse
    
    # set color of quantile
    colorConfig = (c("Q1"="#E69F00", "Q2"="#56B4E9", "Q3"="#009E73",
                     "Q4"="#F0E442", "Q5"="#0072B2", "Q6"="#D55E00",
                     "Q7"="grey", "Q8"="grey"))
    

    dtPlot = mergeTagSNP(tag=tag, snpId=snpId, pheno=pheno, genotype=genotype, cov=cov, regression=regression)
    if(regression){
      dtPlot = dtPlot[, .(stage, tag_regression, snps)] %>%
        dplyr::rename(tag=tag_regression)
    }else{
      dtPlot = dtPlot[, .(stage, tag, snps)]
    }
    dtPlot[, stage:=paste0('Q', stage)]
    
    

    for(i in unique(dtPlot$snps)){
      template1 = data.table(snps=i, stage="Q7", tag=NA)
      template2 = data.table(snps=i, stage="Q8", tag=NA)
      dtPlot = rbind(dtPlot, template1, template2)
    }
    

    dtPlot[, snps:=factor(snps, levels=c(0,1,2), labels=c(paste0(major,major), paste0(major,minor), paste0(minor,minor)))]

    if(reverse){
      dtPlot[, snps:=factor(snps,
                            levels=c(paste0(minor,minor), paste0(major,minor), paste0(major,major)),
                            labels=c(paste0(minor,minor), paste0(minor,major), paste0(major,major)))]}
    

    setorder(dtPlot, snps, stage)
    dtPlot[, group:=paste(snps, stage, sep='-')]
    dtPlot$group = factor(dtPlot$group, levels=unique(dtPlot$group))

    dtPlot = dtPlot[, stage:=factor(stage, levels=unique(stage))]

    outliers = dtPlot[, {low=quantile(tag, 0.25, na.rm=T)-1.5*(quantile(tag, 0.75, na.rm=T)-quantile(tag, 0.25, na.rm=T))
    high=quantile(tag, 0.75, na.rm=T)+1.5*(quantile(tag, 0.75, na.rm=T)-quantile(tag, 0.25, na.rm=T))
    .(tag, outlier=tag<low | tag>high)},
    by=.(snps, stage)]
    
    outliers = outliers[outlier==TRUE]
    outliers = outliers[!duplicated(tag), ]
    outliers[, group:=paste(snps, stage, sep='-')]

    if(TRUE){
      plot = ggplot(dtPlot, aes(x=group, y=tag)) + 
        geom_boxplot(aes(color=stage), width=0.75, alpha=0.9,
                     outlier.shape=NA,
                     position=position_dodge(0.5)) +
        geom_point(data=outliers, aes(color=stage),
                   position=position_dodge(0.5), size=0.3) +
        geom_smooth(method="loess", se=FALSE, aes(group=snps), linewidth=0.5, color="black", span=0.75) +
        scale_color_manual(values=colorConfig) +
        
        labs(title=sprintf("%s(%s)", geneName, rsid), x="", y="", file="") +
        huTheme(subPlot=TRUE) +
        theme(axis.line=element_line(),
        ) +
        theme(axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
              axis.ticks.x=element_blank(),
        ) +
        theme(axis.title.y=element_text(margin=margin(r=15)),
              axis.ticks.y=element_line(),
        ) +
        theme(legend.position="top",
              legend.title=element_blank(),
              legend.direction="horizontal",
              legend.box="horizontal",
        ) +
        guides(color=guide_legend(ncol=6, byrow=TRUE, reverse=F,
                                  label.position="top",
        ),
        ) +
        scale_x_discrete(breaks=c(paste0(major,major,'-',"Q4"), paste0(major,minor,'-',"Q4"), paste0(minor,minor,'-',"Q4")),
                         labels=c(paste0(major,major), paste0(major,minor), paste0(minor,minor)))
    }
    
    return(plot)
  }
}

# Process paraments
if(TRUE){
  plan(sequential)
  setDTthreads(30)
  type = ifelse(leadSNPMarker==TRUE, "lead", "sig")
  options(future.globals.maxSize=2*1024*1024*1024)
}
if(TRUE){
  cat("[Parament]stageMax:", stageMax, '\n')
  cat("[Parament]leadSNPMarker:", leadSNPMarker, '\n')
}

# Loading
if(TRUE){
  cat("[Progress]Loading data...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  config = yaml::read_yaml(CONFIG)
  DISEASE = config[["disease"]] %>%
    unlist()
  # seuratObject
  seuratB = readRDS(SEURATFILE)
  countsMatrix = GetAssayData(seuratB, assay="RNA", slot="counts")
  metadata = seuratB@meta.data
  rm(seuratB);gc()
  # check stages
  if(!sprintf("Q%s", stageMax) %in% colnames(metadata)){
    cat("\t[Info]Adding stageMax info...\n")
    metadata[[sprintf("Q%s", stageMax)]] = Hmisc::cut2(metadata[["pseudoTime"]],
                                                       g=stageMax) %>%
      `levels<-`(paste0("Q",seq_len(stageMax)))
  }
  # load sigSNP
  sigSNP = fread(SIGSNPFILE, sep='\t', quote=FALSE) %>%
    setkey(NULL) %>%
    setorder(-leadSNP)
  linkSNP = fread(LINKSNPFILE, sep='\t', quote=FALSE) %>%
    .[cellType %in% c("B_IN", "B_Mem")] %>%
    .[, .(snps)] %>%
    unique() %>%
    setnames(c("snps"), c("leadSNP")) %>%
    setorder(leadSNP)
  # pheno
  pheno = paste0('Q', seq(1,stageMax)) %>%
    setNames(., .) %>%
    lapply(function(stage){
      cat('+')
      phenoStage = fread(file.path(PHENOPATH, sprintf("%s_phenoImputed.tsv", stage)), sep='\t', quote=FALSE) %>%
        dplyr::rename(sampleId=V1) %>%
        mutate(sampleId=as.character(sampleId))
      return(phenoStage)
    }) %>%
    rbindlist(idcol="stage") %>%
    .[, newId:=paste0(sampleId, '_', stage)] %>%
    select(newId, everything())
  # cov
  cov = paste0('Q', seq(1,stageMax)) %>%
    setNames(., .) %>%
    purrr::map(~ fread(file.path(COVPATH, sprintf("%s_peer_factors.tsv", .x)), sep='\t', quote=FALSE)) %>%
    rbindlist(idcol="stage") %>%
    dplyr::rename(sampleId=sampleid) %>%
    dplyr::mutate(newId=paste0(sampleId, '_', stage)) %>%
    dplyr::mutate(sampleId=as.character(sampleId),
                  Gender=as.character(Gender)) %>%
    dplyr::select(newId, everything())
  # QTL
  snpTemp = sigSNP[["snps"]] %>%
    unique()
  tagTemp = sigSNP[["APATag"]] %>%
    unique()
  matrixQTL = paste0('Q', seq_len(stageMax)) %>%
    setNames(., .) %>%
    lapply(function(stage){
      fread(file.path(MATRIXEQTLPATH, sprintf("%s_cis.tsv", stage)), sep='\t', quote=FALSE) %>%
        .[SNP %in% snpTemp] %>%
        .[gene %in% tagTemp]
    }) %>%
    rbindlist(idcol="stage") %>%
    setnames(c("SNP", "gene"),
             c("snps", "APATag")) %>%
    .[, id:=paste(APATag, snps, sep='@')]
  # geno
  genotype = paste0("chr", seq(1,22)) %>%
    setNames(., .) %>%
    purrr::map(~ list("path"=file.path(GENOPATH, sprintf("%s.tsv", .x)),
                      "chrom"=.x))
  # GWAS
  dtGWAS = fread(GWASFILE, sep='\t', quote=FALSE) %>%
    .[, .SD, .SDcols=c("DISEASE/TRAIT", "MAPPED_TRAIT", "SNP_GENE_IDS", "REPORTED GENE(S)", "SNPS")] %>%
    .[grepl(sprintf("\\b(?:%s)\\b", paste(purrr::map_chr(DISEASE, ~paste0("\\Q", .x, "\\E")), collapse='|')), MAPPED_TRAIT, ignore.case=TRUE)]
}

# main
if(TRUE){
  # Filter genes: keep only those with counts>=3 and expressed in >=3 developmental stages  -->  geneRemained
  if(FALSE){
    cat("[Progress]Filtering gene(counts>=3 and existed in >=3 stage)...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    # Filter genes with counts>=3 in each developmental stage
    geneList = split(metadata, metadata[[paste0('Q', stageMax)]]) %>%
      # Extract barcodes for each developmental stage
      purrr::map(~ rownames(.x)) %>%
      #  Identify genes expressed in each stage
      lapply(function(barcode){
        # Extract count matrix for the current stage
        temp = countsMatrix[, barcode]
        # Consider a gene expressed if counts>=3 in any cell
        gene = rownames(temp)[rowSums(temp>=3) > 0]
        return(gene)
      })
    # Filter genes expressed in at least 3 developmental stages
    geneRemained = geneList %>%
      unlist() %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq>=3) %>%
      dplyr::pull('.') %>%
      as.character()
    geneRemained = sapply(geneRemained, function(g){
      if(grepl(g, pattern="^ENSG[0-9]+\\.[0-9]+$")){return(stringr::str_remove(g, pattern="\\.[0-9]+$"))}
      if(grepl(g, pattern="-ENSG")){return(stringr::str_remove(g, pattern="-ENSG.*"))}
      return(g)
    })
    geneRemained = unique(geneRemained)
  }
  # APATag~snp pairs  --> sigSNP
  if(TRUE){
    cat("[Progress]Filtering sigSNP...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    sigSNP = sigSNP %>%
      .[cellType %in% c("B_IN", "B_Mem")] %>%  # APATag~snp pairs
      .[, cellType:=NULL] %>%
      .[APATag %in% colnames(pheno)]  # APATag
    if(exists("geneRemained")){
      sigSNP = sigSNP%>%
        .[(geneId %in% geneRemained) | (geneName %in% geneRemained)]  # pairs
    }
    cat(sprintf("\t[Info]Sig: there are %s gene-snp pairs, %s APATag and %s snps. geneRemained: %s.",
                nrow(unique(sigSNP[, .(APATag, snps)])), length(unique(sigSNP$APATag)), length(unique(sigSNP$snps)), exists("geneRemained")
    ), '\n')
    cat(sprintf("\t[Info]Lead: there are %s gene-snp pairs, %s APATag and %s snps. geneRemained: %s.",
                nrow(unique(sigSNP[leadSNP==TRUE, .(APATag, snps)])), length(unique(sigSNP[leadSNP==TRUE, APATag])), length(unique(sigSNP[leadSNP==TRUE, snps])), exists("geneRemained")
    ), '\n')
    
    type = ifelse(leadSNPMarker==TRUE, "lead", "sig")
  }

  if(leadSNPMarker==TRUE){
    sigSNP = sigSNP[leadSNP==TRUE] %>%
      unique(by="pairs")
  }else{
    sigSNP = unique(sigSNP, by="pairs")
  }
  
  # genotype
  if(TRUE){
    cat("[Progress]Loading genotype...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\t\t")
    snpV = split(sigSNP, by="chromSNP") %>%
      purrr::map(~ .x[["snps"]] %>% unique())
    # loading
    genotype = lapply(genotype, function(info){
      cat('+')
      chrom = info[["chrom"]]
      path = info[["path"]]
      snpRemained = snpV[[chrom]]
      
      geno = fread(path, select=c("inv", snpRemained)) %>%
        column_to_rownames("inv") %>%
        as.matrix()
      return(geno)
    })
    cat('\n')
    # ref/alt
    snpInfo = fread(SNPREFFILE, sep='\t', quote=FALSE, select=c("rsid", "ref_hg19", "alt_hg19", "chrom_hg19")) %>%
      split(by="chrom_hg19")

    check = all(sapply(genotype, rownames) == rownames(genotype[[1]]))
    if(!check){stop("[Error]The rownames of genotype isn't same.\n")}
    genotype = do.call(cbind, genotype)
  }
  # regression model
  if(TRUE){
    cat("[Progress]Running model...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\t\t")
    # gene-snp pairs
    pairs = sigSNP[, .(snps, APATag)] %>%
      split(., seq_len(nrow(.))) %>%
      purrr::map(~ c("tag"=.x[["APATag"]], "snpId"=.x[["snps"]]))
    for(i in seq_along(pairs)){
      cat(i, '\r')
      names(pairs)[i] = paste0(pairs[[i]][["tag"]], '@', pairs[[i]][["snpId"]])
    }
    
    # linear regression
    progress = 0
    res = future_lapply(pairs, function(x){
      tag = x[["tag"]]
      snpId = x[["snpId"]]
      l = fitLmm(tag=tag, snpId=snpId, #progress=progress,
                 pheno=pheno, genotype=genotype, cov=cov)
      return(l)
    })
    
    saveRDS(res, file.path(CHECKPATH, sprintf("res_%sSNP_Q%s.rds", type, stageMax)))
    res = readRDS(file.path(CHECKPATH, sprintf("res_%sSNP_Q%s.rds", type, stageMax)))

    sigSNP[, anova:=lapply(res, function(x){if(!is.na(x)){x$anova}else{NA}})]
    sigSNP[, singular:=sapply(res, function(x){if(!is.na(x)){x$singular}else{NA}})]
    sigSNP[, pValueModel1:=sapply(anova, function(x){if(is.na(x)){NA}else{x[2, "Pr(>Chisq)"]}})]
    sigSNP[, FDRModel1:=qvalue(pValueModel1)$lfdr]
    
    
    # quadratic regression
    progress = 0
    resSq = future_lapply(pairs, function(x){
      tag = x[["tag"]]
      snpId = x[["snpId"]]
      l = fitLmmSq(tag=tag, snpId=snpId, #progress=progress,
                   pheno=pheno, genotype=genotype, cov=cov)
      return(l)
    })
    saveRDS(resSq, file.path(CHECKPATH, sprintf("resSq_%sSNP_Q%s.rds", type, stageMax)))
    resSq = readRDS(file.path(CHECKPATH, sprintf("resSq_%sSNP_Q%s.rds", type, stageMax)))

    sigSNP[, anovaSq:=lapply(resSq, function(x){if(!is.na(x)){x$anova}else{NA}})]
    sigSNP[, singularSq:=sapply(resSq, function(x){if(!is.na(x)){x$singular}else{NA}})]
    sigSNP[, pValueModel2:=sapply(anovaSq, function(x){if(is.na(x)){NA}else{x[2, "Pr(>Chisq)"]}})]
    sigSNP[, FDRModel2:=qvalue(pValueModel2)$lfdr]
    
    
    # mark model
    sigSNP[, linear_significant:=fifelse(FDRModel1<FDRCUTOFF, TRUE, FALSE)]
    sigSNP[, quadratic_significant:=fifelse(FDRModel2<=FDRCUTOFF, TRUE, FALSE)]
    setorder(sigSNP, FDRModel1, FDRModel2)
    sigSNP = sigSNP[!is.na(FDRModel1) | !is.na(FDRModel2)]
    
    # saving
    filename = file.path(CHECKPATH, sprintf("model_%sSNP_Q%s.rds", type, stageMax))
    cat("\t[Info]Saving file:", filename, '\n')
    saveRDS(sigSNP, filename)
  }
  # filter
  if(TRUE){
    # FDR<0.05 and is not singular
    stageSigSNP = sigSNP[(FDRModel1<FDRCUTOFF & !singular) | (FDRModel2<FDRCUTOFF & !singularSq)]
    stageSigSNP[, method:=fifelse(linear_significant & !quadratic_significant,
                                  "Linear",
                                  fifelse(!linear_significant & quadratic_significant,
                                          "Quadratic",
                                          "Linear & Quadratic"))]
    stageSigSNP[, method:=factor(method, levels=c("Linear", "Linear & Quadratic", "Quadratic"))]
    cat(sprintf("\t[Info]There are %s pairs, %s APATag, %s SNP",
                length(unique(stageSigSNP[["pairs"]])), length(unique(stageSigSNP[["APATag"]])), length(unique(stageSigSNP[["snps"]])) ))
  }
  # combine qtl
  if(TRUE){
    idStageSigSNP = stageSigSNP[["pairs"]]
    matrixQTL = matrixQTL[id %in% idStageSigSNP]

    stageSigSNP = stageSigSNP %>%
      .[, .(pairs,
            APATag, geneId, geneName, chromGene, strand, start, end,
            leadSNP, snps, chromSNP, pos, Ref, Alt, major, minor,
            statistic, pValue, localFDR, estimate,
            method, anova,anovaSq, singular, singularSq,
            pValueModel1, pValueModel2, FDRModel1, FDRModel2, linear_significant, quadratic_significant)]
    

    matrixQTL = matrixQTL[,.(data=list(.SD)), by=.(snps, APATag, id)]
    table(matrixQTL[, sapply(data, nrow)])

    stageSigSNP = merge.data.table(stageSigSNP, matrixQTL[, .(id, data)], by.x="pairs", by.y="id") %>%
      setkey(NULL)
    stageSigSNP[, n_quantiles:=sapply(data, nrow)]
  }
  
  # process beta
  if(TRUE){
    betas = stageSigSNP[, .(pairs, data)] %>%
      tidyr::unnest("data") %>%
      as.data.table() %>%
      dcast.data.table(pairs ~ stage, value.var="beta") %>%
      .[rowSums(!is.na(.[, .SD, .SDcols=patterns("Q")]))>=2] %>%
      .[, lapply(.SD, as.numeric), by=pairs] %>%
      as.data.frame() %>%
      column_to_rownames("pairs")
    # scale
    betaScaled = betas %>%
      t() %>%
      scale() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("id")
  }
  # GWAS
  if(TRUE){
    if(leadSNPMarker==TRUE){
      linkSNP[, disease:=ifelse(leadSNP %in% dtGWAS[["SNPS"]], TRUE, FALSE)]
      leadSNPDisease = linkSNP[disease==TRUE, leadSNP] %>%
        unique()
      stageSigSNP[, disease:=ifelse(snps %in% leadSNPDisease, TRUE, FALSE)]
    }else{
      stageSigSNP[, disease:=ifelse(snps %in% dtGWAS[["SNPS"]], TRUE, FALSE)]
    }
    
  }
  
  # save
  if(TRUE){
    filepath = file.path(CHECKPATH, sprintf("stageSigSNP_%sSNP_Q%s.rds", type, stageMax))
    cat("\t[Info]Saving file:", filepath, '\n')
    saveRDS(stageSigSNP, filepath)
    
    filepath = file.path(CHECKPATH, sprintf("betaScaled_%sSNP_Q%s.rds", type, stageMax))
    cat("\t[Info]Saving file:", filepath, '\n')
    saveRDS(betaScaled, filepath)
    
    filepath = file.path(CHECKPATH, sprintf("genotype_%sSNP_Q%s.rds", type, stageMax))
    cat("\t[Info]Saving file:", filepath, '\n')
    saveRDS(genotype, filepath)
    
    
    temp = copy(stageSigSNP)
    test = list()
    for(i in seq_len(nrow(temp))){
      line = temp[i, ]
      snp = line[, snps]
      apa = line[, APATag]
      data = copy(line[, data][[1]]) %>%
        .[, snps:=snp] %>%
        .[, APATag:=apa]
      test[[i]] = data
    }
    test = rbindlist(test)
    test = merge.data.table(temp, test, by=c("snps", "APATag")) %>%
      .[, data:=NULL] %>%
      .[, anova:=NULL] %>%
      .[, anovaSq:=NULL]
    
    filepath = file.path(CHECKPATH, sprintf("slingshot_summary_%sSNP.tsv", type))
    print(filepath)
    fwrite(test, filepath, sep='\t', quote=FALSE)
    
  }
  
}
