suppressWarnings(suppressMessages(library("tidyverse")))
suppressMessages(library("MatrixEQTL"))
suppressMessages(library("utils"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))

VERSION = "V3.2.0 2024-12-17"

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

# Paraments-debug
if(FALSE){
  cellType = "CD4_CTL"
  
  mode = "pre,formal"
  pvOutputThreshold_cis = 1
  pvOutputThreshold_tra = 1e-4
  cisDist = 1e6  # Distance for local gene-SNP pairs

  # phenotype
  PREPHENOFILENAME = sprintf("%s_format.tsv", cellType)  # Raw phenotype file before processing
  PHENOTEMPFILENAME = sprintf("%s_phenoIntermediate.rds", cellType)  # Temporary file
  PHENOFILENAME = sprintf("%s_format.tsv", cellType)  # Each row represents gene expression, each column represents a sample
  
  # gene location file
  PREGENELOCFILENAME = "geneMetainfo.tsv"  # Raw reference genome file before processing
  GENELOCFILENAME = sprintf("%s_geneloc.tsv", cellType)  # First column: geneid, second column: chr<int>, third column: gene start, fourth column: gene end
  
  # Covariate file
  PRECOVARIATEFILENAME = sprintf("%s_peer_factors.tsv", cellType)  # Each row represents a covariate, each column represents a sample
  COVARIATETEMPFILENAME = sprintf("%s_covariatesIntermediate.rds", cellType)
  COVARIATEFILENAME = sprintf("%s_covariates.tsv", cellType)

  # Genotype file
  PREGENOTYPEPATH = "split"
  PREGENOTYPEREF = "snpTable_uniq.tsv"
  GENOTYPETEMPFILENAME = "genotypeIntermediate.rds"
  GENOTYPEFILENAME = sprintf("%s_genotype.tsv", cellType)  # Each row represents a genotype, each column represents a sample, values are 0, 1, 2
  
  # SNP location file
  SNPLOCFILENAME = "snpsloc.tsv"  # First column: snpid, second column: chr<int>, third column: pos
  ERRCOVFILENAME = NULL  # Error covariance matrix

  # Output files for Matrix eQTL
  OUTPUTCISTABLEFILENAME = sprintf("%s_cis.tsv", cellType)
  OUTPUTTRANSTABLEFILENAME = sprintf("%s_trans.tsv", cellType)
  OUTPUTPLOTFILENAME = sprintf("%s_qqplot.pdf", cellType)
  OUTPUTRDSFILENAME = sprintf("%s.rds", cellType)
}

# Paraments
if(TRUE){
  parser = list()
  parser[["mode"]] = make_option("--mode", dest="mode", type="character", default="pre,formal",
                                 help="mode")
  parser[["p_cis"]] = make_option("--p_cis", dest="p_cis", type="numeric", default=0,
                                  help="=0.01, pValue of cis-eQTL")
  parser[["p_trans"]] = make_option("--p_trans", dest="p_trans", type="numeric", default=0,
                                    help="=0.01, pValue of trans-eQTL")
  parser[["cisDist"]] = make_option("--cisDist", dest="cisDist", type="numeric", default=1e6,
                                    help="=1e6, Distance of cis-eQTL")
  parser[["prePheno"]] = make_option("--prePheno", dest="prePheno", type="character", default=NULL,
                                     help="=<celltype>_phenoImputed.tsv")
  parser[["phenoTemp"]] = make_option("--phenoTemp", dest="phenoTemp", type="character", default=NULL,
                                      help="=<celltype>_phenoIntermediate.rds")
  parser[["pheno"]] = make_option("--pheno", dest="pheno", type="character", default=NULL,
                                  help="<celltype>_phenoImputed.tsv")
  parser[["preGeneloc"]] = make_option("--preGeneloc", dest="preGeneloc", type="character", default="geneMetainfo.tsv",
                                       help="=geneMetainfo.tsv, gene location information extracted from reference genome annotation")
  parser[["geneloc"]] = make_option("--geneloc", dest="geneloc", type="character", default="geneloc.tsv",
                                    help="=geneloc.tsv, formatted gene location file")
  parser[["preCov"]] = make_option("--preCov", dest="preCov", type="character", default=NULL,
                                   help="=<celltype>_peer_factors.tsv")
  parser[["covTemp"]] = make_option("--covTemp", dest="covTemp", type="character", default=NULL,
                                   help="=<celltype>_peer_factors.rds, intermediate covariate file format")
  parser[["cov"]] = make_option("--cov", dest="cov", type="character", default="Mono_covariates.tsv",
                                help="<celltype>_covariates.tsv, formatted covariate file")
  parser[["preGenotype"]] = make_option("--preGenotype", dest="preGenotype", type="character", default="", help="")
  parser[["preGenotypeRef"]] = make_option("--preGenotypeRef", dest="preGenotypeRef", type="character", default="", help="")
  parser[["genotypeTemp"]] = make_option("--genotypeTemp", dest="genotypeTemp", type="character", default="genotypeIntermediate.rds",
                                         help="=<cellType>_genotypeIntermediate.rds, intermediate genotype processing file format")
  parser[["genotype"]] = make_option("--genotype", dest="genotype", type="character", default=NULL,
                                     help="=<celltype>_genotype.tsv")
  parser[["snp"]] = make_option("--snp", dest="snp", type="character", default="snpsloc.tsv",
                                help="=snpsloc.tsv")
  parser[["errcov"]] = make_option("--errcov", dest="errcov", type="character", default=NULL,
                                   help="")
  parser[["out_cis"]] = make_option("--out_cis", dest="out_cis", type="character", default=NULL,
                                    help="=<celltype>_cis.tsv")
  parser[["out_trans"]] = make_option("--out_trans", dest="out_trans", type="character", default=NULL,
                                      help="=<celltype>_trans.tsv")
  parser[["out_plot"]] = make_option("--out_plot", dest="out_plot", type="character", default=NULL,
                                     help="=<celltype>_qqplot.pdf")
  parser[["out_rds"]] = make_option("--out_rds", dest="out_rds", type="character", default=NULL,
                                    help="=<celltype>.rds")
  
  args = parse_args(OptionParser(option_list=parser))
  mode = args$mode
  pvOutputThreshold_cis = args$p_cis
  pvOutputThreshold_tra = args$p_trans
  cisDist = args$cisDist
  PREPHENOFILENAME = args$prePheno
  PHENOTEMPFILENAME = args$phenoTemp
  PHENOFILENAME = args$pheno
  PREGENELOCFILENAME = args$preGeneloc
  GENELOCFILENAME = args$geneloc
  PRECOVARIATEFILENAME = args$preCov
  COVARIATETEMPFILENAME = args$covTemp
  COVARIATEFILENAME = args$cov
  PREGENOTYPEPATH = args$preGenotype
  PREGENOTYPEREF = args$preGenotypeRef
  GENOTYPETEMPFILENAME = args$genotypeTemp
  GENOTYPEFILENAME = args$genotype
  SNPLOCFILENAME = args$snp
  ERRCOVFILENAME = args$errcov
  OUTPUTCISTABLEFILENAME = args$out_cis
  OUTPUTTRANSTABLEFILENAME = args$out_trans
  OUTPUTPLOTFILENAME = args$out_plot
  OUTPUTRDSFILENAME = args$out_rds
}

# Process paraments
if(TRUE){
  setDTthreads(30)
  mode = stringr::str_split(mode, pattern=',')[[1]]
  if("pre" %in% mode){
    phenoLoaded = FALSE
    genotypeLoaded = FALSE
    covariateLoaded = FALSE
  }
  if("formal" %in% mode){
    # Select model - linear model
    useModel = MatrixEQTL::modelLINEAR  # modelANOVA, modelLINEAR, modelLINEAR_CROSS
  }
  
  PREPHENOFILENAME = checkPath(PREPHENOFILENAME, mode="file", exists=TRUE)
  PHENOTEMPFILENAME = checkPath(PHENOTEMPFILENAME, mode="file", exists=FALSE)
  PHENOFILENAME = checkPath(PHENOFILENAME, mode="file", exists=FALSE)
  PREGENELOCFILENAME = checkPath(PREGENELOCFILENAME, mode="file", exists=TRUE)
  GENELOCFILENAME = checkPath(GENELOCFILENAME, mode="file", exists=FALSE)
  PRECOVARIATEFILENAME = checkPath(PRECOVARIATEFILENAME, mode="file", exists=TRUE)
  COVARIATETEMPFILENAME = checkPath(COVARIATETEMPFILENAME, mode="file", exists=FALSE)
  COVARIATEFILENAME = checkPath(COVARIATEFILENAME, mode="file", exists=FALSE)
  PREGENOTYPEPATH = checkPath(PREGENOTYPEPATH, mode="dir", exists=FALSE)
  PREGENOTYPEREF = checkPath(PREGENOTYPEREF, mode="file", exists=TRUE)
  GENOTYPETEMPFILENAME = checkPath(GENOTYPETEMPFILENAME, mode="file", exists=FALSE)
  GENOTYPEFILENAME = checkPath(GENOTYPEFILENAME, mode="file", exists=FALSE)
  SNPLOCFILENAME = checkPath(SNPLOCFILENAME, mode="file", exists=FALSE)
  OUTPUTCISTABLEFILENAME = checkPath(OUTPUTCISTABLEFILENAME, mode="file", exists=FALSE)
  OUTPUTTRANSTABLEFILENAME = checkPath(OUTPUTTRANSTABLEFILENAME, mode="file", exists=FALSE)
  OUTPUTPLOTFILENAME = checkPath(OUTPUTPLOTFILENAME, mode="file", exists=FALSE)
  OUTPUTRDSFILENAME = checkPath(OUTPUTRDSFILENAME, mode="file", exists=FALSE)
}

# Print paraments
if(TRUE){
  cat('\n', rep('=', 47), "CONFIG", rep('=', 47), '\n', sep='')
  cat("[Version]", VERSION, '\n', sep='')
  cat("[Date]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  cat("[Parament]mode:", mode, '\n')
  cat("[Parament]p_cis:", pvOutputThreshold_cis, '\n')
  cat("[Parament]p_trans:", pvOutputThreshold_tra, '\n')
  cat("[Parament]cisDist:", cisDist, '\n')
  cat("[Parament]prePheno:", PREPHENOFILENAME, '\n')
  cat("[Parament]phenoTemp:", PHENOTEMPFILENAME, '\n')
  cat("[Parament]pheno:", PHENOFILENAME, '\n')
  cat("[Parament]preGeneloc:", PREGENELOCFILENAME, '\n')
  cat("[Parament]geneloc:", GENELOCFILENAME, '\n')
  cat("[Parament]preCov:", PRECOVARIATEFILENAME, '\n')
  cat("[Parament]covTemp:", COVARIATETEMPFILENAME, '\n')
  cat("[Parament]cov:", COVARIATEFILENAME, '\n')
  cat("[Parament]preGenotype:", PREGENOTYPEPATH, '\n')
  cat("[Parament]preGenotypeRef:", PREGENOTYPEREF, '\n')
  cat("[Parament]genotypeTemp:", GENOTYPETEMPFILENAME, '\n')
  cat("[Parament]genotype:", GENOTYPEFILENAME, '\n')
  cat("[Parament]snp:", SNPLOCFILENAME, '\n')
  cat("[Parament]errcov:", ERRCOVFILENAME, '\n')
  cat("[Parament]out_cis:", OUTPUTCISTABLEFILENAME, '\n')
  cat("[Parament]out_trans:", OUTPUTTRANSTABLEFILENAME, '\n')
  cat("[Parament]out_plot:", OUTPUTPLOTFILENAME, '\n')
  cat("[Parament]out_rds:", OUTPUTRDSFILENAME, '\n')
  cat(rep('-', 48), "LOG", rep('-', 49), '\n', sep='')
}

# Formatting
if("pre" %in% mode){
  cat("[Progress]Mode pre is on.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  # Formatting phenotype
  cat("[Progress]Checking pheno file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  if(!file.exists(PHENOTEMPFILENAME)){
    cat("\t[Info]Formatting pheno file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    pheno = utils::read.delim(PREPHENOFILENAME) %>%  # loading
      t() %>%
      as.data.frame() %>%  # Convert to a table where columns represent individuals and rows represent tags
      rownames_to_column("geneid")
    phenoLoaded = TRUE  # Flag that phenotype file has been loaded
    cat("\t[Info]Saving pheno intermediate file:\n\t\t", PHENOTEMPFILENAME, '\n')
    saveRDS(pheno, PHENOTEMPFILENAME)
  }else{cat("\t[Info]Existed pheno file:\n\t\t", PHENOTEMPFILENAME, '\n')}
  # Formatting gene position
  cat("[Progress]Checking gene loc file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  if(!file.exists(GENELOCFILENAME)){
    cat("\t[Info]Formatting gene loc file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    # Using gene coordinates
    if(FALSE){
      geneRef = utils::read.delim(PREGENELOCFILENAME, header=FALSE,
                                  col.names=c("chr", "left", "right", "strand", "geneid", "geneVersion", "geneName", "geneBiotype")) %>%
        select(geneid, chr, left, right) %>%
        mutate(chr=paste0("chr", chr))
    }
    # Using PAS coordinates as gene positions
    if(TRUE){
      if(!exists("pheno")){pheno=readRDS(PHENOTEMPFILENAME)}
      geneRef = data.table(geneId=pheno[["geneid"]]) %>%
        .[, chr:=tstrsplit(geneId, '\\.')[2]] %>%
        .[, left:=tstrsplit(geneId, '\\.')[3]] %>%
        .[, left:=as.integer(left)] %>%
        .[, right:=left+1]
    }
    cat("\t[Info]Saving gene loc formatted file:\n\t\t", GENELOCFILENAME, '\n')
    utils::write.table(geneRef, GENELOCFILENAME, sep='\t', row.names=FALSE, quote=FALSE)
  }else{cat("\t[Info]Existed gene loc file:\n\t\t", GENELOCFILENAME, '\n')}
  # Formatting genotype
  cat("[Progress]Checking genotype file..5.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  if(!file.exists(GENOTYPETEMPFILENAME)){
    # Load genotype data
    genotypeFile = paste0("chr", seq(1,22))
    names(genotypeFile) = genotypeFile
    genotype = lapply(genotypeFile, function(chrom){
      filepath = file.path(PREGENOTYPEPATH, sprintf("%s.tsv", chrom))
      cat("\t[Info]Loading file:", filepath, '\n')
      # Convert to matrix format (rows=SNPs, columns=individuals, values=0/1/2)
      data = fread(filepath) %>%
        data.table::transpose(keep.names="snp", make.names="inv")
      return(data)
    }) %>%
      rbindlist() %>%
      column_to_rownames("snp") %>%
      as.matrix()
    # Load SNP position and allele information  
    snpInfo = fread(PREGENOTYPEREF, sep='\t', quote=FALSE) %>%
      .[, .(rsid, chrom_hg19, pos_hg38)] %>%
      .[, chrom_hg19:=fifelse(stringr::str_detect(chrom_hg19, pattern="^chr"), chrom_hg19, paste0("chr",chrom_hg19))] %>%
      setnames(c("snpid", "chr", "pos"))
    cat("\t[Info]Working...\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    # Save SNP position data  
    cat("\t[Info]Saving file:\n\t\t", SNPLOCFILENAME, '\n')
    write.table(snpInfo, SNPLOCFILENAME, sep='\t', quote=FALSE, row.names=FALSE)
    # Save genotype matrix  
    cat("\t[Info]Saving file:\n\t\t", GENOTYPETEMPFILENAME, '\n')
    saveRDS(genotype, GENOTYPETEMPFILENAME, compress=FALSE)
    rm(genotype)
    gc()
  }else{cat("\t[Info]Existed genotype intermediate file:\n\t\t", GENOTYPETEMPFILENAME, '\n')}
  # Formatting cov
  cat("[Progress]Checking covariate file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  if(!file.exists(COVARIATETEMPFILENAME)){
    cov = read.delim(PRECOVARIATEFILENAME)
    covariateLoaded = TRUE
    cov[["sampleid"]] = as.character(cov[["sampleid"]])
    colnames(cov)[which(colnames(cov)=="sampleid")] = "id"
    colnames(cov)[which(colnames(cov)=="Gender")] = "gender"
    colnames(cov)[which(colnames(cov)=="Age")] = "age"
    cov = as.data.frame(t(cov))
    cov[["info"]] = rownames(cov)
    cov = cov[, c("info", colnames(cov)[-which(colnames(cov)=="info")])]
    cov["id",] = as.character(cov["id", ])
    rownames(cov) = NULL
    colnames(cov) = NULL
    cat("\t[Info]Saving covariate intermediate file:\n\t\t", COVARIATETEMPFILENAME, '\n')
    saveRDS(cov, COVARIATETEMPFILENAME)
  }else{cat("\t[Info]Existed covariate file:\n\t\t", COVARIATETEMPFILENAME, '\n')}
  
  cat("[Progress]Formatting pheno and genotype and cov file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  if(!file.exists(PHENOFILENAME) | !file.exists(GENOTYPEFILENAME) | !file.exists(COVARIATEFILENAME)){
    if(!phenoLoaded){cat("\t[Info]Loading file:\n\t\t", PHENOTEMPFILENAME, '\n'); pheno=readRDS(PHENOTEMPFILENAME)}
    if(!genotypeLoaded){cat("\t[Info]Loading file:\n\t\t", GENOTYPETEMPFILENAME, '\n'); genotype=readRDS(GENOTYPETEMPFILENAME)}
    if(!covariateLoaded){cat("\t[Info]Loading file:\n\t\t", COVARIATETEMPFILENAME, '\n'); cov=readRDS(COVARIATETEMPFILENAME)}
    cat("\t[Info]Working...\n")
    pheno = as.data.frame(pheno)
    genotype = as.data.frame(genotype)
    cov = as.data.frame(cov)
    gc()
    sampleidPheno = setdiff(colnames(pheno), "geneid") %>%
      as.character()
    sampleidGenotype = setdiff(colnames(genotype), "sampleid") %>%
      as.character()
    sampleidCov = as.character(cov[1, ][-1])
    sampleid = sampleidPheno %>%
      base::intersect(sampleidGenotype) %>%
      base::intersect(sampleidCov) %>%
      as.character()

    pheno = pheno[, c("geneid", sampleid)]
    genotype = genotype[, sampleid] %>%
      rownames_to_column("snpid")
    colnames(cov) = cov[1, ]
    cov = cov[, c("id", sampleid)]
    colnames(cov) = NULL
    
    cat("\t[Info]Saving phenotype file:\n\t\t", PHENOFILENAME, '\n')
    data.table::fwrite(pheno, PHENOFILENAME, sep='\t', row.names=FALSE, quote=FALSE)
    cat("\t[Info]Saving genotype file:\n\t\t", GENOTYPEFILENAME, '\n')
    data.table::fwrite(genotype, GENOTYPEFILENAME, sep='\t', row.names=FALSE, quote=FALSE)
    cat("\t[Info]Saving covariate file:\n\t\t", COVARIATEFILENAME, '\n')
    data.table::fwrite(cov, COVARIATEFILENAME, sep='\t', row.names=FALSE, quote=FALSE)
    rm(pheno)
    rm(genotype)
    rm(cov)
    gc()
  }else{cat("\t[Info]Existed phenotype file:\n\t\t", PHENOFILENAME, '\n');cat("\t[Info]Existed genotype file:\n\t\t", GENOTYPEFILENAME, '\n');cat("\t[Info]Existed covariate file:\n\t\t", COVARIATEFILENAME, '\n')}
}

# main
if("formal" %in% mode){
  cat("[Progress]Mode formal is on.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  # Loading
  if(TRUE){
    cat("[Progress]Loading file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    # phenotype
    cat("\t[Info]Loading file:", PHENOFILENAME, '\n')
    gene = MatrixEQTL::SlicedData$new()
    gene$fileDelimiter = "\t"      # the TAB character
    gene$fileOmitCharacters = "NA" # denote missing values;
    gene$fileSkipRows = 1          # one row of column labels
    gene$fileSkipColumns = 1       # one column of row labels
    gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
    gene$LoadFile(PHENOFILENAME)
    # gene pos
    cat("\t[Info]Loading file:", GENELOCFILENAME, '\n')
    genepos = read.table(GENELOCFILENAME, header=TRUE, stringsAsFactors=FALSE)
    # genotype
    cat("\t[Info]Loading file:", GENOTYPEFILENAME, '\n')
    snps = MatrixEQTL::SlicedData$new()
    snps$fileDelimiter = "\t"          # tab
    snps$fileOmitCharacters = "NA"    # NA
    snps$fileSkipRows = 1             # 
    snps$fileSkipColumns = 1          # 
    snps$fileSliceSize = 2000         # read file in slices of 2,000 rows
    snps$LoadFile(GENOTYPEFILENAME)
    # genotype pos
    cat("\t[Info]Loading file:", SNPLOCFILENAME, '\n')
    snpspos = read.table(SNPLOCFILENAME, header=TRUE, stringsAsFactors=FALSE)
    # cov
    cat("\t[Info]Loading file:", COVARIATEFILENAME, '\n')
    cvrt = MatrixEQTL::SlicedData$new()
    cvrt$fileDelimiter = "\t"      # the TAB character
    cvrt$fileOmitCharacters = "NA" # denote missing values;
    cvrt$fileSkipRows = 1          # one row of column labels
    cvrt$fileSkipColumns = 1      # one column of row labels
    if(is.null(COVARIATEFILENAME)){NULL}else{cat("\t[Info]Loading file:", COVARIATEFILENAME, '\n');cvrt$LoadFile(COVARIATEFILENAME)}

    if(is.null(ERRCOVFILENAME)){errorCovariance=numeric()}else{errorCovariance=read.table(ERRCOVFILENAME)}
  }
  # Running QTL
  if(TRUE){
    cat("[Progress]Working...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    me = Matrix_eQTL_main(snps=snps,
                          gene=gene,
                          cvrt=cvrt,
                          output_file_name=OUTPUTTRANSTABLEFILENAME,
                          pvOutputThreshold=pvOutputThreshold_tra,
                          useModel=useModel, 
                          errorCovariance=errorCovariance, 
                          verbose=TRUE,
                          #pvalue.hist=TRUE,
                          pvalue.hist="qqplot",
                          min.pv.by.genesnp=FALSE,
                          noFDRsaveMemory=FALSE,
                          cisDist=cisDist,  # maximum distance at which gene-SNP pair is considered local
                          pvOutputThreshold.cis=pvOutputThreshold_cis,  #  p-value threshold for local eQTLs
                          snpspos=snpspos,  #
                          genepos=genepos,  #
                          output_file_name.cis=OUTPUTCISTABLEFILENAME
                          )
  }
  # Saving
  if(TRUE){
    cat("[Progress]Saving file...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
    if(!is.null(OUTPUTRDSFILENAME)){
      cat("\t[Info]Saving file:\n\t\t",OUTPUTRDSFILENAME,'\n')
      saveRDS(me, OUTPUTRDSFILENAME)
    }
    cat("\t[Info]Analysis done in:", me[["time.in.sec"]], "seconds", '\n')
    cat("\t[Info]Num of cis-eQTLs:", me$cis$neqtls, '\n')
    cat("\t[Info]Num of trans-eQTLs:", me$trans$neqtls,'\n')
    cat("\t[Info]Saving file:\n\t\t",OUTPUTPLOTFILENAME,'\n')
    pdf(OUTPUTPLOTFILENAME, width=5, height=5)
    plot(me)
    dev.off()
  }
}

cat(rep('-', 30), paste0('[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ']', "All parts finished"), rep('-', 30), '\n', sep='')
