#' @include ImmuneSpace.R
NULL


#' @importFrom RCurl getCurlHandle curlPerform basicTextGatherer
.ISCon$methods(
  downloadMatrix=function(x, summary = FALSE){
    cache_name <- paste0(x, ifelse(summary, "_sum", ""))
    if(is.null(data_cache[[cache_name]])){
      if(nrow(subset(data_cache[[constants$matrices]], name%in%x)) == 0){
        stop(sprintf("No matrix %s in study\n", x))
      }
      summary <- ifelse(summary, ".summary", "")
      if(config$labkey.url.path == "/Studies/"){
        path <- paste0("/Studies/", data_cache$GE_matrices[name == x, folder], "/")
      } else{
        path <- gsub("/$","",config$labkey.url.path)
      }
      link <- URLdecode(file.path(gsub("http:","https:", gsub("/$","",config$labkey.url.base)),
                                  "_webdav",  path, "@files/analysis/exprs_matrices",
                                  paste0(x, ".tsv", summary)))
      localpath <- .self$.localStudyPath(link)
      if(.self$.isRunningLocally(localpath)){
        fl <- localpath
        message("Reading local matrix")
        data_cache[[cache_name]] <<- fread(fl, header = TRUE, showProgress = FALSE)
      }else{
        opts <- config$curlOptions
        opts$netrc <- 1L
        #opts$httpauth <- 1L
        handle <- getCurlHandle(.opts=opts)
        h <- basicTextGatherer()
        message("Downloading matrix..")
        curlPerform(url = link, curl = handle, writefunction = h$update)
        fl <- tempfile()
        write(h$value(), file = fl)
        EM <- fread(fl, header = TRUE, showProgress = FALSE)
        if(nrow(EM) == 0){
          stop("The downloaded matrix has 0 rows. Something went wrong")
        }
        data_cache[[cache_name]] <<- EM
        file.remove(fl)
      }
      
    }else{
      data_cache[[cache_name]]
    }
  }
)

.ISCon$methods(
  GeneExpressionFeatures=function(matrix_name,summary=FALSE){
    cache_name <- paste0(matrix_name, ifelse(summary, "_sum", ""))
    if(!matrix_name %in% data_cache[[constants$matrices]][, name]){
      stop("Invalid gene expression matrix name");
    }
    annotation_set_id <- .self$.getFeatureId(matrix_name)
    if(is.null(data_cache[[.self$.mungeFeatureId(annotation_set_id)]])){
      if(!summary){
        message("Downloading Features..")
        featureAnnotationSetQuery = sprintf("SELECT * from FeatureAnnotation
                                            where FeatureAnnotationSetId='%s';",
                                            annotation_set_id);
        features <- labkey.executeSql(config$labkey.url.base, 
                                      config$labkey.url.path,
                                      schemaName = "Microarray",
                                      sql = featureAnnotationSetQuery,
                                      colNameOpt = "fieldname")
        setnames(features, "GeneSymbol", "gene_symbol")
      }else{
        features <- data.frame(FeatureId=data_cache[[cache_name]][,gene_symbol],
                             gene_symbol=data_cache[[cache_name]][,gene_symbol])
      }
      data_cache[[.self$.mungeFeatureId(annotation_set_id)]] <<- features
    }
  }
)

.ISCon$methods(
  ConstructExpressionSet=function(matrix_name, summary){
    cache_name <- paste0(matrix_name, ifelse(summary, "_sum", ""))
    #matrix
    message("Constructing ExpressionSet")
    matrix <- data_cache[[cache_name]]
    #features
    features <- data_cache[[.self$.mungeFeatureId(.self$.getFeatureId(matrix_name))]][,c("FeatureId","gene_symbol")]
    
    runID <- data_cache$GE_matrices[name == matrix_name, rowid]
    pheno_filter <- makeFilter(c("Run", "EQUAL", runID), 
                               c("Biosample/biosample_accession", "IN", paste(colnames(matrix), collapse = ";")))
    pheno <- unique(data.table(labkey.selectRows(
      config$labkey.url.base, config$labkey.url.path,
      #"assay.ExpressionMatrix.matrix", "InputSamples", "gene_expression_matrices",
      "study", "HM_InputSamplesQuerySnapshot",
      containerFilter = "CurrentAndSubfolders",
      colNameOpt = "caption", colFilter = pheno_filter)))
    
    setnames(pheno, .self$.munge(colnames(pheno)))

    #pheno <- pheno[, list(biosample_accession, ParticipantId, arm_name,
    pheno <- pheno[, list(biosample_accession, participant_id, cohort,
                          study_time_collected, study_time_collected_unit)]
    
    if(summary){
      fdata <- data.frame(FeatureId = matrix$gene_symbol, gene_symbol = matrix$gene_symbol, row.names = matrix$gene_symbol)
      rownames(fdata) <- fdata$FeatureId
      fdata <- AnnotatedDataFrame(fdata)
    } else{
      try(setnames(matrix, " ", "FeatureId"), silent = TRUE)
      try(setnames(matrix, "V1", "FeatureId"), silent = TRUE)
      fdata <- data.table(FeatureId = as.character(matrix$FeatureId))
      fdata <- merge(fdata, features, by = "FeatureId", all.x = TRUE)
      fdata <- as.data.frame(fdata)
      rownames(fdata) <- fdata$FeatureId
      fdata <- AnnotatedDataFrame(fdata)
    }
    dups <- colnames(matrix)[duplicated(colnames(matrix))]
    if(length(dups) > 0){
      for(dup in dups){
        dupIdx <- grep(dup, colnames(matrix))
        newNames <- paste0(dup, 1:length(dupIdx))
        setnames(matrix, dupIdx, newNames)
        eval(substitute(matrix[, `:=`(dup, rowMeans(matrix[, dupIdx, with = FALSE]))], list(dup = dup)))
        eval(substitute(matrix[, `:=`(newNames, NULL)], list(newNames = newNames)))
      }
      if(.self$config$verbose){
        warning("The matrix contains subjects with multiple measures per timepoint. Averaging expression values.")
      }
    }
    exprs <- as.matrix(matrix[, -1L, with = FALSE])
    exprs <- exprs[, colnames(exprs) %in% pheno$biosample_accession] #At project level, InputSamples may be filtered

    pheno <- data.frame(pheno)
    rownames(pheno) <- pheno$biosample_accession
    pheno <- pheno[colnames(exprs), ]
    ad_pheno <- AnnotatedDataFrame(data = pheno)
    es <- ExpressionSet(assayData = exprs,
                        phenoData = ad_pheno, featureData = fdata)
    data_cache[[cache_name]] <<- es
  }
)

# Downloads a normalized gene expression matrix from ImmuneSpace.
.ISCon$methods(
  getGEMatrix=function(x = NULL, cohort = NULL, summary = FALSE, reload=FALSE){
    "Downloads a normalized gene expression matrix from ImmuneSpace.\n
    `x': A `character'. The name of the gene expression matrix to download.\n
    `cohort': A `character'. The name of a cohort that has an associated gene
    expression matrix. Note that if `cohort' isn't NULL, then `x' is ignored.
    `summary': A `logical'. If set to TRUE, Downloads a matrix with expression
    averaged by gene symbol.
    `reload': A `logical'. If set to TRUE, the matrix will be downloaded again,
    even if a cached cop exist in the ImmuneSpaceConnection object."
    cohort_name <- cohort #can't use cohort = cohort in d.t
    if(!is.null(cohort_name)){
      if(all(cohort_name %in% data_cache$GE_matrices$cohort)){
        x <- data_cache$GE_matrices[cohort %in% cohort_name, name]
      } else{
        validCohorts <- data_cache$GE_matrices[, cohort]
        stop(paste("No expression matrix for the given cohort.",
                   "Valid cohorts:", paste(validCohorts, collapse = ", ")))
      }
    }
    cache_name <- paste0(x, ifelse(summary, "_sum", ""))
    if(length(x) > 1){
      data_cache[cache_name] <<- NULL
      lapply(x, downloadMatrix, summary)
      lapply(x, GeneExpressionFeatures,summary)
      lapply(x, ConstructExpressionSet, summary)
      ret <- .combineEMs(data_cache[cache_name])
      if(dim(ret)[[1]] == 0){
        # No features shared
        warn <- "The returned ExpressionSet has 0 rows. No feature is shared accross the selected runs or cohorts."
        if(!summary){
          warn <- paste(warn, "Try setting summary to TRUE to merge the matrices by gene symbol.")
        }
        warning(warn)
      }
      return(ret)
    } else{
      if (cache_name %in% names(data_cache) && !reload) {
        data_cache[[cache_name]]
      }
      else {
        data_cache[[cache_name]] <<- NULL
        downloadMatrix(x, summary)
        GeneExpressionFeatures(x, summary)
        ConstructExpressionSet(x, summary)
        data_cache[[cache_name]]
      }
    }
  }
)

# Combine EMs and output only genes available in all EMs.
#' @importFrom Biobase fData
.combineEMs <- function(EMlist){
  message("Combining ExpressionSets")
  fds <- lapply(EMlist, function(x){ droplevels(data.table(fData(x)))})
  fd <- Reduce(f = function(x, y){ merge(x, y, by = c("FeatureId", "gene_symbol"))}, fds)
  EMlist <- lapply(EMlist, "[", as.character(fd$FeatureId))
  for(i in 1:length(EMlist)){ fData(EMlist[[i]]) <- fd}
  res <- Reduce(f = combine, EMlist)
}


# Add treatment information to the phenoData of an expression matrix available in the connection object.
.ISCon$methods(
  addTrt=function(x = NULL){
    "Add treatment information to the phenoData of an expression matrix
     available in the connection object.\n
    x: A character. The name of a expression matrix that has been downloaded 
    from the connection."
    if(is.null(x) | !x %in% names(data_cache)){
      stop(paste(x, "is not a valid expression matrix."))
    } 
    bsFilter <- makeFilter(c("biosample_accession", "IN",
                             paste(pData(data_cache[[x]])$biosample_accession, collapse = ";")))
    bs2es <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
                                          "immport", "expsample_2_biosample",
                                          colFilter = bsFilter, colNameOpt = "rname"))
    esFilter <- makeFilter(c("expsample_accession", "IN",
                             paste(bs2es$expsample_accession, collapse = ";")))
    es2trt <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
                                           "immport", "expsample_2_treatment",
                                           colFilter = esFilter, colNameOpt = "rname"))
    trtFilter <- makeFilter(c("treatment_accession", "IN",
                              paste(es2trt$treatment_accession, collapse = ";")))
    trt <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
                                        "immport", "treatment",
                                        colFilter = trtFilter, colNameOpt = "rname"))
    bs2trt <- merge(bs2es, es2trt, by = "expsample_accession")
    bs2trt <- merge(bs2trt, trt, by = "treatment_accession")
    pData(data_cache[[x]])$treatment <<- bs2trt[match(pData(data_cache[[x]])$biosample_accession, biosample_accession), name]
    return(data_cache[[x]])
  }
)

.ISCon$methods(
  .getFeatureId=function(matrix_name){
    subset(data_cache[[constants$matrices]],name%in%matrix_name)[, featureset]
  }
)

.ISCon$methods(
  .mungeFeatureId=function(annotation_set_id){
    return(sprintf("featureset_%s",annotation_set_id))
  }
)


#' @importFrom Biobase pData sampleNames
.ISCon$methods(
  EMNames=function(EM = NULL, colType = "participant_id"){
    "Change the sampleNames of an ExpressionSet fetched by getGEMatrix using the
    information in the phenodData slot.\n
    x: An ExpressionSet, as returned by getGEMatrix.\n
    colType: A character. The type of column names. Valid options are 'expsample_accession'
    and 'participant_id'."
    if(is.null(EM) | !is(EM, "ExpressionSet")){
      stop("EM should be a valid ExpressionSet, as returned by getGEMatrix")
    }
    if(!all(grepl("^BS", sampleNames(EM)))){
      stop("All sampleNames should be biosample_accession, as returned by getGEMatrix")
    }
    pd <- data.table(pData(EM))
    colType <- gsub("_.*$", "", tolower(colType))
    if(colType == "expsample"){
      bsFilter <- makeFilter(c("biosample_accession", "IN",
                                 paste(pd$biosample_accession, collapse = ";")))
      bs2es <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
                                              "immport", "expsample_2_biosample",
                                              colFilter = bsFilter, colNameOpt = "rname"))
      pd <- merge(pd, bs2es[, list(biosample_accession, expsample_accession)], by = "biosample_accession")
      es <- pd[match(sampleNames(EM), pd$biosample_accession), expsample_accession]
      sampleNames(EM) <- pData(EM)$expsample_accession <- es
    } else if(colType %in% c("participant", "subject")){
      pd <- pd[, nID := paste0(participant_id, "_", tolower(substr(study_time_collected_unit, 1, 1)), study_time_collected)]
      sampleNames(EM) <- pd[ match(sampleNames(EM), pd$biosample_accession), nID]
    } else if(colType == "biosample"){
      warning("Nothing done, the column names should already be biosample_accession numbers.")
    } else{
      stop("colType should be one of 'expsample_accession', 'biosample_accession', 'participant_id'.")
    }
    return(EM)
  }
)
