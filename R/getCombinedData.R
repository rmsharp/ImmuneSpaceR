
# Dependencies --------------------------------------------------------

# helper functions ----------------------------------------------------

# make sure arguments are allowed and handle unwanted vals gracefully
.checkArgs <- function(func, args){
  args <- str_
  
  cdata_allowed <- c("hai","neut_ab_titer")
  if( !(clin_data %in% cdata_allowed) ){ 
    stop(paste0("clinical data input must be one of the following: ", cdata_allowed))
  }
  
  csum_allowed <- c(0, 1, 2)
  if( !(clin_sum %in% csum_allowed) ){ 
    stop(paste0("cdata summary input must be one of the following: ", csum_allowed))
  }
  
  if(func = "gCD"){
    merge_allowed <- c("all", "study", "em")
    if( !(merge_by %in% merge_allowed) ){ 
      stop(paste0("merge type must be one of the following: ", merge_allowed))
    }
    
    if(eset == T & clin_sum == 0){
      stop("clinical data must be summarized to output an eset, please use 'clin_sum' = 1 or 2")
    }
    
    if(eset == T & is.null(clin_time)){
      stop("clinical data must be summarized to output an eset, 
           please set 'clin_time' to value in study_time_collected")
    }
  }else if(func = "mM"){
    
  }
}

.getEMData <- function(EM_list){
  suppressMessages(lapply(EM_list, .self$downloadMatrix, summary = T)) 
  init_ems <- data_cache[ -which(names(data_cache) == "GE_matrices")]
}

.getEMSubsPheno <- function(EM_list){
  pheno <- unique(data.table(labkey.selectRows(
    baseUrl = config$labkey.url.base, 
    folderPath = config$labkey.url.path,
    schemaName = sn_study, 
    queryName = qn_InputSmplsShot,
    containerFilter = cf_currandsubs,
    colNameOpt = cn_caption)))
  
  pheno <- pheno[ which(pheno$`Expression Matrix Accession` %in% EM_list), ]
}

.bs2sub <- function(pheno){
  subids <- paste0(pheno$`Participant ID`, "_d", pheno$`Study Time Collected`)
  subids <- gsub("-", "neg", subids, fixed = T)
  bs2id <- hash(pheno$`Biosample Accession`, subids)
}

.getCommonGenes <- function(clean_ems){
  genes_list <- lapply(clean_ems, FUN = function(em){ em$gene_symbol })
  common_genes <- Reduce(intersect, genes_list)
  if(length(common_genes) == 0){ stop("No common genes found") }
  return(common_genes)
}

.filterGenes <- function(em, genes){
  em <- em[ which(em$gene_symbol %in% genes), ]
}

.clean_EM <- function(em, genes, id_hash, time_filt){
  if( !is.null(genes) ){ em <- .filterGenes(em, genes) }
  genes_list <- em$gene_symbol
  em <- em[ , -c("gene_symbol")]
  em <- .map_cols(em, id_hash)
  keep <- grep(paste0("d", time_filt), colnames(em))
  em <- em[ , keep, with = F ]
  colnames(em) <- gsub("_d.*", "", colnames(em))
  em$gene_symbol <- genes_list
  return(em)
}

.getFilteredEMs <- function(EM_list, pheno){
  # Get all GE data into list of dfs and combine
  if( is.null(EM_list) ){ EM_list <- GeneExpressionMatrices()$name }
  init_ems <- .getEMData(EM_list)
  
  # Setup biosample2subject hash
  bs2id <- .bs2sub(pheno)
  clean_ems <- lapply(init_ems, .clean_EM, genes = genes, id_hash = bs2id, time_filt = ge_time)
  
  # Ensure EMs all contain same common genes
  common_genes <- .getCommonGenes(clean_ems)
  filt_ems <- lapply(clean_ems, .filterGenes, genes = common_genes)
}

.cleanHai <- function(cdata, clin_time, clin_sum){
  if(!is.null(clin_time)){ cdata <- cdata[ which(cdata$study_time_collected == clin_time), ] }
  keep <- which(colnames(cdata) != "virus")
  if( clin_sum != 0){
    cdata <- cdata[ , keep, with = F] # avoids use of dplyr::select()
    if( clin_sum == 1 ){
      cdata <- cdata %>%
        group_by(participant_id) %>%
        filter(value_reported == max(value_reported)) %>%
        distinct(.keep_all = TRUE) %>%
        ungroup()
    } else if( clin_sum == 2){
      cdata <- cdata %>%
        group_by(participant_id) %>%
        mutate(value_reported = mean(value_reported)) %>%
        distinct(.keep_all = TRUE) %>%
        ungroup()
    }
  }
  return(cdata)
}

.hashmap <- function(input, hsh){
  values <- unname(unlist(sapply(input, FUN = function(x){
    val <- hsh[[x]]
  })))
}

.map_cols <- function(df, hsh){
  colnames(df) <- .hashmap(hsh = hsh, input = colnames(df))
  return(df)
}

.match_subs <- function(em, cdata){
  gs <- em$gene_symbol
  tmp <- em[ , which(colnames(em) %in% unique(cdata$participant_id)), with = F ]
  tmp$gene_symbol <- gs
  return(tmp)
}

# The method expects a list (elem) with two dataframes representing EM and cData
.prep_eset <- function(elem){
  fData <- data.frame(elem$EM$gene_symbol, stringsAsFactors = F)
  rownames(fData) <- elem$EM$gene_symbol
  colnames(fData) <- "gene_symbol"
  expr_em <- elem$EM[ , -c("gene_symbol") ]
  
  pData <- elem$cData[ order(match(elem$cData$participant_id, colnames(expr_em))), ]
  rownames(pData) <- pData$participant_id
  pData <- droplevels(pData)
  
  eset <- ExpressionSet(assayData = as.matrix(expr_em),
                        featureData = AnnotatedDataFrame(fData),
                        phenoData = AnnotatedDataFrame(pData))
}


tester <- function(a, b, c){
  print(sys.call())
}
# Main Method --------------------------------------------------------
# NOTE: clin_sum controls whether the clinical data is summarized if there is
# more than one row per subject and, if so, how.
# clin_sum options include --
# 0: not summarized, 1: summarized by max value, 2: summarized by average value
# summarizing by average will cause the data to drop the virus column if HAI

#' @importFrom hash hash
#' @importFrom dplyr filter mutate distinct ungroup
#' @importFrom magrittr %>%
#' @importFrom Biobase ExpressionSet
.ISCon$methods(
  getCombinedData = function(clin_data = "hai", 
                             clin_time = NULL,
                             clin_sum = 0,
                             EM_list = NULL,
                             ge_time = 0,
                             genes = NULL, 
                             merge_by = "all",
                             eset = FALSE
                             ){
    # TBD!!!!!
    # .checkArgs(func = "gCD", args = sys.call())
    
    # get GE data from ImmuneSpace and filter so all EMs have common gene set
    pheno <- .getEMSubsPheno(EM_list)
    filt_ems <- .getFilteredEMs(EM_list, pheno)
    
    # pull clinical data with filter for subs and clean
    sub_filter <- makeFilter(c("participant_id", 
                               "IN", 
                               paste0(unique(pheno$`Participant ID`), collapse = ";")))
    
    cdata <- getDataset(clin_data, colFilter = sub_filter)
    cdata <- .cleanHai(cdata, clin_time, clin_sum)
  
    # make sure only cData subs are in EM
    filt_ems <- lapply(clean_ems, .match_subs, cdata)
    
    # for mapping EM names to study ID
    em_runs <- data.table(labkey.selectRows(baseUrl = config$labkey.url.base, 
                                            folderPath = config$labkey.url.path,
                                            schemaName = sn_assayExprMx, 
                                            queryName = qn_Runs))
    
    # return output in either one merged EM, by study, or by EM name
    res <- list()
    if(merge_by == "all"){
      res$EM <- Reduce(merge, filt_ems)
      res$cData <- cdata
      
      message("----Data for single EM and cData set---")
      message(paste0("Common Genes: ", length(common_genes)))
      message(paste0("Number of Subjects: ", length(colnames(res$EM)) - 1))
      if(eset == TRUE){ res <- .prep_eset(res) }
      
    }else if(merge_by == "study"){
      studies <- unique(gsub(".*\\.", "", cdata$participant_id))
      res <- sapply(studies, simplify = F, USE.NAMES = T, FUN = function(sdy){
        tmp <- list()
        tmp$cdata <- cdata[ grep(paste0(".*", sdy), cdata$participant_id), ]
        em_nms <- em_runs$Name[ which(em_runs$Study == paste0("SDY", sdy)) ]
        em_nms <- paste0(em_nms, "_sum")
        em_dfs <- filt_ems[ which(names(filt_ems) %in% em_nms) ]
        tmp$EM <- Reduce(merge, em_dfs)
        
        message(paste0("------- Data for Study ", sdy, "--------"))
        message(paste0("Number of EMs combined: ", length(em_dfs)))
        message(paste0("Number of Subjects: ", length(colnames(tmp$EM)) - 1))
        message(paste0("Number of Genes: ", dim(tmp$EM)[1]))
        return(tmp)
      })
      if(eset == TRUE){ res <- lapply(res, .prep_eset) }
      
    }else{
      em_nms <- em_runs$Name[ which(em_runs$Name %in% gsub("_sum","", names(filt_ems))) ]
      em_nms <- paste0(em_nms, "_sum")
      res <- sapply(em_nms, simplify = F, USE.NAMES = T, FUN = function(em_nm){
        tmp <- list()
        tmp$EM <- filt_ems[[em_nm]]
        tmp$cData <- cdata[ which(cdata$participant_id %in% colnames(tmp$EM)), ]
        
        message(paste0("------- Data for EM: ", em_nm, "--------"))
        message(paste0("Number of Subjects: ", length(colnames(tmp$EM)) - 1))
        message(paste0("Number of Genes: ", dim(tmp$EM)[1]))
        
        return(tmp)
      })
      if(eset == TRUE){ res <- lapply(res, .prep_eset) }
    }
    
    return(res)
  }
)

.ISCon$methods(
  makeMAE = function(cDataList, clin_time, clin_sum, EM_list, ge_time){
    # get GE data from ImmuneSpace and filter so all EMs have common gene set
    pheno <- .getEMSubsPheno(EM_list)
    filt_ems <- .getFilteredEMs(EM_list, pheno)
    
    # Map EM to study ID and separate by study
    em_runs <- data.table(labkey.selectRows(baseUrl = config$labkey.url.base, 
                                            folderPath = config$labkey.url.path,
                                            schemaName = sn_assayExprMx, 
                                            queryName = qn_Runs))
    
    em2sdy <- hash(keys = em_runs$Name, values = em_runs$Study)
    pheno$study <- sapply(pheno$`Expression Matrix Accession`, .hashmap, em2sdy)
    
    studies <- unique(em_runs$Study[ which(em_runs$Name %in% EM_list) ])
    
    EM_list <- sapply(studies, simplify = F, USE.NAMES = T, FUN = function(sdy){
      em_nms <- em_runs$Name[ which(em_runs$Study == paste0("SDY", sdy)) ]
      em_nms <- paste0(em_nms, "_sum")
      em_dfs <- filt_ems[ which(names(filt_ems) %in% em_nms) ]
      EM <- Reduce(merge, em_dfs)
    })
    
    # for all items in cDatalist, go pull using the subfilter
    cdata <- lapply(studies, FUN = function(sdy){
      sdy_con <- CreateConnection(sdy) # using new con for each to limit downloads
      subs_list <- unique(pheno$`Participant ID`[ which(pheno$study == sdy) ])
      sub_filter <- makeFilter(c("participant_id", 
                                 "IN", 
                                 paste0(subs_list, collapse = ";")))
      tmp_list <- lapply(cDataList, FUN = function(data_nm){
        data <- sdy_con$getDataset(data_nm, colFilter = sub_filter)
      })
    }) 
    
    # for each cdata element, clean / prep the data depending on type
    cdata$hai <- .cleanHai(cdata$hai, clin_time, clin_sum)
    
    # prepMAE
    
    # makeMAE
    
    # returnMAE
    
  }
)

