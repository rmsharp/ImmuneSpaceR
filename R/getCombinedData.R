
# Dependencies --------------------------------------------------------
#' @importFrom hash hash

# helper functions ----------------------------------------------------
filterGenes <- function(df, genes){
  df <- df[ which(df[["gene_symbol"]] == genes), ]
}

hashmap <- function(hsh, input){
  values <- unname(unlist(sapply(input, FUN = function(x){
    val <- hsh[[x]]
  })))
}

map_cols <- function(df, hsh){
  colnames(df) <- hashmap(hsh = hsh, input = colnames(df))
  return(df)
}

clean_EM <- function(em, genes, id_hash, time_filt){
  if( !is.null(genes) ){ em <- filterGenes(em, genes) }
  genes_list <- em$gene_symbol
  em <- em[ , -c("gene_symbol")]
  em <- map_cols(em, id_hash)
  keep <- grep(paste0("d", time_filt), colnames(em))
  em <- em[ , keep, with = F ]
  colnames(em) <- gsub("_d.*", "", colnames(em))
  em$gene_symbol <- genes_list
  return(em)
}

filterEM <- function(em, common_genes){
  em[ which(em$gene_symbol %in% common_genes), ]
}

# Main Method --------------------------------------------------------

# expressionSet Maker
#' @importFrom hash hash
.ISCon$methods(
  getCombinedData = function(clin_data, 
                             clin_time = NULL,
                             EM_list = NULL,
                             ge_time = 0,
                             genes = NULL, 
                             output = "all"
                             ){
    
    # make sure clinical data is appropriate
    cdata_allowed <- c("hai")
    if( !(clin_data %in% cdata_allowed) ){ 
      stop(paste0("clinical data type must be one of the following: ", cdata_allowed))
    }
    
    output_allowed <- c("all", "study", "em")
    if( !(output %in% output_allowed) ){ 
      stop(paste0("merge must be one of the following: ", output_allowed))
    }
    
    # Get all GE data into list of dfs and combine
    if( is.null(EM_list) ){ EM_list <- GeneExpressionMatrices()$name }
    suppressMessages(lapply(EM_list, .self$downloadMatrix, summary = T)) 
    init_ems <- data_cache[ -which(names(data_cache) == "GE_matrices")]
    
    # Setup biosample2subject hash
    pheno <- unique(data.table(labkey.selectRows(
      baseUrl = config$labkey.url.base, 
      folderPath = config$labkey.url.path,
      schemaName = sn_study, 
      queryName = qn_InputSmplsShot,
      containerFilter = cf_currandsubs,
      colNameOpt = cn_caption)))
    
    pheno <- pheno[ which(pheno$`Expression Matrix Accession` %in% EM_list), ]
    subids <- paste0(pheno$`Participant ID`, "_d", pheno$`Study Time Collected`)
    subids <- gsub("-", "neg", subids, fixed = T)
    bs2id <- hash(pheno$`Biosample Accession`, subids)
    clean_ems <- lapply(init_ems, clean_EM, genes = genes, id_hash = bs2id, time_filt = ge_time)
    
    # Ensure EMs all contain same common genes
    genes_list <- lapply(clean_ems, FUN = function(em){ em$gene_symbol })
    common_genes <- Reduce(intersect, genes_list)
    if(length(common_genes) == 0){ stop("No common genes found") }
    filt_ems <- lapply(clean_ems, filterEM, common_genes = common_genes)
    
    
    # pull clinical data with filter for subs and filter for values col
    sub_filter <- makeFilter(c("participant_id", 
                               "IN", 
                               paste0(unique(pheno$`Participant ID`), collapse = ";")))
    
    cdata <- getDataset(clin_data, colFilter = sub_filter)
    if(!is.null(clin_time)){ cdata <- cdata[ which(cdata$study_time_collected == clin_time), ] }
    
    # for mapping EM names to study ID
    em_runs <- data.table(labkey.selectRows(
      baseUrl = config$labkey.url.base, 
      folderPath = config$labkey.url.path,
      schemaName = sn_assayExprMx, 
      queryName = qn_Runs))
    
    # return output in either one merged EM, by study, or by EM name
    res <- list()
    if(output == "all"){
      res$EM <- Reduce(merge, filt_ems)
      res$cData <- cdata
      
    }else if(output == "study"){
      studies <- unique(gsub(".*\\.", "", cdata$participant_id))
      res <- sapply(studies, simplify = F, USE.NAMES = T, FUN = function(sdy){
        tmp <- list()
        tmp$cdata <- cdata[ grep(paste0(".*", sdy), cdata$participant_id), ]
        em_nms <- em_runs$Name[ which(em_runs$Study == paste0("SDY", sdy)) ]
        em_nms <- paste0(em_nms, "_sum")
        em_dfs <- filt_ems[ which(names(filt_ems) %in% em_nms) ]
        tmp$EM <- Reduce(merge, em_dfs)
        return(tmp)
      })
      
    }else{
      em_nms <- em_runs$Name[ which(em_runs$Name %in% gsub("_sum","", names(filt_ems))) ]
      em_nms <- paste0(em_nms, "_sum")
      res <- sapply(em_nms, simplify = F, USE.NAMES = T, FUN = function(em_nm){
        tmp <- list()
        tmp$EM <- filt_ems[[em_nm]]
        tmp$cData <- cdata[ which(cdata$participant_id %in% colnames(tmp$EM)), ]
        return(tmp)
      })
    }
    
    return(res)
  }
)

