
source("R/Schema_variables.R")

# helper functions

filterGenes <- function(df, genes){
  df <- df[ which(df[["gene_symbol"]] == genes), ]
}


genesAsrows <- function(df){
  rownames(df) <- df[["gene_symbol"]]
  df <- df[ , -which(colnames(df) == "gene_symbol")]
  return(df)
}

hashmap <- function(hsh, input){
  values <- unname(unlist(sapply(input, FUN = function(x){
    val <- hsh[[x]]
  })))
}

cols_map <- function(df, hsh){
  colnames(df) <- hashmap(hsh = hsh, input = colnames(df))
}



# Clean one study 





# expressionSet Maker
#' @importFrom hash hash
.ISCon$methods(
  get_expressionSet = function(clinical_data, genes = NULL, EM_list = NULL, summary = T, complete = T){
    # get participants from desired expression_matrices
    
    if( is.null(EM_list) ){
      EM_list <- GeneExpressionMatrices()$name
    }
    
    lapply(EM_list, downloadMatrix, summary) 
    
    res_dfs <- data_cache[ -which(names(data_cache) == "GE_matrices")]
    
    if( !is.null(genes) ){
      lapply(res_dfs, filterGenes, genes)
    }
    
    lapply(res_dfs, genesAsrows)
    
    pheno <- unique(data.table(labkey.selectRows(
      baseUrl = con$config$labkey.url.base, 
      folderPath = con$config$labkey.url.path,
      schemaName = sn_study, 
      queryName = qn_InputSmplsShot,
      containerFilter = cf_currandsubs,
      colNameOpt = cn_caption)))
    
    bs2id <- hash(pheno$`Biosample Accession`, pheno$`Participant ID`)
    
    lapply(res_dfs, cols_map, hsh) # now, EM has rownames of Gene_symbols, cols of subids, and only vals
    
      # metaEM <- merge all EM by gene_symbol
      # remove incomplete rows (if complete = T)
    
    # make filter with subs (colnames(metaEM))
    
    # pull clinical data with filter for subs
    
    # create eSet
      # exprs <- metaEM
      # annotatedDataFrame <- rownames(metaEM)
      # pData <- data.frame(clinical_data)
  }
)