########################################################################
# --------------------  HELPER METHODS  -------------------------------#
#        must be declared  before the functions they help              #
########################################################################

.ISCon$methods(
  .munge = function(x){
    new <- tolower(gsub(" ","_",basename(x)))
    idx <- which(duplicated(new) | duplicated(new, fromLast = TRUE))
    if( length(idx) > 0 ){ new[idx] <- .munge(gsub("(.*)/.*$", "\\1", x[idx])) }
    return(new)
  }
)

# Returns TRUE if the connection is at project level ("/Studies")
.ISCon$methods(
  .isProject = function()
    if(config$labkey.url.path == "/Studies/"){
      TRUE
    }else{
      FALSE
    }
)


.ISCon$methods(
  .isRunningLocally = function(path){
    file.exists(path)
  }
)


.ISCon$methods(
  .localStudyPath = function(urlpath){
    LOCALPATH <- "/share/files/"
    PRODUCTION_HOST <- "www.immunespace.org"
    TEST_HOST <- "test.immunespace.org"
    gsub(file.path(gsub("/$","",config$labkey.url.base), "_webdav"), file.path(LOCALPATH), urlpath)
  }
)


#################################################################################
# -----------------------  LIST DATA METHODS -----------------------------------#
#################################################################################

.ISCon$methods(
  listDatasets = function(output = c("datasets", "expression")){
    "List the datasets available in the study or studies of the connection."
    
    if( !(all(output %in% c("datasets","expression"))) ){
      stop("arguments other than 'datasets' and 'expression' not accepted")
    }
    
    if( "datasets" %in% output ){
      cat("Datasets\n")
      for( i in 1:nrow(available_datasets) ){
        cat(sprintf("\t%s\n", available_datasets[i,Name]))
      }
    }
    
    if( "expression" %in% output ){
      if( !is.null(data_cache[[constants$matrices]]) ){
        cat("Expression Matrices\n")
        for( i in 1:nrow(data_cache[[constants$matrices]]) ){
          cat(sprintf("\t%s\n", data_cache[[constants$matrices]][i, name]))
        }
      }else{
        warning("No expression matrices available")
      }
    }
    
  }
)

#################################################################################
# -----------------  GENE EXPRESSION METHODS -----------------------------------#
#################################################################################

.ISCon$methods(
  getGEFiles = function(files, destdir = "."){
    message("Download gene expression raw data files.\n
    files: A character. Filenames as shown on the gene_expression_files dataset.\n
    destdir: A character. The local path to store the downloaded files.")
    links <- paste0(config$labkey.url.base, 
                    "/_webdav",
                    config$labkey.url.path,
                    "/%40files/rawdata/gene_expression/", 
                    files)
    sapply(links, function(x){
      download.file(url = links[1], destfile = file.path(destdir, basename(x)),
                    method = "curl", extra = "-n")
    })
  }
)


.ISCon$methods(
  listGEAnalysis = function(){
    message("List available gene expression analysis for the connection.")
    GEA <- tryCatch(
              data.table(labkey.selectRows(baseUrl = config$labkey.url.base,
                                           folderPath = config$labkey.url.path,
                                           schemaName = sn_geneexpr,
                                           queryName = qn_GEanalysis,
                                           colNameOpt = cn_rname)),
              error = function(e) return(e)
              )
    if( length(GEA$message > 0) ){ stop("Study does not have Gene Expression Analyses") }
    return(GEA)
  }
)


.ISCon$methods(
  getGEAnalysis = function(...){
    message("Downloads data from the gene expression analysis results table.\n
    '...': A list of arguments to be passed to labkey.selectRows.")
    GEAR <- tryCatch(
              data.table(labkey.selectRows(baseUrl = config$labkey.url.base, 
                                           folderPath = config$labkey.url.path,
                                           schemaName = sn_geneexpr, 
                                           queryName = qn_DGEAfilt,  
                                           viewName = vn_DGEAR, 
                                           colNameOpt = cn_caption, 
                                           ...)),
                     error = function(e) return(e)
                    )
    if( length(GEAR$message > 0) ){ stop("Study does not have Gene Expression Analyses") }
    setnames(GEAR, .self$.munge(colnames(GEAR)))
    return(GEAR)
  }
)


.ISCon$methods(
  GeneExpressionInputs = function(){
    if( .isProject() == F ){
      if(!is.null(data_cache[[constants$matrix_inputs]])){
        data_cache[[constants$matrix_inputs]]
      }else{
        ge <- tryCatch(
          data.table(labkey.selectRows(baseUrl = config$labkey.url.base,
                                       folderPath = config$labkey.url.path,
                                       schemaName = sn_assayExprMx,
                                       queryName = qn_InputSmpls,
                                       colNameOpt = cn_fieldname,
                                       viewName = vn_exprmxs,
                                       showHidden = TRUE)),
          error = function(e) stop("GE Inputs not available for study")
        )
        setnames(ge, .self$.munge(colnames(ge)))
        data_cache[[constants$matrix_inputs]] <<- ge
      }
    }else{
      stop("method cannot be run at project level - select individual study")
    }
  }
)


.ISCon$methods(
  GeneExpressionMatrices = function(verbose = FALSE){
    if(!is.null(data_cache[[constants$matrices]])){
      data_cache[[constants$matrices]]
    }else{
      if(verbose){
        ge <- try(data.table(
          labkey.selectRows(baseUrl = config$labkey.url.base,
                            config$labkey.url.path,
                            schemaName = sn_assayExprMx,
                            queryName = qn_Runs,
                            colNameOpt = cn_fieldname,
                            showHidden = TRUE,
                            viewName = vn_EM)),
          silent = TRUE)
      } else {
        suppressWarnings(
          ge <- try(data.table(
            labkey.selectRows(baseUrl = config$labkey.url.base,
                              config$labkey.url.path,
                              schemaName = sn_assayExprMx,
                              queryName = qn_Runs,
                              colNameOpt = cn_fieldname,
                              showHidden = TRUE,
                              viewName = vn_EM)),
            silent = TRUE)
        )
      }
      if(inherits(ge, "try-error") || nrow(ge) == 0){
        #No assay or no runs
        message("No gene expression data")
        data_cache[[constants$matrices]] <<- NULL
      } else{
        setnames(ge,.self$.munge(colnames(ge)))
        data_cache[[constants$matrices]]<<-ge
      }
    }
    return(data_cache[[constants$matrices]])
  }
)


###############################################################################
# -------------- PARTICIPANT FILTERING METHODS -------------------------------#
###############################################################################

.col_lookup <- function(iter, values, keys){
  results <- sapply(iter, FUN = function(x){ res <- values[which(keys == x)] })
}


.getLKtbl <- function(config, schema, query){
  df <- labkey.selectRows(baseUrl = config$labkey.url.base,
                          folderPath = config$labkey.url.path,
                          schemaName = schema,
                          queryName = query,
                          showHidden = TRUE)
}


.ISCon$methods(
  listParticipantGroups = function(){
    if(config$labkey.url.path != "/Studies/"){
      stop("labkey.url.path must be /Studies/. Use CreateConnection with all studies.")
    }
    
    pgrp <- .getLKtbl(config = .self$config, schema = sn_study, query = qn_partGrp)
    pcat <- .getLKtbl(config = .self$config, schema = sn_study, query = qn_partCat)
    pmap <- .getLKtbl(config = .self$config, schema = sn_study, query = qn_partGrpMap)
    user2grp <- .getLKtbl(config = .self$config, schema = sn_core, query = qn_users)
    
    result <- merge(pgrp, pcat, by.x = pgrp_cat_id, by.y = pcat_row_id)
    result <- data.frame(Group_ID = result$`Row Id`, 
                         Label = result$Label.x, 
                         Created = result$Created, 
                         Created_By = result$`Created By`,
                         stringsAsFactors = F)
    
    result$Created_By <- .col_lookup(iter = result$Created_By,
                                     values = user2grp$`Display Name`,
                                     keys = user2grp$`User Id`)
    
    subs <- data.frame(summarize(group_by(pmap, `Group Id`), numsubs = n() ))
    
    result$Subjects <- .col_lookup(iter = result$Group_ID,
                                   values = subs$numsubs,
                                   keys = subs$Group.Id)
    
    return(result)
  }
)

# NOTE: If more than 500 subjects, a con$getDataset() will fail with HTTP400 because 
# the string for the colFilter = filter passed in the GET request is too long.  Need
# to make special version of con$getDataset() to handle this?
.ISCon$methods(
  getParticipantData = function(groupId, dataType){
    allowedData <- c("hai",
                     "neut_ab_titer",
                     "gene_expression_files",
                     "demographics")
    
    if(!(dataType %in% allowedData)){
      stop("DataType must be in ", paste(allowedData, collapse = ", "))
    }
    
    dt <- tolower(dataType)
    sql <- paste0("SELECT ", dt, ".*, pmap.GroupId As groupId ",
                  "FROM ", dt,
                  " JOIN ParticipantGroupMap AS pmap ON ", 
                  dt, ".ParticipantId = pmap.ParticipantId",
                  " WHERE groupId = ", as.character(groupId))
    
    allData <- labkey.executeSql(baseUrl = .self$config$labkey.url.base, 
                                 folderPath = .self$config$labkey.url.path,
                                 schemaName = sn_study,
                                 sql = sql,
                                 colNameOpt = cn_fieldname)
    
    cols2rm <- c("Container", 
                 "lsid",
                 "ParticipantSequenceNum",
                 "sourcelsid",
                 "Created",
                 "CreatedBy",
                 "Modified",
                 "ModifiedBy",
                 "SequenceNum",
                 "workspace_id",
                 "Dataset",
                 "VisitRowId",
                 "Folder",
                 "_key",
                 "filesize",
                 "unique_id",
                 "QCState",
                 "dsrowid",
                 "date")
    
    filtData <- allData[ , !(colnames(allData) %in% cols2rm) ]
    
    return(filtData)
  }
)


###############################################################################
# ----------------  HOUSEKEEPING FUNCTIONS -----------------------------------#
###############################################################################

.ISCon$methods(
  clear_cache = function(){
  "Clear the data_cache. Remove downloaded datasets and expression matrices."
    data_cache[grep("^GE", names(data_cache), invert = TRUE)] <<- NULL
  }
)

.ISCon$methods(
    show = function(){
    "Display information about the object."
      cat(sprintf("Immunespace Connection to study %s\n",study))
      cat(sprintf("URL: %s\n",file.path(gsub("/$","",config$labkey.url.base),gsub("^/","",config$labkey.url.path))))
      cat(sprintf("User: %s\n",config$labkey.user.email))
      cat("Available datasets\n")
      for(i in 1:nrow(available_datasets)){
        cat(sprintf("\t%s\n",available_datasets[i,Name]))
      }
      if(!is.null(data_cache[[constants$matrices]])){
        cat("Expression Matrices\n")
        for(i in 1:nrow(data_cache[[constants$matrices]])){
          cat(sprintf("\t%s\n",data_cache[[constants$matrices]][i, name]))
        }
      }
    }
)


# Returns a list of data frames where TRUE in file_exists column marks files that are accessible.
# This function is used for administrative purposes to check that the flat files
# are properly loaded and accessible to the users.
#' @importFrom RCurl getURL url.exists
#' @importFrom rjson fromJSON
#' @importFrom parallel mclapply detectCores
.ISCon$methods(
  .test_files = function(what = c("gene_expression_files", "fcs_sample_files", "protocols")){
    list_files <- function(link){
      response <- NULL
      if (url.exists(url = link, netrc = TRUE)){
        response_raw <- getURL(url = link, netrc = TRUE)
        response_json <- fromJSON(response_raw)
        response <- unlist(lapply(response_json$files, function(x) x$text))
      }
      response
    }
    
    check_links <- function (name1, name2){
      res <- data.frame(file_info_name = NULL, study_accession = NULL, 
                        file_link = NULL, file_exists = NULL, 
                        stringsAsFactors = FALSE)
      
      if (name1 %in% studies$available_datasets$Name){
        temp <- .self$getDataset(name1, original_view = TRUE)
        temp <- temp[!is.na(file_info_name)]
        temp <- unique(temp[, list(study_accession, file_info_name)])
        
        file_link <- paste0(.self$config$labkey.url.base, "/_webdav/Studies/", 
                            temp$study_accession, "/%40files/rawdata/", name2, "/", 
                            sapply(temp$file_info_name, URLencode))
        
        studies <- unique(temp$study_accession)
        folder_link <- paste0(.self$config$labkey.url.base, "/_webdav/Studies/", 
                              studies, "/%40files/rawdata/", name2, "?method=JSON")
        
        file_list <- unlist(mclapply(folder_link, list_files, mc.cores = detectCores()))
        
        file_exists <- temp$file_info_name %in% file_list
        res <- data.frame(study = temp$study_accession, file_link = file_link, file_exists = file_exists, 
                          stringsAsFactors = FALSE) 
        print(paste0(sum(res$file_exists), "/", nrow(res), " ", name1, " with valid links."))
      }
      
      res
    }
    
    ret <- list()
    what <- tolower(what)
    
    if("gene_expression_files" %in% what){
      ret$gene_expression_files <- check_links("gene_expression_files", "gene_expression")
    }
    if("fcs_sample_files" %in% what){
      ret$fcs_sample_files <- check_links("fcs_sample_files", "flow_cytometry")
    }
    if("protocols" %in% what){
      if(.self$.isProject()){
        folders_list <- labkey.getFolders(baseUrl = .self$config$labkey.url.base, 
                                          folderPath = "/Studies/")
        folders <- folders_list[, 1]
        folders <- folders[!folders %in% c("SDY_template","Studies")]
      } else{
        folders <- basename(.self$config$labkey.url.path)
      }
      
      file_link <- paste0(.self$config$labkey.url.base, "/_webdav/Studies/", folders, 
                              "/%40files/protocols/", folders, "_protocol.zip")
      file_exists <- unlist(mclapply(file_link, url.exists, netrc = TRUE, mc.cores = detectCores()))
      
      res <- data.frame(study = folders, file_link = file_link, file_exists = file_exists, 
                        stringsAsFactors = FALSE)      
      print(paste0(sum(res$file_exists), "/", nrow(res), " protocols with valid links."))
      ret$protocols <- res
    }
    return(ret)
  }
)

#################################################################################
# -----------------------  INITIALIZE METHODS ----------------------------------#
#################################################################################

#' @importFrom gtools mixedsort
#' @importFrom dplyr summarize group_by
.ISCon$methods(
  checkStudy = function(verbose = FALSE){
    validStudies <- mixedsort(grep("^SDY", 
                                   basename(lsFolders(getSession(config$labkey.url.base, "Studies"))), 
                                   value = TRUE))
    req_study <- basename(config$labkey.url.path)
    if( !req_study %in% c("", validStudies) ){
      if( !verbose ){
        stop(paste0(req_study, " is not a valid study"))
      } else{
        stop(paste0(req_study, " is not a valid study\nValid studies: ",
                    paste(validStudies, collapse = ", ")))
      }
    }
  }
)

.ISCon$methods(
  getAvailableDataSets = function(){
    if( length(available_datasets) == 0 ){
      df <- labkey.selectRows(baseUrl = config$labkey.url.base,
                              folderPath = config$labkey.url.path,
                              schemaName = sn_study,
                              queryName =  qn_ISC)
      available_datasets <<- data.table(df) # [,list(Label,Name,Description,`Key Property Name`)]
    }
  }
)


.ISCon$methods(
  initialize = function(..., config = NULL){
    
    #invoke the default init routine in case it needs to be invoked 
    #(e.g. when using $new(object) to construct the new object based on the exiting object)
    callSuper(...)
    
    constants <<- list(matrices = "GE_matrices", matrix_inputs = "GE_inputs")
    
    if(!is.null(config))
      config <<- config
    
    study <<- basename(config$labkey.url.path)
    if(config$verbose){
      checkStudy(config$verbose)
    }
    
    getAvailableDataSets()
    
    gematrices_success <- .self$GeneExpressionMatrices(verbose = FALSE)
    
  }
)


