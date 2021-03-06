.ISCon$methods(
  .munge=function(x){
    new <- tolower(gsub(" ","_",basename(x)))
    idx <- which(duplicated(new) | duplicated(new, fromLast = TRUE))
    if(length(idx)>0)
      new[idx] <- .munge(gsub("(.*)/.*$", "\\1", x[idx]))
    return(new)
  }
)

# Returns TRUE if the connection is at project level ("/Studies")
.ISCon$methods(
  .isProject=function()
    if(config$labkey.url.path == "/Studies/"){
      TRUE
    } else{
      FALSE
    }
)
.ISCon$methods(
  GeneExpressionInputs=function(){
    if(!is.null(data_cache[[constants$matrix_inputs]])){
      data_cache[[constants$matrix_inputs]]
    }else{
      ge<-data.table(labkey.selectRows(baseUrl = config$labkey.url.base,config$labkey.url.path,schemaName = "assay.ExpressionMatrix.matrix",queryName = "InputSamples",colNameOpt = "fieldname",viewName = "gene_expression_matrices",showHidden=TRUE))
      setnames(ge,.self$.munge(colnames(ge)))
      data_cache[[constants$matrix_inputs]]<<-ge
    }
  }
)


.ISCon$methods(
  .isRunningLocally=function(path){
    file.exists(path)
  }
)
.ISCon$methods(
  .localStudyPath=function(urlpath){
    LOCALPATH <- "/share/files/"
    PRODUCTION_HOST <- "www.immunespace.org"
    TEST_HOST <- "test.immunespace.org"
    gsub(file.path(gsub("/$","",config$labkey.url.base), "_webdav"), file.path(LOCALPATH), urlpath)
  }
)

.ISCon$methods(
    listDatasets=function(which = c("datasets", "expression")){
      "List the datasets available in the study or studies of the connection."
      
      if("datasets" %in% which){
        cat("datasets\n")
        for(i in 1:nrow(available_datasets)){
          cat(sprintf("\t%s\n",available_datasets[i,Name]))
        }
      }
      if("expression" %in% which){
        if(!is.null(data_cache[[constants$matrices]])){
          cat("Expression Matrices\n")
          for(i in 1:nrow(data_cache[[constants$matrices]])){
            cat(sprintf("\t%s\n",data_cache[[constants$matrices]][i, name]))
          }
        }
      }
    })


.ISCon$methods(
    listGEAnalysis = function(){
      "List available gene expression analysis for the connection."
      GEA <- data.table(labkey.selectRows(config$labkey.url.base,
                                          config$labkey.url.path,
                                          "gene_expression",
                                          "gene_expression_analysis",
                                          colNameOpt = "rname"))
      return(GEA)
    })

.ISCon$methods(
  getGEAnalysis = function(...){
    "Downloads data from the gene expression analysis results table.\n
    '...': A list of arguments to be passed to labkey.selectRows."
    GEAR <- data.table(labkey.selectRows(config$labkey.url.base, config$labkey.url.path,
        "gene_expression", "DGEA_filteredGEAR",  "DGEAR", colNameOpt = "caption", ...))
    setnames(GEAR, .self$.munge(colnames(GEAR)))
    return(GEAR)
  }
)

.ISCon$methods(
  clear_cache = function(){
  "Clear the data_cache. Remove downloaded datasets and expression matrices."
    data_cache[grep("^GE", names(data_cache), invert = TRUE)] <<- NULL
  }
)
.ISCon$methods(
    show=function(){
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

.ISCon$methods(
  getGEFiles=function(files, destdir = "."){
    "Download gene expression raw data files.\n
    files: A character. Filenames as shown on the gene_expression_files dataset.\n
    destdir: A character. The loacal path to store the downloaded files."
    links <- paste0(config$labkey.url.base, "/_webdav/",
                    config$labkey.url.path,
                    "/%40files/rawdata/gene_expression/", files)
    sapply(links, function(x){
      download.file(url = links[1], destfile = file.path(destdir, basename(x)),
                    method = "curl", extra = "-n")
    })
  }
)

# Returns a list of data frames where TRUE in file_exists column marks files that are accessible.
# This function is used for administrative purposes to check that the flat files
# are properly loaded and accessible to the users.
#' @importFrom RCurl getURL url.exists
#' @importFrom rjson fromJSON
#' @importFrom parallel mclapply detectCores
.ISCon$methods(
  .test_files=function(what = c("gene_expression_files", "fcs_sample_files", "protocols", "ge_matrices")){
    list_files <- function(link){
      response <- NULL
      if (url.exists(url = link, netrc = TRUE)){
        response_raw <- getURL(url = link, netrc = TRUE)
        response_json <- fromJSON(response_raw)
        response <- unlist(lapply(response_json$files, function(x) x$text))
      }
      response
    }
    
    check_links <- function (dataset, folder){
      res <- data.frame(file_info_name = NULL, study_accession = NULL, 
                        file_link = NULL, file_exists = NULL, 
                        stringsAsFactors = FALSE)
      
      if (dataset %in% .self$available_datasets$Name){
        temp <- .self$getDataset(dataset, original_view = TRUE)
        temp <- temp[!is.na(file_info_name)]
        temp <- unique(temp[, list(study_accession, file_info_name)])
        
        file_link <- paste0(config$labkey.url.base, "/_webdav/Studies/", 
                            temp$study_accession, "/%40files/rawdata/", folder, "/", 
                            sapply(temp$file_info_name, URLencode))
        
        studies <- unique(temp$study_accession)
        folder_link <- paste0(config$labkey.url.base, "/_webdav/Studies/", 
                              studies, "/%40files/rawdata/", folder, "?method=JSON")
        
        file_list <- unlist(mclapply(folder_link, list_files, mc.cores = detectCores()))
        
        file_exists <- temp$file_info_name %in% file_list
        res <- data.frame(study = temp$study_accession, file_link = file_link, file_exists = file_exists, 
                          stringsAsFactors = FALSE) 
        print(paste0(sum(res$file_exists), "/", nrow(res), " ", dataset, " with valid links."))
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
        folders_list <- labkey.getFolders(baseUrl = config$labkey.url.base, 
                                          folderPath = "/Studies/")
        folders <- folders_list[, 1]
        folders <- folders[!folders %in% c("SDY_template","Studies")]
      } else{
        folders <- basename(config$labkey.url.path)
      }
      
      file_link <- paste0(config$labkey.url.base, "/_webdav/Studies/", folders, 
                              "/%40files/protocols/", folders, "_protocol.zip")
      file_exists <- unlist(mclapply(file_link, url.exists, netrc = TRUE, mc.cores = detectCores()))
      
      res <- data.frame(study = folders, file_link = file_link, file_exists = file_exists, 
                        stringsAsFactors = FALSE)      
      print(paste0(sum(res$file_exists), "/", nrow(res), " protocols with valid links."))
      ret$protocols <- res
    }
    if ("ge_matrices" %in% what){
      matrix_queries <- labkey.getQueries(baseUrl = config$labkey.url.base, 
                                          folderPath = config$labkey.url.path, 
                                          schemaName = "assay.ExpressionMatrix.matrix")
      
      if ("OutputDatas" %in% matrix_queries$queryName) {
        ge<-data.frame(labkey.selectRows( baseUrl = config$labkey.url.base, folderPath = config$labkey.url.path,  
                        schemaName = "assay.ExpressionMatrix.matrix", queryName = "OutputDatas", colNameOpt = "rname", viewName = "links"))
        output <- lapply(ge[4], function(x) gsub("@", "%40", gsub("file:/share/files", 
                        paste0(config$labkey.url.base, "/_webdav"), x)))
        file_exists <- unlist(mclapply(output$data_datafileurl, url.exists, netrc = TRUE, mc.cores = detectCores()))
        res <- data.frame(file_link = output$data_datafileurl, file_exists = file_exists, 
                        stringsAsFactors = FALSE)
        print(paste0(sum(res$file_exists), "/", nrow(res), " ge_matrices with valid links."))
      } else {
        res <- data.frame(file_link = NULL, file_exists = NULL, stringsAsFactors = FALSE)
      }
      
      ret$ge_matrices <- res
    }
    return(ret)
  }
)


