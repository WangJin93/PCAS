#' @title Get CPTAC expression data
#' @description
#' Get the mRNA/protein expression data in CPTAC database.
#' @import tibble reshape2
#' @param datasets Dataset names, you can input one or multiple datasets. Use 'dataset$Abbre' to get all datasets.
#' @param genes Gene symbols, you can input one or multiple symbols.
#' @examples
#' \dontrun{
#' results <- get_expr_data(datasets = "LUAD_CPTAC_mRNA", genes = c("GAPDH","TNS1"))
#' results <- get_expr_data(datasets = c("LUAD_CPTAC_protein","LSCC_CPTAC_protein"), genes = "GAPDH")
#' results <- get_expr_data(datasets = c("CCRCC_CPTAC_mRNA","GBM_CPTAC_mRNA","HNSCC_CPTAC_mRNA","LSCC_CPTAC_mRNA","LUAD_CPTAC_mRNA","PDAC_CPTAC_mRNA","UCEC_CPTAC2_mRNA","UCEC_CPTAC1_mRNA"), genes = c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA"))
#' }
#' @export
#'
get_expr_data <- function(datasets=c("LUAD_CPTAC_protein","LSCC_CPTAC_protein"), genes= c("TP53","TNS1")) {
  # 创建缓存目录（如果不存在）
  base_dir <- "data_files"
  action_dir <- file.path(base_dir, "data_temp")
  if (!dir.exists(base_dir)) {
    dir.create(base_dir)
  }
  if (!dir.exists(action_dir)) {
    dir.create(action_dir)
  }
  if (length(genes)==0){
    return(NULL)
  }else{
    message("Querying data of identifier ", paste0(genes,collapse = ", "), " from datasets ", paste0(datasets,collapse = ", "))
    cptac_data<- data.frame()
    for (x in datasets) {
      if (stringr::str_detect(x,"mRNA")) {
        id <- idmap_RNA[which(idmap_RNA$Symbol %in% genes),]
        ids <- id$row_names
      }else{
        ids <- genes
      }
      # 生成基因的哈希值作为文件名的一部分
      genes_hash <- digest::digest(ids, algo = "md5")
      cache_file <- file.path(action_dir, paste0(x, "_", genes_hash, ".RData"))

      # 检查缓存文件是否存在
      if (file.exists(cache_file)) {
        message("Loading cached data from ", cache_file)
        load(cache_file)
        cptac_data <- plyr::rbind.fill(cptac_data, cached_data)
        next()
      }
      data <- get_data(x,
                       "expression",
                       ids)

      if (is.null(nrow(data))){
        message(paste0(genes,collapse = ", "), " retrive no results in ", x)
        next()
      }else{
        if (stringr::str_detect(x,"mRNA")) {
          data <- merge(id[2:3],data,by="row_names")[-1]
          colnames(data)[1] <- "row_names"
        }
        row.names(data) <- NULL
        cached_data <- data %>% tibble::column_to_rownames(var = "row_names") %>% t() %>% as.data.frame()

        cached_data$type <- lapply(row.names(cached_data), function(x) strsplit(x,"_")[[1]][length(strsplit(x,"_")[[1]])])%>% as.character()
        cached_data$ID <- stringr::str_remove(row.names(cached_data),"_Tumor|_Normal|_Other")
        # data <- tibble::rownames_to_column(data,var = "ID")
        # data <- reshape2::melt(data, measure.vars = genes)
        # data$value <- as.numeric(data$value)
        cached_data <- cached_data[c("ID","type",colnames(cached_data)[1:(ncol(cached_data)-2)])]
        cached_data$dataset <- x
        # 保存到缓存
        save(cached_data, file = cache_file)
        cptac_data <- plyr::rbind.fill(cptac_data,cached_data)
        }
    }
    if (nrow(cptac_data) == 0){
      message("Retrive no data.")
      return(NULL)
    }
    if (length(genes)==1){
      cptac_data[,3] <- as.numeric(cptac_data[,3])
    }else{
      cptac_data[intersect(colnames(cptac_data),genes)] <- apply(cptac_data[intersect(colnames(cptac_data),genes)], 2, as.numeric)
    }
    cptac_data <- cptac_data[c("ID","type","dataset",intersect(colnames(cptac_data),genes))]
    return(cptac_data)
  }
}

