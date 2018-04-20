#' @title Vectorized function to enrich GO exact terms with a list of gene sets against a common background.
#' @return a \code{data.frame} for GO exact enrichment result.
#' @details \code{GOEA} conduct regular GOEA analysis with Fisher's exact test / permutation test.
#' @param gene_set a list of character vectors contains gene IDs of the query gene sets.
#' @param back_ground a character vector contains the IDs of the background genes.
#' @param orgDb an \code{OrgDb} object defined by AnnotationDbi.
#' @param category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param ID_type a character specifying the type of the ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param min_bg_count term minimum number of count in background gene set; default 1.
#' @param max_bg_count term maximum number of count in background gene set; default 1000.
#' @param EASE_Score whether or not use EASE score method. (a more conservative hypergeomatric test calculation used in DAVID)
#'for more details please refer to \url{https://david.ncifcrf.gov/helps/functional_annotation.html#fisher}, default FALSE.
#' @param pvalue_correction method used for multiple hypothesis adjustment, can be one in "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", and "none".
#' @param interpret_term whether to let the GO term readable, default FALSE.
#' @param show_gene_name whether to let the gene name readable, default FALSE.
#' @param attach_slim whether to include additional collumn indicating the slim terms for each GO, default FALSE.
#' @param exclude_root_slim whether to exclude the GO slim terms that labeling "BP", "CC", or "MF", default TRUE.
#' @export
GOEA <- function(gene_set,
                       back_ground,
                       orgDb,
                       category="BP",
                       ID_type = "ENTREZID",
                       min_bg_count = 1,
                       max_bg_count = 1000,
                       EASE_Score= F,
                       pvalue_correction = "BH",
                       interpret_term = F,
                       show_gene_name = F,
                       attach_slim = F,
                       exclude_root_slim = F,
                       exclude_zero_slim = T
                       )
{
  Mapp_result <- suppressMessages( select(orgDb, keys = as.character(back_ground), columns ="GO", keytype=ID_type) )
  Mapp_result <- na.omit(Mapp_result)
  idx_cat <- Mapp_result$ONTOLOGY == category
  Mapp_result <- Mapp_result[idx_cat,]
  Mapp_result <- Mapp_result[!duplicated( paste0(Mapp_result$SYMBOL,Mapp_result$GO) ),]
  GO_tb <- table(Mapp_result$GO)
  idx_bg <- GO_tb >= min_bg_count & GO_tb <= max_bg_count
  Mapp_result <- Mapp_result[Mapp_result$GO%in%(names(GO_tb)[idx_bg]),]


  Pop_bg <- table(Mapp_result$GO)
  Pop_ls <- table(Mapp_result$GO[Mapp_result[[ID_type]]%in%gene_set])

  result <- data.frame(term = names(Pop_ls),
                       count_ls = as.numeric(Pop_ls),
                       count_bg = as.numeric(Pop_bg[match(names(Pop_ls),names(Pop_bg))]))
  bg_number <- length(back_ground)
  gl_number <- length(gene_set)
  p_hyper <- phyper(result$count_ls - ifelse(EASE_Score,2,1),
                    result$count_bg,
                    bg_number-result$count_bg,
                    gl_number,
                    lower.tail = FALSE,
                    log.p = FALSE)
  result$p <- p_hyper
  result$FDR <- p.adjust(p_hyper,method = "BH")
  result$OR <- round((result$count_ls/gl_number)/(result$count_bg/bg_number),2)
  if(show_gene_name) {
    Map_sub <- Mapp_result[Mapp_result[[ID_type]]%in%gene_set,]
    gene_names <- suppressWarnings( select(orgDb, keys = Map_sub$SYMBOL, columns = c("GENENAME"), keytype="SYMBOL")$GENENAME )
    result$genes_onlist = split(gene_names,Map_sub$GO)[result$term]
  }
  if(interpret_term) { result$term = Term(GOTERM)[result$term] }
  rownames(result) = NULL
  result[order(result$p),]
}
