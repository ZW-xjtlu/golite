#' @title Vectorized function to enrich GO slim terms with a list of gene sets against a common background.
#' @return a \code{list} of \code{data.frame} for GO slim enrichment result.
#'
#' @param gene_set a list of character vectors contains gene IDs of the query gene sets.
#' @param back_ground a character vector contains the IDs of the background genes.
#' @param orgDb an \code{OrgDb} object defined by AnnotationDbi.
#' @param category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param gene_key a character specifying the type of the gene ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param min_gs_count GO slim term minimum number of occurence in the query gene set; default 10.
#' @param max_gs_count GO slim term maximum number of occurence in the query gene set; default 500.
#' @param exclude_zero_slim whether to exclude the involvement of the genes that belong to no GO slim terms at all, default TRUE (increase power).
#' @param exclude_root_slim whether to drop the GO slim terms that are roots of "BP", "CC", or "MF" at the very first, default TRUE.
#' @param EASE_Score whether or not use EASE score method. (a more conservative hypergeomatric test calculation used in DAVID)
#'for more details please refer to \url{https://david.ncifcrf.gov/helps/functional_annotation.html#fisher}, default FALSE.
#' @param pvalue_correction method used for multiple hypothesis adjustment, can be one in "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", and "none".
#' @param interpret_term whether to let the GO term readable, default FALSE.
#' @param show_gene_name whether to attach readable gene names for each GO term, default FALSE.
#'
#' @details \code{GOslimEA} conduct GO slim EA analysis using Fisher's exact test / permutation test on GO slim terms. GO slim terms refer to the reduced GO terms that are roots in GO hierachical structures.
#'
#' @importFrom AnnotationDbi select Term
#' @importFrom GO.db GOTERM
#' @export
GOslimEA <- function(gene_set,
                       back_ground,
                       orgDb,
                       category ="BP",
                       gene_key = "ENTREZID",
                       min_gs_count = 10,
                       max_gs_count = 500,
                       exclude_root_slim = T,
                       exclude_zero_slim = T,
                       EASE_Score= F,
                       pvalue_correction = "BH",
                       interpret_term = F,
                       show_gene_name = F
                      ) {

  stopifnot(category %in% c("BP","CC","MF"))

  if(any(duplicated(back_ground))) {
    warning("back ground gene IDs contains duplicated terms, the duplicates are removed", call. = TRUE)
    back_ground = unique(back_ground)
  }

  if(class( gene_set ) == "character") { gene_set = list(gene_set) }

  if(any(sapply(gene_set, function(x) any(duplicated(x)) ))) {
    warning("gene set gene IDs contains duplicated terms, the duplicates are removed", call. = TRUE)
    gene_set = lapply(gene_set, function(x) unique(x))
  }

  GO_indx <- gene2goslim(back_ground,gene_key,orgDb,category,
                         Exclude_self_slim = exclude_root_slim,
                         Drop_zero_match = exclude_zero_slim)

  GO_tb <- table(GO_indx$GO)

  Freq_bg <- table( as.character( GO_indx$GO_slim ) )

  result_lst <- list()

  bg_genes_num <- length( unique(GO_indx[[gene_key]]) )

  for(i in 1:length(gene_set) ) {

    indx_match <- GO_indx[[gene_key]] %in% gene_set[[i]]

    Freq_gs <- table( as.character( GO_indx$GO_slim [ indx_match ] ) )

    Freq_gs <- Freq_gs[ Freq_gs >= min_gs_count & Freq_gs <= max_gs_count ]

    result_lst[[i]] <- gsea(
      freq_gs = Freq_gs,
      freq_bg = Freq_bg,
      gs_total_gene = length(unique(GO_indx[[gene_key]] [indx_match])),
      bg_total_gene = bg_genes_num,
      adj_method = pvalue_correction,
      ease = EASE_Score
    )

  }

  if(show_gene_name) {
    gene_names <- suppressWarnings( select(orgDb, keys = GO_indx[[gene_key]], columns = c("GENENAME"), keytype = gene_key) )
    gene_names <- gene_names [!duplicated( gene_names$ENTREZID ),"GENENAME"]
    gene_names <- split(gene_names,GO_indx$GO_slim)
    result_lst <- Map(function(x,y) {
      y$genes_names = sapply( gene_names[as.character( y$term )], paste0, collapse = ", ")
      return(y)
    }, gene_set,
    result_lst)
  }

  if(interpret_term) {
    term_defs <- Term(GOTERM)
    result_lst <- lapply(result_lst, function(x) {
      x$definition = term_defs[x$term]
      return(x[,c(1,7,2,3,4,5,6)])
    })
  }


  if(length(result_lst) == 1){
    return(result_lst[[1]])
  } else {
    names(result_lst) = names(gene_set)
    return( result_lst )
  }
}
