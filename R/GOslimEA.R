#' @title Vectorized function to enrich GO slim terms with a list of gene sets against a common background.
#' @return a \code{data.frame} for GO slim enrichment result.
#'
#' @param gene_set a list of character vectors contains gene IDs of the query gene sets.
#' @param back_ground a character vector contains the IDs of the background genes.
#' @param orgDb an \code{OrgDb} object defined by AnnotationDbi.
#' @param category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param ID_type a character specifying the type of the ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param min_bg_count term minimum number of count in background gene set; default 1.
#' @param max_bg_count term maximum number of count in background gene set; default 1000.
#' @param exclude_zero_slim whether to exclude the involvement of the genes that belong to no GO slim terms at all, default TRUE.
#' @param exclude_root_slim whether to drop the GO slim terms that are roots of "BP", "CC", or "MF" at the very first, default TRUE.
#' @param EASE_Score whether or not use EASE score method. (a more conservative hypergeomatric test calculation used in DAVID)
#'for more details please refer to \url{https://david.ncifcrf.gov/helps/functional_annotation.html#fisher}, default FALSE.
#' @param pvalue_correction method used for multiple hypothesis adjustment, can be one in "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", and "none".
#' @param interpret_term whether to let the GO term readable, default FALSE.
#' @param show_gene_name whether to let the gene name readable, default FALSE.
#' @param slim_collection an \code{OBOCollection} object define the mapping between GO exact and GO slim, if not provided, it will use the generic collection.
#'
#' @details \code{GOslimEA} conduct GO slim EA analysis using Fisher's exact test / permutation test on GO slim terms. GO slim terms refer to the reduced GO terms that are roots in GO hierachical structures.
#'
#' @export
GOslimEA <- function(gene_set,
                       back_ground,
                       orgDb,
                       category ="BP",
                       ID_type = "ENTREZID",
                       min_bg_count = 1,
                       max_bg_count = 1000,
                       exclude_root_slim = T,
                       exclude_zero_slim = T,
                       EASE_Score= F,
                       interpret_term = F,
                       show_gene_name = F,
                       slim_collection = NULL
                      ) {

}
