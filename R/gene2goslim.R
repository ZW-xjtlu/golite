#' @title A function that map Gene IDs to GO slim terms.
#' @return a \code{data.frame} that define the mapping relationship.
#' @param Gene_ID gene ids to query.
#' @param Gene_key_type a character specifying the type of the ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param OrgDB an \code{OrgDb} object defined by AnnotationDbi.
#' @param Category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param Drop_zero_match whether to drop the genes that mapped to no terms at all, default FALSE.
#' @param Exclude_self_slim whether the GO slim terms of its own category i.e. remove terms of c("GO:0008150","GO:0005575","GO:0003674").
#' @importFrom AnnotationDbi select
#' @importFrom dplyr left_join
#' @export
gene2goslim <- function(Gene_ID,
                    Gene_key_type = "ENTREZID",
                    OrgDB,
                    Category = "BP",
                    Drop_zero_match = F,
                    Exclude_self_slim = F) {

  stopifnot(Category %in% c("BP","CC","MF"))

  if(species( OrgDB ) != "Homo sapiens") {
    stop("Currently don't support species other than Homo sapiens. ",call. = F)
  }

  converted_id  <- suppressMessages( AnnotationDbi::select(OrgDB, keys = as.character(Gene_ID), columns = "ENTREZID", keytype = Gene_key_type) )

  etable_goslim <- eval( parse(text = paste0("GOslim_gene_lst_",Category)) )

  Map_result <- suppressWarnings( left_join( converted_id, etable_goslim, by = "ENTREZID") )

  Map_result <- Map_result[,colnames(Map_result) %in% c(Gene_key_type,"GO_slim")]

  Map_result <- Map_result[!duplicated( paste0(Map_result[[Gene_key_type]],Map_result$GO_slim) ),] #drop the potential duplicates in GO slim ids given gene.

  if(Exclude_self_slim) {
    Map_result$GO_slim[Map_result$GO_slim %in% c("GO:0008150","GO:0005575","GO:0003674")] = NA
  }

  if(Drop_zero_match) {
    Map_result <- na.omit(Map_result)
  }

  rownames(Map_result) = NULL

  return(Map_result)
}
