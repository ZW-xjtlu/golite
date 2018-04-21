#' @title A function that map Gene IDs to GO terms.
#' @return a \code{data.frame} that define the mapping relationship.
#' @param Gene_ID gene ids to query.
#' @param Gene_key_type a character specifying the type of the ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param OrgDB an \code{OrgDb} object defined by AnnotationDbi.
#' @param Category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param Drop_zero_match whether to drop the genes that mapped to no terms at all, default FALSE.
#' @importFrom select AnnotationDbi
#' @export
gene2go <- function(Gene_ID,
                   Gene_key_type = "ENTREZID",
                   OrgDB,
                   Category = "BP",
                   Drop_zero_match = F) {

  Map_result <- suppressMessages( select(OrgDB, keys = as.character(Gene_ID), columns = "GO", keytype = Gene_key_type) )

  Map_result <- na.omit(Map_result) #The Genes that map to no GO ids are dropped

  Map_result <- Map_result[Map_result$ONTOLOGY == Category,]

  Map_result <- Map_result[!duplicated( paste0(Map_result[[Gene_key_type]],Map_result$GO) ),] #remove the duplicated GO terms for each gene ID

  Map_result <- Map_result[,colnames(Map_result) %in% c(Gene_key_type,"GO")]

  if(!Drop_zero_match) {

   df_attach <- data.frame( ENTREZID = Gene_ID[!Gene_ID %in% Map_result$ENTREZID],
                GO = NA)

   Map_result <- rbind(Map_result,df_attach)

  }

  rownames(Map_result) = NULL

  return(Map_result)
}
