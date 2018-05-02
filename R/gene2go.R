#' @title A function that map Gene IDs to GO terms (including the ancestors).
#' @return a \code{data.frame} that define the mapping relationship between genes and GO terms.
#' @param Gene_ID gene ids to query.
#' @param Gene_key_type a character specifying the type of the ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param OrgDB an \code{OrgDb} object defined by AnnotationDbi.
#' @param Category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param Drop_zero_match whether to drop the genes that mapped to no terms at all, default FALSE.
#' @param Slim Wheather to return a mapping to GO slim terms (a certain subset to GO terms), default FALSE.
#' @param Slim_subset A character sting of GO terms that define the scope of GO Slim. if not provided, the GO slim would be the generic subset defined in : \url{http://geneontology.org/ontology/subsets/goslim_generic.obo}
#' @param Exclude_self_slim whether the GO slim terms of its own category i.e. remove terms of c("GO:0008150","GO:0005575","GO:0003674").
#'
#'
#' @importFrom AnnotationDbi select
#' @importFrom GO.db GOBPANCESTOR GOMFANCESTOR GOCCANCESTOR
#' @export
gene2go <- function(Gene_ID,
                   Gene_key_type = "ENTREZID",
                   OrgDB,
                   Category = "BP",
                   Drop_zero_match = F,
                   Slim = F,
                   Slim_subset = NULL,
                   Exclude_self_slim = T) {

  stopifnot(Category %in% c("BP","CC","MF"))

  Map_result <- suppressMessages( select(OrgDB, keys = as.character(Gene_ID), columns = "GO", keytype = Gene_key_type) )

  Map_result <- na.omit(Map_result) #The Genes that map to no GO ids are dropped

  Map_result <- Map_result[Map_result$ONTOLOGY == Category,]

  Map_result <- Map_result[!duplicated( paste0(Map_result[[Gene_key_type]],Map_result$GO) ),] #remove the duplicated GO terms for each gene ID

  Map_result <- Map_result[,colnames(Map_result) %in% c(Gene_key_type,"GO")]

  #Retrieve the ancestors

  GO_pergene <- split(Map_result$GO,Map_result$ENTREZID)

  ancest_lst <- as.list(eval( parse(text = paste0("GO",Category,"ANCESTOR")) ))

  goancest_pergene <- lapply(GO_pergene, function(x) {
   unique_ancest <- unique(unlist(ancest_lst[x],use.names = F))
   return(unique_ancest[unique_ancest!="all"])
  } )

  if(Slim) {

   if(is.null(Slim_subset)) {
     Slim_subset <- slim_generic
    }
   if(Exclude_self_slim) {
     Slim_subset <- Slim_subset[!Slim_subset %in% c("GO:0008150","GO:0005575","GO:0003674")]
    }
    goancest_pergene <- lapply(goancest_pergene, function(x) x[x%in%Slim_subset])
}

  Map_result <- data.frame(
                   ENTREZID = rep(names(goancest_pergene),elementNROWS(goancest_pergene)),
                   GO = unlist(goancest_pergene)
                  )


  if(!Drop_zero_match) {

  indx_missing = !Gene_ID %in% Map_result$ENTREZID

  if(sum(indx_missing) != 0) {

   df_attach <- data.frame( ENTREZID = Gene_ID[indx_missing],
                GO = NA)

   Map_result <- rbind(Map_result,df_attach)

  }

  }

  rownames(Map_result) = NULL

  return(Map_result)
}
