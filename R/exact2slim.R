#' @title Convert a set of GO exact terms to GO slim terms.
#' @param Exact_terms vector of query GO exact terms.
#' @param Slim_collection the \code{OBOCollection} object, if not provided, it will use the generic GO slim collection downloaded from \url{"http://geneontology.org/ontology/subsets/goslim_generic.obo"}.
#' @param Category a character specifying the gene ontology category, can be one in "BP", "CC", and "MF", default "BP".
#' @param Exclude_self_slim whether to remove the slim terms of the Category it self, default FALSE.
#' @param Unique wheather to unique the mapped GO slim term of each GO exact term.
#' @return a \code{data.frame} indicating the mapping index between query GO exact terms and its GO slim terms.
#'
#' @importFrom goSlim GSEABase
#' @importFrom GOCollection GSEABase
#' @export
exact2slim <- function(Exact_terms,
                       Slim_collection = NULL,
                       Category = "BP",
                       Exclude_self_slim = F,
                       Unique = T) {

  if(is.na(Exact_terms)) return(NA)

  if(is.null(Slim_collection)) Slim_collection = slim_generic

  goSlim_all <- goSlim(GOCollection(Exact_terms[!is.na(Exact_terms)]),
                       Slim_collection,
                       Category)

  if(Exclude_self_slim) {
  goSlim_all <- goSlim_all[!rownames(  goSlim_all ) %in% c("GO:0008150",
                                                           "GO:0005575",
                                                           "GO:0003674"),]
  }

  if(!Unique){
   Freq_vec  <- goSlim_all$Count
   names(Freq_vec) = rownames(goSlim_all)
   return(Freq_vec)
  } else {
    return(rownames(goSlim_all)[goSlim_all$Count > 0])
  }
}
