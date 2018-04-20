#' @importFrom goSlim GSEABase
#' @export
exact2slim <- function(exact_terms = "",
                       slim_collection = NULL,
                       category = "BP",
                       Exclude_root_slim = F) {

  goSlim_all <- goSlim(exact_terms, slim_generic, category)

}
