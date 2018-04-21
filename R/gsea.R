#'@title permutation test based gene set enrichment analysis given count of each terms and total gene numbers for the gene set and background.
#'@param freq_gs The named vector of term frequencies in gene set.
#'@param freq_bg The named vector of term frequencies in background.
#'@param gs_total_gene Total gene number in gene set.
#'@param bg_total_gene Total gene number in back ground.
#'@param adj_method the p adjust method.
#'@param ease To calculate EASE score or not.

gsea <- function(freq_gs,
                 freq_bg,
                 gs_total_gene,
                 bg_total_gene,
                 adj_method,
                 ease) {

  result <- data.frame(term = names(freq_gs),
                       freq_gs = as.numeric(freq_gs),
                       freq_bg = as.numeric(freq_bg[match(names(freq_gs),names(freq_bg))]))

  p_hyper <- phyper(result$freq_gs - ifelse(ease,2,1),
                    result$freq_bg,
                    bg_total_gene-result$freq_bg,
                    gs_total_gene,
                    lower.tail = FALSE,
                    log.p = FALSE)

  result$p <- p_hyper

  result[[ paste0("adj_", adj_method) ]] <- p.adjust(p_hyper,method = adj_method)

  result$OR <- round((result$freq_gs/gs_total_gene)/(result$freq_bg/bg_total_gene),2)

  result <- result[order(result$p),]

  row.names(result) = NULL

  return(result)
}
