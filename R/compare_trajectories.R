#' compare_populations
#'
#' @description Compare different populations along the same pseudotime trajectory
#'
#' @param pops1.pt A named vector of pseudotime values forming the first group to test
#' @param pops2.pt A named vector of pseudotime values forming the second group to test
#' @param escores The gene set by cell matrix of enrichment z-scores produced by `run_enrichment()`
#' @param gene_set_name The name of the gene set to test
#' @param region.start Optional - start of pseudotime region to compare (the smallest pseudotime value)
#' @param region.end Optional - end of pseudotime region to compare (the largest pseudotime value)
#' @param window.size Optional - if region.start==region.end==NA, then test consecutive windows of this
#' size and report statistics and corrected pvalues
#' @param step.size Optional - If scanning for significance, fraction of pseudotime to move the window
#' on each test
compare_populations <- function(pops1.pt, pops2.pt, escores, gene_set_name,
                                region.start=NA, region.end=NA,
                                window.size=NA, step.size=NA) {
  if (is.na(window.size)) {
    pops1.pt.subs <- which(pops1.pt >= region.start & pops1.pt <= region.end)
    pops2.pt.subs <- which(pops2.pt >= region.start & pops2.pt <= region.end)

    pops1.subs <- escores[gene_set_name, names(pops1.pt.subs)]
    pops2.subs <- escores[gene_set_name, names(pops2.pt.subs)]
    print(c(length(pops1.subs), length(pops2.subs)))
    return(wilcox.test(pops1.subs, pops2.subs)[c("statistic", "p.value")])
  } else {
    ## TODO
  }
}
