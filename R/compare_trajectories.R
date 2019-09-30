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
    return(test_one_window(pops1.pt, pops2.pt, escores, gene_set_name,
                           region.start, region.end))
  } else {
    ## Define windows for scanning
    windows <- data.frame(
      window.starts = seq(0.0, 1.0-window.size, step.size),
      window.ends = seq(window.size, 1.0, step.size)
    )
    ## Function for testing within a window
    tests_res <- apply(windows, MARGIN=1,
                       function(row)
                         test_one_window(pops1.pt, pops2.pt, escores, gene_set_name,
                                         row[1], row[2]))
    res.df <- as.data.frame(t(matrix(unlist(tests_res), nrow=2)))
    colnames(res.df) <- c("W", "p.value")
    return(res.df)
  }
}

test_one_window <- function(pops1.pt, pops2.pt, escores, gene_set_name,
                            region.start, region.end) {
  ## Subset pseudotime within queried region
  pops1.pt.subset <- which(pops1.pt >= region.start & pops1.pt <= region.end)
  pops2.pt.subset <- which(pops2.pt >= region.start & pops2.pt <= region.end)

  ## Extract the z scores for cells within the queried region for
  ## the queried gene set
  pops1.subset <- escores[gene_set_name, names(pops1.pt.subset)]
  pops2.subset <- escores[gene_set_name, names(pops2.pt.subset)]

  ## Return the MW-U test result
  return(wilcox.test(pops1.subset, pops2.subset)[c("statistic", "p.value")])
}
