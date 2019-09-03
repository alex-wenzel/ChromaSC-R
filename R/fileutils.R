#' load_gene_sets
#'
#' @description This function takes a path to a GMT format file and returns a named list of
#' character vectors corresponding to each gene set.
#' For more on the GMT format, see the
#' \href{http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT}{GenePattern documentation}
#'
#' @param path A path to a .gmt
load_gene_sets <- function(path) {
  gene_sets <- read.table(path, sep='\t', stringsAsFactors=F, fill=T, header=F)
  gene_set_names <- gene_sets[,1]
  gene_sets <- gene_sets[,3:dim(gene_sets)[2]]
  gene_sets <- apply(gene_sets, MARGIN=1, FUN=function(row) as.character(row[row!=""]))
  names(gene_sets) <- gene_set_names
  return(gene_sets)
}
