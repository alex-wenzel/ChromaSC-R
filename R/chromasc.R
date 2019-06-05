#' ChromaSC
#'
#' @param seurat_obj A Seurat object
#' @param traj_vec A vector of continuous values named by a subset of seurat_obj@cell.names
#' @param window_size A double in (0.0, 1.0) giving the fraction of the trajectory to use for each window
#' @param step_size A double in (0.0, 1.0) giving the fraction the window should be shifted each iteration
#' @param gmt_path A path to a GMT file with one or more gene sets
#'
#' @return A gene-sets-by-windows data.frame with the enrichment values for each gene set across all windows
ChromaSC <- function(seurat_obj, traj_vec, window_size, step_size, gmt_path) {
  ## Load gene sets
  gene_sets_full <- GSA.read.gmt(gmt_path)
  gene_sets <- gene.sets.full$genesets
  names(gene_sets) <- gene_sets_full$geneset.names
}
