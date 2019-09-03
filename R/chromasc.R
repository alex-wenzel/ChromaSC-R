#' run_enrichment
#'
#' @description Hi
#'
#' @param seurat_obj A Seurat object
#' @param traj_vecs A list of vectors of continuous values named by a subset of seurat_obj@cell.names
#' @param gmt_path A path to a GMT file with one or more gene sets
#' @param window_size The fraction of pseudotime to use for the sliding window
#' @param step_size The fraction of pseudotime values to move the sliding window on each iteration

run_enrichment <- function(seurat_obj, traj_vecs, gmt_path, window_size, step_size) {
  ## Load gene sets
  gene_sets_full <- GSA::GSA.read.gmt(gmt_path)
  gene_sets <- gene_sets_full$genesets
  names(gene_sets) <- gene_sets_full$geneset.names

  ## Make sure trajectory vectors are sorted
  traj_vecs <- lapply(traj_vecs, function(l) l[order(l)])

  ## Build a list of empty results data.frames
  all_res <- replicate(length(traj_vecs), data.frame(row.names=names(gene_sets)))

  ## Calculate the null distribution (means and vars for gs and cell sizes)

}
