get.random.cells <- function(n_cells, n_genes_per_cell, frac_expressed) {
  replicate(
    n_cells,
    sample(names(frac_expressed), size=n_genes_per_cell)
  )
}

get.es.2 <- function(cell, transcriptome, gene_set) {
  gs_genes_in_cell <- base::intersect(gene_set, cell)
  gs_genes_frac <- length(gs_genes_in_cell)/length(cell)
  return(gs_genes_frac)
}

#' compute_null_params
#'
#' @description This function computes the null distribution parameters for the input dataset. Use
#' `save_null_params()` and `load_null_params()` to skip this step for future analyses with the same
#' gene sets and data.
#'
#' @param seurat_obj A seurat object
#' @param gene_sets Gene sets loaded with `load_gene_sets()`
#' @param n_random_cells The number of random cells to generate to calculate the null parameters
#' for each combination of cell size and gene set
compute_null_params <- function(seurat_obj, gene_sets, n_random_cells=50) {
  ## Get a vector of what fraction of cells express each gene
  frac_expressed <- apply(seurat_obj@data, 1, function(r)sum(r!=0))/seurat_obj@raw.data@Dim[2]

  ## Get all unique sizes of all cells
  cell_sizes = sort(unique(seurat_obj@meta.data$nGene))

  ## Get all unique sizes of gene sets
  gene_set_sizes = sort(unique(sapply(gene_sets, length)))

  ## Initialize result matrices
  null_means <- matrix(data=NA, nrow=length(cell_sizes), ncol=length(gene_set_sizes))
  null_vars <- matrix(data=NA, nrow=length(cell_sizes), ncol=length(gene_set_sizes))

  ## Build a progress bar for the null calculations specifically
  prog_bar <- txtProgressBar(min=0, max=length(cell_sizes), style=3)

  ## Compute null means and variance
  for (i  in 1:length(cell_sizes)) {
    rand_cells <- get.random.cells(n_random_cells, cell_sizes[i], frac_expressed)
    for (j in 1:length(gene_set_sizes)) {
      #rand_gene_set <- get.random.cell(names(frac_expressed), gene_set_sizes[j], frac_expressed)
      rand_gene_set <- get.random.cells(1, gene_set_sizes[j], frac_expressed)[,1]
      escores <- apply(rand_cells, MARGIN=2, FUN=get.es.2,
                       transcriptome=names(frac_expressed),
                       gene_set=rand_gene_set)
      null_means[i, j] <- mean(escores)
      null_vars[i, j] <- var(escores)
    }
    setTxtProgressBar(prog_bar, i)
  }
  rownames(null_means) <- cell_sizes
  rownames(null_vars) <- cell_sizes
  colnames(null_means) <- gene_set_sizes
  colnames(null_vars) <- gene_set_sizes

  return(list("null_means"=null_means, "null_vars"=null_vars))
}

#' save_null_params
#'
#' @description This function saves the result of `compute_null_params()`
save_null_params <- function(null_params, path) {

}

#' load_null_params
#'
load_null_params <- function(path) {

}
