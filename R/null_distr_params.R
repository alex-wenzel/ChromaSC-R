#' compute_null_params
#'
#' @description This function computes the null distribution parameters for the input dataset. Use
#' `save_null_params()` and `load_null_params()` to skip this step for future analyses with the same
#' gene sets and data.
#'
#' @param seurat_obj A seurat object
#' @param gene_sets Gene sets loaded with `load_gene_sets()`
compute_null_params <- function(seurat_obj, gene_sets) {

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
