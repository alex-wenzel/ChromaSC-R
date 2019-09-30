get_es <- function(cell, gene_set) {
  transcriptome <- names(cell)
  all_genes_in_cell <- names(cell[cell>0])
  gs_genes_in_cell <- base::intersect(gene_set, all_genes_in_cell)
  gs_genes_frac <- length(gs_genes_in_cell)/length(all_genes_in_cell)
  return(gs_genes_frac)
}

get_z_es <- function(cell_name, seurat_obj, gene_set, null_params) {
  cell <- seurat_obj@data[,cell_name]
  cell_size <- as.character(length(names(cell[cell>0])))
  gene_set_size <- as.character(length(gene_set))

  escore <- get_es(cell, gene_set)

  if(!cell_size %in% row.names(null_params$null_means)) {
    return(0)
  }

  nm <- null_params$null_means[cell_size, gene_set_size]
  nv <- null_params$null_vars[cell_size, gene_set_size]

  z_escore <- (escore-nm)/sqrt(nv)
  return(z_escore)
}

#' run_enrichment
#'
#' @description Hi
#'
#' @param seurat_obj A Seurat object
#' @param gene_sets A gene set object loaded by `load_gene_sets()`

run_enrichment <- function(seurat_obj, gene_sets, null_params) {
  return(
    pbsapply(
      seurat_obj@cell.names,
      function(cell.name) sapply(
        names(gene_sets),
        function(gs_name) get_z_es(cell.name, gene_set=gene_sets[[gs_name]],
                                   seurat_obj=seurat_obj, null_params=null_params)
      ), simplify=T
    )
  )
}
