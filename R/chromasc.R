get_es <- function(cell, gene_set) {
  transcriptome <- names(cell)
  all_genes_in_cell <- names(cell[cell>0])
  gs_genes_in_cell <- base::intersect(gene_set, all_genes_in_cell)
  gs_genes_frac <- length(gs_genes_in_cell)/length(all_genes_in_cell)
  return(gs_genes_frac)
}

get_z_es <- function(col_num, expr_mat, gene_set, null_params) {
  cell <- expr_mat[,col_num]
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
#' @param expr_mat An expression matrix
#' @param gene_sets A gene set object loaded by `load_gene_sets()`
run_enrichment <- function(expr_mat, gene_sets, null_params) {
  unsparse_expr <- as.matrix(expr_mat)
  res <- pbsapply(
        seq(expr_mat@Dim[2]),
        function(col_num) sapply(
          names(gene_sets),
          function(gs_name) get_z_es(col_num, gene_set=gene_sets[[gs_name]],
                                    expr_mat=unsparse_expr, null_params=null_params)
        ), simplify=T
    )
  colnames(res) <- expr_mat@Dimnames[[2]]
  return(res)
}
