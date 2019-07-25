# Utilities/helper functions (not exposed) ---------------

get_window_cells <- function(traj, interval_start, interval_end) {
  filt_traj <- traj[traj>interval_start & traj<interval_end]
  return(filt_traj[order(filt_traj)])
}

get_n_cells_window <- function(traj_vecs, window_starts, window_size) {
  n_cells_window <- c()
  for (ws in window_starts) {
    we <- ws + window_size
    trajs_windows <- lapply(traj_vecs, function(v) get_window_cells(v, ws, we))
    trajs_window_sizes <- lapply(trajs_windows, length)
    n_cells_window <- c(n_cells_window, min(unlist(trajs_window_sizes)))
  }
  return(n_cells_window)
}

get_downsampled_cells <- function(cells, min_cells) {
  select_vec <- c(rep(TRUE, min_cells), rep(FALSE, length(cells)-min_cells))
  select_vec <- sample(select_vec)
  chosen_cells <- cells[select_vec]
  return(chosen_cells[order(chosen_cells)])
}

get_effect_from_window <- function(seur_obj, traj_cells, gs, n_avail_genes) {
  ## Identify gene set genes overlapping with expression matrix genes
  avail_gs_genes <- intersect(rownames(seur_obj@data), unlist(gs))

  ## Check if enough genes available to do enrichment
  if (length(avail_gs_genes) < 10) { ## TODO: Parameterize this cutoff
    return(0)
  }

  ## Extract gene expression for gene set genes
  cells_mat <- seur_obj@data[avail_gs_genes, names(traj_cells)]

  ## Count total # of genes that are expressed in this window
  all_ngene_obj <- seur_obj@meta.data$nGene
  names(all_ngene_obj) <- seur_obj@cell.names
  n_genes_traj <- sum(all_ngene_obj[names(traj_cells)])

  ## Calculate expected fraction of gene set genes expressed
  exp_frac <- length(avail_gs_genes)/n_avail_genes
  traj_n_exp <- n_genes_traj * exp_frac

  traj_rowsums <- rowSums(cells_mat>0)
  traj_gs_count <- sum(traj_rowsums)

  traj_ratio <- traj_gs_count/traj_n_exp

  return(log2(traj_ratio))
}

# Driver/exposed function --------------------------------

#' ChromaSC
#'
#' @param seurat_obj A Seurat object
#' @param traj_vecs A list of vectors of continuous values named by a subset of seurat_obj@cell.names
#' @param window_size A double in (0.0, 1.0) giving the fraction of the trajectory to use for each window
#' @param step_size A double in (0.0, 1.0) giving the fraction the window should be shifted each iteration
#' @param gmt_path A path to a GMT file with one or more gene sets
#'
#' @return A gene-sets-by-windows data.frame with the enrichment values for each gene set across all windows
ChromaSC <- function(seurat_obj, traj_vecs, window_size, step_size, gmt_path) {
  ## Load gene sets
  gene_sets_full <- GSA::GSA.read.gmt(gmt_path)
  gene_sets <- gene_sets_full$genesets
  names(gene_sets) <- gene_sets_full$geneset.names

  ## Make sure trajectory vectors are sorted
  traj_vecs <- lapply(traj_vecs, function(l) l[order(l)])

  ## Build a list of empty results data.frames
  all_res <- replicate(length(traj_vecs), data.frame(row.names=names(gene_sets)))

  ## Identify starting locations of windows
  window_starts <- seq(0.0, 1.0-window_size, step_size)
  print(paste("Beginning", length(window_starts), "window enrichment calculations..."))

  ## Calculate the number of cells in each window
  n_cells_windows <- get_n_cells_window(traj_vecs, window_starts, window_size)

  ## Build progress bar
  prog_bar <- txtProgressBar(min=0, max=length(window_starts), style=3)

  ## Calculate available genes
  n_avail_genes <- length(rownames(seurat_obj@data))

  ## Loop over window starts
  for (i in 1:length(window_starts)) {
    ## Set window start
    ws <- window_starts[i]

    ## Set window end
    we <- ws + window_size

    ## Get cells in each population within the window
    traj_cells <- lapply(traj_vecs, function(v) get_window_cells(v, ws, we))
    traj_cells_lens <- lapply(traj_cells, length)

    ## Skip window if too few cells (TODO: add min num of cells parameter)
    if (min(unlist(traj_cells_lens)) < 2) {
      ## Populate enrichment result data.frames with 0s
      for (j in 1:length(traj_vecs)) {
        all_res[[j]] <- cbind(all_res[[j]], rep(0, length(gene_sets)))
      }
      setTxtProgressBar(prog_bar, i)
      next()
    }

    ## Downsample to smallest # of cells in a population in this window
    trajs_ds <- lapply(traj_cells, function(v) get_downsampled_cells(v, n_cells_windows[i]))

    ## Make empty results lists
    window_result <- rep(list(list()), length(trajs_ds))

    ## Calculate enrichment scores in window
    for (gs in gene_sets) {
      eff_sizes <- lapply(trajs_ds, function(p) get_effect_from_window(seurat_obj, p, gs, n_avail_genes))
      for (j in 1:length(trajs_ds)) {
        window_result[[j]] <- c(window_result[[j]], eff_sizes[[j]])
      }
    }

    ## Add to all results
    for (j in 1:length(trajs_ds)) {
      all_res[[j]] <- cbind(all_res[[j]], unlist(window_result[[j]]))
    }

    ## Update progress bar
    setTxtProgressBar(prog_bar, i)
  }

  ## Post-process results
  for (j in 1:length(traj_vecs)) {
    colnames(all_res[[j]]) <- paste("eff", 1:length(window_starts), sep="_")
  }

  #return(all_res)
  return(
    list(
      enrichment=all_res,
      winsizes=n_cells_windows
    )
  )
}
