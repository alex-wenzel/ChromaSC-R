#' show_populations_enrichment
#'
#' @description Plots the enrichment of a multiple populations along pseudotime
#'
#' @param escores The gene set by cell enrichment score matrix calculated by `run_enrichment()`
#' @param trajs A list of vectors of pseudotime values named by cells found in the columns of `escores`
#' @param gene_sets The names of gene sets to visualize (subset of the `rownames` of `escores`)
show_populations_enrichment <- function(escores, trajs, gene_sets) {
  flat_trajs <- as.data.frame(unlist(trajs))
  colnames(flat_trajs) <- c("pseudotime")
  lineages <- sapply(rownames(flat_trajs), function(s)strsplit(s, ".", fixed=T)[[1]][1])
  cell_names <- sapply(rownames(flat_trajs), function(s)strsplit(s, ".", fixed=T)[[1]][2])
  names(lineages) <- cell_names
  names(cell_names) <- cell_names
  rownames(flat_trajs) <- cell_names

  for (gene_set_name in gene_sets) {
    print(gene_set_name)
    plot_gs_enrichment(escores, gene_set_name, flat_trajs, lineages)
  }
}

plot_gs_enrichment <- function(escores, gene_set_name, flat_trajs, lineages) {
  enrichment_melted <- melt(data = escores[gene_set_name,])
  enrichment_melted$pseudotime <- sapply(rownames(enrichment_melted), function(cid)flat_trajs[cid, "pseudotime"])
  enrichment_melted$lineage <- sapply(rownames(enrichment_melted), function(cid)lineages[cid])
  enrichment_melted <- enrichment_melted[complete.cases(enrichment_melted),]

  print(ggplot(data = enrichment_melted, aes(x = pseudotime, y = value, group = lineage,
                                             colour = lineage, fill = lineage)) +
          geom_smooth(method = 'loess', level = 0.95) +
          ggtitle(gene_set_name) +
          ylab("z-Score") +
          theme(
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            plot.title = element_text(size = 20, hjust = 0.5),
            legend.title = element_text(size=16),
            legend.text = element_text(size=14)
          )
  )
}
