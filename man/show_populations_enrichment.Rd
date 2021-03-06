% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{show_populations_enrichment}
\alias{show_populations_enrichment}
\title{show_populations_enrichment}
\usage{
show_populations_enrichment(escores, trajs, gene_sets)
}
\arguments{
\item{escores}{The gene set by cell enrichment score matrix calculated by `run_enrichment()`}

\item{trajs}{A list of vectors of pseudotime values named by cells found in the columns of `escores`}

\item{gene_sets}{The names of gene sets to visualize (subset of the `rownames` of `escores`)}
}
\description{
Plots the enrichment of a multiple populations along pseudotime
}
