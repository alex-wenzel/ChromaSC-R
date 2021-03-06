% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_trajectories.R
\name{compare_populations}
\alias{compare_populations}
\title{compare_populations}
\usage{
compare_populations(pops1.pt, pops2.pt, escores, gene_set_name,
  region.start = NA, region.end = NA, window.size = NA,
  step.size = NA)
}
\arguments{
\item{pops1.pt}{A named vector of pseudotime values forming the first group to test}

\item{pops2.pt}{A named vector of pseudotime values forming the second group to test}

\item{escores}{The gene set by cell matrix of enrichment z-scores produced by `run_enrichment()`}

\item{gene_set_name}{The name of the gene set to test}

\item{region.start}{Optional - start of pseudotime region to compare (the smallest pseudotime value)}

\item{region.end}{Optional - end of pseudotime region to compare (the largest pseudotime value)}

\item{window.size}{Optional - if region.start==region.end==NA, then test consecutive windows of this
size and report statistics and corrected pvalues}

\item{step.size}{Optional - If scanning for significance, fraction of pseudotime to move the window
on each test}
}
\description{
Compare different populations along the same pseudotime trajectory
}
