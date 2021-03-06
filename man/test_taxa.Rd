% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_taxa.R
\name{test_taxa}
\alias{test_taxa}
\title{Test a taxonomical group of features}
\usage{
test_taxa(
  ps,
  rank = "Phylum",
  time_variable = "TimePoint",
  treatment_variable = "Treatment",
  sample_variable = "Sample",
  taxa_to_test = c("Actinobacteria", "Verrucomicrobia"),
  test_by = "time",
  strat_by = NULL,
  comparisons_list = NULL,
  show_test = TRUE,
  adj = FALSE,
  y_position = NULL
)
}
\arguments{
\item{ps}{A phyloseq object.}

\item{rank}{The taxonomical rank, found in the tax_table columns, for taxa agglomeration.}

\item{time_variable}{Name of the column, in sample_data slot of the phyloseq object, where repeated measures/time variable is stored.}

\item{treatment_variable}{Name of the column, in sample_data slot of the phyloseq object, where treatment/condition variable is stored.}

\item{sample_variable}{Name of the column, in sample_data slot of the phyloseq object, where sample names variable is stored.}

\item{taxa_to_test}{A vector of character containings taxonomical group names of interest. If \code{NULL}, all the taxonomical features of the agglomerated rank are tested.}

\item{test_by}{\code{"time"} character string if the difference between repeated measures/time variable is to be tested or \code{"treatment"} character string if the difference between cndition/treatment variable is to be tested.}

\item{strat_by}{Name of the column, in sample_data slot of the phyloseq object, for the stratification variable. Default set to \code{NULL}.}

\item{comparisons_list}{A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the groups of interest, to be compared. If empty and \code{test_by = "time"} each level of repeated measures/time variable is compared with his baseline (the first level). If empty and \code{test_by = "treatment"} each level of condition/treatment variable is compared with his baseline (the first level).}

\item{show_test}{If TRUE, significant p-values or adjusted p-values are shown.}

\item{adj}{If TRUE (default = FALSE) and \code{"show_test = TRUE"}, significant adjusted (FDR) p-values are shown.}

\item{y_position}{If specified (default = NULL) set the y-axis position of the tests.}
}
\value{
a list with two slots: \code{plot} which contains a \code{ggplot} object with taxa_to_test faceted rows and repeated measures/time or condition/treatment levels as faceted columns. And \code{data} slot, which contains the summary mean and sd for relative abundances.
}
\description{
This function takes a phyloseq object with otu_table, sam_data and tax_table slots and compute the relative abundances for the features in a specified taxonomical rank.
If a repeated measures/time variable and a treatment/condition variables are specified it performs test for the specified taxa.
The graphical representation can be of two types: by the repeated measures/time variable or by the treatment/condition variable.
Non parametric Mann-Whitney or Wilcox test are performed to study differences in the relative abundances.
}
