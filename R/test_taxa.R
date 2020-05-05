#' Test a taxonomical group of features
#'
#' This function takes a phyloseq object with otu_table, sam_data and tax_table slots and compute the relative abundances for the features in a specified taxonomical rank.
#' If a repeated measures/time variable and a treatment/condition variables are specified it performs test for the specified taxa.
#' The graphical representation can be of two types: by the repeated measures/time variable or by the treatment/condition variable.
#' Non parametric Mann-Whitney or Wilcox test are performed to study differences in the relative abundances.
#'
#' @param ps A phyloseq object.
#' @param rank The taxonomical rank, found in the tax_table columns, for taxa agglomeration.
#' @param time_variable Name of the column, in sample_data slot of the phyloseq object, where repeated measures/time variable is stored.
#' @param treatment_variable Name of the column, in sample_data slot of the phyloseq object, where treatment/condition variable is stored.
#' @param sample_variable Name of the column, in sample_data slot of the phyloseq object, where sample names variable is stored.
#' @param taxa_to_test A vector of character containings taxonomical group names of interest. If \code{NULL}, all the taxonomical features of the agglomerated rank are tested.
#' @param test_by \code{"time"} character string if the difference between repeated measures/time variable is to be tested or \code{"treatment"} character string if the difference between cndition/treatment variable is to be tested.
#' @param strat_by Name of the column, in sample_data slot of the phyloseq object, for the stratification variable. Default set to \code{NULL}.
#' @param comparisons_list A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the groups of interest, to be compared. If empty and \code{test_by = "time"} each level of repeated measures/time variable is compared with his baseline (the first level). If empty and \code{test_by = "treatment"} each level of condition/treatment variable is compared with his baseline (the first level).
#' @return a \code{ggplot} object with taxa_to_test faceted columns and repeated measures/time or condition/treatment levels as faceted rows.
#' @export

test_taxa <- function(ps,
                      rank = "Phylum",
                      time_variable = "TimePoint",
                      treatment_variable = "Treatment",
                      sample_variable = "Sample",
                      taxa_to_test = c("Actinobacteria","Verrucomicrobia"),
                      test_by = "time",
                      strat_by = NULL,
                      comparisons_list = NULL) {
  # Glom by rank
  physeq <- tax_glom(ps, rank)
  # Get the relative abundances
  TSS <- apply(physeq@otu_table@.Data, 2, function(x) x/sum(x))
  # Get the taxa names from the tax_table
  name_indexes <- match(rownames(TSS), rownames(tax_table(physeq)))
  rownames(TSS) = as.vector(tax_table(physeq)[name_indexes, rank])
  # Convert to data.frame the SAMPLE x TAXA matrix
  df <- data.frame(t(TSS))
  # Add metadata about time, treatment and sample ID
  df[, c(time_variable, treatment_variable, sample_variable)] <- physeq@sam_data[, c(time_variable, treatment_variable, sample_variable)]

  # A recursive fast method to order df by specified variables
  iterative_ordering <- function(df, var_names, i = 1) {
    ord <- order(df[, var_names[i]])
    df_ordered <- df[ord, ]
    if (i < length(var_names))
      iterative_ordering(df_ordered, var_names, i + 1) else return(df_ordered)
  }

  # Order the df
  df_ordered <- iterative_ordering(df = df, var_names = c(sample_variable, time_variable, treatment_variable))
  # Melt the df
  df_melted <- reshape2::melt(df_ordered[, -ncol(df_ordered)])
  # Check taxa to test
  if (length(taxa_to_test) > 0) {
    df_to_plot <- df_melted[df_melted$variable %in% taxa_to_test, ]
  } else df_to_plot <- df_melted
  # Get mean and sd for samples grouped by specified variables
  df_summary <- ddply(df_to_plot, .variables = c(time_variable, treatment_variable, "variable"), function(x) return(data.frame(sd = sd(x[, "value"]),
                                                                                                                               value = mean(x[, "value"]))))

  if (test_by == "treatment") {
    if(is.null(comparisons_list)){
      comparisons_list = list()
      levels_treatment_var <- unique(df_to_plot[, treatment_variable])
      combinations <- expand.grid(levels_treatment_var[1], levels_treatment_var[-1])
      for (i in 1:nrow(combinations)) {
        comparisons_list[[i]] <- as.character(unname(unlist(c(combinations[i, ]))))
      }
    }
    # Plot the data
    ggplot(data = df_to_plot, aes(x = eval(parse(text = treatment_variable)), y = value, fill = eval(parse(text = treatment_variable)))) +
      geom_col(data = df_summary,position = position_dodge(), color = "black") +
      geom_errorbar(data = df_summary, aes(ymin = value, ymax = value + sd), width = 0.2, position = position_dodge(0.9)) +
      facet_grid(variable ~ eval(parse(text = time_variable)) + eval(parse(text = strat_by)), scales = "free_y") +
      stat_compare_means(data = df_to_plot, method = "wilcox", comparisons = comparisons_list,label.y.npc = 0.6) +
      # scale_y_continuous(limits = c(0,1)) +
      labs(fill = treatment_variable, x = treatment_variable, y = "Relative Abundance")
  } else if (test_by == "time") {
    # Get the levels for the time variable. If they are more than 2, each level is compared with the baseline (first level).
    if(is.null(comparisons_list)){
      comparisons_list = list()
      levels_time_var <- unique(df_to_plot[, time_variable])
      combinations <- expand.grid(levels_time_var[1], levels_time_var[-1])
      for (i in 1:nrow(combinations)) {
        comparisons_list[[i]] <- as.character(unname(unlist(c(combinations[i, ]))))
      }
    }
    # Plot the data
    ggplot(data = df_to_plot, aes(x = eval(parse(text = time_variable)), y = value, fill = eval(parse(text = time_variable)))) +
      geom_col(data = df_summary, position = position_dodge(), color = "black") +
      geom_errorbar(data = df_summary, aes(ymin = value, ymax = value + sd), width = 0.2, position = position_dodge(0.9)) +
      facet_grid(variable ~ eval(parse(text = treatment_variable)) + eval(parse(text = strat_by)), scales = "free_y") +
      stat_compare_means(data = df_to_plot, method = "wilcox", comparisons = comparisons_list, paired = TRUE, label.y.npc = 0.6) +
      # scale_y_continuous(limits = c(0,0.1+max(df_to_plot$value))) +
      labs(fill = time_variable, x = time_variable, y = "Relative Abundance")
  }
}
