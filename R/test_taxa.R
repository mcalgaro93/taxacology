#' Test a taxonomical group of features
#'
#' This function takes a phyloseq object with otu_table, sam_data and tax_table slots and compute the relative abundances for the features in a specified taxonomical rank.
#' If a repeated measures/time variable and a treatment/condition variables are specified it performs test for the specified taxa.
#' The graphical representation can be of two types: by the repeated measures/time variable or by the treatment/condition variable.
#' Non parametric Mann-Whitney or Wilcox test are performed to study differences in the relative abundances.
#'
#' @import rstatix
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
#' @param show_test If TRUE, significant p-values or adjusted p-values are shown.
#' @param adj If TRUE (default = FALSE) and \code{"show_test = TRUE"}, significant adjusted (FDR) p-values are shown.
#' @param y_position If specified (default = NULL) set the y-axis position of the tests.
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
                      comparisons_list = NULL,
                      show_test = TRUE,
                      adj = FALSE,
                      y_position = NULL) {
  # Glom by rank
  physeq <- tax_glom(ps, rank)
  # Get the relative abundances
  TSS <- apply(physeq@otu_table@.Data, 2, function(x) x/sum(x))
  # Get the taxa names from the tax_table
  name_indexes <- match(rownames(TSS), rownames(tax_table(physeq)))
  rownames(TSS) = as.vector(tax_table(physeq)[name_indexes, rank])
  # Convert to data.frame the SAMPLE x TAXA matrix
  df <- data.frame(t(TSS))

  # A recursive fast method to order df by specified variables
  iterative_ordering <- function(df, var_names, i = 1) {
    ord <- order(df[, var_names[i]])
    df_ordered <- df[ord, ]
    if (i < length(var_names))
      iterative_ordering(df_ordered, var_names, i + 1) else return(df_ordered)
  }

  # Add metadata about time, treatment and sample ID + stratification variable if present
  if(!is.null(strat_by)){
    df[, c(time_variable, treatment_variable, strat_by, sample_variable)] <- physeq@sam_data[, c(time_variable, treatment_variable, strat_by, sample_variable)]
    # Order the df
    df_ordered <- iterative_ordering(df = df, var_names = c(sample_variable, strat_by, time_variable, treatment_variable))
  } else {
    df[, c(time_variable, treatment_variable, sample_variable)] <- physeq@sam_data[, c(time_variable, treatment_variable, sample_variable)]
    # Order the df
    df_ordered <- iterative_ordering(df = df, var_names = c(sample_variable, time_variable, treatment_variable))
  }

  # Melt the df
  df_melted <- reshape2::melt(df_ordered[, -ncol(df_ordered)])
  if(is.null(strat_by)){
      colnames(df_melted) <- c("time","treatment","feature","relative_abundance")
  } else {
    colnames(df_melted) <- c("time","treatment","strat","feature","relative_abundance")
    df_melted[,"time_strat"] <- paste0(df_melted$time,"\n",df_melted$strat)
    df_melted[,"treatment_strat"] <- paste0(df_melted$treatment,"\n",df_melted$strat)
  }
  # Check taxa to test
  if (length(taxa_to_test) > 0) {
    df_to_plot <- df_melted[df_melted$feature %in% taxa_to_test, ]
  } else df_to_plot <- df_melted
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
    if(!is.null(strat_by)){ # With an additional stratification variable
      # Get mean and sd for samples grouped by specified variables
      df_summary <- ddply(df_to_plot, .variables = ~ time_strat + treatment + feature, function(x) return(data.frame(sd = sd(x[, "relative_abundance"]), relative_abundance = mean(x[, "relative_abundance"]))))
      df_to_plot[,"time_strat"] <- paste0(df_to_plot$time,"\n",df_to_plot$strat)
      if(show_test){
        stat.test <- df_to_plot %>%
          group_by(feature, time_strat) %>%
          wilcox_test(relative_abundance ~ treatment) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p") %>%
          add_xy_position(x = "treatment")
      }
      g <- ggbarplot(df_to_plot, x = "treatment", y = "relative_abundance",
        add = "mean_se",
        fill = "treatment",
        facet.by = c("feature","time_strat"),
        error.plot = "upper_errorbar", ggtheme = theme_grey(),scales = "free_y") +
        geom_text(data = df_summary, aes(label = round(relative_abundance*100,2)), vjust = 1) +
        labs(fill = time_variable, x = time_variable, y = "Relative Abundance")
    } else { # or without stratification variable
      # Get mean and sd for samples grouped by specified variables
      df_summary <- ddply(df_to_plot, .variables = ~ time + treatment + feature, function(x) return(data.frame(sd = sd(x[, "relative_abundance"]), relative_abundance = mean(x[, "relative_abundance"]))))
      if(show_test){
        stat.test <- df_to_plot %>%
          group_by(feature, time) %>%
          wilcox_test(relative_abundance ~ treatment) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p") %>%
          add_xy_position(x = "treatment")
      }
      g <- ggbarplot(df_to_plot, x = "treatment", y = "relative_abundance",
        add = "mean_se",
        fill = "treatment",
        facet.by = c("feature", "time"),
        error.plot = "upper_errorbar", ggtheme = theme_grey(), scales = "free_y") +
        geom_text(data = df_summary, aes(label = round(relative_abundance*100,2)), vjust = 1) +
        labs(fill = time_variable, x = time_variable, y = "Relative Abundance")
    }

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
    if(!is.null(strat_by)){ # With an additional stratification variable
      # Get mean and sd for samples grouped by specified variables
      df_summary <- ddply(df_to_plot, .variables = ~ treatment_strat + time + feature, function(x) return(data.frame(sd = sd(x[, "relative_abundance"]), relative_abundance = mean(x[, "relative_abundance"]))))
      df_to_plot[,"treatment_strat"] <- paste0(df_to_plot$treatment,"\n",df_to_plot$strat)
      if(show_test){
        stat.test <- df_to_plot %>%
          group_by(feature, treatment_strat) %>%
          wilcox_test(relative_abundance ~ time, paired = TRUE) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p") %>%
          add_xy_position(x = "time")
      }
      g <- ggbarplot(df_to_plot, x = "time", y = "relative_abundance",
        add = "mean_se",
        fill = "time",
        facet.by = c("feature","treatment_strat"),
        error.plot = "upper_errorbar", ggtheme = theme_grey(),scales = "free_y") +
        geom_text(data = df_summary, aes(label = round(relative_abundance*100,2)), vjust = 1) +
        labs(fill = time_variable, x = time_variable, y = "Relative Abundance")
    } else { # or without stratification variable
      # Get mean and sd for samples grouped by specified variables
      df_summary <- ddply(df_to_plot, .variables = ~ time + treatment + feature, function(x) return(data.frame(sd = sd(x[, "relative_abundance"]), relative_abundance = mean(x[, "relative_abundance"]))))
      if(show_test){
        stat.test <- df_to_plot %>%
          group_by(feature, treatment) %>%
          wilcox_test(relative_abundance ~ time) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p") %>%
          add_xy_position(x = "time")
      }
      g <- ggbarplot(df_to_plot, x = "time", y = "relative_abundance",
        add = "mean_se",
        fill = "time",
        facet.by = c("feature", "treatment"),
        error.plot = "upper_errorbar", ggtheme = theme_grey(), scales = "free_y") +
        geom_text(data = df_summary, aes(label = round(relative_abundance*100,2)), vjust = 1) +
        labs(fill = treatment_variable, x = treatment_variable, y = "Relative Abundance")
    }
  }
  # If show_test is TRUE, add p-value or adjusted p-values to plot
  if(show_test){
    if(!adj){
      stat.test <- stat.test[,!grepl(pattern = "adj", x = colnames(stat.test))]
      g + stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p", y.position = ifelse(is.null(y_position), "y.position", y_position))
    } else g + stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.adj", y.position = ifelse(is.null(y_position), "y.position", y_position))
  } else g
}
