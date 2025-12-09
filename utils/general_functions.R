library(multcompView)
library(FSA)
library(latex2exp)
library(kableExtra)



# Compute Compact Letter Display
compute_cld <- function(dunn_result, short_label_order, label_order, threshold = 0.05) {
  # Computes Compact Letter Display

  # Diversity
  # Make a zero matrix with rownames and colnames
  mat <- matrix(
    0,
    nrow = length(short_label_order),
    ncol = length(short_label_order),
    dimnames = list(short_label_order, short_label_order)
  )


  # Create zero matrix
  mat <- matrix(0,
    nrow = length(short_label_order), ncol = length(short_label_order),
    dimnames = list(short_label_order, short_label_order)
  )

  # Fill matrix
  for (i in 1:nrow(dunn_result$res)) {
    comps <- strsplit(dunn_result$res$Comparison[i], " - ")[[1]] # split the comparison
    mat[comps[1], comps[2]] <- dunn_result$res$P.adj[i] # assign P.adj
    mat[comps[2], comps[1]] <- dunn_result$res$P.adj[i] # symmetric assignment
  }

  rownames(mat) <- label_order
  colnames(mat) <- label_order

  mat_letters_div <- multcompLetters(mat, threshold = threshold)
}


# Format p-values with asterisks
format_pval <- function(p) {
  case_when(
    p < 0.001 ~ paste0(sprintf("%.4f", p), " ***"),
    p < 0.01 ~ paste0(sprintf("%.4f", p), " ** "),
    p < 0.05 ~ paste0(sprintf("%.4f", p), " *  "),
    TRUE ~ sprintf("%.4f", p)
  )
}





plot_scatter_and_box_plot <- function(
    df,
    multiplicity_col,
    diversity_col,
    group_col,
    label_map,
    short_label_map,
    custom_colors,
    title,
    multiplicity_name,
    diversity_name,
    group_name,
    plot_location,
    cld_div,
    cld_mul,
    width = 14,
    height = 12,
    sig_scatter_factor = 0) {
  # Method start

  df_plot <- data.frame(SampleID = df$SampleID, div = df[[diversity_col]], mul = df[[multiplicity_col]], Group = df[[group_col]])


  df_plot$Group <- factor(df_plot$Group, levels = label_order)

  custom_colors <- setNames(custom_colors, label_order)

  # Random noise
  df_scatter_plot <- df_plot
  if (sig_scatter_factor > 0) {
    df_scatter_plot$mul <- df_scatter_plot$mul + rnorm(length(df_scatter_plot$mul), 0, max(df_scatter_plot$mul) / sig_scatter_factor)
    df_scatter_plot$div <- df_scatter_plot$div + rnorm(length(df_scatter_plot$div), 0, max(df_scatter_plot$div) / sig_scatter_factor)
  }

  df_scatter_plot$Group <- factor(df_scatter_plot$Group,
    levels = label_order,
    labels = label_map[label_order]
  )


  renamed_custom_colors <- setNames(custom_colors, label_map[label_order])

  # Create scatter plot
  p_main <- ggplot(df_scatter_plot, aes(x = mul, y = div, color = Group)) +
    geom_point(size = 3, alpha = 0.9) + # Customize points
    scale_color_manual(values = renamed_custom_colors) +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      color = group_name
    ) + # Add labels
    theme(
      axis.title = element_text(size = 18), # Increase axis title font size
      axis.text = element_text(size = 16), # Increase axis label font size
      plot.title = element_text(size = 20, hjust = 0), # Increase plot title font size
      legend.title = element_text(size = 16), # Increase legend title font size
      legend.text = element_text(size = 14) # Increase legend text font size
    ) # Use a clean theme


  p_legend <- suppressWarnings(cowplot::get_legend(p_main))

  p_main <- p_main + theme(legend.position = "none")



  # Labels
  label_positions <- df_plot %>%
    group_by(Group) %>%
    summarize(y_pos = quantile(div, 0.75)) # Get the max y-value for each Group

  label_positions$labels <- short_label_map[label_positions$Group]

  # Compact Letter Display
  cld_positions <- df_plot %>%
    group_by(Group) %>%
    summarize(y_pos = quantile(div, 0.25)) # Get the max y-value for each Group

  cld_positions$labels <- cld_div$Letters[cld_positions$Group]




  # Diversity
  # Create the boxplot
  p_left <- ggplot(df_plot, aes(x = Group, y = div, color = Group)) +
    geom_boxplot() +
    labs(
      title = NULL,
      x = group_name,
      y = diversity_name
    ) + # Names
    geom_text(
      data = label_positions,
      aes(x = Group, y = y_pos, label = labels),
      vjust = -0.5,
      hjust = 0.5,
      size = 5.5, # Adjust text size
      inherit.aes = FALSE
    ) + # Compact Letter displays
    geom_text(
      data = cld_positions,
      aes(x = Group, y = y_pos, label = labels),
      vjust = 1.2,
      hjust = 0.5, # center horizontally
      size = 5, # Adjust text size
      inherit.aes = FALSE
    ) +
    scale_x_discrete(labels = label_map) +
    scale_color_manual(values = custom_colors) + # Add labels
    theme(
      axis.title = element_text(size = 18), # Increase axis title font size
      axis.text = element_text(size = 16), # Increase axis label font size
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.text.x = element_blank(), # Remove x-axis text
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    theme(legend.position = "none") # Hides the legend


  # Labels
  label_positions <- df_plot %>%
    group_by(Group) %>%
    summarize(y_pos = quantile(mul, 0.75)) # Get the max y-value for each Group

  label_positions$labels <- short_label_map[label_positions$Group]


  #  CLD
  cld_positions <- df_plot %>%
    group_by(Group) %>%
    summarize(y_pos = quantile(mul, 0.25)) # Get the max y-value for each Group

  cld_positions$labels <- cld_mul$Letters[cld_positions$Group]

  # cld_positions$hjust <- ifelse(nchar(cld_positions$labels) > 1, 1.5, 1.5)


  # Create the boxplot
  p_down <- ggplot(df_plot, aes(x = Group, y = mul, color = Group)) +
    geom_boxplot() +
    labs(
      # title = TeX("Distance Multiplicity"),
      title = NULL,
      x = group_name,
      y = multiplicity_name
    ) +
    geom_text(
      data = label_positions,
      aes(x = Group, y = y_pos, label = labels),
      vjust = -0.2,
      hjust = -0.5,
      size = 5.5, # Adjust text size
      inherit.aes = FALSE
    ) +
    # Compact Letter displays
    geom_text(
      data = cld_positions,
      aes(x = Group, y = y_pos, label = labels, hjust = hjust),
      vjust = -0.2,
      hjust = 1.5,
      size = 5, # Adjust text size
      inherit.aes = FALSE
    ) +
    scale_x_discrete(labels = label_map) +
    scale_color_manual(values = custom_colors) + # Add labels
    theme(
      axis.title = element_text(size = 18), # Increase axis title font size
      axis.text = element_text(size = 16), # Increase axis label font size
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.text.y = element_blank(), # Remove x-axis text
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    theme(legend.position = "none") +
    coord_flip() # Hides the legend


  # Combine left and main (side by side)
  top_row <- p_left + p_main + plot_layout(widths = c(3, 9))

  # Add the bottom row with p_down and the legend
  p_legend_plot <- ggdraw(p_legend)
  bottom_row <- p_legend_plot + p_down + plot_layout(widths = c(3, 8.3))


  # Stack the two rows
  final_plot <- top_row / bottom_row +
    plot_layout(heights = c(9, 3)) # Approximate height ratio

  # Save the full layout
  ggsave(paste0(plot_location), plot = final_plot, width = width, height = height, dpi = 300)

  return(final_plot)
}


kw_and_dunn_test <- function(df_test,
                             order,
                             labels,
                             multiplicity_col,
                             diversity_col,
                             group_col,
                             diversity_name,
                             group_name,
                             multiplicity_name,
                             kw_caption,
                             dun_caption,
                             kw_label,
                             dunn_label) {
  set.seed(1)


  df_test[group_col] <- factor(df_test[[group_col]], levels = order, labels = labels)

  # Run Kruskal-Wallis tests

  kw_div <- kruskal.test(reformulate(group_col, diversity_col), data = df_test)
  kw_mul <- kruskal.test(reformulate(group_col, multiplicity_col), data = df_test)

  kw_div$p.value[kw_div$p.value <= 0.0001] <- 0.0001
  kw_mul$p.value[kw_mul$p.value <= 0.0001] <- 0.0001


  # Create a data frame for the Kruskal-Wallis results
  kw_df <- data.frame(
    Metric = c("DiversityTempName", "MultiplicityTempName"),
    Chi_squared = c(kw_div$statistic, kw_mul$statistic),
    df = c(kw_div$parameter, kw_mul$parameter),
    p_value = c(format_pval(kw_div$p.value), format_pval(kw_mul$p.value))
  )

  # Output Kruskal-Wallis table for LaTeX with [!h] float placement
  kw_table <- kable(kw_df,
    format = "latex", booktabs = TRUE,
    caption = kw_caption
  ) %>%
    kable_styling(latex_options = "hold_position")

  kw_table <- gsub("end{tabular}", paste0("end{tabular}\n\\label{", kw_label, "}"), kw_table, fixed = TRUE)
  kw_table <- gsub("Chi\\_squared", "$\\chi^2$", kw_table, fixed = TRUE)
  kw_table <- gsub("p\\_value", "p-val", kw_table, fixed = TRUE)
  kw_table <- gsub("DiversityTempName", diversity_name, kw_table, fixed = TRUE)
  kw_table <- gsub("MultiplicityTempName", multiplicity_name, kw_table, fixed = TRUE)
  kw_table <- gsub(" \\, ", " ", kw_table, fixed = TRUE)


  # Run Dunn's tests
  dunn_div <- dunnTest(reformulate(group_col, diversity_col), data = df_test, method = "bh")
  dunn_mul <- dunnTest(reformulate(group_col, multiplicity_col), data = df_test, method = "bh")

  dunn_div$res$P.unadj[dunn_div$res$P.unadj <= 0.0001] <- 0.0001
  dunn_div$res$P.adj[dunn_div$res$P.adj <= 0.0001] <- 0.0001

  dunn_mul$res$P.unadj[dunn_mul$res$P.unadj <= 0.0001] <- 0.0001
  dunn_mul$res$P.adj[dunn_mul$res$P.adj <= 0.0001] <- 0.0001

  # Extract and label results
  dunn_div_table <- dunn_div$res %>%
    mutate(
      Z_div = sprintf("%.4f", Z),
      P.unadj_div = format_pval(P.unadj),
      P.adj_div = format_pval(P.adj)
    ) %>%
    select(Comparison, Z_div, P.unadj_div, P.adj_div)

  dunn_mul_table <- dunn_mul$res %>%
    mutate(
      Z_mul = sprintf("%.4f", Z),
      P.unadj_mul = format_pval(P.unadj),
      P.adj_mul = format_pval(P.adj)
    ) %>%
    select(Comparison, Z_mul, P.unadj_mul, P.adj_mul)

  # Join the two tables by Comparison
  dunn_df <- left_join(dunn_div_table, dunn_mul_table, by = "Comparison")

  # Output Dunn's test table for LaTeX with grouped headers
  dunn_table <- kable(dunn_df,
    format = "latex", booktabs = TRUE,
    caption = dun_caption
  ) %>%
    add_header_above(setNames(
      c(1, 3, 3),
      c(" ", "DiversityTempName", "MultiplicityTempName")
    )) %>%
    kable_styling(latex_options = "hold_position")

  # Final Fixes
  dunn_table <- gsub("end{tabular}", paste0("end{tabular}\n\\label{", dunn_label, "}"), dunn_table, fixed = TRUE)
  dunn_table <- gsub("Comparison & Z\\_div & P.unadj\\_div & P.adj\\_div & Z\\_mul & P.unadj\\_mul & P.adj\\_mul", "Comparison & Z & p-val (unadj) & p-val (adj) & Z & p-val (unadj) & p-val (adj)", dunn_table, fixed = TRUE)
  dunn_table <- gsub("DiversityTempName", diversity_name, dunn_table, fixed = TRUE)
  dunn_table <- gsub("MultiplicityTempName", multiplicity_name, dunn_table, fixed = TRUE)
  dunn_table <- gsub(" \\, ", " ", dunn_table, fixed = TRUE)


  resp <- list(
    kw_df = kw_df,
    kw_table = kw_table,
    kw_div = kw_div,
    kw_mul = kw_mul,
    dunn_df = dunn_df,
    dunn_table = dunn_table,
    dunn_div = dunn_div,
    dunn_mul = dunn_mul
  )

  return(resp)
}


plot_taxa_boxplot <- function(df,
                              taxa_col,
                              taxa_name,
                              group_col,
                              value_col,
                              df_cld = NULL,
                              label_order = NULL,
                              custom_colors = NULL,
                              value_name = NULL,
                              group_name = NULL,
                              label_map = NULL,
                              title = NULL) {
  # Convert column names to symbols
  taxa_col <- sym(taxa_col)
  group_col <- sym(group_col)
  value_col <- sym(value_col)

  # Reorder group factor if label_order provided
  if (!is.null(label_order)) {
    df <- df %>%
      mutate(!!group_col := factor(!!group_col, levels = label_order))
  }

  # Y-axis label with LaTeX
  y_label <- if (!is.null(value_name)) latex2exp::TeX(value_name) else rlang::as_name(value_col)

  # Boxplot dodge width
  dodge_width <- 0.75

  # Base plot
  p <- ggplot(df, aes(x = !!group_col, y = !!value_col, fill = !!group_col)) +
    geom_boxplot(width = 0.7, position = position_dodge(width = dodge_width)) +
    facet_wrap(vars(!!taxa_col), scales = "free_x", nrow = 1) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(size = 20, hjust = 0.5), # Increase plot title font size
    ) +
    labs(
      title = title,
      x = taxa_name,
      y = y_label,
      fill = group_name %||% rlang::as_name(group_col)
    )

  # Apply custom legend labels and colors
  if (!is.null(label_map)) {
    p <- p + scale_fill_manual(values = custom_colors, labels = label_map)
  } else if (!is.null(custom_colors)) {
    p <- p + scale_fill_manual(values = custom_colors)
  }

  # Add letters from df_cld if provided
  if (!is.null(df_cld)) {
    # Ensure group_col factor levels match
    df_cld <- df_cld %>%
      mutate(!!group_col := factor(!!group_col, levels = label_order))

    # Compute 75th percentile for label positions
    label_positions <- df %>%
      group_by(!!taxa_col, !!group_col) %>%
      summarize(y_pos = quantile(!!value_col, 0.75, na.rm = TRUE), .groups = "drop") %>%
      left_join(df_cld %>% rename(label = !!value_col),
        by = c(rlang::as_name(taxa_col), rlang::as_name(group_col))
      )

    # Add letters aligned with boxplots
    p <- p + geom_text(
      data = label_positions,
      aes(x = !!group_col, y = y_pos, label = label),
      inherit.aes = FALSE,
      # position = position_dodge(width = dodge_width),
      hjust = 1.3, # Shifts to the left
      vjust = -0.5
    )
  }

  return(p)
}




single_kw_and_dunn_test <- function(df_test,
                                    order,
                                    labels,
                                    value_col,
                                    group_col,
                                    value_name,
                                    group_name,
                                    kw_caption,
                                    dun_caption,
                                    kw_label,
                                    dunn_label) {
  set.seed(1)


  df_test[group_col] <- factor(df_test[[group_col]], levels = order, labels = labels)

  # Run Kruskal-Wallis tests

  kw_val <- kruskal.test(reformulate(group_col, value_col), data = df_test)
  kw_val$p.value[kw_val$p.value <= 0.0001] <- 0.0001



  # Create a data frame for the Kruskal-Wallis results
  kw_df <- data.frame(
    Metric = "ValueTempName",
    Chi_squared = kw_val$statistic,
    df = kw_val$parameter,
    p_value = format_pval(kw_val$p.value)
  )


  # Output Kruskal-Wallis table for LaTeX with [!h] float placement
  kw_table <- kable(kw_df,
    format = "latex", booktabs = TRUE,
    caption = kw_caption
  ) %>%
    kable_styling(latex_options = "hold_position")

  kw_table <- gsub("end{tabular}", paste0("end{tabular}\n\\label{", kw_label, "}"), kw_table, fixed = TRUE)
  kw_table <- gsub("Chi\\_squared", "$\\chi^2$", kw_table, fixed = TRUE)
  kw_table <- gsub("p\\_value", "p-val", kw_table, fixed = TRUE)
  kw_table <- gsub("ValueTempName", value_name, kw_table, fixed = TRUE)
  kw_table <- gsub(" \\, ", " ", kw_table, fixed = TRUE)


  # Run Dunn's tests
  dunn_val <- dunnTest(reformulate(group_col, value_col), data = df_test, method = "bh")

  dunn_val$res$P.unadj[dunn_val$res$P.unadj <= 0.0001] <- 0.0001
  dunn_val$res$P.adj[dunn_val$res$P.adj <= 0.0001] <- 0.0001

  # Extract and label results
  dunn_val_table <- dunn_val$res %>%
    mutate(
      Z_div = sprintf("%.4f", Z),
      P.unadj_div = format_pval(P.unadj),
      P.adj_div = format_pval(P.adj)
    ) %>%
    select(Comparison, Z_div, P.unadj_div, P.adj_div)


  # Join the two tables by Comparison
  dunn_df <- dunn_val_table

  # Output Dunn's test table for LaTeX with grouped headers
  dunn_table <- kable(dunn_df,
    format = "latex", booktabs = TRUE,
    caption = dun_caption
  ) %>%
    add_header_above(setNames(
      c(1, 3),
      c(" ", "ValueTempName")
    )) %>%
    kable_styling(latex_options = "hold_position")

  # Final Fixes
  dunn_table <- gsub("end{tabular}", paste0("end{tabular}\n\\label{", dunn_label, "}"), dunn_table, fixed = TRUE)
  dunn_table <- gsub("Comparison & Z\\_div & P.unadj\\_div & P.adj\\_div", "Comparison & Z & p-val (unadj) & p-val (adj)", dunn_table, fixed = TRUE)
  dunn_table <- gsub("ValueTempName", value_name, dunn_table, fixed = TRUE)
  dunn_table <- gsub(" \\, ", " ", dunn_table, fixed = TRUE)


  resp <- list(
    kw_df = kw_df,
    kw_table = kw_table,
    kw_val = kw_val,
    dunn_df = dunn_df,
    dunn_table = dunn_table,
    dunn_val = dunn_val
  )

  return(resp)
}
