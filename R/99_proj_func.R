# For cleaning
pivot_selected_variables <- function(data, ID_col, variables, tissue){
  data |>
    select(ID_col, all_of(ends_with(variables))) |>
    pivot_longer(cols = ends_with(variables),
                 names_to = c("cell_type", ".value"),
                 names_pattern = "(.*)_(.*)",
                 values_drop_na = TRUE) |>
    mutate(tissue = tissue,
           .before = ID_col)
}

# For augmenting
add_lineage <- function(data){
  data |>
    mutate(lineage = case_when(tissue == "thymus" ~ "Thymocytes",
                               cell_type == "Lymphs" ~ "Lymphocytes",
                               str_detect(cell_type, "CD4") ~ "CD4 T cells",
                               str_detect(cell_type, "CD8") ~ "CD8 T cells",
                               str_starts(cell_type, "T") ~ "T cells",
                               .default = "B cells"),
           .before = cell_type)
}

# General functions
cell_type_order <- function(data){
  data |>
    mutate(cell_type = factor(cell_type, 
                              levels = c("DN34p",
                                         "DN34m1ap",
                                         "CD4ISP",
                                         "DP3m",
                                         "DP3p",
                                         "CD4SP1ap",
                                         "CD4SP1am",
                                         "CD8SP1ap",
                                         "CD8SP1am",
                                         "BnaiveTo",
                                         "CC",
                                         "CB",
                                         "UnswtMem",
                                         "SwtMem",
                                         "PC",
                                         "CD138negPC",
                                         "CD138posPC",
                                         "Lymphs",
                                         "T",
                                         "Tgd",
                                         "TCD4",
                                         "TCD4naive",
                                         "TCD4CM",
                                         "TCD4EM",
                                         "TCD4TEMRA",
                                         "TCD8",
                                         "TCD8naive",
                                         "TCD8CM",
                                         "TCD8EM",
                                         "TCD8TEMRA",
                                         "TCD8RAdim",
                                         "B",
                                         "Bnaive",
                                         "BnatEff",
                                         "BIgM",
                                         "BswMem",
                                         "Bdn",
                                         "B27high")))
}

CD_order <- function(data, reverse = FALSE){
  if (!reverse) {
    data |>
      mutate(CD_num = parse_number(CD),
             CD = fct_reorder(CD, CD_num)) |>
      select(-CD_num)
  }
  
  else {
    data |>
      mutate(CD_num = parse_number(CD),
             CD = fct_reorder(CD, CD_num),
             CD = fct_rev(CD)) |>
      select(-CD_num)
  }

}

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(here::here(paste0("results/", filename)), 
         plot, 
         dpi = 300, 
         width = width, 
         height = height)
}

# Analysis 1 - plots
#Scatter plots
PCA_plot <- function(fit, data, subtitle){
  
  augment(fit, data = data) |>
    ggplot(aes(
      .fittedPC1,
      .fittedPC2,
      shape = tissue,
      color = lineage
    )) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      subtitle = subtitle,
      x = "PC1",
      y = "PC2",
      color = "Lineage"
    ) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 13) +
    theme(
      legend.key.size  = unit(0.4, "lines"),
      legend.text      = element_text(size = 7),
      legend.spacing.y = unit(0.1, "lines"),
      legend.spacing.x = unit(0.2, "lines")
    )
  
}

#Rotation plots
PCA_rotation <- function(data, prefix, amount, xlimits, ylimits, subtitle) {
  arrow_style <- arrow(angle  = 30, 
                       ends = "first", 
                       type  = "closed", 
                       length = grid::unit(5, "pt"))
  
  data |>
    tidy(matrix = "rotation") |>
    pivot_wider(names_from = "PC", 
                names_prefix = "PC", 
                values_from  = "value") |>
    mutate(contrib = abs(PC1) + abs(PC2),
           arrow_length = sqrt(PC1^2 + PC2^2)) |>
    slice_max(contrib, n = amount) |>
    ggplot(aes(PC1, PC2)) +
    geom_segment(aes(xend = 0, 
                     yend = 0),
                 arrow = arrow_style,
                 color = "hotpink") +
    geom_text_repel(aes(label = str_remove(column, prefix)),
                    size = 3,
                    color = "black",
                    min.segment.length = 0,
                    direction = "both",
                    max.overlaps = Inf) +
    xlim(xlimits) +
    ylim(ylimits) +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 13) +
    labs(subtitle = subtitle)
}

# Analysis 2 - plots
heat_dot_plot <- function(data){
  data |>
    CD_order(reverse = TRUE) |>
    ggplot(aes(x = cell_type, 
               y = CD)) +
    geom_point(aes(size = PEpos, 
                   color = log10(MedQb))) +
    geom_tile(aes(x = cell_type, 
                  y = tube_bar, 
                  fill = tissue), 
              height = 0.5) +
    scale_size_continuous(breaks = c(5, 20, 60, 100), 
                          range = c(1, 6)) +
    scale_color_gradientn(colors = c("lightblue", 
                                     "orange", 
                                     "darkred"), 
                          values = scales::rescale(c(2, 3, 4, 5))) +
    scale_fill_manual(values = c("blood" = "#93aa9b", 
                                 "tonsil" = "#c8b145", 
                                 "thymus" = "#bb6e36")) +
    scale_x_discrete(position = "top") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 0),
          axis.text.y = element_text(),
          axis.title = element_text(), 
          legend.position = "right",
          plot.margin = margin(5, 5, 5, 5),
          panel.grid = element_blank()) +
    labs(x = "Cell Type",
         y = "CD",
         size = "PEpos (%)",
         color = "log10(MedQb)")
}

# Analysis 3 - plots
asc_CDs_boxplot <- function(data, lineage, hierarchy, tissue, measure){
  measure <- enquo(measure)
  
  data|>
    filter(lineage == !!lineage,
           hierarchy == !!hierarchy,
           tissue == !!tissue) |>
    drop_na(!!measure) |>
    group_by(CD) |>
    mutate(MedQb_grp = median(MedQb),
           measure_grp = median(!!measure)) |>
    ungroup() |>
    ggplot(aes(x = fct_reorder(CD, MedQb_grp),
               y = !!measure,
               fill = log10(measure_grp))) +
    geom_boxplot() +
    scale_y_log10(labels = trans_format("log10", 
                                        math_format(10^.x))) +
    scale_fill_gradient2(low = "lightblue", 
                         mid = "orange", 
                         high = "darkred", 
                         midpoint = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6,
                                     angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom") +
    labs(x = "CD",
         y = measure,
         fill = paste0("log10(", quo_name(measure), ")"),
         title = paste("Fluorescence of", lineage, "from the", tissue))
}

plot_scatter <- function(data, params, measure, lineage, tissue, color){
  measure <- enquo(measure)
  measure_str <- as_name(measure)
  
  measure_intercept <- params |> 
    filter(lineage == !!lineage,
           tissue == !!tissue) |>
    pull(paste0(measure_str, "_intercept"))
  
  measure_slope <- params |> 
    filter(lineage == !!lineage,
           tissue == !!tissue) |>
    pull(paste0(measure_str, "_slope"))
  
  measure_cor <- params |> 
    filter(lineage == !!lineage,
           tissue == !!tissue) |>
    pull(paste0(measure_str, "_cor"))
  
  data |>
    filter(lineage == !!lineage,
           tissue == !!tissue,
           !hierarchy %in% c(1, 2),
           lineage != "T cells") |>
    drop_na(MedQb, !!measure) |>
    ggplot(aes(x = log(!!measure),
               y = log(MedQb))) +
    geom_point(color = "black",
               alpha = 0.5) +
    geom_smooth(method = "lm",
                color = color) +
    theme_bw() + 
    labs(title = paste0(lineage, " from ", tissue),
         subtitle = paste0("y = ", round(measure_slope, 2), 
                           " * x + ", 
                           round(measure_intercept, 2), 
                           ", cor = ", 
                           round(measure_cor, 2)))
}

# Analysis 4 - plots
plot_cd_expression <- function(data, lineage_filter, cds, cell_type_exclude = NULL,
                               title = "Expression of Selected Markers") {
  
  df <- data |>
    filter(lineage %in% lineage_filter,
           CD %in% cds)
  
  if (!is.null(cell_type_exclude)) {
    df <- df |> filter(!str_detect(cell_type, 
                                   cell_type_exclude))
  }
  
  df |>
    cell_type_order() |>
    ggplot(aes(x = cell_type,
               y = estimate,
               color = factor(is_significant, 
                              levels = c("reference", 
                                         "yes", 
                                         "no")),
               shape = tissue)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = conf.low, 
                      ymax = conf.high), 
                  width = 0.4) +
    facet_wrap(~ CD, ncol = 2, 
               scales = "free_y") +
    scale_color_manual(
      name = "Significance",
      values = c(
        "reference" = "black",
        "yes"       = "#29cf8f",
        "no"        = "#fa2840"
      )
    ) +
    scale_shape_discrete(name = "Tissue") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, 
                                 vjust = 0.5, 
                                 hjust = 1),
      strip.text = element_text(size = 10, 
                                face = "bold"),
      axis.line = element_line(color = "black")
    ) +
    labs(
      title = title,
      x = "Cell Type",
      y = "log(MedQb)"
    )
}

# Analysis 5 - plots
plot_CD4vsCD8 <- function(data, pair, legend_position) {
  data |>
    filter(pair == !!pair) |>
    ggplot(aes(x = CD,
               y = estimate,
               color = lineage)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width = 0.4) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1),
          legend.position = legend_position) +
    ylim(c(0, 13)) +
    scale_color_manual(values = c("CD4 T cells" = "navy",
                                  "CD8 T cells" = "hotpink")) +
    labs(subtitle = paste0("CD4", pair, " vs CD8", pair),
         x = "Significant CDs",
         y = "log(MedQb)",
         color = "Lineage")
}

plot_sig_CDs <- function(data, group_selection, fill_color) {
  
  sig_CD_count <- data |>
    filter(is_significant == "yes") |>
    group_by({{ group_selection }}) |>
    summarise(count_CD = n_distinct(CD), .groups = "drop") |>
    mutate({{ group_selection }} := fct_rev(fct_reorder({{ group_selection }}, count_CD))) |>
    ggplot(aes(x = {{ group_selection }}, y = count_CD)) +
    geom_col(fill = fill_color) +
    theme_bw(base_size = 13) +
    labs(
      subtitle = paste("Count of significant CDs per", as_label(enquo(group_selection))),
      x = as_label(enquo(group_selection)),
      y = "Count"
    )
}

fit_linear_model <- function(data, formula) {
  data |>
    nest() |> 
    mutate(
      model_object = map(.x = data,
                         .f = ~lm(formula = formula, data = .x))
    ) |>
    mutate(
      model_object_tidy = map(.x = model_object,
                              .f = ~tidy(.x,
                                         conf.int = TRUE,
                                         conf.level = 0.95))
    ) |>
    unnest(model_object_tidy)
}

