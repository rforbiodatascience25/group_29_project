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

asc_CDs_boxplot <- function(data, lineage, hierarchy_level, tissue, measure){
  measure <- enquo(measure)
  
  data|>
    filter(lineage == !!lineage,
           hierarchy_level == !!hierarchy_level,
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
    scale_color_manual(values = c("CD4 T cells" = "navy",
                                  "CD8 T cells" = "hotpink")) +
    labs(subtitle = paste0("CD4", pair, " vs CD8", pair),
         x = "Significant CDs",
         y = "log(MedQb)",
         color = "Lineage")
}

PCA_rotation <- function(filename, data, prefix, amount, xlimits, ylimits, title) {
  arrow_style <- arrow(
    angle  = 30, 
    ends   = "first", 
    type   = "closed", 
    length = grid::unit(5, "pt")
  )
  
  rotation_data <- data |>
    tidy(matrix = "rotation") |>
    pivot_wider(
      names_from   = "PC", 
      names_prefix = "PC", 
      values_from  = "value"
    )
  
  plot_data <- rotation_data |>
    filter(str_starts(column, prefix)) |>
    mutate(
      contrib      = abs(PC1) + abs(PC2),
      arrow_length = sqrt(PC1^2 + PC2^2)
    ) |>
    slice_max(contrib, n = amount)    # <- no arrow_length filter
  
  ggplot(plot_data, aes(PC1, PC2)) +
    geom_segment(
      aes(xend = 0, yend = 0),
      arrow = arrow_style
    ) +
    geom_text_repel(
      aes(label = str_remove(column, prefix)),
      size = 3,
      color = "hotpink",
      min.segment.length = 0,
      direction = "both",
      max.overlaps = Inf
    ) +
    xlim(xlimits) +
    ylim(ylimits) +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 16) +
    labs(title = title)
}


heat_dot_plot <- function(data, bar){
  data |>
  CD_order(reverse = TRUE) |>
  ggplot(aes(x = cell_type, 
             y = CD)) +
  geom_point(aes(size = PEpos, 
                 color = log10(MedQb))) +
  geom_tile(data = bar, 
            aes(x = cell_type, 
                y = y, 
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
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 0, 
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20), 
        legend.position = "right",
        plot.margin = margin(5, 5, 5, 5),
        panel.grid = element_blank()) +
  labs(x = "Cell Type",
       y = "CD",
       size = "Frequency of Positive Cells (%)",
       color = "Fluorescence [log10(Median ABC)]")
}


save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(here::here(paste0("results/", filename)), 
         plot, 
         dpi = 300, 
         width = width, 
         height = height)
}
