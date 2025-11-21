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

asc_CDs_boxplot <- function(data, lineage, hierarchy_level, tissue, measure){
  measure <- enquo(measure)
  
  data|>
    cell_type_order() |>
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
    df <- df |> filter(!str_detect(cell_type, cell_type_exclude))
  }
  
  df |>
    cell_type_order() |>
    ggplot(aes(x = cell_type,
               y = estimate,
               color = factor(is_significant, levels = c("reference", "yes", "no")),
               shape = tissue)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.4) +
    facet_wrap(~ CD, ncol = 2, scales = "free_y") +
    scale_color_manual(
      name = "Significance",
      values = c(
        "reference" = "black",
        "yes"       = "#29cf8f",
        "no"        = "#fa2840"
      )
    ) +
    scale_shape_discrete(name = "Tissue") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = element_text(size = 14, face = "bold"),
      axis.line = element_line(color = "black")
    ) +
    labs(
      title = title,
      x = "Cell Type",
      y = "Estimate"
    )
}
