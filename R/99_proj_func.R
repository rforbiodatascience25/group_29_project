pivot_selected_variables <- function(data, ID_col, variables, tissue){
  data |>
    select(ID_col, ends_with(variables)) |>
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