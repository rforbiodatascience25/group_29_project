pivot_selected_variables <- function(data, ID_col, variables, tissue){
  data |>
    select(ID_col, ends_with(variables)) |>
    pivot_longer(cols = ends_with(variables),
                 names_to = c("cell_type", ".value"),
                 names_pattern = "(.*)_(.*)",
                 values_frop_na = TRUE) |>
    mutate(tissue = tissue,
           .before = ID_col)
}
