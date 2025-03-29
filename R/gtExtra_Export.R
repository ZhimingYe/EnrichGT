# Function export from gtExtra with litte modifying
# copyright Thomas Mock
# MIT license
# https://jthomasmock.github.io/gtExtras/index.html
gt_merge_stack2 <- function(
  gt_object,
  col1,
  col2,
  palette = c("black", "grey"),
  ...,
  small_cap = TRUE,
  font_size = c("14px", "10px"),
  font_weight = c("bold", "bold")
) {
  colors <- scales::col2hcl(palette, ...)
  col1_bare <- rlang::enexpr(col1) |> rlang::as_string()
  row_name_var <- gt_object[["_boxhead"]][["var"]][which(
    gt_object[["_boxhead"]][["type"]] == "stub"
  )]
  data_in <- gt::gt_index2(gt_object, column = {{ col2 }})
  gt_object |>
    gt::text_transform(
      locations = if (isTRUE(row_name_var == col1_bare)) {
        gt::cells_stub(rows = gt::everything())
      } else {
        gt::cells_body(columns = {{ col1 }})
      },
      fn = function(x) {
        if (small_cap) {
          font_variant <- "small-caps"
        } else {
          font_variant <- "normal"
        }
        glue::glue(
          "<div style='line-height:{font_size[1]}'><span style='font-weight:{font_weight[1]};font-variant:{font_variant};color:{colors[1]};font-size:{font_size[1]}'>{x}</span></div>\n        <div style='line-height:{font_size[2]}'><span style ='font-weight:{font_weight[2]};color:{colors[2]};font-size:{font_size[2]}'>{data_in}</span></div>"
        )
      }
    ) |>
    gt::cols_hide(columns = {{ col2 }})
}

gt_index2 <- function(gt_object, column, as_vector = TRUE) {
  if (length(gt_object[["_row_groups"]]) >= 1) {
    gt_row_grps <- gt_object[["_row_groups"]]
    grp_vec_ord <- gt_object[["_stub_df"]] |>
      dplyr::mutate(group_id = factor(group_id, levels = gt_row_grps)) |>
      dplyr::arrange(group_id) |>
      dplyr::pull(rownum_i)
    df_ordered <- gt_object[["_data"]] |> dplyr::slice(grp_vec_ord)
  } else {
    df_ordered <- gt_object[["_data"]]
  }
  if (isTRUE(as_vector)) {
    df_ordered |> dplyr::pull({{ column }})
  } else {
    df_ordered
  }
}
gt_hulk_col_numeric2 <- function(
  gt_object,
  columns = NULL,
  pal,
  domain = NULL,
  ...,
  trim = FALSE
) {
  pal_hex <- pal
  if (isTRUE(trim)) pal_hex <- pal_hex[2:6]
  hulk_pal <- function(x) {
    (scales::col_numeric(pal_hex, domain = domain, ...))(x)
  }
  gt::data_color(gt_object, columns = {{ columns }}, fn = hulk_pal)
}
