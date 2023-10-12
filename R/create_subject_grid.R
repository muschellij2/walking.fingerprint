
#' Create Grid Data for single subject from raw accelerometry
#'
#' @param data a `data.frame` with the columns `time` and `vm`
#' @param max_signal maximum signal of vector magnitude to create grid over,
#' grid will be `0` to `max_signal` by `bin_size`
#' @param bin_size Size of the bins. `max_signal` must be
#' evenly divisible by `bin_size`.
#'
#' @return A `data.frame` of subject level lags and cut-values
#' @export
#'
#' @examples
create_subject_predictor = function(
    data,
    max_signal = 3,
    bin_size = 0.25,
    n_lags = 3
) {
  cut_lagsig = cut_sig = vm = lag_vm = time = second = NULL
  rm(list = c("time", "second", "vm", "lag_vm",
              "cut_lagsig", "cut_sig"))

  assertthat::assert_that(
    is.data.frame(data),
    assertthat::has_name(data, "time"),
    assertthat::has_name(data, "vm")
  )
  assertthat::assert_that(
    assertthat::is.count(max_signal/bin_size),
    is.finite(max_signal/bin_size)
  )

  data = data %>%
    dplyr::mutate(second = lubridate::floor_date(
      dplyr::.data[["time"]], "second"))
  check = data %>%
    dplyr::count(dplyr::.data[["second"]])
  assertthat::assert_that(
    any(check$n > 1)
  )

  cell_breaks = seq(0, max_signal, by = bin_size)
  lags = seq(n_lags)

  result = purrr::map_dfr(lags, .f = function(lag) {
    out = data %>%
      dplyr::group_by(second) %>%
      # for each second, calculate signal and lagged signal
      dplyr::mutate(lag_vm = dplyr::lag(vm, n = lag)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        cut_sig = cut(
          vm,
          breaks = cell_breaks,
          include.lowest = TRUE
        ),
        cut_lagsig = cut(
          lag_vm,
          breaks = cell_breaks,
          include.lowest = TRUE
        )
      ) %>%
      dplyr::filter(!is.na(cut_sig) | !is.na(cut_lagsig))
    out = out %>%
      dplyr::count(second, cut_sig, cut_lagsig, .drop = FALSE) %>%
      dplyr::mutate(
        lag = lag
      )
  })

  return(result)
}
