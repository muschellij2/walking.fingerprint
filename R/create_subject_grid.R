
#' Create Grid Data for single subject from raw accelerometry
#'
#' @param data a `data.frame` with the columns `time` and `vm`
#' @param max_signal maximum signal of vector magnitude to create grid over,
#' grid will be `0` to `max_signal` by `bin_size`
#' @param bin_size Size of the bins. `max_signal` must be
#' evenly divisible by `bin_size`.
#' @param lags the lag (in samples, not seconds) to compute bin inclusion.
#'
#' @return A `data.frame` of subject level lags and cut-values
#' @export
#'
#' @examples
#' df = example_walking %>%
#'    dplyr::mutate(time = lubridate::floor_date(Sys.time(), "seconds") + time)
#' create_subject_predictor(df)
#' create_subject_predictor(df, lags = c(1, 15))
create_subject_predictor = function(
    data,
    max_signal = 3,
    bin_size = 0.25,
    lags = c(15, 30, 45)
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
    is.finite(max_signal/bin_size),
    is.integer(lags) ||
      (is.numeric(lags) && all(lags == trunc(lags))),
    !any(is.na(lags))
  )

  data = data %>%
    dplyr::mutate(second = lubridate::floor_date(
      .data[["time"]], "second"))
  check = data %>%
    dplyr::count(.data[["second"]])
  assertthat::assert_that(
    any(check$n > 1)
  )

  n_second = data %>%
    dplyr::count(second)
  max_records = max(n_second$n)
  if (max(lags) >= max_records) {
    warning(
      paste0("There are ", max_records, " max samples in a second,",
             " but max lag is ", max(lags), ", likely an error!")
      )
  }

  cell_breaks = seq(0, max_signal, by = bin_size)

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
      dplyr::filter(!is.na(cut_sig) & !is.na(cut_lagsig))
    if (nrow(out) == 0) {
      return(NULL)
    }
    out = out %>%
      dplyr::count(second, cut_sig, cut_lagsig, .drop = FALSE) %>%
      dplyr::mutate(
        lag = paste0("lag_", lag)
      )
  })

  result = result %>%
    tidyr::pivot_wider(
      id_cols = c(second),
      names_from = c(lag, cut_sig, cut_lagsig),
      names_sep = "_",
      values_from = n
    )

  return(result)
}
