
#' Create Grid Data for single subject from raw accelerometry
#'
#' @param data a `data.frame` with the columns `time` and `vm`
#' @param max_signal maximum signal of vector magnitude to create grid over,
#' grid will be `0` to `max_signal` by `bin_size`
#' @param bin_size Size of the bins. `max_signal` must be
#' evenly divisible by `bin_size`.
#'
#' @return
#' @export
#'
#' @examples
create_subject_predictor = function(
    data,
    max_signal = 3,
    bin_size = 0.25
    ) {
  # time = second = NULL
  # rm(list = c("time", "second"))

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
    dplyr::mutate(second = lubridate::floor_date(.data[["time"]], "second"))
  check = data %>%
    dplyr::count(.data[["second"]])
  assertthat::assert_that(
    any(check$n > 1)
  )


}
