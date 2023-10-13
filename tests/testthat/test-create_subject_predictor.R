testthat::test_that(
  "create_subject_predictor gives correct output",
  {

    testthat::expect_error(
      create_subject_predictor(example_walking),
      regexp = "Unsupported date-time class"
    )
    df = example_walking %>%
      mutate(time = lubridate::floor_date(Sys.time(), "seconds") + time)

    out = create_subject_predictor(df, lags = 1:3)
    testthat::expect_true(
      all(c("lag_3_(2.75,3]_(2.75,3]", "second") %in% colnames(out))
    )
    testthat::expect_true(!anyNA(out$second))

    testthat::expect_warning(
      {
        out = create_subject_predictor(df, lags = c(1, 101))
      },
      regexp = "max samples in a second"
    )
    # ensure all "bad" lags are dropped
    testthat::expect_true(
      !any(grepl("lag_101", colnames(out)))
    )
    testthat::expect_true(!anyNA(out$second))

  })
