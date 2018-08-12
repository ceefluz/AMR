context("portion.R")

test_that("portions works", {
  # amox resistance in `septic_patients`
  expect_equal(portion_R(septic_patients$amox), 0.6603, tolerance = 0.0001)
  expect_equal(portion_I(septic_patients$amox), 0.0030, tolerance = 0.0001)
  expect_equal(1 - portion_R(septic_patients$amox) - portion_I(septic_patients$amox),
               portion_S(septic_patients$amox))
  expect_equal(portion_R(septic_patients$amox) + portion_I(septic_patients$amox),
               portion_IR(septic_patients$amox))
  expect_equal(portion_S(septic_patients$amox) + portion_I(septic_patients$amox),
               portion_SI(septic_patients$amox))

  # pita+genta susceptibility around 98.09%
  expect_equal(suppressWarnings(rsi(septic_patients$pita,
                                    septic_patients$gent,
                                    interpretation = "S")),
               0.9535,
               tolerance = 0.0001)

  # percentages
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(R = portion_R(cipr, as_percent = TRUE),
                           I = portion_I(cipr, as_percent = TRUE),
                           S = portion_S(cipr, as_percent = TRUE),
                           n = n_rsi(cipr),
                           total = n()) %>%
                 pull(n) %>%
                 sum(),
               1404)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_p = portion_S(cipr, as_percent = TRUE),
                           cipro_n = n_rsi(cipr),
                           genta_p = portion_S(gent, as_percent = TRUE),
                           genta_n = n_rsi(gent),
                           combination_p = portion_S(cipr, gent, as_percent = TRUE),
                           combination_n = n_rsi(cipr, gent)) %>%
                 pull(combination_n),
               c(202, 482, 201, 499))

  expect_warning(portion_R(as.character(septic_patients$amcl)))
  expect_warning(portion_S(as.character(septic_patients$amcl)))
  expect_warning(portion_S(as.character(septic_patients$amcl,
                                             septic_patients$gent)))
  expect_equal(n_rsi(as.character(septic_patients$amcl,
                                  septic_patients$gent)),
               1570)


  # check for errors
  expect_error(portion_IR(septic_patients %>% select(amox, amcl)))
  expect_error(portion_IR("test", minimum = "test"))
  expect_error(portion_IR("test", as_percent = "test"))
  expect_error(portion_I(septic_patients %>% select(amox, amcl)))
  expect_error(portion_I("test", minimum = "test"))
  expect_error(portion_I("test", as_percent = "test"))
  expect_error(portion_S("test", minimum = "test"))
  expect_error(portion_S("test", as_percent = "test"))
  expect_error(portion_S(septic_patients %>% select(amox, amcl)))
  expect_error(portion_S("R", septic_patients %>% select(amox, amcl)))
  expect_error(n_rsi(septic_patients %>% select(amox, amcl)))
  expect_error(n_rsi(septic_patients$amox, septic_patients %>% select(amox, amcl)))


  # check too low amount of isolates
  expect_identical(portion_R(septic_patients$amox, minimum = nrow(septic_patients) + 1),
                   NA)
  expect_identical(portion_I(septic_patients$amox, minimum = nrow(septic_patients) + 1),
                   NA)
  expect_identical(portion_S(septic_patients$amox, minimum = nrow(septic_patients) + 1),
                   NA)

  # warning for speed loss
  expect_warning(portion_R(as.character(septic_patients$gent)))
  expect_warning(portion_I(as.character(septic_patients$gent)))
  expect_warning(portion_S(septic_patients$amcl, as.character(septic_patients$gent)))

})

test_that("old rsi works", {
  # amox resistance in `septic_patients` should be around 66.33%
  expect_equal(suppressWarnings(rsi(septic_patients$amox)), 0.6633, tolerance = 0.0001)
  expect_equal(suppressWarnings(rsi(septic_patients$amox, interpretation = "S")), 1 - 0.6633, tolerance = 0.0001)

  # pita+genta susceptibility around 98.09%
  expect_equal(suppressWarnings(rsi(septic_patients$pita,
                                    septic_patients$gent,
                                    interpretation = "S",
                                    info = TRUE)),
               0.9535,
               tolerance = 0.0001)

  # count of cases
  expect_equal(septic_patients %>%
                 group_by(hospital_id) %>%
                 summarise(cipro_S = suppressWarnings(rsi(cipr, interpretation = "S",
                                                          as_percent = TRUE, warning = FALSE)),
                           cipro_n = n_rsi(cipr),
                           genta_S = suppressWarnings(rsi(gent, interpretation = "S",
                                                          as_percent = TRUE, warning = FALSE)),
                           genta_n = n_rsi(gent),
                           combination_S = suppressWarnings(rsi(cipr, gent, interpretation = "S",
                                                                as_percent = TRUE, warning = FALSE)),
                           combination_n = n_rsi(cipr, gent)) %>%
                 pull(combination_n),
               c(202, 482, 201, 499))

  # portion_df
  expect_equal(
    septic_patients %>% select(amox) %>% portion_df(TRUE) %>% pull(Percentage),
    c(septic_patients$amox %>% portion_S(),
      septic_patients$amox %>% portion_I(),
      septic_patients$amox %>% portion_R())
  )

})

test_that("prediction of rsi works", {
  amox_R <- septic_patients %>%
    filter(bactid == "ESCCOL") %>%
    rsi_predict(col_ab = "amox",
                col_date = "date",
                minimum = 10,
                info = TRUE) %>%
    pull("value")
  # amox resistance will increase according to data set `septic_patients`
  expect_true(amox_R[3] < amox_R[20])

  expect_output(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                            model = "binomial",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                            model = "loglin",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))
  expect_output(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                            model = "lin",
                            col_ab = "amox",
                            col_date = "date",
                            info = TRUE))

  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           model = "INVALID MODEL",
                           col_ab = "amox",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           col_ab = "NOT EXISTING COLUMN",
                           col_date = "date",
                           info = TRUE))
  expect_error(rsi_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                           col_ab = "amox",
                           col_date = "NOT EXISTING COLUMN",
                           info = TRUE))
  # almost all E. coli are mero S in the Netherlands :)
  expect_error(resistance_predict(tbl = filter(septic_patients, bactid == "ESCCOL"),
                                  col_ab = "mero",
                                  col_date = "date",
                                  info = TRUE))
})