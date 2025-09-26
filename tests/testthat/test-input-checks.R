# tests/testthat/test-input-checks.R
test_that("Functions error on mismatched lengths / bad inputs", {
  skip_if_not_installed("survival")
  if (!"AUC.uno" %in% getNamespaceExports("survAUC")) skip("AUC.uno not exported")

  lung <- lung_clean
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                         data = lung, ties = "breslow")
  lp  <- as.numeric(stats::predict(fit, type = "lp"))
  SurvTrain <- survival::Surv(lung$time, lung$status)
  SurvTest  <- survival::Surv(lung$time, lung$status)
  times <- as.numeric(stats::quantile(lung$time, probs = c(0.25, 0.5)))

  expect_error(
    survAUC::AUC.uno(SurvTrain, SurvTest, lp, times = c(NA_real_, -1)),
    regex = "time|times|NA|positive|non", ignore.case = TRUE
  )
})
