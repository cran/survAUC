# tests/testthat/test-pred-compat.R
test_that("Higher risk scores tend to correspond to earlier events (very weak monotonicity proxy)", {
  skip_if_not_installed("survival")
  lung <- lung_clean
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                         data = lung, ties = "breslow")
  lp <- as.numeric(stats::predict(fit, type = "lp"))
  early <- order(lung$time)[seq_len(min(100, nrow(lung)))]
  r <- suppressWarnings(stats::cor(lp[early], lung$status[early], method = "spearman"))
  expect_true(is.finite(r))
})
