# tests/testthat/test-auc-uno.R
test_that("AUC.uno computes sensible values on a small Cox model", {
  skip_if_not_installed("survival")
  skip_if_not("AUC.uno" %in% getNamespaceExports("survAUC"),
              "AUC.uno not exported")

  lung <- lung_clean
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                         data = lung, ties = "breslow")
  lp  <- as.numeric(stats::predict(fit, type = "lp"))
  times <- as.numeric(stats::quantile(lung$time, probs = c(0.25, 0.5, 0.75)))

  out <- survAUC::AUC.uno(
    survival::Surv(lung$time, lung$status),
    survival::Surv(lung$time, lung$status),
    lp,
    times = times
  )

  nm <- names(out)
  expect_true(any(nm %in% c("AUC", "auc")))
  auc <- if ("AUC" %in% nm) out$AUC else out$auc
  expect_length(auc, length(times))
  expect_true(all(is.finite(auc)))
  expect_true(all(auc >= 0 & auc <= 1))
  tms <- if ("times" %in% nm) out$times else times
  expect_equal(as.numeric(tms), as.numeric(times))
})
