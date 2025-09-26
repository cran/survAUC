# tests/testthat/helper-data.R
lung_clean <- local({
  if (!requireNamespace("survival", quietly = TRUE)) return(NULL)
  set.seed(1)
  lung <- survival::lung
  lung <- subset(lung, complete.cases(time, status, age, sex, ph.ecog))
  lung$status <- as.integer(lung$status == 2)  # 1=event, 0=censor
  lung
})
