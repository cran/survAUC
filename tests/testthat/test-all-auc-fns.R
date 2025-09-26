# tests/testthat/test-all-auc-fns.R
test_that("All exported AUC.* functions accept the common (Surv, Surv, lp, times) calling convention", {
  skip_if_not_installed("survival")
  
  ex <- getNamespaceExports("survAUC") 
  auc_fns <- grep("^AUC\\.", ex, value = TRUE) 
  
  skip_if_not(length(auc_fns) > 0, "No AUC.* exports to test") 
  
  lung <- lung_clean 
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog, data = lung, ties = "breslow") 
  lp <- as.numeric(stats::predict(fit, type = "lp")) 
  times <- as.numeric(stats::quantile(lung$time, probs = c(0.25, 0.5)))

  SurvTrain <- survival::Surv(lung$time, lung$status)
  SurvTest  <- survival::Surv(lung$time, lung$status)

  for (fn in auc_fns) {
    f <- get(fn, envir = asNamespace("survAUC"))
    fmls <- names(formals(f))
  
    args <- list()
  
    # First survival arg
    if ("Surv.rsp" %in% fmls) args$Surv.rsp <- SurvTrain else args[[fmls[1]]] <- SurvTrain
    # Optional test survival
    if ("Surv.rsp.new" %in% fmls) args$Surv.rsp.new <- SurvTest
  
    # Linear predictors (train + optional test)
    if ("lp" %in% fmls) args$lp <- lp
    if ("lpnew" %in% fmls) args$lpnew <- lp
    if ("new.lp" %in% fmls && !"lpnew" %in% fmls) args$`new.lp` <- lp  # alt naming in some APIs
    if ("marker" %in% fmls && !"lp" %in% fmls) args$marker <- lp
    if ("marker.new" %in% fmls && !"lpnew" %in% fmls) args$`marker.new` <- lp
  
    # Times if supported
    if ("times" %in% fmls) args$times <- times
  
    expect_no_error({
      out <- do.call(f, args)
      expect_true(is.list(out) || is.numeric(out) || is.data.frame(out))
    }, message = paste("Function:", fn, "| formals:", paste(fmls, collapse = ", ")))
  }
  ex <- getNamespaceExports("survAUC")
  auc_fns <- grep("^AUC[.]", ex, value = TRUE)  # literal dot

  SurvTrain <- survival::Surv(lung$time, lung$status)
  SurvTest  <- survival::Surv(lung$time, lung$status)

  for (fn in auc_fns) {
    f <- get(fn, envir = asNamespace("survAUC"))
    fmls <- names(formals(f))
  
    args <- list()
  
    # First survival arg
    if ("Surv.rsp" %in% fmls) args$Surv.rsp <- SurvTrain else args[[fmls[1]]] <- SurvTrain
    # Optional test survival
    if ("Surv.rsp.new" %in% fmls) args$Surv.rsp.new <- SurvTest
  
    # Linear predictors (train + optional test)
    if ("lp" %in% fmls) args$lp <- lp
    if ("lpnew" %in% fmls) args$lpnew <- lp
    if ("new.lp" %in% fmls && !"lpnew" %in% fmls) args$`new.lp` <- lp  # alt naming in some APIs
    if ("marker" %in% fmls && !"lp" %in% fmls) args$marker <- lp
    if ("marker.new" %in% fmls && !"lpnew" %in% fmls) args$`marker.new` <- lp
  
    # Times if supported
    if ("times" %in% fmls) args$times <- times
  
    expect_no_error({
      out <- do.call(f, args)
      expect_true(is.list(out) || is.numeric(out) || is.data.frame(out))
    }, message = paste("Function:", fn, "| formals:", paste(fmls, collapse = ", ")))
  }
})
