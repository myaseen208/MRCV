# genloglin.fit()
# A function that estimates the model of interest
# limit.output argument needed when function is used with apply
# model.vars argument is needed for estimation involving bootstrap resamples

genloglin.fit <- function(data, model, nvars, limit.output = FALSE, model.vars = NULL) {
  op <- options()
  on.exit(options(op))
  model.data <- data
  if (!is.null(model.vars)) {
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter <- MRCV_globals$pb.counter + MRCV_globals$B / MRCV_globals$B.use
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
    model.data <- data.frame(model.vars, data)
    colnames(model.data)[ncol(model.data)] <- "count"
  }
  if (limit.output & is.null(model.vars)) {
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter <- MRCV_globals$pb.counter + 1
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
  }
  ##################
  # Replaced with suppressWarnings() whenever glm() is used
  #   Need because of adjustments when 0 count occurs (replaced with 0.5)
  # options(warn = -1)
  if (model == "spmi") {
    mod.fit <- suppressWarnings(glm(
      formula = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y,
      data = model.data, family = poisson(link = log)
    ))
  }
  if (model == "homogeneous") {
    mod.fit <- suppressWarnings(glm(
      formula = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y + wi:yj,
      data = model.data, family = poisson(link = log)
    ))
  }
  if (model == "w.main") {
    mod.fit <- suppressWarnings(glm(formula = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y + wi:yj
      + wi:yj %in% W, data = model.data, family = poisson(link = log)))
  }
  if (model == "y.main") {
    mod.fit <- suppressWarnings(glm(formula = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y + wi:yj
      + wi:yj %in% Y, data = model.data, family = poisson(link = log)))
  }
  if (model == "wy.main") {
    mod.fit <- suppressWarnings(glm(
      formula = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y + wi:yj
        + wi:yj %in% W + wi:yj %in% Y, data = model.data,
      family = poisson(link = log)
    ))
  }
  if (model == "saturated") {
    mod.fit <- suppressWarnings(glm(
      formula = count ~ -1 + W:Y + wi %in% W:Y + yj %in% W:Y + wi:yj
        + wi:yj %in% W + wi:yj %in% Y + wi:yj %in% W:Y,
      data = model.data, family = poisson(link = log)
    ))
  }

  # if (class(model)=="formula")
  if (inherits(model, "formula")) {
    mod.fit <- suppressWarnings(glm(formula = model, data = model.data, family = poisson(link = log)))
    formula.char <- Reduce(paste, deparse(model, width.cutoff = 499L))
    mod.fit$call <- paste("glm(formula =", formula.char, ", family = poisson(link = log), data = model.data)",
      collapse = " "
    )
  }
  ###################
  # options(warn = 0)
  output <- mod.fit
  if (limit.output) {
    output <- c(mod.fit$fitted.values, mod.fit$deviance)
    if (is.null(model.vars)) {
      output <- mod.fit$fitted.values
    }
  }
  output
}
