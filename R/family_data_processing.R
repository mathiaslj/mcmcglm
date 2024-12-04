#' Converts string or function to a family object
#' @noRd
check_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  return(family)
}

#' Uses formula to retrieve response variable and model matrix
#' @noRd
extract_model_data <- function(formula, data) {
  mf <- stats::model.frame(formula, data = data)
  mt <- attr(mf, "terms")

  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    as.matrix(model.matrix(mt, mf))
  else matrix(, NROW(Y), 0L)

  return(list(Y = Y, X = X))
}
