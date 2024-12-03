check_family_add_class <- function(family) {
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

  family_with_dist_as_class <- structure(family, class = c(family$family, class(family)))

  return(family_with_dist_as_class)
}

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

# extract_modelframe_from_call <- function(call) {
#
#   m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
#              names(call),
#              0L)
#   call <- call[c(1L, m)]
#   call$drop.unused.levels <- TRUE
#   call[[1L]] <- quote(stats::model.frame)
#   mf <- eval(call, parent.frame())
#
#   return(mf)
# }
