# Parallel functions

# from plyr


#' Set up parallel environment
#'
#' @return Nothing
setup_parallel <- function() {
  if (!requireNamespace("foreach", quietly = TRUE)) {
    # EXCLUDE COVERAGE START
    stop("foreach package required for parallel plyr operation",
         call. = FALSE)
    # EXCLUDE COVERAGE END
  }
  if (foreach::getDoParWorkers() == 1) {
    # EXCLUDE COVERAGE START
    warning("No parallel backend registered", call. = TRUE)
    # EXCLUDE COVERAGE END
  }
}

# if (.parallel) {
#   setup_parallel()
#
#   i <- seq_len(n)
#   fe_call <- as.call(c(list(quote(foreach::foreach), i = i), .paropts))
#   fe <- eval(fe_call)
#
#   result <- foreach::`%dopar%`(fe, do.ply(i))
# } else {
#   result <- loop_apply(n, do.ply)
# }

