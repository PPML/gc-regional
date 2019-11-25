# Helper Objects and Functions

# Convert a Vector from Real-Space to the log/logit Transformed Space
transform_theta <- function(theta) {
  stopifnot(length(theta) == length(gc_env$theta))

  vapply(1:length(theta), function(x) {
    gc_env$transforms[[strsplit(names(gc_env$theta)[[x]], "\\.")[[1]][[1]]]](
      theta[[x]]
    )
  }, 0)
}

# Convert a Vector from log/logit transformed space to Real space
untransform_theta <- function(theta) {

  out <- vapply(1:length(theta), function(x) {
    gc_env$inverse_transforms[[strsplit(names(theta)[[x]], "\\.")[[1]][[1]]]](
      theta[[x]]
    )
  }, 0)
  names(out) <- sapply(names(theta), function(x) paste0(strsplit(x, "\\.")[[1]][-1], collapse="."))
  return(out)
}


# Sitename object
sitenames <- c(SF = 'San Francisco', BA = 'Baltimore')
