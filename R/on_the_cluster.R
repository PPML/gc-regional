#' A Function to Test if we're on the Harvard cluster.
on_the_cluster <- function() {
  e <- Sys.info()
  return(
    grepl("Linux", e[['sysname']]) &&
      grepl("harvard.edu", e[['nodename']]))
}
