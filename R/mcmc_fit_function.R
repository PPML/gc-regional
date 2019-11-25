###########################################################
###     function to sample from mcmc parameter sets     ###
###########################################################

#to compare model outputs to calibration targets
#trace = burned and thinned mcmc trace
#sample size = number of samples to draw (want 1000+ for actual analysis)
#' @export
model_fits<- function(trace, sample.size, use_trace_without_sampling = F){
  if (use_trace_without_sampling) {
    ind <- 1:nrow(trace)
  } else {
    sample.size <- min(c(sample.size, nrow(trace)))
    ind <- sample(1:nrow(trace), sample.size, replace = TRUE)
    names(ind) <- ind
  }
  outputs = NULL
  theta.list = NULL
  for(i.iter in 1:length(ind)) {
    theta.list <- rbind(theta.list, trace[ind[i.iter],colnames(trace)])
  }
  for(i.iter in 1:nrow(theta.list)) {
    res<-model_epi_loglik(unlist(theta.list[i.iter,]))
    outputs <- rbind(outputs, c(i.iter,res))
  }
  return(list(theta.list=theta.list, outputs=outputs))
}
