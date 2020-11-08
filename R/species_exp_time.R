#' Get expansion time for each species
#'
#' A short function for get expansion time for each species
#'
#' @param species.assignment a list showing the species-to-event assignment
#' @param exp.time The time of the co-expansion events
#' @return A numeric vector containing the expansion time for each species.
#' @seealso \code{\link[CoExpansionValidation]{assign_species_to_events}}
#' @export
species_exp_time<-function(species.assignment,exp.time){
  nsp<-length(unlist(species.assignment))
  sptime<-c()
  for(l in 1:length(species.assignment)){
    sptime<-c(sptime,rep(exp.time[l],length(species.assignment[[l]])))
  }
  names(sptime)<-unlist(species.assignment)
  sptime<-sptime[order(names(sptime))]
  return(sptime)
}
