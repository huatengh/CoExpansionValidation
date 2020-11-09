#' Generate the number of events according to the dirichlet process
#'
#' This function puts species into coexpansion events according to dirichlet process.
#' It will go through the list of species one by one, and determine whether each one should
#' be assigned to a existing event or to a new event, with probabilities determined by the alpha
#' parameter.
#'
#' @param species A vector of species names
#' @param alpha A numeric. The parameter alpha for the dirichlet process
#' @param maxevent An integer for the maximum number of event allowed. The default is NA-- no limits on the maximum number of events.
#' @return A list, the length of which is equal to the number of events. Each element of the list is a vector of species names assigned to the event
#'
#' @export
generate_coevent_number_dirichlet<-function(species,alpha,maxevent=NA){
  #this is a function for generating the number of divergence events
  #and the assignment of species to divergence events
  #according to the dirichlet process given alpha
  nspecies<-length(species)
  if(is.na(maxevent)){maxevent=nspecies}
  divevent<-vector(mode = "list", length = maxevent+1)
  while(length(divevent)>maxevent){
    divevent<-list()
    divevent[[1]]<-c(species[1])
    n.event<-rep(1,nspecies)
    for(v in 2:nspecies){
      prob<-unlist(lapply(divevent, function(x){length(x)/(alpha+v-1)}))
      prob<-c(prob,alpha/(alpha+v-1))
      cumprob<-cumsum(prob)
      draw<-runif(1,min=0,max=cumprob[length(cumprob)])
      eventnum<-min(which(cumprob>draw))
      if(eventnum>length(divevent)){
        divevent[[eventnum]]<-c(species[v])
      }else{
        divevent[[eventnum]]<-c(divevent[[eventnum]],species[v])
      }
      n.event[v]<-eventnum
    }
  }
  return(divevent)
}
