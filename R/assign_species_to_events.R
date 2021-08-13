#' Assign species to events
#'
#' This function puts species into co-expansion events randomly or evenly
#' The number of co-expansion events is specified by the user
#'
#' @param species A vector of species names, integers are treated as characters
#' @param nco.events Integer. The number of co-expansion events
#' @param even \code{TRUE} or \code{FALSE}(Default). The default is to assign species randomly. IF \code{TRUE}, will spread species evenly acrosss events
#' @return A list, the length of which is equal to the number of events. Each element of the list is a vector of species names assigned to the event.
#'
#' @export
assign_species_to_events<-function(species,nco.events,even=F){
  #a function for assigning species to events
  divevent<-list()
  nspecies<-length(species)
  species_to_events<-rep(0,nspecies)
  #if Even=TRUE, then spread species among divergence events evenly
  if(even==T){
    x<-nspecies/nco.events
    for(i in 1:nco.events){
      divevent[[i]]<-species[(floor(x*(i-1))+1):floor(x*i)]
      if(i==nco.events)divevent[[i]]<-species[(floor(x*(i-1))+1):nspecies]
      #print(i)
      #print(divevent)
      #scan()
    }

    return(divevent)
  }

  #make sure each events have at least one species
  l<-sample(1:nspecies,nco.events,replace = F)
  species_to_events[l]<-1:nco.events
  #assign the rest of species randomly
  species_to_events[-l]<-sample(1:nco.events,size = nspecies-nco.events,replace = T)
  for(i in 1:nco.events){
    divevent[[i]]<-species[which(species_to_events==i)]
  }
  return(divevent)
}
