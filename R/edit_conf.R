#' Sampling the configuration table for each species
#'
#' A function for randomly draw values from prior distributions in the configuration table. BayeSSC can take prior distribution for parameters as well as fixed values. However, loci of the same species need to share the same generation time, effective population size and expansion history. So these values need to be fixed for each species instead of letting BayeSSC to draw a different sample for each locus.
#'
#' @param conf A data frame with configurations
#' @param col A character vector containing the name of the columns, for which values need to be sampled from prior distributions. The default contains three columns in the \code{conf} data frame: "gen"(generation time), "Ne"(effective population size), and "popratio"(the population expansion ratio)
#' @return A configuration data frame with randomly sampled values from prior distribution
#' @export
edit_conf<-function(conf,col=c("gen","Ne","popratio")){
  x<-conf[,c("species",col)]
  #print(x)
  #check to see if there is distributions
  x<-x[!duplicated(x$species),]
  if(length(grep("\\{",x))==0)return(conf)#if no distribution is specified in the columns, return the original
  edited_conf<-conf
  for(i in 1:dim(x)[1]){
    for(j in col){
      y<-0
      while(y<=0) y<-distribution_prior(text = x[i,j],sample = 1)
      edited_conf[edited_conf$species==x$species[i],j]<-y
    }
  }
  for(i in 1:dim(edited_conf)[2]){
    if(colnames(edited_conf)[i]=="gamma")next;
    edited_conf[,i]<-as.numeric(edited_conf[,i])
  }

  return(edited_conf)
}
