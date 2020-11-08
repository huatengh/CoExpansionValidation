#' Sampling the configuration table for each species
#'
#' A function for randomly draw values from distributions included in the configuration table. BayeSSC can take prior distribution for parameters instead of fixed values. However, loci of the same species need to share the same generation time, effective population size and expansion history. So these values need to be fixed for each species instead of letting BayeSSC to draw a difference sample for each locus.
#'
#' @param conf A data frame with configurations
#' @param col A character vector containing the columns for which the value need to be sampled from distributions. The default contains three columns in the \code{conf} data frame: "gen"(generation time), "Ne"(effective population size), and "popratio"(the population expansion ratio)
#' @return A configuration data frame with randomly sampled values
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

  return(edited_conf)
}
