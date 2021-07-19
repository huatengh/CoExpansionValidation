#' Calculate hyperstats from observation stats
#'
#' A function for calculating hyperstats from BayeSSC simulated stats. Here we calculate four hyperstats for each summary statistics: mean, variance, skewness and kurtosis as explained on https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
#'
#' @param obs A data frame containing the BayeSSC simulated summary statistics
#' @param hbayessclike \code{TRUE} (Default) or \code{FALSE}. If \code{TRUE}, directly calculate the hyperstat across all locus in all species. If \code{FALSE}, still need develop
#' @param col A character vector providing the names of columns for hyperstat calculation. Currently this parameter is not in use.
#'
#' @return A data frame containing the calculated hyperstats
#' @export
calculate_hyperstat<-function(obs,hbayessclike=T,col=NULL){
  if(hbayessclike==T){
    col<-c("Haptypes","HapDiver","nucdiv","TajimasD","F(\\*|\\.)$","pairdiffs","SegSites")
    colnum<-sapply(col,function(i)grep(i,colnames(obs),ignore.case = T))
    obs<-as.data.frame(obs)
    for(i in colnum) obs[,i]<-as.numeric(obs[,i])
    newcol<-c("haptypes","hapdiv","nucdiv","tajimasD","fusf","pairdiffs","segSites")
    x<-sapply(1:length(col),function(i)mean_var_skewness_kurtosis(obs[,colnum[i]]))
    x<-as.vector((x))
    statnames<-c("Mean","Variance","Skewness","Kurtosis")
    names(x)<-paste(rep(newcol,each=4),rep(statnames,7),sep='_')
    #x<-c(rep('nan',16),x)
    return(x)
  }
}
