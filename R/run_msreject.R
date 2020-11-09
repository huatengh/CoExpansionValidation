#' Run MsReject
#'
#' This is wrapper function for running the msReject in msbayes.For installing/compiling msReject see [msbayes webpage](http://msbayes.sourceforge.net/).
#'
#' @param hyperstat A vector containing the "observed" hyperstats from step 1
#' @param reference.table A data frame containing the hyperstats from ABC simulation(step 2). The column name of this data frame must match the names of \code{hyperstat} vector
#' @param MsRejectallocation A character string. The location of the msReject (including the executable itself)
#' @param samplingtolerance Numeric. The sampling fraction for msReject
#' @param prefix A character string for writing the "observed" hyperstat and reference table to files, \code{prefix}_hyperstat_file and \code{prefix}_reference_table, respectively
#' @param hbayeSSC.style \code{TRUE}(Default) or \code{FALSE}.In hBayeSSC, running msReject involves using the hyperstats of haptypes, hapdiv, nucdiv, tajimasD.
#' @param write.posterior.file \code{TRUE} or \code{FALSE}(Default). If \code{TRUE}, the retained abc simulation replicates will be write to a file named \code{prefix}_Posterior. If \code{False}, they will be returned as a data frame
#'
#' @return If \code{write.posterior.file} is \code{FALSE}, this function returns a data frame containing the retained rows from reference table after rejection. The first column is 'uid', the replicates number, the second column is the number of events, and the third column is the total number of species. Following are columns for the "real" expansion time for each species in the simulation, and then, the columns of hyperstats.If \code{write.posterior.file} is \code{TRUE}, this function returns the file name (including the path) of the posterior file.
#'
#' @export
run_msreject<-function(hyperstat,reference.table,MsRejectallocation,samplingtolerance,prefix="temp",hbayeSSC.style=T,write.posterior.file=F){

  hsfile<-paste0(prefix,"_hyperstat_file")
  cat(paste(hyperstat,sep='',collapse = "\t"),"\n",sep='',file = hsfile)
  reffile<-paste0(prefix,"_reference_table")
  write.table(reference.table,file=reffile,sep="\t",col.names = F,row.names = F,quote=F)
  outposteriorname<-sub("hyperstat_file","Posterior",hsfile)
  if(hbayeSSC.style==T){
    colnum<-grep("haptypes_Mean",names(hyperstat)):grep("tajimasD_Kurtosis",names(hyperstat))
  }
  command<-paste(MsRejectallocation,hsfile,reffile, samplingtolerance,paste(colnum,sep='',collapse = ' '),">", outposteriorname, sep=' ')
  system(command)
  if(write.posterior.file==F){
    posteriortable<-read.table(outposteriorname,sep="\t")
    colnames(posteriortable)<-colnames(reference.table)
    return(posteriortable)
    unlink(outposteriorname)
  }else{
    return(outposteriorname)
  }
}
