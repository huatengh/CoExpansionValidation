#' Do ABC simulation with a configuration table
#'
#' A function for doing abc simulation with a configuration table. see Readme on github for required columns in the table
#'
#' @param npod An integer, the number of simulation replicates
#' @param conf A data frame with configurations
#' @param time.range A vector containing two numeric elements specifying the minimal and maximum event time.
#' @param buffer The minimal amount of time separating two events. The default is 0, when events' time are sampled randomly
#' @param concentrationscale A numeric. The \code{scale} for the gamma distribution. See 'Details'
#' @param concentrationShape A numeric. The \code{shape} for the gamma distribution. See 'Details'
#' @param prefix BayeSSC will generate many intermediate files, this function will create temporary folders for storing these files. One folder for one row in the configuration table, and the folders are named as \code{prefix}\code{simulation replicates}_row number.If folders with the same name already exist, user will be prompt for providing another \code{prefix} or overwriting the folders.
#' @param BayeSSCallocation A character string providing the location of the BayeSSC executable (including the file name)
#' @param do.parallel A integer for the number of parallel threads to use. if equals \code{1},no parallel execution. The parallel R package required
#' @param write.reference.file \code{TRUE} or \code{FALSE}(Default). If \code{TRUE}, the simulated hyperstat will be append to a file \code{prefix}_reference_table. Useful for running large number of replicates.If \code{False}, the simulated hyperstat will be collected and returned as a data frame
#'
#' @return If \code{write.reference.file} is \code{FALSE}, this function returns a data frame containing the simulated hyperstats. The first column is 'uid' (\code{prefix}+the replicates number), the second column is the number of events, and the third column is the total number of species. Following are columns of the "real" expansion time for each species in the simulation, and then, the columns of hyperstats.If \code{write.reference.file} is \code{TRUE}, this function returns a character string providing the path to the reference file.
#'
#' @details For the prior on the number of co-expansion events, hBayeSSC used a flat distribution. In the paper, we adopt the [PyMsbayes-style](http://joaks1.github.io/PyMsBayes/) prior: a gamma distribution for the alpha parameter of the dirichelete process. User need to select the two parameters for gamma distribution: concentrationShape and consentrationScale. Check [prior selection in PyMsbayes](http://joaks1.github.io/PyMsBayes/tutorials/selecting-priors.html) for how to select these two parameters.
#'
#' @export
ABC_simulation_with_conf<-function(npod,conf,time.range,buffer=0,concentrationscale,concentrationShape,prefix='temp',BayeSSCallocation,do.parallel=1,write.reference.file=F){

  fx<-function(i,conf,time.range,buffer,concentrationscale,concentrationShape,prefix,BayeSSCallocation,write.reference.file){
    species<-as.character(unique(conf$species))
    nspecies<-length(species)
    alpha<-rgamma(1,shape = concentrationShape,scale = concentrationscale)
    event<-generate_coevent_number_dirichlet(species = species,alpha = alpha)
    exp.time<-generate_cotime_with_buffer(time.range = time.range,nco.events = length(event),buffer = buffer)
    simulatedobs<-runbayeSSC_with_conf(BayeSSCallocation = BayeSSCallocation,conf = conf,prefix = paste0(prefix,i),species.assignment = event,exp.time = exp.time)
    simulatehyperstat<-calculate_hyperstat(simulatedobs)
    a<-rep(0,3)
    a[1]<-paste0(prefix,i)
    a[2]<-length(exp.time)
    a[3]<-nspecies
    names(a)<-c("uid","nevent","nspecies")
    species.time<-species_exp_time(species.assignment = event,exp.time = exp.time)
    a<-c(a,species.time)
    simulatehyperstat<-c(a,simulatehyperstat)
    if(write.reference.file==T){
      hsfile<-paste0(prefix,"_reference_table")
      cat(paste(simulatehyperstat,sep='',collapse = "\t"),"\n",sep='',file = hsfile,append = T)
      return(hsfile)
    }else{
      return(simulatehyperstat)
    }
  }

  if(write.reference.file==T){
    species<-as.character(unique(conf$species))
    nspecies<-length(species)
    data("reference.table")
    colheads<-c("uid","nevent","nspecies")
    colheads<-c(colheads,paste0("species",1:nspecies),colnames(reference.table)[grep("haptypes_Mean",colnames(reference.table)):dim(reference.table)[2]])
    hsfile<-paste0(prefix,"_reference_table")
    cat(paste(colheads,sep='',collapse = "\t"),"\n",sep='',file = hsfile,append = F)
  }

  if(do.parallel==1){
    x<-sapply(1:npod,function(i)fx(i,conf,time.range,buffer,concentrationscale,concentrationShape,prefix,BayeSSCallocation,write.reference.file))
  }else{
    cluster <- parallel::makeCluster(do.parallel)
    parallel::clusterExport(cl = cluster, varlist=ls("package:CoExpansionValidation"))
    x<-parallel::parSapply(cluster, X = 1:npod, fx, conf,time.range,buffer,concentrationscale,concentrationShape,prefix,BayeSSCallocation,write.reference.file)
    parallel::stopCluster(cluster)

  }
  if(write.reference.file==F){
    x<-t(x)
    return(x)
  }else{x[1]}

}
