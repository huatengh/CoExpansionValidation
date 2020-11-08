#' Run BayeSCC with a configuration table
#'
#' A wrapper function for running BayeSCC with a configuration table. See the example in the vignettes for the configurations needed.
#'
#' @param BayeSSCallocation A character string providing the location of the BayeSSC executable (including the file name)
#' @param conf A data.frame with configurations
#' @param prefix BayeSSC will generate some intermediate files, this function will create temporary folders for storing these files. One folder for one row in the configuration table, and the folder is named as \code{prefix}_row number.If folders with the same name already exist, user will be prompt for choosing another \code{prefix} or overwrite the folders.
#' @param species.assignment A list showing species-to-event assignment. This can be \code{NULL},if \code{eventtime.generation} is included as a column in the \code{conf} data frame.
#' @param exp.time A numeric vector with the events time.This can be \code{NULL},if \code{eventtime.generation} is included as a column in the \code{conf} data frame.
#' @param intern Default is \code{TRUE}: the screen output of BayeSSC will be suppressed. If \code{FALSE}, all screen output of BayeSSC will be shown
#' @param deletetempfolder \code{TRUE}(Default) or \code{FALSE}, whether the temporary folders for holding BayeSCC's intermediate files will be deleted after running
#' @param write.conf.to.file \code{TRUE} or \code{FALSE}(Default), whether to write the sampled configuration to a file. If \code{TRUE}, the file name will be \code{prefix}_conf_file
#' @param write.obs.to.file \code{TRUE} or \code{FALSE}(Default), whether to write the observation data simulated by BayeSSC to a file. If \code{TRUE}, the file name will be \code{prefix}_obs_file
#'
#' @return A data frame storing the the observation summary statistics simulated BayeSSC
#'
#' @export
runbayeSSC_with_conf<-function(BayeSSCallocation,conf,prefix='temp',species.assignment=NULL,exp.time=NULL,intern=TRUE,deletetempfolder=T,write.conf.to.file=F,write.obs.to.file=F){

  checkfolder=F
  while (checkfolder==F){
   currentsubfolders<-dir(path = ".",pattern = paste(prefix,"_\\d+$",sep=''))
   conflictingsubfolder<-currentsubfolders[as.integer(sub(prefix,"",currentsubfolders)) %in% 1:dim(conf)[1]]
   if(length(conflictingsubfolder)==0)break
   cat("Folder:\n")
   cat(paste(conflictingsubfolder,sep = '',collapse = "\n"))
   cat("already exist\n")
   print("enter C to clean up these folders")
   print("or another folder name:")
   x<-scan(what = "character",n = 1)
   if(x=='C'){
     unlink(conflictingsubfolder, recursive=TRUE)
   }else{
     prefix<-x
   }
  }
  #sample the prior distribution in conf
  #make sure one species only have one Ne, expansion ratio and generation time
  conf<-edit_conf(conf = conf)
  #get the event time in generation
  if(! "eventtime.generation" %in% colnames(conf)){
    if(is.null(species.assignment)|is.null(exp.time)){
      stop("eventtime not provided in configuration table\n need species-to-event assignment and time\n")
    }
    conf$eventtime.generation<-0
    sptime<-species_exp_time(species.assignment,exp.time)
    conf$eventtime.generation<-sptime[as.character(conf$species)]/conf$gen
  }
  if(write.conf.to.file==T){
  write.table(conf,sep="\t",row.names = F,quote=F,file=paste(prefix,"_conf_file",sep=''))
  }
  #run bayescc for every row
  csvfiles<-rep("",dim(conf)[1])
  for(i in 1:dim(conf)[1]){
    csvfiles[i]<-runbayeSCC_with_confrow(BayeSSCallocation,conf[i,],tempfolder=paste(prefix,i,sep='_'),intern=TRUE)
  }

  #collecting data together
  datalist <- lapply(csvfiles, function(x){read.csv(file=x,header=T,as.is = T)})
  coln<-table(unlist(lapply(datalist,colnames)))
  select.col<-names(coln)[coln==length(datalist)]
  datalist<-lapply(datalist,function(x)x[,select.col])
  myfulldata <- do.call("rbind", datalist)
  temp_matrix<-  matrix(0,nrow = dim(myfulldata)[1], ncol = 17)
  colnames(temp_matrix)<- c("species", "nsam", "nsites", "tstv", "gamma", "gen", "locuslow","locushigh","Nelow","Nehigh", "SegSites", "Nucdiv","Haptypes", "HapDiver", "Pairdiffs", "TajimasD", "F*" )

  sp<-c()
  for(j in 1:dim(conf)[1]){
    if(j==1 || conf$species[j] !=conf$species[j-1]){
      l<-1
    }
    sp<-c(sp,paste("sp",conf$species[j],"loci",l:(l+conf$locinum[j]-1),sep=''))
    l<-l+conf$locinum[j]
    #print(sp)
  }
  temp_matrix[1:dim(myfulldata)[1],"species"]<-sp
  temp_matrix[1:dim(myfulldata)[1],"nsam"] <- rep(conf$nsam,conf$locinum)
  temp_matrix[1:dim(myfulldata)[1],"nsites"] <- rep(conf$nsite,conf$locinum)
  temp_matrix[1: dim(myfulldata)[1],"tstv"] <- rep(conf$tstv,conf$locinum)
  temp_matrix[1: dim(myfulldata)[1],"gamma"] <- sub(" \\w+$","",rep(conf$gamma,conf$locinum))
  temp_matrix[1:dim(myfulldata)[1],"gen"] <- rep(sapply(conf$gen,function(i)distribution_prior(i,sample = 0,prob = 0.01)),conf$locinum)
  temp_matrix[1:dim(myfulldata)[1],"locuslow"] <- rep(sapply(conf$mutationrate,function(i)distribution_prior(i,sample = 0,prob = 0.01)),conf$locinum)
  temp_matrix[1:dim(myfulldata)[1],"locushigh"] <- rep(sapply(conf$mutationrate,function(i)distribution_prior(i,sample = 0,prob = 0.99)),conf$locinum)
  temp_matrix[1:dim(myfulldata)[1],"Nelow"] <-as.integer( rep(sapply(conf$Ne,function(i)distribution_prior(i,sample = 0,prob = 0.01)),conf$locinum))
  temp_matrix[1:dim(myfulldata)[1],"Nehigh"] <- as.integer(rep(sapply(conf$Ne,function(i)distribution_prior(i,sample = 0,prob = 0.99)),conf$locinum))
  temp_matrix[1:dim(myfulldata)[1],"SegSites"] <-myfulldata$SegSites
  temp_matrix[1:dim(myfulldata)[1],"Nucdiv"] <- myfulldata$NucltdDiv
  temp_matrix[1:dim(myfulldata)[1],"Haptypes"] <- myfulldata$Haptypes
  temp_matrix[1:dim(myfulldata)[1],"HapDiver"] <- myfulldata$HapDiver
  temp_matrix[1:dim(myfulldata)[1],"Pairdiffs"] <- myfulldata$PairDiffs
  temp_matrix[1:dim(myfulldata)[1],"TajimasD"] <- myfulldata$TajimasD
  temp_matrix[1:dim(myfulldata)[1],"F*"] <- myfulldata$F.
  if(write.obs.to.file==T){
  write.table(temp_matrix, file=paste(prefix,"_obs_file",sep=''), quote= FALSE, sep = '\t',row.names = F)
  }
  if(deletetempfolder==T){
    unlink(paste(prefix,1:dim(conf)[1],sep='_'),recursive = T)
  }
  return(temp_matrix)
}
