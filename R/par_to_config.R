#' Extract information from a par file
#'
#' A function extracting information from a par file and put it into a configuration data frame. Note that because we only consider one historical expansion event, the par file format is semi-fixed, and a template par file is included in the package for using as a template.
#'
#' @param bayessc_par_file A character string providing the path to the .par file
#' @param species A character vector containing the species names
#' @param nloci A integer. The number of loci to simulate for each species
#' @param gen A vector of integers containing the generation time for each species. If it is of length 1, function will assume the same generation time for all species. If it is \code{NA}(Default), assuming one generation per time unit.
#' @return A configuration data frame
#' @export
par_to_config<-function(bayessc_par_file,species,nloci,gen=NA){
  parfile<-scan(bayessc_par_file,what="character",sep="\n")
  if(length(parfile)!=21){stop("the number of lines in the template file can not be altered!\n")}

  if(is.na(gen))gen<-1
  if(length(gen)==1)gen<-rep(gen,length(species))
  conf<-data.frame(species=species,locinum=rep(nloci,length(species)),gen=gen,stringsAsFactors = F)
  conf$Ne<-parfile[4]
  conf$nsam<-parfile[6]
  conf$nsite<-parfile[17]
  conf$tstv<-gsub("^\\D+\\s","",parfile[19])
  conf$gamma<-parfile[21]
  l<-strsplit(parfile[13],split = ' ')
  conf$popratio<-l[[1]][5]
  conf$mutationrate<-parfile[15]
  return(conf)
}
