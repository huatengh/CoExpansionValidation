write_parfile<-function(confrow,tempfolder=NULL){
  parfile<-templateparfile
  parfile[4]<-format(confrow$Ne,scientific = F)
  parfile[6]<-confrow$nsam
  l<-unlist(strsplit(parfile[13],split = ' '))
  l[5]<-confrow$popratio
  l[1]<-floor(confrow$eventtime.generation)
  parfile[13]<-paste(l,sep='',collapse = ' ')
  parfile[15]<-confrow$mutationrate

  parfile[17]<-confrow$nsite
  parfile[19]<-paste0("DNA ",confrow$tstv)
  parfile[21]<-confrow$gamma
  if(is.null(tempfolder)){
    outfile="row.par"
  }  else{outfile=paste(tempfolder,"row.par",sep='/')}
  sink(outfile,append = F)
  cat(paste(parfile,sep = '',collapse = '\n'))
  cat("\n",sep='')
  sink()
  return(outfile)


}
