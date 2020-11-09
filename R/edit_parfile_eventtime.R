
edit_parfile_eventtime<-function(templateparfilename,outparfilename,eventtime,overwrite=T){
  templateparfile<- scan(file = templateparfilename,what = 'character',sep = "\n")
  lline<-grep("1 historical event",templateparfile)
  templateparfile[lline+1]<-sub("^\\d+(.*)$","\\1",templateparfile[lline+1])
  templateparfile[lline+1]<-paste(eventtime,templateparfile[lline+1],sep='')
  sink(outparfilename)
  for(n in templateparfile){
    cat(n,"\n",sep='')
  }
  sink()
}
