runbayeSCC_with_confrow<-function(BayeSSCallocation,confrow,tempfolder='temp',intern=TRUE){

# a function for running  BayeSSC with only one row of the configuation rable
# return the path to the simulated stat file
  if(dir.exists(tempfolder)){
    print("output folder for BayeSSC already exists")
    print("enter C to clean up this folder")
    print("or another folder name:")
    x<-scan(what = "character",n = 1)
    if(x=='C'){
      unlink(tempfolder, recursive=TRUE)
    }else{
      tempfolder<-x
    }
  }
  dir.create(tempfolder)
  file.copy(BayeSSCallocation,tempfolder)
  old.dir<-unlist(stringr::str_split(BayeSSCallocation,"\\/"))
  if(old.dir[1]=='.'){
    new.dir<-paste(old.dir[1],old.dir[length(old.dir)],sep="/")
  }else{
    new.dir<-old.dir[length(old.dir)]
  }
  currentdir<-getwd()
  setwd(tempfolder)
  parfile<-write_parfile(confrow)

  runbayeSCC(BayeSSCallocation = new.dir,parfilename = parfile,nloci = confrow$locinum,intern=TRUE)
  setwd(currentdir)
#  unlink(paste(tempfolder,"*.distr",sep="/"))
#  unlink(paste(tempfolder,"*.trees",sep="/"))
#  unlink(paste(tempfolder,"*.gen",sep="/"))
#  unlink(paste(tempfolder,"*.sum",sep="/"))
#  unlink(paste(tempfolder,"*.par",sep="/"))
  unlink(paste(tempfolder,old.dir[length(old.dir)],sep="/"))

  return(paste(tempfolder,sub(".par","_stat.csv",parfile),sep='/'))
}
