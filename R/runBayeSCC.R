runbayeSCC<-function(BayeSSCallocation, parfilename,nloci,intern=TRUE){
  command<-paste(BayeSSCallocation,"-f",parfilename,nloci, sep=' ')
  system(command,intern=TRUE)
  #unlink("*.distr")
  #unlink("*.sum")
  #unlink("*.trees")
  #unlink("*.gen")
}
