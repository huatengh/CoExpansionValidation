#' Generate co-expansion time
#'
#' This function randomly draw events' time from a time range.
#' User can specify the minimal time length (buffer) between events
#'
#' @param time.range A vector containing two numeric elements specifying the minimal and maximum event time.
#' @param nco.events The number of coexpansion events. See 'Details'.
#' @param buffer The minimal amount of time separating two events. The default is 0, when events' time are sampled randomly
#' @details Note that the \code{buffer} is the time length on each side of a expansion time. Also, if the \code{buffer} is non-zero,the maximum number of events (max \code{nco.event})=\code{time.range}/\code{buffer} +2. However, if \code{nco.events} is close to the it maximum (even if not exceeding), this function might take a while to find a set of events time that are at least buffer apart from each other.
#' @return A numeric vector of the sampled expansion time
#'
#' @export
generate_cotime_with_buffer<-function(time.range,nco.events,buffer=0){
  #this is a function for generating divergence time for testing data
  #note that buffer is the buffer time on *each* side of a divergence event
  #so the maximum number of nco.events=time.range/buffer+2
  #but may need to retry many times with large nco.events even if it does not exceed the maximum
  #if the nco.events is close to its maximum value, this scripts will take a long time to run.
  if(nco.events>(time.range[2]-time.range[1])/buffer+2){
    stop("too many nco.events\n")
  }
  #if buffer is 0, just randomly sample the time
  if(buffer==0){
    divtime<-runif(min=time.range[1],max=time.range[2],n = nco.events)
    return(divtime)
  }

  n<-0
  while (n<nco.events){
    divtime<-rep(0,nco.events)
    divtime[1]<-runif(min=time.range[1],max=time.range[2],n = 1)
    bufferzone<-matrix(c(max(time.range[1],divtime[1]-buffer),min(time.range[2],divtime[1]+buffer)),nrow = 1)
    n<-1
    for(i in 2:nco.events){
      timelength<-time.range[2]-time.range[1]-sum(bufferzone[,2]-bufferzone[,1])
      if(timelength<=0){break;}
      dt<-runif(min=time.range[1], max=(time.range[1]+timelength),n = 1)
      j<-1
      while(j<=dim(bufferzone)[1] && dt>bufferzone[j,1]){
        dt<-dt+bufferzone[j,2]-bufferzone[j,1]
        j<-j+1
      }
      divtime[i]<-dt
      bufferzone<-rbind(bufferzone,c(max(dt-buffer,time.range[1]),min(time.range[2],dt+buffer)))
      bufferzone<-bufferzone[order(bufferzone[,1],bufferzone[,2]),]
      newbufferzone<-bufferzone[1,,drop=F]
      for(j in 2:dim(bufferzone)[1]){
        k<-dim(newbufferzone)[1]
        if((newbufferzone[k,2]-bufferzone[j,1])*(bufferzone[j,2]-newbufferzone[k,1])>0){
          newbufferzone[k,1]<-min(newbufferzone[k,1],bufferzone[j,1])
          newbufferzone[k,2]<-max(newbufferzone[k,2],bufferzone[j,2])
        }else{
          newbufferzone<-rbind(newbufferzone,bufferzone[j,])
        }
        j<-j+1
      }
      bufferzone<-newbufferzone
      n<-n+1
    }
    if(n<nco.events)print("sampling failed, trying another time\n")
  }
  return(divtime)
}
