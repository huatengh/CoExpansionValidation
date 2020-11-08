# a function for dealing with the prior distribution in the configuration or the par file
distribution_prior<-function(text,sample=1,prob=NA){
  if(length(grep("\\{",text))==0) return(text)
  distribution<-sub("\\{(\\w):.*","\\1",text)
  if(distribution=="U"){#uniform distritution
    min<-as.numeric(sub(".*:(\\d+\\S*),\\d+\\S*.*","\\1",text))
    max<-as.numeric(sub(".*:\\d+\\S*,(\\d+.*)\\s*\\}","\\1",text))
    if(sample>0){
      return(runif(n=sample,min=min,max=max))
    }else if(sample==0 && !(is.na(prob[1]))){
      return(qunif(p = prob,min=min,max=max))
    }
  }else if(distribution=='G'){#gamma distribution
    shape<-as.numeric(sub(".*:(\\d+\\S*),\\d+\\S*.*","\\1",text))
    scale<-as.numeric(sub(".*:\\d+\\S*,(\\d+.*)\\s*\\}","\\1",text))
    if(sample>0){
      return(rgamma(n=sample,shape = shape,scale = scale))
    }else if(sample==0 && !(is.na(prob[1]))){
      return(qgamma(p=prob,shape = shape,scale = scale))
    }

  }else if(distribution=="M"){#geom distribution
    p<-as.numeric(sub(".*:(\\d+.*)\\s*\\}","\\1",text))
    if(sample>0){
      return(rgeom(n=sample,prob = p))
    }else if(sample==0 && !(is.na(prob[1]))){
      return(qgeom(p=prob,prob = p))
    }

  }else if(distribution=='N'){#normal distribution
    mean<-as.numeric(sub(".*:(\\d+\\S*),\\d+\\S*.*","\\1",text))
    sd<-as.numeric(sub(".*:\\d+\\S*,(\\d+.*)\\s*\\}","\\1",text))
    if(sample>0){
      return(rnorm(n=sample,mean=mean,sd=sd))
    }else if(sample==0 && !(is.na(prob[1]))){
      return(qnorm(p=prob,mean=mean,sd=sd))
    }

  }else if(distribution=='E'){#exponential distribution
    scale<-as.numeric(sub(".*:(\\d+.*)\\s*\\}","\\1",text))
    if(sample>0){
      return(rexp(n=sample,rate = 1/scale))
    }else if(sample==0 && !(is.na(prob[1]))){
      return(qexp(p=prob,rate = 1/scale))
    }

  }else {
    stop("Distribution not supported\n")
  }
}
