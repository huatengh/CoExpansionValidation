mean_var_skewness_kurtosis<-function(v){
#https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
  n<-M1<-M2<-M3<-M4<-0

  for (x in v){
    n1 = n
    n = n + 1
    delta = x - M1
    delta_n = delta / n
    delta_n2 = delta_n * delta_n
    term1 = delta * delta_n * n1
    M1 = M1 + delta_n
    M4 = M4 + term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3
    M3 = M3 + term1 * delta_n * (n - 2) - 3 * delta_n * M2
    M2 = M2 + term1
  }
  mean=M1
  variance=M2/(n-1)
  skewness=sqrt(n) * M3/ M2^1.5
  kurtosis = (n * M4) / (M2 * M2) - 3
  return(c(mean,variance,skewness,kurtosis))
}
