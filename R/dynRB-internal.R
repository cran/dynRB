.trapz <- 
function(x,y){
  n <- length(x)
  return(sum((y[1:(n-1)]+y[2:n])/2*(x[2:n]-x[1:(n-1)])))
}
.intersection_onecolumn <-
function(alpha, SSN1, SSN2, nrow_S1, nrow_S2, Ncol, .quantile_intersection){
  SP_inters <- .quantile_intersection(vec1 = SSN1, vec2 = SSN2, a = alpha/2, b = 1 - alpha/2, n = nrow_S1)
  vol_intersect <- (SP_inters[6])
  vol_B <- (SP_inters[3])
  r <- 1
  if (min(SP_inters[4] == SP_inters[1]) == 0 | min(SP_inters[5] == SP_inters[2]) == 0){
    if (vol_intersect < vol_B) {r <- vol_intersect/vol_B}
    if (vol_B == 0) {r <- 0}
  }
  return(r)
}
.intersection_severalcol <-
function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol, .quantile_intersection){
  SP_inters<-mapply(.quantile_intersection, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=((1-(1-alpha)^(1/ncol(SSN1)))/2),  b=(1-(1-(1-alpha)^(1/ncol(SSN1)))/2),n=nrow_S1))  
  product<-prod(SP_inters[6,])
  vol_intersect<-c(product,mean(SP_inters[6,]),product^{1/length(SP_inters[6,])})  
  product<-prod(SP_inters[3,])
  vol_B<-c(product,mean(SP_inters[3,]),product^{1/length(SP_inters[3,])})  
  r<-c(1,1,1)                                                                                             
  if(min(SP_inters[4,]==SP_inters[1,])==0 | min(SP_inters[5,]==SP_inters[2,])==0){                  
    r[vol_intersect<vol_B]<-vol_intersect[vol_intersect<vol_B]/vol_B[vol_intersect<vol_B]
    r[vol_B==0]<-0                                                                             
  }
  return(r)                                                                                        
}
.portionAinB_coordinates_full <-
  function(S1,S2,steps=101 ){ 
    nt<-ncol(S1)
    integral_coord<-rep(0,length=nt)
    portionAinB_function<-function(x1,x2,steps){
      r<-.portionAinB_full_onecolumn(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
      return(r$integral_approx)
    }
    portionAinB_function2<-function(x1,x2,steps){
      r<-.portionAinB_full_onecolumn(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
      return(r$overlap)
    }
    integral_coord<-mapply(portionAinB_function, x1=S1, x2=S2, MoreArgs =list(steps=steps))
    plot_data_overlap<-mapply(portionAinB_function2, x1=S1, x2=S2, MoreArgs =list(steps=steps))
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_data_overlap=plot_data_overlap)
    return(erg)
  }
.portionAinB_coordinates_full_dim1 <-
function(S1,S2,steps=101 ){ 
  nt<-1
  integral_coord<-rep(0,length=nt)
  r<-.portionAinB_full_onecolumn(data.frame(v1=S1) ,data.frame(v1=S2),steps=steps  )
  integral_coord<-r$integral_approx
  plot_data_overlap<-r$overlap
  alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
  erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_data_overlap=plot_data_overlap)
  return(erg)
}
.portionAinB_full_onecolumn <-
function (S1, S2, steps = 101) {
  S <- rbind(S1, S2)
  Ncol <- ncol(S1)
  alpha_grid <- seq(0, 1, length = steps)
  Spans <- data.frame(trait_nr = 1:Ncol, min = min(S, na.rm = TRUE), max = max(S, na.rm = TRUE))
  ab <- Spans$max - Spans$min
  SSN <- (S - Spans$min)/(ab)
  SSN[, ab == 0] <- 0.5
  nrow_S1 <- nrow(S1)
  nrow_S2 <- nrow(S2)
  z <- mapply(.intersection_onecolumn, alpha = alpha_grid,  MoreArgs = list(SSN1 = SSN[1:nrow_S1, ], SSN2 = SSN[(nrow_S1 + 1):(nrow_S1 + nrow_S2), ], nrow_S1 = nrow_S1, nrow_S2 = nrow_S2, Ncol = Ncol, .quantile_intersection = .quantile_intersection))
  integral_approx <- c(prod = .trapz(alpha_grid, z))
  erg <- list(alpha_grid = alpha_grid, overlap = z, integral_approx = integral_approx)
  return(erg)
}
.portionAinB2_full <-
function(S1,S2,steps=101,alpha_grid){
  steps0<-steps
  alpha_grid0<-alpha_grid
  .portionAinB2_full_severalcol(S1,S2,steps=steps0,alpha_grid=alpha_grid0)
}
.portionAinB2_full_dim1 <-
  function(S1,S2,steps=101,alpha_grid){
    steps0<-steps
    alpha_grid0<-alpha_grid
    .portionAinB2_full_severalcol_dim1(S1,S2,steps=steps0,alpha_grid=alpha_grid0)
}
.portionAinB2_full_severalcol_dim1 <- 
  function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)[1:(steps-1)]){
    S<-c(S1,S2)
    alpha_grid<-seq(0,1,length=steps) 
    Ncol<-1
    ab<-max(S)-min(S)         
    SSN<-t((t(S)-min(S))/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    nrow_S1<-length(S1)
    nrow_S2<-length(S2)
    z<-mapply(.intersection_onecolumn, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,.quantile_intersection=.quantile_intersection)) 
    integral_approx<-c(prod=.trapz(seq(0,1,length=steps), z),mean=.trapz(seq(0,1,length=steps),z),gmean=.trapz(seq(0,1,length=steps),z)) 
    plot_data_prod<-z
    erg<-list(alpha_grid=seq(0,1,length=steps)[1:(steps-1)],overlap=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod)                                                                                             
    return(erg)
  }
.portionAinB2_full_severalcol <-
function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)[1:(steps-1)]){
  S<-rbind(S1,S2)
  alpha_grid<-seq(0,1,length=steps) 
  Ncol<-ncol(S1)
  Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
  ab<-Spans$max-Spans$min          
  SSN<-t((t(S)-Spans$min)/(ab))    
  SSN[,ab==0]<-0.5                 
  SSN<-as.data.frame(SSN)          
  nrow_S1<-nrow(S1)
  nrow_S2<-nrow(S2)
  z<-mapply(.intersection_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,.quantile_intersection=.quantile_intersection)) 
  integral_approx<-c(prod=.trapz(seq(0,1,length=steps), z[1,]),mean=.trapz(seq(0,1,length=steps),z[2,]),gmean=.trapz(seq(0,1,length=steps),z[3,])) 
  plot_data_prod<-z[1,]
  erg<-list(alpha_grid=seq(0,1,length=steps)[1:(steps-1)],overlap=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod)                                                                                             
  return(erg)
}
.quantile_intersection <-
function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs = c(a,b), na.rm = TRUE),  quantile(vec2, probs = c(a,b), na.rm = TRUE) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                        
  return(r)                                                            
}
.volume_onecol <-
function(alpha,SSN1){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2),na.rm = TRUE)  
  vol<-(SP1[2]-SP1[1])  
  return(vol)                                                                                      
}
.volume_severalcol <-
function(alpha,SSN1){
  alpha<-(1-(1-alpha)^(1/ncol(SSN1))) 
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN1, MoreArgs = list(probs=c(a,b),na.rm=TRUE))   
  product<-prod(SP1[2,]-SP1[1,])
  vol<-c(prod=product,mean=mean(SP1[2,]-SP1[1,]),gmean=product^{1/length(SP1[2,]-SP1[1,])})   
  return(vol)                                                                                        
}
.volumeA_coordinates_full_dim1 <-
  function(S1,S2,steps=101 ){ 
    nt<-1
    names(S2)<-"v1"
    r<-.volumeA_full_onecol(data.frame(v1=S1), data.frame(v1=S2), steps=steps)
    integral_coord<-r$integral_approx
    plot_volume<-r$volume
    alpha_grid<-seq(0,1,length=steps)    
    erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_volume=plot_volume)
    return(erg)
}
.volumeA_coordinates_full <-
function(S1,S2,steps=101 ){ 
  nt<-ncol(S1)
  volumeA_function<-function(x1,x2,steps){
    r<-.volumeA_full_onecol(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
    return(r$integral_approx)
  }
  volumeA_function2<-function(x1,x2,steps){
    r<-.volumeA_full_onecol(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
    return(r$volume)
  }
  integral_coord<-mapply(volumeA_function, x1=S1, x2=S2, MoreArgs =list(steps=steps))
  plot_volume<-mapply(volumeA_function2, x1=S1, x2=S2, MoreArgs =list(steps=steps))
  alpha_grid<-seq(0,1,length=steps)    
  erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_volume=plot_volume)
  return(erg)
}
.volumeA_full_onecol <-
function(S1,S2,steps=101){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  alpha_grid<-seq(0,1,length=steps) 
  Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))      
  ab<-Spans$max-Spans$min          
  SSN<-(S-Spans$min)/(ab)          
  SSN[,ab==0]<-0.5                 
  z<-mapply(.volume_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),])) 
  z<-1/(1-alpha_grid[-length(alpha_grid)])*z[-length(alpha_grid)]
  z <- ifelse(z <= 1, z, 1)
  integral_approx<-c(prod=.trapz(alpha_grid[-length(alpha_grid)], z))
  erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx) 
  return(erg)
}
.volumeA2_full <-
  function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
    steps0<-steps
    alpha0_grid<-alpha_grid
    .volumeA2_full_severalcol(S1,S2,steps=steps0,alpha_grid=alpha0_grid  )
  }
.volumeA2_full_dim1 <-
  function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
    steps0<-steps
    alpha0_grid<-alpha_grid
    .volumeA2_full_severalcol_dim1(S1,S2,steps=steps0,alpha_grid=alpha0_grid  )
  }
.volumeA2_full_severalcol_dim1 <-
  function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
    S<-c(S1,S2)
    Ncol<-1
    ab<-max(S)-min(S)          
    SSN<-t((t(S)-min(S))/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN) 
    z<-mapply(.volume_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN))
    z<-rbind(1/(1-alpha_grid[-length(alpha_grid)])*z[-length(alpha_grid)])
    z <- ifelse(z <= 1, z, 1)
    integral_approx<-c(prod=.trapz(alpha_grid[-length(alpha_grid)], z),mean=.trapz(alpha_grid[-length(alpha_grid)],z),gmean=.trapz(alpha_grid[-length(alpha_grid)],z))
    plot_data_prod<-z
    erg<-list(alpha_grid=seq(0,1,length=steps),volume=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod) 
    return(erg)
  }
.volumeA2_full_severalcol <-
  function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    z<-mapply(.volume_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    z<-rbind(1/(1-alpha_grid[-length(alpha_grid)])*z[1,][-length(alpha_grid)],1/(1-alpha_grid[-length(alpha_grid)])*z[2,][-length(alpha_grid)],1/(1-alpha_grid[-length(alpha_grid)])*z[3,][-length(alpha_grid)])
    z <- ifelse(z <= 1, z, 1)
    integral_approx<-c(prod=.trapz(alpha_grid[-length(alpha_grid)], z[1,]),mean=.trapz(alpha_grid[-length(alpha_grid)],z[2,]),gmean=.trapz(alpha_grid[-length(alpha_grid)],z[3,]))
    plot_data_prod<-z[1,]
    erg<-list(alpha_grid=seq(0,1,length=steps),volume=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod) 
    return(erg)
  }
.trpca <- function(data, va){  # data = dataset as for dynRB, va = how much variance (0-1) should included axes explain
  PCA <- prcomp(data[,-1])
    vars <- apply(PCA$x, 2, var)
    prop <- cumsum(vars / sum(vars))
    k <- which(prop >= va)[1]
    k <- ifelse(k==1, 2, k)
    data1 <- data.frame(data[,1], PCA$x[,1:k])
    colnames(data1)[1] <- "Species"
  return(data1)
}
.getVarianceEstimate <-function(x,y){
  # X
  x_smallerM = x[1:round(length(x)/2 + 0.1)] #x[x < median(x)]
  x_greaterM = x[round(length(x)/2+1.1):length(x)]  #x[x >= median(x)]
  n = as.numeric(length(x))
  k = as.numeric(length(x_smallerM))
  x_smallerM_ranks = rank(x_smallerM)
  x_greaterM_ranks = rank(x_greaterM)
  # Y
  y_smallerM =  y[1:round(length(y)/2 +0.1)] # y[y < median(y)]
  y_greaterM =  y[round(length(y)/2+1.1):length(y)]# y[y >= median(y)]
  m = as.numeric(length(y))
  l = as.numeric(length(y_smallerM))
  y_smallerM_ranks = rank(y_smallerM)
  y_greaterM_ranks = rank(y_greaterM)
  xy_smallerM_ranks = rank(c(x_smallerM,y_smallerM))[1:k]
  yx_smallerM_ranks = rank(c(y_smallerM,x_smallerM))[1:l]
  xy_greaterM_ranks = rank(c(x_greaterM,y_greaterM))[1:(n-k)]
  yx_greaterM_ranks = rank(c(y_greaterM,x_greaterM))[1:(m-l)]
  s2_X1 = 1/(l^2*(k-1))       * sum((xy_smallerM_ranks - x_smallerM_ranks - mean(xy_smallerM_ranks) + (k+1)/2)^2)
  s2_X2 = 1/((m-l)^2*(n-k-1)) * sum((xy_greaterM_ranks - x_greaterM_ranks - mean(xy_greaterM_ranks) + (n-k+1)/2)^2)
  s2_Y1 = 1/(k^2*(l-1))       * sum((yx_smallerM_ranks - y_smallerM_ranks - mean(yx_smallerM_ranks) + (l+1)/2)^2)
  s2_Y2 = 1/((n-k)^2*(m-l-1)) * sum((yx_greaterM_ranks - y_greaterM_ranks - mean(yx_greaterM_ranks) + (m-l+1)/2)^2)
  s2_2 = (l+k)*(s2_X1/k + s2_Y1/l) + (n+m-l-k)*(s2_X2/(n-k) + s2_Y2/(m-l))
  return(s2_2)
}
.getRankedSingleOverlapIndex <- function(x,y){
  x_length = as.numeric(length(x))
  if(round(x_length/2)==x_length/2){
    k=x_length/2
  }else{
    k=(x_length+1)/2
  }
  y_length = as.numeric(length(y))
  x = sort(x)
  xy_ranks = rank(c(x,y))
  x_rank = xy_ranks[1:x_length]
  y_rank = xy_ranks[(x_length+1):(x_length+y_length)]
  rank_gtm = sum(x_rank[(k+1) : x_length])
  rank_stm = sum(x_rank[1:k])
  K= sum(rank(x)[1:k])
  c = 2/((x_length)*(y_length)) * ( - x_length*(x_length+1) + K*4)  
  overlap_ranked = 2*(rank_gtm - rank_stm)/((x_length)*(y_length)) + c/2
  return(overlap_ranked)
}
.geo.mean <- function (x){
  n<-length(x) 
  mittelwert<-prod(x)^(1/n) 
  return (mittelwert)
}