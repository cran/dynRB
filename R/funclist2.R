
funclist2<-list(
  volumeA2=function(S1,S2,alpha=0){
    alpha0<-alpha
    if(ncol(S1)>1){
      funclist1$volumeA2_severalcol(S1,S2,alpha=alpha0)
    }else{
      funclist1$volumeA2_onecol(S1,S2,alpha=alpha0)
    }
  },
  portionAinB2_full_severalcol_graphics_off=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func2_severalcolumns_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc1_NA=funclist1$helpingfunc1_NA)) # (6) bzw. (3)+(4)+(5)
    }else{
      z<-mapply(funclist1$func2_severalcolumns, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc1=funclist1$helpingfunc1)) # (6) bzw. (3)+(4)+(5)
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                             
    return(erg)
  },
  portionAinB2_full_severalcol_graphics=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func2_severalcolumns_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc1_NA=funclist1$helpingfunc1_NA)) # (6) bzw. (3)+(4)+(5)
    }else{
      z<-mapply(funclist1$func2_severalcolumns, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc1=funclist1$helpingfunc1)) # (6) bzw. (3)+(4)+(5)
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                             
    plot(erg$alpha_grid,erg$overlap[1,],type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx[1],4)),xlab="alpha","ylab"="overlap prod")    
    return(erg)
  },
  portionAinB2_full_onecolumn_graphics_off=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))    
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)         
    SSN[,ab==0]<-0.5                
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func2_onecolumn_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc2_NA=funclist1$helpingfunc2_NA)) 
    }else{
      z<-mapply(funclist1$func2_onecolumn, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc2=funclist1$helpingfunc2)) 
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                               
    return(erg)
  },
  portionAinB2_full_onecolumn_graphics=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))     
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)          
    SSN[,ab==0]<-0.5                 
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func2_onecolumn_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc2_NA=funclist1$helpingfunc2_NA)) 
    }else{
      z<-mapply(funclist1$func2_onecolumn, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helpingfunc2=funclist1$helpingfunc2)) 
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                               
    plot(erg$alpha_grid,erg$overlap[1,],type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx[1],4)),xlab="alpha","ylab"="overlap_prod")                                
    return(erg)
  },
  volumeA2_full_severalcol_graphics=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)  
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$helpingfunc3_severalcol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    }else{
      z<-mapply(funclist1$helpingfunc3_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)  
    plot(erg$alpha_grid,erg$volume[1,],type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx[1],4)),xlab="alpha","ylab"="volume prod")
    return(erg)
  },
  volumeA2_full_severalcol_graphics_off=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)   
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$helpingfunc3_severalcol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    }else{
      z<-mapply(funclist1$helpingfunc3_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)  
    return(erg)
  },
  volumeA2_full_onecol_graphics=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps) 
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))      
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)          
    SSN[,ab==0]<-0.5                 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$helpingfunc3_onecol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),])) 
    }else{
      z<-mapply(funclist1$helpingfunc3_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)  
    plot(erg$alpha_grid,erg$volume[1,],type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx[1],4)),xlab="alpha","ylab"="volume prod")
    return(erg)
  },
  volumeA2_full_onecol_graphics_off=function(S1,S2,steps=101){
    S<-rbind(S1,S2)
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps) 
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))     
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)         
    SSN[,ab==0]<-0.5                 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$helpingfunc3_onecol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),])) 
    }else{
      z<-mapply(funclist1$helpingfunc3_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
    }
    integral_approx<-c(prod=mean(z[1,]),mean=mean(z[2,]) ,gmean=mean(z[3,]) )                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)
    return(erg)
  },
  portionAinB_full_severalcol_graphics_off=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    #calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func_severalcolumns_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g_NA=funclist1$g_NA)) # (6) bzw. (3)+(4)+(5)
    }else{
      z<-mapply(funclist1$func_severalcolumns, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g=funclist1$g)) # (6) bzw. (3)+(4)+(5)
    }
    integral_approx<-mean(z)                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                                 
    return(erg)
  },
  portionAinB_full_severalcol_graphics=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func_severalcolumns_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g_NA=funclist1$g_NA)) # (6) bzw. (3)+(4)+(5)
    }else{
      z<-mapply(funclist1$func_severalcolumns, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g=funclist1$g)) # (6) bzw. (3)+(4)+(5)
    }
    integral_approx<-mean(z)                                                                                                                                                   # (6)
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                                 # (6)
    plot(erg$alpha_grid,erg$overlap,type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx,4)),xlab="alpha","ylab"="overlap")    # (6)
    return(erg)
  },
  portionAinB_full_onecolumn_graphics_off=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))    
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)         
    SSN[,ab==0]<-0.5                
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func_onecolumn_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g1_NA=funclist1$g1_NA)) 
    }else{
      z<-mapply(funclist1$func_onecolumn, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g1=funclist1$g1)) 
    }
    integral_approx<-mean(z)                                                                                                                                                  
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                               
    return(erg)
  },
  portionAinB_full_onecolumn_graphics=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))     
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)          
    SSN[,ab==0]<-0.5                 
    # calculate quantile intervals depending on alpha and calculate intersection 
    nrow_S1<-nrow(S1)
    nrow_S2<-nrow(S2)
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$func_onecolumn_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g1_NA=funclist1$g1_NA)) 
    }else{
      z<-mapply(funclist1$func_onecolumn, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,agg=agg,g1=funclist1$g1)) 
    }
    integral_approx<-mean(z)                                                                                                                                                   
    erg<-list(alpha_grid=alpha_grid,overlap=z,integral_approx=integral_approx)                                                                                                
    plot(erg$alpha_grid,erg$overlap,type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx,4)),xlab="alpha","ylab"="overlap")                                
    return(erg)
  },
  volumeA=function(S1,S2,alpha=0,aggregation="product"){
    alpha0<-alpha
    agg<-aggregation
    if(ncol(S1)>1){
      funclist1$volumeA_severalcol(S1,S2,alpha=alpha0,aggregation=agg)
    }else{
      funclist1$volumeA_onecol(S1,S2,alpha=alpha0,aggregation=agg)
    }
  },
  volumeA_full_severalcol_graphics=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)  
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$fun_severalcol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg))
    }else{
      z<-mapply(funclist1$fun_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg))
    }
    integral_approx<-mean(z)                                                                                                                                                  
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)
    plot(erg$alpha_grid,erg$volume,type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx,4)),xlab="alpha","ylab"="volume")
    return(erg)
  },
  volumeA_full_severalcol_graphics_off=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps)   
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
    ab<-Spans$max-Spans$min          
    SSN<-t((t(S)-Spans$min)/(ab))    
    SSN[,ab==0]<-0.5                 
    SSN<-as.data.frame(SSN)          
    # calculate quantile intervals depending on alpha and calculate intersection 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$fun_severalcol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg))
    }else{
      z<-mapply(funclist1$fun_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg))
    }
    integral_approx<-mean(z)                                                                                                                                                   
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)
    return(erg)
  },
  volumeA_full_onecol_graphics=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps) 
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))      
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)          
    SSN[,ab==0]<-0.5                 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$fun_onecol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg)) 
    }else{
      z<-mapply(funclist1$fun_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg))
    }
    integral_approx<-mean(z)                                                                                                                                                 
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)
    plot(erg$alpha_grid,erg$volume,type="l",ylim=c(0,1.1),main=paste("Integral=", round(erg$integral_approx,4)),xlab="alpha","ylab"="volume")
    return(erg)
  },
  volumeA_full_onecol_graphics_off=function(S1,S2,steps=101,aggregation="product"){
    S<-rbind(S1,S2)
    agg<-function(x){prod(x)}
    if(aggregation=="mean"){agg<-function(x){mean(x)}}
    if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
    Ncol<-ncol(S1)
    alpha_grid<-seq(0,1,length=steps) 
    # normalize data (same transformation for ALL insects)
    Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))     
    ab<-Spans$max-Spans$min          
    SSN<-(S-Spans$min)/(ab)         
    SSN[,ab==0]<-0.5                 
    if(length(S1[is.na(S1)])+length(S2[is.na(S2)])>0){
      z<-mapply(funclist1$fun_onecol_NA, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg)) 
    }else{
      z<-mapply(funclist1$fun_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),],agg=agg))
    }
    integral_approx<-mean(z)                                                                                                                                                  
    erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx)
    return(erg)
  }
)

