funclist2<-list(portionAinB2_full_severalcol=function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)[1:(steps-1)]){
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
  z<-mapply(funclist1$func2_severalcolumns, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,helping_func1=funclist1$helping_func1)) 
  integral_approx<-c(prod=trapz(seq(0,1,length=steps), z[1,]),mean=trapz(seq(0,1,length=steps),z[2,]),gmean=trapz(seq(0,1,length=steps),z[3,])) 
  plot_data_prod<-z[1,]
  erg<-list(alpha_grid=seq(0,1,length=steps)[1:(steps-1)],overlap=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod)                                                                                             
  return(erg)
},
volumeA2_full_severalcol=function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
  ab<-Spans$max-Spans$min          
  SSN<-t((t(S)-Spans$min)/(ab))    
  SSN[,ab==0]<-0.5                 
  SSN<-as.data.frame(SSN)          
  z<-mapply(funclist1$helping_func3_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
  z<-rbind(1/(1-alpha_grid[-length(alpha_grid)])*z[1,][-length(alpha_grid)],1/(1-alpha_grid[-length(alpha_grid)])*z[2,][-length(alpha_grid)],1/(1-alpha_grid[-length(alpha_grid)])*z[3,][-length(alpha_grid)])
  z <- ifelse(z <= 1, z, 1)
  integral_approx<-c(prod=trapz(alpha_grid[-length(alpha_grid)], z[1,]),mean=trapz(alpha_grid[-length(alpha_grid)],z[2,]),gmean=trapz(alpha_grid[-length(alpha_grid)],z[3,]))
  plot_data_prod<-z[1,]
  erg<-list(alpha_grid=seq(0,1,length=steps),volume=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod) 
  return(erg)
},
volumeA_full_onecol=function(S1,S2,steps=101){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  alpha_grid<-seq(0,1,length=steps) 
  Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))      
  ab<-Spans$max-Spans$min          
  SSN<-(S-Spans$min)/(ab)          
  SSN[,ab==0]<-0.5                 
  z<-mapply(funclist1$fun_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),])) 
  z<-1/(1-alpha_grid[-length(alpha_grid)])*z[-length(alpha_grid)]
  z <- ifelse(z <= 1, z, 1)
  integral_approx<-c(prod=trapz(alpha_grid[-length(alpha_grid)], z))
  erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx) 
  return(erg)
},
portionAinB_full_onecolumn=function (S1, S2, steps = 101) {
  S <- rbind(S1, S2)
  Ncol <- ncol(S1)
  alpha_grid <- seq(0, 1, length = steps)
  Spans <- data.frame(trait_nr = 1:Ncol, min = min(S, na.rm = TRUE), max = max(S, na.rm = TRUE))
  ab <- Spans$max - Spans$min
  SSN <- (S - Spans$min)/(ab)
  SSN[, ab == 0] <- 0.5
  nrow_S1 <- nrow(S1)
  nrow_S2 <- nrow(S2)
  z <- mapply(funclist1$func_onecolumn, alpha = alpha_grid,  MoreArgs = list(SSN1 = SSN[1:nrow_S1, ], SSN2 = SSN[(nrow_S1 + 1):(nrow_S1 + nrow_S2), ], nrow_S1 = nrow_S1, nrow_S2 = nrow_S2, Ncol = Ncol, helping_func2 = funclist1$helping_func2))
  integral_approx <- c(prod = trapz(alpha_grid, z))
  erg <- list(alpha_grid = alpha_grid, overlap = z, integral_approx = integral_approx)
  return(erg)
}
)