dynRB_VPac <-
function(A=A,steps=201,graphic=FALSE){
  #aggregation has to be one of "product", "mean", "gmean"
  #tranSLate aggregation
  names(A)[1]<-"Species"
  #agg0<-aggregation
  steps0<-steps
  
  #expand grid  
  dims<-ncol(A)-1
  insects<-levels(factor(A$Species))
  G<-expand.grid(insects,insects)
  names(G)[1:2]<-c("V1","V2")
  G$overlap_prod<-0
  G$overlap_mean<-0
  G$overlap_gmean<-0
  G$niche_size_V1_prod<-0
  G$niche_size_V1_mean<-0
  G$niche_size_V1_gmean<-0
  G$niche_size_V2_prod<-0
  G$niche_size_V2_mean<-0
  G$niche_size_V2_gmean<-0
  
  if(graphic){
    for(i in 1:nrow(G)){
      S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
      S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
      SH<-A[,2:(dims+1)]
      SH<-as.data.frame(SH)
      
      SL<-cbind(seq(0,1,0.001))
      for(j in 1:(dims-1)){
        SL<-cbind(SL,seq(0,1,0.001))
      }
      SL<-as.data.frame(SL)
      VlHI <- funclist3$volumeA2_full_graphics(SL,SL,steps=steps0)
      V <- funclist3$volumeA2_full_graphics(S1,SH,steps=steps0)      ## das 2. S kann auch das Hyperinsect sein # das 1. S ist das 2. S aus portionAinB_full
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- apply(quotVl, 1, funclist1$func_quotVl_cor)
      NSL <-  mapply(mean, data.frame(quotVl_cor),MoreArgs = list(na.rm=TRUE))
      G$niche_size_V1_prod[i]<-NSL[1]
      G$niche_size_V1_mean[i]<-NSL[2]
      G$niche_size_V1_gmean[i]<-NSL[3]
      
      V <- funclist3$volumeA2_full_graphics(S2,SH,steps=steps0)      ## das 2. S kann auch das Hyperinsect sein # das 1. S ist das 2. S aus portionAinB_full
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- apply(quotVl, 1, funclist1$func_quotVl_cor)
      NSL <-  mapply(mean, data.frame(quotVl_cor),MoreArgs = list(na.rm=TRUE))
      G$niche_size_V2_prod[i]<-NSL[1]
      G$niche_size_V2_mean[i]<-NSL[2]
      G$niche_size_V2_gmean[i]<-NSL[3]
      
      EE<-funclist3$portionAinB2_full_graphics(S1,S2,steps=steps0)
      RO <- t(EE$overlap) * quotVl_cor[1:(steps0-1),]             ## hier Fehlermeldung!
      #ROm <- c(mean(RO[1,], na.rm=TRUE),mean(RO[2,], na.rm=TRUE),mean(RO[3,], na.rm=TRUE)) 
      ROm <- mapply(mean, data.frame(RO),MoreArgs = list(na.rm=TRUE))
      r<-ROm/NSL
      G$overlap_prod[i]<-r[1]
      G$overlap_mean[i]<-r[2]
      G$overlap_gmean[i]<-r[3]
      print(i)
    }
  }else{
    for(i in 1:nrow(G)){
      S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
      S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
      SH<-A[,2:(dims+1)]
      SH<-as.data.frame(SH)
      
      SL<-cbind(seq(0,1,0.001))
      for(j in 1:(dims-1)){
        SL<-cbind(SL,seq(0,1,0.001))
      }
      SL<-as.data.frame(SL)
      VlHI <- funclist3$volumeA2_full_graphics_off(SL,SL,steps=steps0)
      V <- funclist3$volumeA2_full_graphics_off(S1,SH,steps=steps0)      ## das 2. S kann auch das Hyperinsect sein # das 1. S ist das 2. S aus portionAinB_full
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- apply(quotVl, 1, funclist1$func_quotVl_cor)
      NSL <-  mapply(mean, data.frame(quotVl_cor),MoreArgs = list(na.rm=TRUE))
      G$niche_size_V1_prod[i]<-NSL[1]
      G$niche_size_V1_mean[i]<-NSL[2]
      G$niche_size_V1_gmean[i]<-NSL[3]
      
      V <- funclist3$volumeA2_full_graphics_off(S2,SH,steps=steps0)      ## das 2. S kann auch das Hyperinsect sein # das 1. S ist das 2. S aus portionAinB_full
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- apply(quotVl, 1, funclist1$func_quotVl_cor)
      NSL <-  mapply(mean, data.frame(quotVl_cor),MoreArgs = list(na.rm=TRUE))
      G$niche_size_V2_prod[i]<-NSL[1]
      G$niche_size_V2_mean[i]<-NSL[2]
      G$niche_size_V2_gmean[i]<-NSL[3]
      
      EE<-funclist3$portionAinB2_full_graphics_off(S1,S2,steps=steps0)
      RO <- t(EE$overlap) * quotVl_cor[1:(steps0-1),]             ## hier Fehlermeldung!
      #ROm <- c(mean(RO[1,], na.rm=TRUE),mean(RO[2,], na.rm=TRUE),mean(RO[3,], na.rm=TRUE)) 
      ROm <- mapply(mean, data.frame(RO),MoreArgs = list(na.rm=TRUE))
      r<-ROm/NSL
      G$overlap_prod[i]<-r[1]
      G$overlap_mean[i]<-r[2]
      G$overlap_gmean[i]<-r[3]
      print(i)
    }
  }
  print(G[,c(1,2,3,6,9)])
  return(G)
}
