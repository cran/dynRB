
funclist4<-list(
  portionAinB_coordinates_full_adjusted_graphics_off=function(S1,S2,steps=101){
    nt<-ncol(S1)
    steps0<-steps
    integral_coord<-rep(0,length=nt)
    o_coord_adj<-rep(0,nt)
    for(k in 1:nt){
      T1<-subset(S1,select=k)
      T2<-subset(S2,select=k)
      EE<-funclist3$portionAinB_full_graphics_off(T1,T2,steps=steps0)
      #--
      SH<-rbind(T1,T2)
      SL<-cbind(seq(0,1,0.001))
      
      VlHI <- funclist3$volumeA_full_graphics_off(SL,SL,steps=steps0)
      
      V <- funclist3$volumeA_full_graphics_off(T2,SH,steps=steps0)      
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- ifelse(quotVl>quantile(quotVl, p=0.9, na.rm=TRUE) | quotVl<quantile(quotVl, p=0.1, na.rm=TRUE), median(quotVl, na.rm=TRUE), quotVl)
      NSL <-  mean(quotVl_cor, na.rm=T)
      #G$niche_size_V2[i]<-NSL
      
      EE<-funclist3$portionAinB_full_graphics_off(T1,T2,steps=steps0)
      RO <- EE$overlap * quotVl_cor[1:(steps0-1)]
      ROm <- mean(RO, na.rm=TRUE)
      o_coord_adj[k]<-ROm/NSL
    }
    
    erg<-list(integral_coord=o_coord_adj)
    return(erg)
  },
  portionAinB_coordinates_full_adjusted_graphics=function(S1,S2,steps=101){
    nt<-ncol(S1)
    steps0<-steps
    integral_coord<-rep(0,length=nt)
    o_coord_adj<-rep(0,nt)
    for(k in 1:nt){
      T1<-subset(S1,select=k)
      T2<-subset(S2,select=k)
      EE<-funclist3$portionAinB_full_graphics(T1,T2,steps=steps0)
      #--
      SH<-rbind(T1,T2)
      SL<-cbind(seq(0,1,0.001))
      
      VlHI <- funclist3$volumeA_full_graphics(SL,SL,steps=steps0)
      
      V <- funclist3$volumeA_full_graphics(T2,SH,steps=steps0)      
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- ifelse(quotVl>quantile(quotVl, p=0.9, na.rm=TRUE) | quotVl<quantile(quotVl, p=0.1, na.rm=TRUE), median(quotVl, na.rm=TRUE), quotVl)
      NSL <-  mean(quotVl_cor, na.rm=T)
      #G$niche_size_V2[i]<-NSL
      
      EE<-funclist3$portionAinB_full_graphics(T1,T2,steps=steps0)
      RO <- EE$overlap * quotVl_cor[1:(steps0-1)]
      ROm <- mean(RO, na.rm=TRUE)
      o_coord_adj[k]<-ROm/NSL
    }
    
    erg<-list(integral_coord=o_coord_adj)
    return(erg)
  },
  volumeA_coordinates_full_adjusted_graphics_off=function(S1,S2,steps=101){
    nt<-ncol(S1)
    steps0<-steps
    integral_coord<-rep(0,length=nt)
    s_coord_adj<-rep(0,nt)
    for(k in 1:nt){
      T1<-subset(S1,select=k)
      T2<-subset(S2,select=k)
      #--
      SH<-rbind(T1,T2)
      SL<-cbind(seq(0,1,0.001))
      
      VlHI <- funclist3$volumeA_full_graphics_off(SL,SL,steps=steps0)
      
      V <- funclist3$volumeA_full_graphics_off(T1,SH,steps=steps0)      
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- ifelse(quotVl>quantile(quotVl, p=0.9, na.rm=TRUE) | quotVl<quantile(quotVl, p=0.1, na.rm=TRUE), median(quotVl, na.rm=TRUE), quotVl)
      NSL <-  mean(quotVl_cor, na.rm=T)
      
      s_coord_adj[k]<-NSL
    }
    
    erg<-list(integral_coord=s_coord_adj)
    return(erg)
  },
  volumeA_coordinates_full_adjusted_graphics=function(S1,S2,steps=101){
    nt<-ncol(S1)
    steps0<-steps
    integral_coord<-rep(0,length=nt)
    s_coord_adj<-rep(0,nt)
    for(k in 1:nt){
      T1<-subset(S1,select=k)
      T2<-subset(S2,select=k)
      #--
      SH<-rbind(T1,T2)
      SL<-cbind(seq(0,1,0.001))
      
      VlHI <- funclist3$volumeA_full_graphics(SL,SL,steps=steps0)
      
      V <- funclist3$volumeA_full_graphics(T1,SH,steps=steps0)      
      difl <- VlHI$volume-V$volume
      
      quotVl <- ifelse(difl>=0, V$volume/VlHI$volume, 1)
      quotVl_cor <- ifelse(quotVl>quantile(quotVl, p=0.9, na.rm=TRUE) | quotVl<quantile(quotVl, p=0.1, na.rm=TRUE), median(quotVl, na.rm=TRUE), quotVl)
      NSL <-  mean(quotVl_cor, na.rm=T)
      
      s_coord_adj[k]<-NSL
    }
    
    erg<-list(integral_coord=s_coord_adj)
    return(erg)
  }
)
