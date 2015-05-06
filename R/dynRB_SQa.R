dynRB_SQa <-
function(A=A,steps=201,graphic=FALSE){
  #aggregation has to be one of "product", "mean", "gmean"
  #tranSLate aggregation
  names(A)[1]<-"Species"
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
      
      V <- funclist3$volumeA2_full_graphics(S1,SH,steps=steps0)     
      G$niche_size_V1_prod[i]<-V$integral_approx[1]
      G$niche_size_V1_mean[i]<-V$integral_approx[2]
      G$niche_size_V1_gmean[i]<-V$integral_approx[3]
      
      V <- funclist3$volumeA2_full_graphics(S2,SH,steps=steps0)      
      G$niche_size_V2_prod[i]<-V$integral_approx[1]
      G$niche_size_V2_mean[i]<-V$integral_approx[2]
      G$niche_size_V2_gmean[i]<-V$integral_approx[3]
      
      EE<-funclist3$portionAinB2_full_graphics(S1,S2,steps=steps0)
      G$overlap_prod[i]<-EE$integral_approx[1]
      G$overlap_mean[i]<-EE$integral_approx[2]
      G$overlap_gmean[i]<-EE$integral_approx[3]
      print(i)
    }
  }else{
    for(i in 1:nrow(G)){
      S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
      S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
      SH<-A[,2:(dims+1)]
      
      V <- funclist3$volumeA2_full_graphics_off(S1,SH,steps=steps0)     
      G$niche_size_V1_prod[i]<-V$integral_approx[1]
      G$niche_size_V1_mean[i]<-V$integral_approx[2]
      G$niche_size_V1_gmean[i]<-V$integral_approx[3]
      
      V <- funclist3$volumeA2_full_graphics_off(S2,SH,steps=steps0)      
      G$niche_size_V2_prod[i]<-V$integral_approx[1]
      G$niche_size_V2_mean[i]<-V$integral_approx[2]
      G$niche_size_V2_gmean[i]<-V$integral_approx[3]
      
      EE<-funclist3$portionAinB2_full_graphics_off(S1,S2,steps=steps0)
      G$overlap_prod[i]<-EE$integral_approx[1]
      G$overlap_mean[i]<-EE$integral_approx[2]
      G$overlap_gmean[i]<-EE$integral_approx[3]
      print(i)
    }
  }
  print(G[,c(1,2,3,6,9)])
  return(G)
}
