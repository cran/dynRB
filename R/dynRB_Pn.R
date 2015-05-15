dynRB_Pn <-
function(A=A,steps=201,graphic=FALSE){
  #aggregation has to be one of "product", "mean", "gmean"
  #tranSLate aggregation
  steps0<-steps
  
  #expand grid  
  dims<-ncol(A)-1
  insects<-levels(factor(A$Species))
  G<-expand.grid(insects,insects)
  names(G)[1:2]<-c("V1","V2")
  for(i in 1:dims){
    G<-cbind(G,data.frame(x=0))
    names(G)[i+2]<-names(A)[i+1]
  }
  
  if(graphic){
    for(i in 1:nrow(G)){
      S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
      S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
      SH<-A[,2:(dims+1)]
      
      EE<-funclist3$portionAinB_coordinates_full_graphics(S1,S2,steps=steps0)
      G[i,3:(dims+2)]<-EE$integral_coord
      print(i)
    }
  }else{
    for(i in 1:nrow(G)){
      S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
      S2<-subset(A,A$Species==G$V2[i])[,2:(dims+1)]
      SH<-A[,2:(dims+1)]
      
      EE<-funclist3$portionAinB_coordinates_full_graphics_off(S1,S2,steps=steps0)
      G[i,3:(dims+2)]<-EE$integral_coord
      print(i)
    }
  }
  
  
  
  return(G)
}
