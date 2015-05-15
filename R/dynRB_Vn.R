dynRB_Vn <-
function(A=A,steps=201){
  #aggregation has to be one of "product", "mean", "gmean"
  #tranSLate aggregation
  steps0<-steps
  
  #expand grid  
  dims<-ncol(A)-1
  S2<-A[2:(dims+1)]
  insects<-as.character(levels(factor(A$Species)))
  G<-data.frame(V1=insects)
  for(i in 1:dims){
    G<-cbind(G,data.frame(x=0))
    names(G)[i+1]<-names(A)[i+1]
  }
  
  for(i in 1:nrow(G)){
    S1<-subset(A,A$Species==G$V1[i])[,2:(dims+1)]
    
    EE<-funclist3$volumeA_coordinates_full_graphics_off(S1,S2,steps=steps0)
    G[i,2:(dims+1)]<-EE$integral_coord
    print(i)
  }
  return(G)
}
