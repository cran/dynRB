funclist1<-list(helping_func1=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b),na.rm=TRUE),  quantile(vec2, probs=c(a,b),na.rm=TRUE) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                        
  return(r)                                                            
},
helping_func2=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b),na.rm = TRUE),  quantile(vec2, probs=c(a,b),na.rm = TRUE) )   
  x<-c(x,x[4]-x[3])  
  x[5]<-ifelse(x[5]<=0,0,x[5]) 
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                     
  return(r)                                                         
},
helping_func3_severalcol=function(alpha,SSN1){
  alpha<-(1-(1-alpha)^(1/ncol(SSN1))) 
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN1, MoreArgs = list(probs=c(a,b),na.rm=TRUE))   
  product<-prod(SP1[2,]-SP1[1,])
  vol<-c(prod=product,mean=mean(SP1[2,]-SP1[1,]),gmean=product^{1/length(SP1[2,]-SP1[1,])})   
  return(vol)                                                                                        
},
fun_onecol=function(alpha,SSN1){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2),na.rm = TRUE)  
  vol<-(SP1[2]-SP1[1])  
  return(vol)                                                                                      
},
func2_severalcolumns=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,helping_func1){
  SP_inters<-mapply(helping_func1, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=((1-(1-alpha)^(1/ncol(SSN1)))/2),  b=(1-(1-(1-alpha)^(1/ncol(SSN1)))/2),n=nrow_S1))  
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
},
func_onecolumn=function (alpha, SSN1, SSN2, nrow_S1, nrow_S2, Ncol, helping_func2) {
  SP_inters <- helping_func2(vec1 = SSN1, vec2 = SSN2, a = alpha/2, b = 1 - alpha/2, n = nrow_S1)
  vol_intersect <- (SP_inters[6])
  vol_B <- (SP_inters[3])
  r <- 1
  if (min(SP_inters[4] == SP_inters[1]) == 0 | min(SP_inters[5] == SP_inters[2]) == 0) {
    if (vol_intersect < vol_B) {r <- vol_intersect/vol_B}
    if (vol_B == 0) {r <- 0}
  }
  return(r)
}
)