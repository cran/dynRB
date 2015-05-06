
funclist1<-list(helpingfunc1_NA=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b),na.rm=TRUE),  quantile(vec2, probs=c(a,b),na.rm=TRUE) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                        
  return(r)                                                            
},
helpingfunc1=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b)),  quantile(vec2, probs=c(a,b)) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))  
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                          
  return(r)                                                             
},
helpingfunc2_NA=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b),na.rm = TRUE),  quantile(vec2, probs=c(a,b),na.rm = TRUE) )   
  x<-c(x,x[4]-x[3])  
  x[5]<-ifelse(x[5]<=0,0,x[5]) 
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                     
  return(r)                                                         
},
helpingfunc2=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b)),  quantile(vec2, probs=c(a,b)) )   
  x<-c(x,x[4]-x[3])                                       
  x[5]<-ifelse(x[5]<=0,0,x[5]) 
  r<-c(max(x[c(1,3)]),min(x[c(2,4)])) 
  r<-c(x[3:5],r,r[2]-r[1])  
  r[6]<-ifelse(r[6]<=0,0,r[6])                                          
  return(r)                                                             
},
volumeA2_severalcol=function(S1,S2,alpha=0){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  # step 1: normalize data (same transformation for ALL insects)
  Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
  ab<-Spans$max-Spans$min         
  SSN<-t((t(S)-Spans$min)/(ab))    
  SSN[,ab==0]<-0.5                 
  SSN<-as.data.frame(SSN)          
  # step 2: calculate quantile intervals depending on alpha 
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN[1:nrow(S1),], MoreArgs = list(probs=c(a,b),na.rm=TRUE)) 
  #step 3: calculate volume
  product<-prod(SP1[2,]-SP1[1,])
  vol<-c(prod=product,mean=mean(SP1[2,]-SP1[1,]),gmean=product^{1/length(SP1[2,]-SP1[1,])}) 
  erg<-list(volume=vol)
  return(erg)
},
volumeA2_onecol=function(S1,S2,alpha=0){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  #step 1: normalize data (same transformation for ALL insects)
  Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))    
  ab<-Spans$max-Spans$min          
  SSN<-(S-Spans$min)/(ab)          
  SSN[,ab==0]<-0.5                 
  #step 2: calculate quantile intervals depending on alpha 
  SP1<-quantile(SSN[1:nrow(S1),],probs=c(alpha/2,1-alpha/2),na.rm = TRUE)
  #step 3: calculate volume
  product<-prod(SP1[2]-SP1[1])
  vol<-c(prod=product,mean=mean(SP1[2]-SP1[1]),gmean=product^{1/length(SP1[2]-SP1[1])}) 
  erg<-list(volume=vol)
  return(erg)
},
helpingfunc3_severalcol_NA=function(alpha,SSN1){
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN1, MoreArgs = list(probs=c(a,b),na.rm=TRUE))   # calculate quantile intervals depending on alpha 
  product<-prod(SP1[2,]-SP1[1,])
  vol<-c(prod=product,mean=mean(SP1[2,]-SP1[1,]),gmean=product^{1/length(SP1[2,]-SP1[1,])})   # calculate volume
  return(vol)                                                                                        
},
helpingfunc3_severalcol=function(alpha,SSN1){
  a<-alpha/2
  b<-1-alpha/2
  SP1<-vapply(SSN1, quantile, probs=c(a,b), name=FALSE,FUN.VALUE=double(2))   # calculate quantile intervals depending on alpha 
  product<-prod(SP1[2,]-SP1[1,])
  vol<-c(prod=product,mean=mean(SP1[2,]-SP1[1,]),gmean=product^{1/length(SP1[2,]-SP1[1,])})   # calculate volume
  return(vol)                                                                                        
},
helpingfunc3_onecol_NA=function(alpha,SSN1){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2),na.rm = TRUE)  # calculate quantile intervals depending on alpha
  product<-prod(SP1[2]-SP1[1])
  vol<-c(prod=product,mean=mean(SP1[2]-SP1[1]),gmean=product^{1/length(SP1[2]-SP1[1])})   # calculate volume
  return(vol)                                                                                      
},
helpingfunc3_onecol=function(alpha,SSN1){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2))  # calculate quantile intervals depending on alpha
  product<-prod(SP1[2]-SP1[1])
  vol<-c(prod=product,mean=mean(SP1[2]-SP1[1]),gmean=product^{1/length(SP1[2]-SP1[1])})   # calculate volume
  return(vol)                                                                                      
},
helpingfunc4_NA=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b),na.rm=TRUE),  quantile(vec2, probs=c(a,b),na.rm=TRUE) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                        
  return(r)                                                            
},
helpingfunc4=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b)),  quantile(vec2, probs=c(a,b)) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))  
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                          
  return(r)                                                             
},
g1_NA=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b),na.rm = TRUE),  quantile(vec2, probs=c(a,b),na.rm = TRUE) )   
  x<-c(x,x[4]-x[3])  
  x[5]<-ifelse(x[5]<=0,0,x[5]) 
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                     
  return(r)                                                         
},
g1=function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs=c(a,b)),  quantile(vec2, probs=c(a,b)) )   
  x<-c(x,x[4]-x[3])                                       
  x[5]<-ifelse(x[5]<=0,0,x[5]) 
  r<-c(max(x[c(1,3)]),min(x[c(2,4)])) 
  r<-c(x[3:5],r,r[2]-r[1])  
  r[6]<-ifelse(r[6]<=0,0,r[6])                                          
  return(r)                                                             
},
volumeA_severalcol=function(S1,S2,alpha=0,aggregation="product"){
  S<-rbind(S1,S2)
  agg<-function(x){prod(x)}
  if(aggregation=="mean"){agg<-function(x){mean(x)}}
  if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
  Ncol<-ncol(S1)
  # step 1: normalize data (same transformation for ALL insects)
  Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
  ab<-Spans$max-Spans$min         
  SSN<-t((t(S)-Spans$min)/(ab))    
  SSN[,ab==0]<-0.5                 
  SSN<-as.data.frame(SSN)          
  # step 2: calculate quantile intervals depending on alpha 
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN[1:nrow(S1),], MoreArgs = list(probs=c(a,b),na.rm=TRUE)) 
  #step 3: calculate volume
  vol<-agg(SP1[2,]-SP1[1,])
  erg<-list(volume=vol)
  return(erg)
},
volumeA_onecol=function(S1,S2,alpha=0,aggregation="product"){
  S<-rbind(S1,S2)
  agg<-function(x){prod(x)}
  if(aggregation=="mean"){agg<-function(x){mean(x)}}
  if(aggregation=="gmean"){agg<-function(x){prod(x)^{1/length(x)}}}
  Ncol<-ncol(S1)
  #step 1: normalize data (same transformation for ALL insects)
  Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))    
  ab<-Spans$max-Spans$min          
  SSN<-(S-Spans$min)/(ab)          
  SSN[,ab==0]<-0.5                 
  #step 2: calculate quantile intervals depending on alpha 
  SP1<-quantile(SSN[1:nrow(S1),],probs=c(alpha/2,1-alpha/2),na.rm = TRUE)
  #step 3: calculate volume
  vol<-agg(SP1[2]-SP1[1])
  erg<-list(volume=vol)
  return(erg)
},
fun_severalcol_NA=function(alpha,SSN1,agg){
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN1, MoreArgs = list(probs=c(a,b),na.rm=TRUE))   # calculate quantile intervals depending on alpha 
  vol<-agg(SP1[2,]-SP1[1,])   # calculate volume
  return(vol)                                                                                        
},
fun_severalcol=function(alpha,SSN1,agg){
  a<-alpha/2
  b<-1-alpha/2
  SP1<-vapply(SSN1, quantile, probs=c(a,b), name=FALSE,FUN.VALUE=double(2))   # calculate quantile intervals depending on alpha 
  vol<-agg(SP1[2,]-SP1[1,])  # calculate volume
  return(vol)                                                                                        
},
fun_onecol_NA=function(alpha,SSN1,agg){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2),na.rm = TRUE)  # calculate quantile intervals depending on alpha
  vol<-agg(SP1[2]-SP1[1])  # calculate volume
  return(vol)                                                                                      
},
fun_onecol=function(alpha,SSN1,agg){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2))  # calculate quantile intervals depending on alpha
  vol<-agg(SP1[2]-SP1[1]) # calculate volume
  return(vol)                                                                                      
},
func2_severalcolumns=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,helpingfunc1){
  #calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-mapply(helpingfunc1, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=alpha/2,b=1-alpha/2,n=nrow_S1))  
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
func2_severalcolumns_NA=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,helpingfunc1_NA){
  #calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-mapply(helpingfunc1_NA, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=alpha/2,b=1-alpha/2,n=nrow_S1))  
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
func2_onecolumn_NA=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,helpingfunc2_NA){
  # calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-helpingfunc2_NA(vec1 = SSN1, vec2 = SSN2,a=alpha/2,b=1-alpha/2,n=nrow_S1)
  product<-prod(SP_inters[6])
  vol_intersect<-c(product,mean(SP_inters[6]),product^{1/length(SP_inters[6])})  
  product<-prod(SP_inters[3])
  vol_B<-c(product,mean(SP_inters[3]),product^{1/length(SP_inters[3])})  
  r<-c(1,1,1)                                                                                             
  if(min(SP_inters[4]==SP_inters[1])==0 | min(SP_inters[5]==SP_inters[2])==0){                  
    r[vol_intersect<vol_B]<-vol_intersect[vol_intersect<vol_B]/vol_B[vol_intersect<vol_B]
    r[vol_B==0]<-0                                                                             
  }
  return(r)                                                                                        
},
func2_onecolumn=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,helpingfunc2){
  # calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-helpingfunc2(vec1 = SSN1, vec2 = SSN2,a=alpha/2,b=1-alpha/2,n=nrow_S1)
  product<-prod(SP_inters[6])
  vol_intersect<-c(product,mean(SP_inters[6]),product^{1/length(SP_inters[6])})  
  product<-prod(SP_inters[3])
  vol_B<-c(product,mean(SP_inters[3]),product^{1/length(SP_inters[3])})  
  r<-c(1,1,1)                                                                                             
  if(min(SP_inters[4]==SP_inters[1])==0 | min(SP_inters[5]==SP_inters[2])==0){                  
    r[vol_intersect<vol_B]<-vol_intersect[vol_intersect<vol_B]/vol_B[vol_intersect<vol_B]
    r[vol_B==0]<-0                                                                             
  }
  return(r)                                                                                        
},
func_severalcolumns_NA=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,agg,g_NA){
  #calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-mapply(g_NA, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=alpha/2,b=1-alpha/2,n=nrow_S1))  
  vol_intersect<-agg(SP_inters[6,])                                                                 
  vol_B<-agg(SP_inters[3,])                                                                                     
  r<-1                                                                                              
  if(min(SP_inters[4,]==SP_inters[1,])==0 | min(SP_inters[5,]==SP_inters[2,])==0){                  
    if(vol_intersect<vol_B){r<-vol_intersect/vol_B}                                                 
    if(vol_B==0){r<-0}                                                                              
  }
  return(r)                                                                                        
},
func_severalcolumns=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,agg,g){
  #calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-mapply(g, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=alpha/2,b=1-alpha/2,n=nrow_S1))  
  vol_intersect<-agg(SP_inters[6,])                                                                 
  vol_B<-agg(SP_inters[3,])                                                                                       
  r<-1                                                                                              
  if(min(SP_inters[4,]==SP_inters[1,])==0 | min(SP_inters[5,]==SP_inters[2,])==0){                  
    if(vol_intersect<vol_B){r<-vol_intersect/vol_B}                                                 
    if(vol_B==0){r<-0}                                                                            
  }
  return(r)                                                                                        
},
func_onecolumn_NA=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,agg,g1_NA){
  # calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-g1_NA(vec1 = SSN1, vec2 = SSN2,a=alpha/2,b=1-alpha/2,n=nrow_S1)
  vol_intersect<-agg(SP_inters[6])                                                                
  vol_B<-agg(SP_inters[3])                                                                                       
  r<-1                                                                                            
  if(min(SP_inters[4]==SP_inters[1])==0 | min(SP_inters[5]==SP_inters[2])==0){                      
    if(vol_intersect<vol_B){r<-vol_intersect/vol_B}                                                
    if(vol_B==0){r<-0}                                                                             
  }
  return(r)                                                                                        
},
func_onecolumn=function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol,agg,g1){
  # calculate quantile intervals depending on alpha and calculate intersection 
  SP_inters<-g1(vec1 = SSN1, vec2 = SSN2,a=alpha/2,b=1-alpha/2,n=nrow_S1)
  vol_intersect<-agg(SP_inters[6])                                                                
  vol_B<-agg(SP_inters[3])                                                                                      
  r<-1                                                                                            
  if(min(SP_inters[4]==SP_inters[1])==0 | min(SP_inters[5]==SP_inters[2])==0){                      
    if(vol_intersect<vol_B){r<-vol_intersect/vol_B}                                                 
    if(vol_B==0){r<-0}                                                                              
  }
  return(r)                                                                                        
},
func_quotVl_cor=function(x){
  r<-ifelse(x>quantile(x, p=0.9, na.rm=TRUE) | x<quantile(x, p=0.1, na.rm=TRUE), median(x, na.rm=TRUE), x)
  return(r)
}
)
