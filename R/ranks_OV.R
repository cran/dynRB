ranks_OV <- function(A=A, alpha=0.05, reps4boot = 1000, digit=3){
  insects<-levels(factor(A[,1]))
  if(length(insects) != 2) stop ("Need two groups")
  A[,1][A[,1]==insects[1]] <- 0
  A[,1][A[,1]==insects[2]] <- 1
  Species_1 = filter(A, A[,1]==0)
  Species_2 = filter(A, A[,1]==1)
  
  n = nrow(Species_1)
  m = nrow(Species_2)
  if(n <= 2) stop ("Need more observation of Data Set 1")
  if(m <= 2) stop ("Need more observation of Data Set 2")
  if(ncol(Species_1)!=ncol(Species_2)) stop("Data Sets must have same dimension")
  
  z = qnorm(1- alpha/2)
  dim = ncol(Species_1)-1
  if(dim == 0) stop ("No data contained in the Data Sets")
  
  for (i in 1:dim) {
    if(is.numeric(A[,i+1]) !=TRUE ) stop("Data Set contains non-numeric entries")
  }
  
  ranks_s2 = rep(0,dim)
  ranks_I2_hat = rep(0,dim)
  
  ranks_lowerBoundary = rep(0,dim)
  ranks_upperBoundary = rep(0,dim)
  ranks_CI_length = rep(0,dim)
  
  vect_CI = matrix(rep(0,2*dim),ncol=2,nrow=dim)
  
  ranks_s1 = rep(0,dim)
  ranks_I1_hat = rep(0,dim)
  
  ranks_lowerBoundary_1 = rep(0,dim)
  ranks_upperBoundary_1 = rep(0,dim)
  ranks_CI_length_1 = rep(0,dim)
  
  vect_CI_1 = matrix(rep(0,2*dim),ncol=2,nrow=dim)
  
  ranks_NO = rep(0,dim)
  ranks_NO_lB =rep(0,dim)
  ranks_NO_uB = rep(0,dim)
  vect_NO = matrix(rep(0,2*dim),ncol=2,nrow=dim)
  
  for(i in 2 : (dim+1)){
    
    x1 = Species_1[,i]
    x2 = Species_2[,i]
    
    ############## I2 Ranks #################################
    ranks_s2[i-1] = .getVarianceEstimate(x1,x2)
    ranks_I2_hat[i-1] = .getRankedSingleOverlapIndex(x1,x2)
    
    ranks_lowerBoundary[i-1] = pmax(ranks_I2_hat[i-1] - z*(1/sqrt(n))*sqrt(ranks_s2[i-1])/2,0)
    ranks_upperBoundary[i-1] = pmin(ranks_I2_hat[i-1] + z*(1/sqrt(n))*sqrt(ranks_s2[i-1])/2,1)
    ranks_CI_length[i-1] = abs(ranks_upperBoundary[i-1] - ranks_lowerBoundary[i-1])
    
    ############## I2 Bootstrap #################################
    vect_I2 = rep(0,length(reps4boot))
    for(j in 1:reps4boot){
      x1_bt = sample(x1,n,replace = TRUE)
      x2_bt = sample(x2,m,replace = TRUE)
      vect_I2[j] = .getRankedSingleOverlapIndex(x1_bt,x2_bt)
    }
    #Lower and Upper Boundary of the CI
    vect_CI[i-1,1] = max(as.numeric(quantile(vect_I2,0.025)),0)
    vect_CI[i-1,2] = min(as.numeric(quantile(vect_I2,0.975)),1)
    
    ############## I1 Ranks #################################
    ranks_s1[i-1] = .getVarianceEstimate(x2,x1)
    ranks_I1_hat[i-1] = .getRankedSingleOverlapIndex(x2,x1)
    
    ranks_lowerBoundary_1[i-1] = pmax(ranks_I1_hat[i-1] - z*(1/sqrt(n))*sqrt(ranks_s1[i-1])/2,0)
    ranks_upperBoundary_1[i-1] = pmin(ranks_I1_hat[i-1] + z*(1/sqrt(n))*sqrt(ranks_s1[i-1])/2,1)
    ranks_CI_length_1[i-1] = abs(ranks_upperBoundary_1[i-1] - ranks_lowerBoundary_1[i-1])
    
    ############## I1 Bootstrap #################################
    vect_I1 = rep(0,length(reps4boot))
    
    for(j in 1:reps4boot){
      x1_bt = sample(x1,n,replace = TRUE)
      x2_bt = sample(x2,m,replace = TRUE)
      vect_I1[j] = .getRankedSingleOverlapIndex(x2_bt,x1_bt)
    }
    
    #Lower and Upper Boundary of the CI
    vect_CI_1[i-1,1] = max(as.numeric(quantile(vect_I1,0.025)),0)
    vect_CI_1[i-1,2] = min(as.numeric(quantile(vect_I1,0.975)),1)
    
    ######### NO #####
    ranks_NO[i-1] = ranks_I1_hat[i-1] * ranks_I2_hat[i-1] * 4
    ranks_NO_lB[i-1] = ranks_lowerBoundary_1[i-1] * ranks_lowerBoundary[i-1] * 4
    ranks_NO_uB[i-1]= ranks_upperBoundary_1[i-1] * ranks_upperBoundary[i-1] * 4
    
    vect_NO[i-1,1] = vect_CI_1[i-1,1] * vect_CI[i-1,1] * 4
    vect_NO[i-1,2] = min(vect_CI_1[i-1,2] * vect_CI[i-1,2] * 4 ,1)
  }

  I1_d = .geo.mean(ranks_I1_hat)
  I2_d = .geo.mean(ranks_I2_hat)
  Values = data.frame(ranks_I1_hat, ranks_lowerBoundary_1, ranks_upperBoundary_1,vect_CI_1,
                      ranks_I2_hat, ranks_lowerBoundary, ranks_upperBoundary, vect_CI)
  Values <- round(Values, digit)
  colnames(Values) <- c("I1", "CI1_lower", "CI1_upper", "Boot_CI1_lower", "Boot_CI1_upper",
                        "I2", "CI2_lower", "CI2_upper", "Boot_CI2_lower", "Boot_CI2_upper")
  Values <- rbind(Values, data.frame(I1=round(I1_d, digit), CI1_lower=" ", CI1_upper=" ", Boot_CI1_lower=" ", Boot_CI1_upper = " ",
                                     I2 = round(I2_d, digit), CI2_lower=" ", CI2_upper= " ", Boot_CI2_lower= " ", Boot_CI2_upper=" "))
  
  return(Values)
  
}

