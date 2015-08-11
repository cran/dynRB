funclist3<-list(
  portionAinB2_full=function(S1,S2,steps=101,alpha_grid){
    steps0<-steps
    alpha_grid0<-alpha_grid
    funclist2$portionAinB2_full_severalcol(S1,S2,steps=steps0,alpha_grid=alpha_grid0)
  },
  volumeA2_full=function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
    steps0<-steps
    alpha0_grid<-alpha_grid
    funclist2$volumeA2_full_severalcol(S1,S2,steps=steps0,alpha_grid=alpha0_grid  )
  },
  portionAinB_coordinates_full=function(S1,S2,steps=101 ){ 
    nt<-ncol(S1)
    integral_coord<-rep(0,length=nt)
    hilf<-function(x1,x2,steps){
      r<-funclist2$portionAinB_full_onecolumn(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
      return(r$integral_approx)
    }
    hilf2<-function(x1,x2,steps){
      r<-funclist2$portionAinB_full_onecolumn(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
      return(r$overlap)
    }
    integral_coord<-mapply(hilf, x1=S1, x2=S2, MoreArgs =list(steps=steps))
    plot_data_overlap<-mapply(hilf2, x1=S1, x2=S2, MoreArgs =list(steps=steps))
    alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
    erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_data_overlap=plot_data_overlap)
    return(erg)
  },
  volumeA_coordinates_full=function(S1,S2,steps=101 ){ 
    nt<-ncol(S1)
    hilf<-function(x1,x2,steps){
      r<-funclist2$volumeA_full_onecol(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
      return(r$integral_approx)
    }
    hilf2<-function(x1,x2,steps){
      r<-funclist2$volumeA_full_onecol(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
      return(r$volume)
    }
    integral_coord<-mapply(hilf, x1=S1, x2=S2, MoreArgs =list(steps=steps))
    plot_volume<-mapply(hilf2, x1=S1, x2=S2, MoreArgs =list(steps=steps))
    alpha_grid<-seq(0,1,length=steps)    
    erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_volume=plot_volume)
    return(erg)
  }
)