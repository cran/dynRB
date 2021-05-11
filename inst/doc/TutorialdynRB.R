## ----library, message=FALSE---------------------------------------------------
library(dynRB)
library(ggplot2)
library(reshape2)
library(vegan)
library(RColorBrewer)

## ----data---------------------------------------------------------------------
data(finch)
head(finch)[,1:5]

## ----finch, results = "hide"--------------------------------------------------
r <- dynRB_VPa(finch)

## ----result-------------------------------------------------------------------
res <- r$result
res[1:6, 1:5]

## ----result2------------------------------------------------------------------
res[1:6, c(1:2, 6:11)]

## ----result3, results = "hide"------------------------------------------------
r <- dynRB_Pn(finch)

## ----result31-----------------------------------------------------------------
head(r$result)

## ----result4, results = "hide"------------------------------------------------
r <- dynRB_Vn(finch)

## ----result41-----------------------------------------------------------------
head(r$result)

## ----result5, results = "hide"------------------------------------------------
r <- dynRB_VPa(finch)
# aggregation method product  
om1 <- reshape(r$result[,1:3], direction="wide", idvar="V1", timevar="V2") 
#aggregation method mean  
om2 <- reshape(r$result[,c(1,2,4)], direction="wide", idvar="V1", timevar="V2") 
#aggregation method gmean  
om3 <- reshape(r$result[,c(1,2,5)], direction="wide", idvar="V1", timevar="V2") 


## ----result51-----------------------------------------------------------------
om1

## ----plot1, results = "hide", fig.width=7, fig.height=5, message=FALSE, warning=FALSE----
r <- dynRB_VPa(finch)
#aggregation method mean  
om <- reshape(r$result[,c(1,2,4)], direction="wide", idvar="V1", timevar="V2") 
mantel(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)])), permutations = 1000)
plot(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)]))) 

## ----plot2, results = "hide", fig.width=6, fig.height=4-----------------------
r <- dynRB_VPa(finch)  #main function of dynRB to calculate size and overlap of hypervolumes
result <- r$result
Overlap <- as.numeric(ifelse(result$V1 == result$V2, "NA", result$port_prod))  
# 'result$port_prod' may be changed to 'result$port_mean' or 'result$port_gmean'
is.numeric(Overlap)
Result2 <- cbind(result, Overlap)
ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black")

## ----plot3, results = "hide", fig.width=6, fig.height=4-----------------------
r <- dynRB_VPa(finch)  
theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(colour="black", size = rel(1.5), angle=35, hjust = 1),
  axis.text.y = element_text(colour="black", size = rel(1.5)),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)
result <- r$result
Overlap <- as.numeric(ifelse(result$V1 == result$V2, "NA", result$port_prod))  
# 'result$port_prod' may be changed to 'result$port_mean' or 'result$port_gmean'
is.numeric(Overlap)
Result2<-cbind(result, Overlap)
breaks <- seq(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE),  
              by=round(max(Overlap, na.rm=TRUE)/10, digits=3))
col1 <- colorRampPalette(c("white", "navyblue")) #define color gradient
ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black") +
  scale_fill_gradientn(colours=col1(8), breaks=breaks, guide="colorbar",  
                       limits=c(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE))) +
  theme_change

