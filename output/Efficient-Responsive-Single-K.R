#################
# R code for analyzing output and plot figures
# v1.0 (organized on 2013-04-17)
#################

#NEED TO FIRST SET R WORKING DIRECTORY TO WHERE THE FILES ARE LOCATED!!!
    setwd("/Users/buxx/Desktop/test/")

#read the output file
    data <- read.table("Efficient-Responsive-Single-Nash-K.txt", header=TRUE)



#prepare data for Figure 1
    data.K1 <- data[data$K==1,]
    if (nrow(data.K1[data.K1$TS01==0,]) >0 ) data.K1[data.K1$TS01==0,]$TS01 <- NA
    if (nrow(data.K1[data.K1$TS12==20,]) >0 ) data.K1[data.K1$TS12==20,]$TS12 <- NA

    data.K3 <- data[data$K==3,]
    if (nrow(data.K3[data.K3$NashE==20,]) >0 ) data.K3[data.K3$NashE==20,]$NashE <- NA
    if (nrow(data.K3[data.K3$NashR==20,]) >0 ) data.K3[data.K3$NashR==20,]$NashR <- NA
    if (nrow(data.K3[data.K3$TS01==0,]) >0 ) data.K3[data.K3$TS01==0,]$TS01 <- NA
    if (nrow(data.K3[data.K3$TS12==20,]) >0 ) data.K3[data.K3$TS12==20,]$TS12 <- NA


#plot figure for K=1
    pdf('Figure-K1.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K1$b, data.K1$NashE, type="l", lwd=4, col="red")
    lines(data.K1$b, data.K1$NashR, type="l", lwd=4, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K1$b, data.K1$TS01, type="l", lwd=4, col="darkolivegreen")
    lines(data.K1$b, data.K1$TS12, type="l", lwd=4, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main="K = 1", outer=T)

    dev.off()


#plot figure for K=3
    pdf('Figure-K3.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K3$b, data.K3$NashE, type="l", lwd=4, col="red")
    lines(data.K3$b, data.K3$NashR, type="l", lwd=4, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K3$b, data.K3$TS01, type="l", lwd=4, col="darkolivegreen")
    lines(data.K3$b, data.K3$TS12, type="l", lwd=4, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main="K = 3", outer=T)

    dev.off()
