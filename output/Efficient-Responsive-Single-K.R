#################
# R code for analyzing output and plot figures
# v1.0 (organized on 2013-04-17)
# for Figure 10 and 11 in the paper
#################

#NEED TO FIRST SET R WORKING DIRECTORY TO WHERE THE FILES ARE LOCATED!!!
    setwd("/Users/buxx/Desktop/test/")

#read the output file
    data <- read.table("Efficient-Responsive-Single-Nash-K.txt", header=TRUE)



#prepare data for Figure 10 and 11

    data.K1 <- data[data$K==1,]
    if (nrow(data.K1[data.K1$TS12==20,]) >0 ) data.K1[data.K1$TS12==20,]$TS12 <- NA

    data.K3 <- data[data$K==3,]
    if (nrow(data.K3[data.K3$NashE==20,]) >0 ) data.K3[data.K3$NashE==20,]$NashE <- NA
    if (nrow(data.K3[data.K3$NashR==20,]) >0 ) data.K3[data.K3$NashR==20,]$NashR <- NA
    if (nrow(data.K3[data.K3$TS12==20,]) >0 ) data.K3[data.K3$TS12==20,]$TS12 <- NA


#plot figure 10 (K=1)

    pdf('Figure10-K1.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K1$b, data.K1$NashE, lty=1, lwd=3, col="red")
    lines(data.K1$b, data.K1$NashR, lty=1, lwd=3, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K1$b, data.K1$TS01, lty=1, lwd=3, col="darkolivegreen")
    lines(data.K1$b, data.K1$TS12, lty=1, lwd=3, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main="Figure 10. K = 1", outer=T)

    dev.off()


#plot figure 11 (K=3)
    pdf('Figure11-K3.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K3$b, data.K3$NashE, lty=1, lwd=3, col="red")
    lines(data.K3$b, data.K3$NashR, lty=1, lwd=3, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.K3$b, data.K3$TS01, lty=1, lwd=3, col="darkolivegreen")
    lines(data.K3$b, data.K3$TS12, lty=1, lwd=3, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main="Figure 11. K = 3", outer=T)

    dev.off()
