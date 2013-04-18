#################
# R code for analyzing output and plot figures
# v1.0 (organized on 2013-04-17)
#################

#NEED TO FIRST SET R WORKING DIRECTORY TO WHERE THE FILES ARE LOCATED!!!
    setwd("/Users/buxx/Desktop/test/")

#read the output file
    data <- read.table("Efficient-Responsive-Single-Nash-delta.txt", header=TRUE)



#prepare data for Figure 1

    data.delta2 <- data[data$delta==1.2,]
    if (nrow(data.delta2[data.delta2$TS12==20,]) >0 ) data.delta2[data.delta2$TS12==20,]$TS12 <- NA

    data.delta5 <- data[data$delta==1.5,]
    if (nrow(data.delta5[data.delta5 $NashE==20,]) >0 ) data.delta5[data.delta5 $NashE==20,]$NashE <- NA
    if (nrow(data.delta5[data.delta5 $NashR==20,]) >0 ) data.delta5[data.delta5 $NashR==20,]$NashR <- NA
    if (nrow(data.delta5[data.delta5 $TS12==20,]) >0 ) data.delta5[data.delta5 $TS12==20,]$TS12 <- NA


#plot figure for delta=1.2

    pdf('Figure-delta2.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.delta2 $b, data.delta2 $NashE, lty=1, lwd=3, col="red")
    lines(data.delta2 $b, data.delta2 $NashR, lty=1, lwd=3, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.delta2 $b, data.delta2 $TS01, lty=1, lwd=3, col="darkolivegreen")
    lines(data.delta2 $b, data.delta2 $TS12, lty=1, lwd=3, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main=expression(paste(Delta, " = 1.2")), outer=T)

    dev.off()


#plot figure for delta=1.5
    pdf('Figure-delta5.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.delta5 $b, data.delta5 $NashE, lty=1, lwd=3, col="red")
    lines(data.delta5 $b, data.delta5 $NashR, lty=1, lwd=3, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.delta5 $b, data.delta5 $TS01, lty=1, lwd=3, col="darkolivegreen")
    lines(data.delta5 $b, data.delta5 $TS12, lty=1, lwd=3, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main=expression(paste(Delta, " = 1.5")), outer=T)

    dev.off()
