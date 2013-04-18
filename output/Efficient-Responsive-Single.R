#################
# R code for analyzing output and plot figures
# v1.0 (organized on 2013-04-17)
# for Figure 4 in the paper
#################

#NEED TO FIRST SET R WORKING DIRECTORY TO WHERE THE FILES ARE LOCATED!!!
    setwd("/Users/buxx/Desktop/test/")

#read the output file
    data <- read.table("Efficient-Responsive-Single-Nash.txt", header=TRUE)


#prepare data for figure 4
#calculate the analytical approximation
    u <- 10
    c <- 1
    data$VOCE = data$b^3 * (16 - 8*data$b^2 - data$b^3) * (u - c)^2 / (4*(2-data$b^2)^2*(2+data$b)^2)
    data$VOCR = data$b^4 * (u-c)^2 / ( 8*(2-data$b^2))

#plot figure 4

    pdf('Figure4-Single.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data$b, data$NashE, lty=1, lwd=3, col="red")
    lines(data$b, data$NashR, lty=1, lwd=3, col="blue")
    lines(data$b, data$VOCE, lty=2, lwd=2, col="red")
    lines(data$b, data$VOCR, lty=2, lwd=2, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data$b, data$TS01, lty=1, lwd=3, col="darkolivegreen")
    lines(data$b, data$TS12, lty=1, lwd=3, col="darkolivegreen")
    lines(data$b, data$NashE, lty=1, lwd=2, col="red")
    lines(data$b, data$NashR, lty=1, lwd=2, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main="Figure 4. Single", outer=T)

    dev.off()
