#################
# R code for analyzing output and plot figures
# v1.0 (organized on 2013-04-17)
#################

#NEED TO FIRST SET R WORKING DIRECTORY TO WHERE THE FILES ARE LOCATED!!!
    setwd("/Users/buxx/Desktop/test/")

#read the output file
    data <- read.table("Efficient-Responsive-Holdback-Nash.txt", header=TRUE)



#prepare data for Figure 1

    data.s0 <- data[data$s==0,]
    data.s05 <- data[data$s==0.5,]
    data.s09 <- data[data$s==0.9,]
    data.s1 <- data[data$s==-1,]
    data.s20 <- data[data$s==-20,]



#plot figure

    pdf('Figure-Holdback.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)



    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.s0$b, data.s0$NashEh, lty=5, lwd=3, col="red")
    lines(data.s0$b, data.s0$NashRh, lty=5, lwd=3, col="blue")
    lines(data.s09$b, data.s09$NashEh, lty=1, lwd=3, col="red")
    lines(data.s09$b, data.s09$NashRh, lty=1, lwd=3, col="blue")
    lines(data.s20$b, data.s20$NashEh, lty=1, lwd=3, col="red")
    lines(data.s20$b, data.s20$NashRh, lty=1, lwd=3, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")


    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data.s09$b, data.s09$TS12h, lty=1, lwd=3, col="darkolivegreen")
    lines(data.s05$b, data.s05$TS12h, lty=1, lwd=3, col="darkolivegreen")
    lines(data.s0$b, data.s0$TS12h, lty=1, lwd=3, col="darkolivegreen")
    lines(data.s1$b, data.s1$TS12h, lty=1, lwd=3, col="darkolivegreen")
    lines(data.s20$b, data.s20$TS12h, lty=1, lwd=3, col="darkolivegreen")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

 
    title(main="Holdback", outer=T)

    dev.off()

