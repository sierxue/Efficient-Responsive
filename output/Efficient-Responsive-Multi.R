#################
# R code for analyzing output and plot figures
# v1.0 (organized on 2013-04-17)
# for Figure 6 in the paper
#################

#NEED TO FIRST SET R WORKING DIRECTORY TO WHERE THE FILES ARE LOCATED!!!
    setwd("/Users/buxx/Desktop/test/")

#read the output file
    data <- read.table("Efficient-Responsive-Multi-Nash.txt", header=TRUE)



#prepare data for Figure 6
    data.K1 <- data[data$K==1,]
    if (nrow(data.K1[data.K1$TS12==20,]) >0 ) data.K1[data.K1$TS12==20,]$TS12 <- NA

    data.K3 <- data[data$K==3,]
    if (nrow(data.K3[data.K3$NashE==20,]) >0 ) data.K3[data.K3$NashE==20,]$NashE <- NA
    if (nrow(data.K3[data.K3$NashR==20,]) >0 ) data.K3[data.K3$NashR==20,]$NashR <- NA
    if (nrow(data.K3[data.K3$TS12==20,]) >0 ) data.K3[data.K3$TS12==20,]$TS12 <- NA


#plot figure 6
    pdf('Figure6-Multi.pdf', width = 12, height = 7)
    par(oma=c(0,0,2,0))
    par(mfrow=c(1,2))


    xrange = c(0, 1)
    yrange = c(0, 20)

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data$b, data$NashEm, lty=1, lwd=3, col="red")
    lines(data$b, data$NashRm, lty=1, lwd=3, col="blue")
    lines(data$b, data$NashE, lty=2, lwd=2, col="red")
    lines(data$b, data$NashR, lty=2, lwd=2, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(a) Firms' Strategic Choices")

    plot(xrange, yrange, type="n", xlab="b", ylab=expression(paste("VOI (", sigma^2, ")")) , xaxt="n", yaxt="n")
    lines(data$b, data$TS01m, lty=1, lwd=3, col="darkolivegreen")
    lines(data$b, data$TS12m, lty=1, lwd=3, col="darkolivegreen")
    lines(data$b, data$NashEm, lty=1, lwd=2, col="red")
    lines(data$b, data$NashRm, lty=1, lwd=2, col="blue")
    axis(side=1, at=seq(0,1,0.1), labels=seq(0,1,0.1))
    axis(side=2, at=seq(0,20,2), labels=seq(0,20,2))
    title(main="(b) Socially Efficient Strategies")

    title(main="Figure 6. Multi", outer=T)

    dev.off()

