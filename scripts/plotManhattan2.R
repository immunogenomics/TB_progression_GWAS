plotManhattan2<-function (bedfile, chrom = NULL, chromstart = NULL, chromend = NULL, 
                         pvalues, ld,genome = NULL, col = SushiColors(5), space = 0.01, 
                         ymax = 1.05,ldColors = c("lightskyblue","navy","green","orange","red"), ...) 
{
  if (is.null(genome) == FALSE) {
    chromoffsets = chromOffsets(genome, space)
    if (class(col) == "function") {
      col = col(nrow(chromoffsets))
    }
    else {
      col = rep(col, ceiling(nrow(chromoffsets)/length(col)))
    }
    bedfile = bedfile[bedfile[, 1] %in% chromoffsets[, 1], 
                      ]
    columber = 0
    columbers = rep(0, nrow(bedfile))
    for (i in (1:nrow(chromoffsets))) {
      columber = columber + 1
      if (columber > length(col)) {
        columber = 1
      }
      rowsofinterest = which(bedfile[, 1] == chromoffsets[i, 
                                                          1])
      bedfile[rowsofinterest, 2] = bedfile[rowsofinterest, 
                                           2] + chromoffsets[i, 3]
      columbers[rowsofinterest] = columber
    }
    cumsums = cumsum(as.numeric(genome[, 2]))
    spacer = cumsums[length(cumsum(as.numeric(genome[, 2])))] * 
      space
    yrange = c(min(-log10(bedfile[, 5])), max(-log10(bedfile[, 
                                                             5])) * ymax)
    #LD Colors from locusZoom
    ldColors = c("navy","lightskyblue","green","orange","red")
    
    #This adds a column of color values
    # based on the y values
    #ld<-bedfile$ld
    #ldCol <- ldColors[as.numeric(cut(bedfile[,6],breaks = 5))]
    
    plot(bedfile[, 2], -log10(bedfile[, 5]), col = ldCol, 
         xlim = c(0, max(chromoffsets[, 4]) + spacer), ylim = yrange, 
         pch = 20, xaxt = "n", xlab = "", yaxt = "n", ylab = "", 
         xaxs = "i", yaxs = "i",cex.axis=2, ...)
  }
  if (is.null(genome) == TRUE) {
    bedfile = bedfile[which(bedfile[, 1] == chrom), ]
    plot(bedfile[, 2], -log10(bedfile[, 5]), col = ldColors[as.numeric(cut(bedfile[,6],breaks = 5))], xlim = c(chromstart, 
                                                                 chromend), pch = 19, xaxt = "n", xlab = "", yaxt = "n", 
         ylab = "", xaxs = "i", yaxs = "i", cex.axis=2,...)
    #customized gnomewide significant
    abline(h=-log10(1.83e-7),col=1,lty=2)
    legend.text<-c("0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2")
    legend("topright", legend=legend.text,
           fill=rev(ldColors), cex=2,title=expression(r^2),y.intersp=.5)
  }
}

zoomsregion2<-function (region, chrom = NULL, genome = NULL, space = 0.01, 
          padding = 0.00, col = NA, zoomborder = "red", lty = 2, 
          lwd = 1, extend = 0, wideextend = 0.1, offsets = c(0, 0), 
          highlight = FALSE) 
{
  if (is.null(genome) == FALSE) {
    chromoffsets = chromOffsets(genome, space)
    region[1] = region[1] + chromoffsets[which(chromoffsets[, 
                                                            1] == chrom), 3]
    region[2] = region[2] + chromoffsets[which(chromoffsets[, 
                                                            1] == chrom), 3]
    if (abs(region[2] - region[1]) < 2 * padding * (max(chromoffsets[, 
                                                                     4]))) {
      center = mean(region)
      region[1] = center - padding * max(chromoffsets[, 
                                                      4])
      region[2] = center + padding * max(chromoffsets[, 
                                                      4])
    }
  }
  orgimal_xpd = par()$xpd
  par(xpd = TRUE)
  plotlef = par("usr")[1]
  plotrig = par("usr")[2]
  plotbot = par("usr")[3]
  plottop = par("usr")[4]
  xrange = abs(plotrig - plotlef)
  yrange = abs(plotbot - plottop)
  if (length(extend) == 1) {
    extend.upper = extend * yrange
    extend.lower = extend * yrange
  }
  if (length(extend) > 1) {
    extend.upper = extend[1] * yrange
    extend.lower = extend[2] * yrange
  }
  plotlef = plotlef + xrange * offsets[1]
  plotrig = plotrig - xrange * offsets[2]
  miny = plotbot - extend.lower - (wideextend * yrange)
  lowy = plotbot - extend.lower - (5 * yrange)
  if (is.null(region) == FALSE) {
    currentxpd = par()$xpd
    par(xpd = TRUE)
    if ((abs(region[1] - region[2])/xrange) < (2 * padding)) {
      center = mean(region)
      region[1] = center - padding * xrange
      region[2] = center + padding * xrange
    }
    if (highlight == TRUE) {
      plotrig = region[2]
      plotlef = region[1]
      miny = plotbot - extend.lower
      lowy = plotbot - extend.lower
    }
    polygon(x = c(region[1], region[1], region[2], region[2], 
                  plotrig, plotrig, plotlef, plotlef, region[1]), y = c(plotbot - 
                                                                          extend.lower, plottop + extend.upper, plottop + extend.upper, 
                                                                        plotbot - extend.lower, miny, lowy, lowy, miny, plotbot - 
                                                                          extend.lower), col = col, border = zoomborder, 
            lty = lty, lwd = lwd)
    par(xpd = currentxpd)
  }
  par(xpd = orgimal_xpd)
}

zoomsregion3<-function (region, chrom = NULL, top=.1,genome = NULL, space = 0.01, 
                        padding = 0.00, col = NA, zoomborder = "red", lty = 2, 
                        lwd = 1, extend = 0, wideextend = 0.1, offsets = c(0, 0), 
                        highlight = FALSE) 
{
  if (is.null(genome) == FALSE) {
    chromoffsets = chromOffsets(genome, space)
    region[1] = region[1] + chromoffsets[which(chromoffsets[, 
                                                            1] == chrom), 3]
    region[2] = region[2] + chromoffsets[which(chromoffsets[, 
                                                            1] == chrom), 3]
    if (abs(region[2] - region[1]) < 2 * padding * (max(chromoffsets[, 
                                                                     4]))) {
      center = mean(region)
      region[1] = center - padding * max(chromoffsets[, 
                                                      4])
      region[2] = center + padding * max(chromoffsets[, 
                                                      4])
    }
  }
  orgimal_xpd = par()$xpd
  par(xpd = TRUE)
  plotlef = par("usr")[1]
  plotrig = par("usr")[2]
  plotbot = par("usr")[3]
  plottop2 = top
  plottop = par("usr")[4]
  xrange = abs(plotrig - plotlef)
  yrange = abs(plotbot - plottop)
  if (length(extend) == 1) {
    extend.upper = extend * yrange
    extend.lower = extend * yrange
  }
  if (length(extend) > 1) {
    extend.upper = extend[1] * yrange
    extend.lower = extend[2] * yrange
  }
  plotlef = plotlef + xrange * offsets[1]
  plotrig = plotrig - xrange * offsets[2]
  miny = plotbot - extend.lower - (wideextend * yrange)
  lowy = plotbot - extend.lower - (5 * yrange)
  if (is.null(region) == FALSE) {
    currentxpd = par()$xpd
    par(xpd = TRUE)
    if ((abs(region[1] - region[2])/xrange) < (2 * padding)) {
      center = mean(region)
      region[1] = center - padding * xrange
      region[2] = center + padding * xrange
    }
    if (highlight == TRUE) {
      plotrig = region[2]
      plotlef = region[1]
      miny = plotbot - extend.lower
      lowy = plotbot - extend.lower
    }
    polygon(x = c(region[1], region[1], region[2], region[2], 
                  plotrig, plotrig, plotlef, plotlef, region[1]), y = c(plotbot*0.1 - 
                                                                          extend.lower, plottop2 + extend.upper, plottop2 + extend.upper, 
                                                                        plotbot - extend.lower, miny, lowy, lowy, miny, plotbot - 
                                                                          extend.lower), col = col, border = zoomborder, 
            lty = lty, lwd = lwd)
    par(xpd = currentxpd)
  }
  par(xpd = orgimal_xpd)
}