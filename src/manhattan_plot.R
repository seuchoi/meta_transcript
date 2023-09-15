#######
####### manhattan plot
#### manhattan plot for Jun 02 2018
man<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
    genomewideline = -log10(5e-08), highlight = NULL, highlight2=NULL, logp = TRUE,
    ...)
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x)))
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x)))
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x)))
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x)))
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]]))
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]]))
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]]))
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]]))
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        options(scipen = 999)
        d$pos = d$BP/1e+06
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs = ticks
    }else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index ==
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP +
                  lastbase
            }
            ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.001)
    xmin = floor(max(d$pos) * -0.001)
#    def_args <- list(xaxt = "n",yaxt = "n", bty = "l",  xaxs = "i", yaxs = "i",
#        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
#            ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10]("p-value")),mgp=c(2,.5,0))
    def_args <- list(xaxt = "n",yaxt = "n", bty = "l",  xaxs = "i", yaxs = "i",
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
            ceiling(max(d$logp)+0.5)), xlab = "", ylab = "",mgp=c(2,.5,0))


    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1,cex.axis=0.8,tck=-0.03, line = 0, ...)
    }
    else {
        axis(1,cex.axis=0.8,tck=-0.03, line = 0, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
         for (i in unique(d$index)) {

                 with(d[d$index == unique(d$index)[i], ], points(pos,
                    logp, col = col[icol], pch = 20, ...))

     	#d1<-subset(d,logp<=20)
      #d2<-subset(d,logp>=45)
      #      with(d1[d1$index == unique(d1$index)[i], ], points(pos,
      #          logp, col = col[icol], pch = 20, ...))
      #      with(d2[d2$index == unique(d2$index)[i], ], points(pos,
      #          logp-25, col = col[icol], pch = 20, ...))
            icol = icol + 1

        }
    }
    if (suggestiveline)
        abline(h = suggestiveline, col = "blue",lty=3)
    if (genomewideline)
        abline(h = genomewideline, col = "black",lty=3)
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")

      d.highlight = d[which(d$SNP %in% highlight), ]
      with(d.highlight, points(pos, logp, col = "#377eb8", pch = 20, ...))

      #  d1<-subset(d,logp<=20)
      #  d2<-subset(d,logp>=45)
      #  d1.highlight = d1[which(d1$SNP %in% highlight), ]
      #  with(d1.highlight, points(pos, logp, col = "#377eb8", pch = 20, ...))
      #  d2.highlight = d2[which(d2$SNP %in% highlight), ]
      #  with(d2.highlight, points(pos, logp-25, col = "#377eb8", pch = 20,...))
    }
       if (!is.null(highlight2)) {
		if (any(!(highlight2 %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")

            d.highlight2 = d[which(d$SNP %in% highlight2), ]
            with(d.highlight2, points(pos, logp, col = "#e41a1c", pch = 20, ...))

      #  d1<-subset(d,logp<20)
		  # d1.highlight2 = d1[which(d1$SNP %in% highlight2), ]
      #  with(d1.highlight2, points(pos, logp, col="#e41a1c", pch = 20, ...))

    }

}
