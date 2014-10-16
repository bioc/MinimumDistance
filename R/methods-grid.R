#' @aliases plotDenovo,MinDistExperiment,MDRanges-method
#' @rdname plotDenovo
setMethod("plotDenovo",  c("MinDistExperiment", "MDRanges"),
          function(object, g, param){
            ##if(missing(maxgap)) maxgap <- 50*width(g)
            FUN <- VanillaICE:::expandFun(param)
            ylimits <- VanillaICE:::yLimits(param)
            maxgap <- FUN(g)

            index <- subjectHits(findOverlaps(g, object, maxgap=maxgap))
            me <- object
            r <- threshold(as.numeric(lrr(me)[index, ]), lim=ylimits[[2]])

            df <- data.frame(r=r,
                             b=as.numeric(baf(me)[index, ]),
                             pos=rep(start(me)[index], ncol(me)),
                             id=as.character(matrix(colnames(me), length(index), ncol(me), byrow=TRUE)))

            locs <- pretty(range(df$pos), n=5)
            labels <- prettyNum(locs/1000, big.mark=",")

            ylim <- ylimits[[2]]
            exponent <- pretty(seq(ylim[1], ylim[2]), n=6)
            exponent <- exponent[exponent %% 1 == 0]
            is_notpositive <- exponent <= 0
            denom <- 2^(-exponent[is_notpositive])
            ylabels <- paste0("1/", denom)
            if(any(!is_notpositive)){
              ylabels <- c(ylabels, 2^exponent[!is_notpositive])
            }

            rfig <- xyplot(r~pos | id, df, pch=20, cex=0.3, col="gray", layout=c(1,ncol(me)),
                           panel=function(x,y, pch, col, ..., denovo_region,  id, acfs, subscripts){
                             id <- unique(id[subscripts])
                             ##local_acf <- round(acfs[id], 2)
                             panel.xyplot(x, y, pch=pch, col=col, ...)
                             panel.abline(h=0, col="gray")
                             panel.abline(v=c(start(denovo_region), end(denovo_region)))
                             ##panel.text(min(x), 1.3, paste0("  lag10 acf:", local_acf), col="black", cex=0.8)
                             if(id == gsub("md_", "", denovo_region$sample)){
                               panel.points(x, y, col="black", pch=20,...)
                             }
                           }, denovo_region=g, ylim=ylim, id=df$id,
                           scales=list(cex=0.6, x=list(at=locs, labels=labels, axs="i"),
                             y=list(at=exponent, labels=ylabels)),
                           xlab="kb", ylab="fold change")
            bfig <- xyplot(b~pos | id, df, pch=20, cex=0.3, col="gray", layout=c(1,ncol(me)),
                           panel=function(x,y, col, pch, ...,
                             denovo_region, id, subscripts){
                             id <- unique(id[subscripts])
                             panel.xyplot(x, y, col=col, pch=pch, ...)
                             panel.abline(v=c(start(denovo_region), end(denovo_region)))
                             if(id == gsub("md_", "", denovo_region$sample)){
                               panel.points(x, y, col="black", pch=20,...)
                             }
                           }, denovo_region=g, id=df$id,
                           scales=list(cex=0.6, x=list(at=locs, labels=labels, axs="i")),
                           xlab="kb")
            list(rfig, bfig)
          })

#' Default viewports for plotting log R ratios, BAFs, chromosome
#' idiogram, and a legend for a case-parent trio
#'
#' @seealso \code{\link{plotDenovo}} \code{\link{pedigreeGrid}}
#' @examples
#' vps <- pedigreeViewports()
#' @export
pedigreeViewports <- function(){
  lrr <- viewport(x=0,
                  y=0.01,
                  width=unit(0.5, "npc"),
                  height=unit(0.88, "npc"), just=c("left", "bottom"),
                  name="lvp1")
  baf <- viewport(x=0.5,
                  y=0.01,
                  width=unit(0.5, "npc"),
                  height=unit(0.88, "npc"), just=c("left", "bottom"),
                  name="lvp2")
  idiogram <- viewport(x=0.5, y=0.87,
                       width=unit(0.75, "npc"),
                       height=unit(0.15, "npc"),
                       just=c("center", "bottom"))

  legend <- viewport(x=0.88, y=0.9, width=unit(0.15, "npc"),
                     height=unit(0.1, "npc"),
                     just=c("left", "bottom"))

  list(lrr=lrr, baf=baf, idiogram=idiogram, legend=legend)
}

#' Plot the log R ratios and BAFs on a grid given by precomputed viewports
#'
#' @param g a \code{MDRanges} object
#' @param vps a list of viewports.  See \code{\link{pedigreeViewports}}.
#' @param figs a list of trellis objects created by the function \code{\link{plotDenovo}}.
#' @examples
#' library(GenomicRanges)
#' library(VanillaICE)
#' require(grid)
#' ##marker-level summaries
#' data(md_exp)
#' seqlevels(md_exp, force=TRUE) <- "chr22"
#' ## segmentation results
#' data(md_gr)
#' posteriorCalls <- MAP2(md_exp, md_gr, MinDistParam())
#' g <- denovoHemizygous(posteriorCalls)
#' g
#' vps <- MinimumDistance:::pedigreeViewports()
#' param <- HmmTrellisParam()
#' p <- plotDenovo(md_exp, g[1], param)
#' p <- pedigreeGrid(g=g[1], vps=vps, figs=p)
#' leg <- mdLegend(g[1])
#' upViewport(0)
#' pushViewport(vps[["legend"]])
#' grid.text(leg, x=unit(0.02, "npc"), y=unit(0.95, "npc"), just=c("left", "top"),
#'           gp=gpar(cex=0.6, fontfamily="mono"))
#' ##
#' ## combine adjacent denovo hemizygous
#' ##
#' g2 <- reduce(denovoHemizygous(posteriorCalls), min.gapwidth=500e3)
#' post <- MAP2(md_exp, g2)
#' g2 <- denovoHemizygous(post)
#' p <- plotDenovo(md_exp, g2, param)
#' p <- pedigreeGrid(g=g2, vps=vps, figs=p)
#' leg <- mdLegend(g2)
#' upViewport(0)
#' pushViewport(vps[["legend"]])
#' grid.text(leg, x=unit(0.02, "npc"), y=unit(0.95, "npc"), just=c("left", "top"),
#'           gp=gpar(cex=0.6, fontfamily="mono"))
#' @seealso \code{\link{plotDenovo}} \code{\link{pedigreeViewports}}
#' @export
pedigreeGrid <- function(g, vps, figs){
  xlim <- VanillaICE:::.find_xlim_percent(g, 0.05)
  chr <- chromosome(g)
  sl <- seqlengths(g)[chr]
  iparams <- IdiogramParams(seqnames=chr,
                            genome="hg19",
                            seqlengths=sl,
                            box=list(xlim=xlim, color="blue"))
  idiogram <- VanillaICE::plot(iparams)
  grid.newpage()
  vp <- vps[["idiogram"]]
  pushViewport(vp)
  print(idiogram, vp=vp, newpage=FALSE)
  upViewport(0)
  vp <- vps[["lrr"]]
  pushViewport(vp)
  print(figs[[1]], vp=vp, newpage=FALSE)
  upViewport(0)
  vp <- vps[["baf"]]
  pushViewport(vp)
  print(figs[[2]], vp=vp, newpage=FALSE)
  upViewport(0)
  id <- g[1]$sample
  grid.text(paste0(id, "\n", chr), x=unit(0.01, "npc"), y=unit(0.98, "npc"),
            just=c("left", "top"))
}

#' Text summary of information encapculated in a MDRanges object for a particular interval
#'
#' @param g a \code{MDRanges} object
#' @export
mdLegend <- function(g){
  size <- prettyNum(round(width(g)/1000, 1), big.mark=",")
  fold_change <- round(2^g$seg.mean, 1)
  md_map <- g$calls
  pr.map <- round(posteriorMAP(g), 2)
  podds <- round(posteriorOdds(g), 1)
  rr <- pretty(round(posteriorRR(g), 2))
  labels <- paste0("size (kb)     :", size, "\n",
                   "fold change   :", fold_change, "\n",
                   "MAP estimate  :", md_map, "\n",
                   "posterior Pr  :", pr.map, "\n",
                   "posterior odds\n",
                   "221 or 220    :", podds, "\n",
                   "posterior RR  :", rr)
  labels
}
