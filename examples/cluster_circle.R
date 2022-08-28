##
draw.line <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
  w  <- (w/360)*2*pi;
  x1 <- xc+l1*cos(w);
  y1 <- yc-l1*sin(w);
  x2 <- xc+l2*cos(w);
  y2 <- yc-l2*sin(w);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
}
##
draw.arc.s <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  # s = simple
  # r = radius
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 20){
    pix.n <- 20;
  }
  
  ang.seq <- rev(seq(w1, w2, length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  ## lend=0(round); lend=1(butt); lend=2(square)
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
  #points(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
}

### the options should be the same as plot.phylo {ape}
cluster.circle <- function(cx, cy, r1, r2, data=data, angle.start=0,  angle.end=360, tip=TRUE, tip.cex=1,
                            edge.color = NULL, tip.color=NULL, edge.width = 1, lab.d=1, cex=1){
  angle.start <- angle.start + 270;
  angle.end   <- angle.end + 270;
  dat.h  <- data;
  num    <- length(dat.h$labels);
  if (missing(edge.color) | is.null(edge.color)){
    col.e <- rep("blue", num);
  } else {
    col.e <- rep(edge.color, num, length.out=num);
  }
  if (missing(tip.color) | is.null(tip.color)){
    col.t <- rep("blue", num);
  } else {
    col.t <- rep(tip.color, num, length.out=num);
  }
  lwd    <- edge.width;

  cols <- col.t;
  ##  radius  
  r.l <- r2 - r1;

  max.h        <- max(dat.h$height);
  #dat.h$height <- max.h - dat.h$height;
  max.n        <- length(dat.h$labels);
  
  lab   <- (angle.end - angle.start)/max.n;
  ratio <- r.l/max.h;
  
  ####################  
  sub2h   <- c();
  sub2po  <- c();
  sub2col <- c();

  for (i in 1:nrow(dat.h$merge)){
    le <- dat.h$merge[i,1];
    ri <- dat.h$merge[i,2];
    if (le < 0 && ri < 0){
      col.le <- cols[abs(le)];
      col.ri <- cols[abs(ri)];
      sub2col <- rbind(sub2col, c(i, col.le, col.ri));
      po1 <- which(dat.h$order==abs(le));
      po2 <- which(dat.h$order==abs(ri));
      r.l1  <- r2;
      r.l2  <- r2;
    } else if (le < 0){    
      col.le <- cols[abs(le)];
      col.ri <- sub2col[ri,3];
      sub2col <- rbind(sub2col, c(i, col.le, col.ri));
      po1 <- which(dat.h$order==abs(le));
      po2 <- sub2po[ri,2];
      r.l1  <- r2;
      r.l2  <- r2 - ratio * sub2h[ri,2];
    } else if (ri < 0){
      col.le <- sub2col[le,2];
      col.ri <- cols[abs(ri)];
      sub2col <- rbind(sub2col, c(i, col.le, col.ri));
      po1 <- sub2po[le, 2];
      po2 <- which(dat.h$order==abs(ri));
      r.l1  <- r2 - ratio * sub2h[le,2];
      r.l2  <- r2;
    } else {
      col.le <- sub2col[le,2];
      col.ri <- sub2col[ri,3];
      sub2col <- rbind(sub2col, c(i, col.le, col.ri));
      po1   <- sub2po[le, 2];
      po2   <- sub2po[ri, 2];
      r.l1  <- r2 - ratio * sub2h[le,2];
      r.l2  <- r2 - ratio * sub2h[ri,2];
    }
    sub2po <- rbind(sub2po, c(i, (po1+po2)/2));
    sub2h  <- rbind(sub2h,  c(i, dat.h$height[i]));
    ang1   <- angle.start + lab*po1;
    ang2   <- angle.start + lab*po2;
    r.l3   <- r2 - ratio*dat.h$height[i];
    # arc
    if (sub2col[i,2]==sub2col[i,3]){
      draw.arc.s(cx, cy, r.l3, ang1, ang2, col=col.le, lwd=lwd, lend=1);
    } else {
      col1 <- sub2col[i,2];
      col2 <- sub2col[i,3];
      ang1.5 <- (ang1+ang2)/2;
      draw.arc.s(cx, cy, r.l3, ang1.5, ang2, col=col2, lwd=lwd, lend=1);
      draw.arc.s(cx, cy, r.l3, ang1, ang1.5, col=col1, lwd=lwd, lend=1);
    }
    # left line (small angle)
    draw.line(cx, cy, ang1, r.l1, r.l3, col=col.le, lwd=lwd, lend=1);
    # right line (large angle)
    draw.line(cx, cy, ang2, r.l2, r.l3, col=col.ri, lwd=lwd, lend=1);
  }
  ## tip names
  if (tip){
    angle.s <- seq(angle.start, angle.end, length.out=max.n);
    cols    <- tip.color[dat.h$order];
    for (i in 1:length(dat.h$order)){
       n  <- dat.h$label[dat.h$order[i]];
       w3 <- angle.s[i];
       radialtext(n, center = c(cx, cy), start = r2 + ratio + 1, middle = 1, end = r2 + ratio, 
                  angle = pi/180 * (359.5 - w3), deg = NA, expand = 0, stretch = 1, nice = TRUE, 
                  cex = tip.cex, col=cols[i]) 
    }
  }
}
### end cluster.circle

### library(plotrix);
arctext <- function (x, center = c(0, 0), radius = 1, start = NA, middle = pi/2, 
                     end = NA, stretch = 1, clockwise = TRUE, cex = NA, ...) 
{
  oldcex <- par("cex")
  if (is.na(cex)) 
    cex <- oldcex
  par(cex = cex)
  xvec <- strsplit(x, "")[[1]]
  lenx <- length(xvec)
  xwidths <- stretch * strwidth(xvec)
  charangles <- xwidths/radius
  changrang <- range(charangles)
  charangles[charangles < changrang[2]/2] <- changrang[2]/2
  if (!is.na(end)) 
    start <- end + sum(charangles)
  if (is.na(start)) {
    if (clockwise) 
      start <- middle + sum(charangles)/2
    else start <- middle - sum(charangles)/2
  }
  if (clockwise) {
    charstart <- c(start, start - cumsum(charangles)[-lenx])
    charpos <- charstart - charangles/2
  }
  else {
    charstart <- c(start, start + cumsum(charangles)[-lenx])
    charpos <- charstart + charangles/2
  }
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
  for (xchar in 1:lenx) {
    srt <- 180 * charpos[xchar]/pi - 90
    text(center[1] + radius * cos(charpos[xchar]), center[2] + 
           radius * sin(charpos[xchar]) * ymult, xvec[xchar], 
         adj = c(0.5, 0.5), srt = srt + 180 * (!clockwise), 
         ...)
  }
  par(cex = oldcex)
}

### library(plotrix);
radialtext <- function (x, center = c(0, 0), start = NA, middle = 1, end = NA, 
                        angle = 0, deg = NA, expand = 0, stretch = 1, nice = TRUE, 
                        cex = NA, ...) {
  oldcex <- par("cex")
  if (is.na(cex)) 
    cex <- oldcex
  par(cex = cex)
  if (is.na(deg)) 
    deg <- angle * 180/pi
  deg <- deg%%360
  chardeg <- deg
  if (nice && (deg > 90 && deg < 270)) {
    chardeg <- (deg + 180)%%360
    x <- paste(rev(unlist(strsplit(x, ""))), collapse = "")
  }
  angle <- deg * pi/180
  xvec <- strsplit(x, "")[[1]]
  lenx <- length(xvec)
  xwidths <- stretch * strwidth(xvec)
  xwrange <- range(xwidths)
  xwidths[xwidths < xwrange[2]/2] <- xwrange[2]/2
  chexp <- rep(1, lenx)
  if (!is.na(start) && expand != 0) {
    expfactor <- expand/start
    deltar <- 0
    for (xchar in 2:lenx) {
      deltar <- deltar + xwidths[xchar - 1]
      expansion <- 1 + deltar * expfactor
      if (expansion < 0.1) 
        expansion <- 0.1
      chexp[xchar] <- expansion
      xwidths[xchar] <- xwidths[xchar] * expansion
    }
  }
  if (is.na(start)) {
    if (is.na(end)) 
      start <- middle - sum(xwidths)/2
    else start <- end - sum(xwidths)
  }
  cosang <- cos(angle)
  sinang <- sin(angle)
  charstart <- c(start, start + cumsum(xwidths)[-lenx])
  charpos <- charstart + xwidths/2
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
  for (xchar in 1:lenx) {
    par(cex = cex * chexp[xchar])
    text(center[1] + charpos[xchar] * cosang, center[2] + 
           charpos[xchar] * ymult * sinang, xvec[xchar], adj = c(0.5, 
                                                                 0.5), srt = chardeg, ...)
  }
  par(cex = oldcex)
}

