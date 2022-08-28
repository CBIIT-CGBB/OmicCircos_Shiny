### heatmap.cluster
### the options should be the same as plot.phylo {ape}
heatmap.cluster.h <- function(x1, y1, x2, y2, data=data, 
                            edge.color = NULL, tip.color=NULL, edge.width = 1, lab.d=1, cex=1){

    dat.h  <- data;
    cols   <- edge.color;
    cols   <- tip.color;
    lwd    <- edge.width;
    ####################
    c.x1   <- x1;
    c.x2   <- x2;
    c.y1   <- y1;
    c.y2   <- y2;
    max.h  <- max(dat.h$height);
    max.n  <- length(dat.h$labels);

    lab   = (c.x2-c.x1)/max.n;
    ratio = (c.y2-c.y1)/max.h;

    ####################  
    sub2h   <- c();
    sub2po  <- c();
    sub2col <- c();
    y0 <- c.y1;
    for (i in 1:nrow(dat.h$merge)){
      le <- dat.h$merge[i,1];
      ri <- dat.h$merge[i,2];
      if (le < 0 && ri < 0){
        col.le <- cols[abs(le)];
        col.ri <- cols[abs(ri)];
        sub2col <- rbind(sub2col, c(i, col.le, col.ri));
        po1 <- which(dat.h$order==abs(le));
        po2 <- which(dat.h$order==abs(ri));
        y1  <- y0;
        y2  <- y0;
      } else if (le < 0){ 
        col.le <- cols[abs(le)];
        col.ri <- sub2col[ri,3];
        sub2col <- rbind(sub2col, c(i, col.le, col.ri));
        po1 <- which(dat.h$order==abs(le));
	      po2 <- sub2po[ri,2];
	      y1  <- y0;
	      y2  <- c.y1 + ratio * sub2h[ri,2];
      } else if (ri < 0){
        col.le <- sub2col[le,2];
        col.ri <- cols[abs(ri)];
        sub2col <- rbind(sub2col, c(i, col.le, col.ri));
	      po1 <- sub2po[le, 2];
	      po2 <- which(dat.h$order==abs(ri));
	      y1  <- c.y1 + ratio * sub2h[le,2];
	      y2  <- y0;
      } else {
        col.le <- sub2col[le,2];
        col.ri <- sub2col[ri,3];
        sub2col <- rbind(sub2col, c(i, col.le, col.ri));
	      po1 = sub2po[le, 2];
	      po2 = sub2po[ri, 2];
	      y1  = c.y1 + ratio * sub2h[le,2];
	      y2  = c.y1 + ratio * sub2h[ri,2];
      }
      sub2po <- rbind(sub2po, c(i, (po1+po2)/2));
      sub2h  <- rbind(sub2h,  c(i, dat.h$height[i]));
      x1 = c.x1 + lab*po1;
      x2 = c.x1 + lab*po2;
      y3 = c.y1 + ratio*dat.h$height[i];
      # horizontal line
      if (sub2col[i,2]==sub2col[i,3]){
        segments(x1,y3,x2,y3, lwd=lwd, col=col.le);
      } else {
        col1 <- sub2col[i,2];
        col2 <- sub2col[i,3];
        x1.5 <- (x1+x2)/2;
        segments(x1,y3,x1.5,y3, lwd=lwd, col=col1);
        segments(x1.5,y3,x2,y3, lwd=lwd, col=col2);
      }
      # left vertical line
      segments(x1,y3,x1,y1, lwd=lwd, col=col.le);
      # right vertical line
      segments(x2,y3,x2,y2, lwd=lwd, col=col.ri);
    }
    ## tip names
    #y0 <- c.y1 - lab.d;
    #for (i in 1:length(dat.h$order)){
    #   n  <- dat.h$label[dat.h$order[i]];
    #   x0 <- c.x1 + lab * i;
    #   text(x0, y0, n, srt=90, col=cols[dat.h$order[i]], cex=cex, adj=1, offset=0);
    #}
    ####################
}
### end heatmap.cluster
