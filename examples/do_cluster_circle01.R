rm(list=ls());

library(OmicCircos2);

source("make_heatmap_dat.R");
source("cluster_circle.R");
source("heatmap_cluster_h.R");
set.seed(1234);
numc  <- 80;
dat   <- make.heatmap.dat(cnum=numc, rnum=10, clu.num=3, clu.mean=c(0.8,1.2,2.4));
colnames(dat)  <- paste0("name", 1:ncol(dat));
row.names(dat) <- paste0("gene", 1:nrow(dat));

dat.d <- dist(t(dat));
dat.h <- hclust(dat.d);

order.n <- dat.h$labels[dat.h$order];
col.i   <- match(order.n, colnames(dat));
dat     <- dat[,col.i];

seg.End   <- 1:ncol(dat);
seg.Start <- seg.End - 1;
nodat     <- rep("NA", ncol(dat));

seg.d      <- data.frame(seg.name=rep("seg1",ncol(dat)), seg.Start=seg.Start, 
                    seg.End=seg.End, the.v=nodat, NO=nodat);
cir        <- segAnglePo(seg.d, seg="seg1", angle.start=90, angle.end=360);
heatmap.d  <- data.frame(chr=rep("seg1", ncol(dat)), po=seg.End, t(dat));
lab.d      <- data.frame(heatmap.d[,1:2], names=dat.h$labels[dat.h$order]);

cols <- rainbow(10, alpha=0.8);
col1 <- c(rep(cols[1],numc),rep(cols[7],numc),rep(cols[9],numc));

pdf("do_cluster_circle01.pdf", 6, 6);
par(mar=c(2,2,2,2));
plot(c(1:800), c(1:800), type="n", axes=F, xlab="", ylab="", main="", asp=1);
cluster.circle(400, 400, 20, 200, data=dat.h, angle.start=90,  angle.end=358, tip=T, tip.cex=0.01,
               edge.color = col1, tip.color=col1, edge.width = 1, lab.d=1, cex=1);
circos(mapping=heatmap.d, xc=400, yc=400, R=210, W=150,
       cir=cir, type="heatmap2", lwd=0.2)
circos(mapping=lab.d, xc=400, yc=400, R=350, W=10,
       cir=cir, type="label2", side="out", cex=0.3, col="black")

y <- 800;
x <- 550;
for (i in 1:20){
  t.s <- paste0("annotation x=", x, " y=", y);
  text(x, y, t.s, cex=0.6);
  y <- y - 20;
}
dev.off();

