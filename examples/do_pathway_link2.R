rm(list=ls());

options(stringsAsFactors = F);
library(OmicCircos2);
library(OmicPath);

data(msigdb.h)
pathway_table <- GSEA.db$pathway.table;
pathway_table <- pathway_table[c(1,2,9,12,13),];
pathway_num   <- nrow(pathway_table);

pathway_gene <- NULL;
segment_data <- NULL;
for (i in 1:pathway_num){
  gene  <- unlist(strsplit(pathway_table[i,2], " "));
  seg.s <- paste0("seg", i);
  df    <- data.frame(pathway.name=rep(seg.s, length(gene)), gene=gene, gene.id=1:length(gene));
  pathway_gene <- rbind(pathway_gene, df);
  segment_data <- rbind(segment_data, c(seg.s, 1,  length(gene), pathway_table[i,3], "HALLMARK"))
}
colnames(segment_data) <- c("seg.name", "seg.Start", "seg.End", "the.v", "note");
cir <- segAnglePo(seg.dat=segment_data, seg=segment_data[,1], angle.start=1, angle.end=360);

ppi_file <- system.file("extdata", "PPI_hs.txt", package="OmicPath");
ppi      <- read.table(ppi_file);
ppi      <- ppi[!duplicated(ppi[,c(1,2)]),]
ppi.i    <- which(ppi[,1]==ppi[,2]);
ppi      <- ppi[-ppi.i,];

link_data <- NULL;
for (i in 1:nrow(pathway_gene)){

  g1.i <- which(ppi[,1]==pathway_gene[i,2]);
  g2.i <- which(ppi[,2]==pathway_gene[i,2]);
  g.s  <- NULL;
  if (length(g1.i) > 0){
    g.s <- c(g.s, ppi[g1.i,2]);  
  }
  if (length(g2.i) > 0){
    g.s <- c(g.s, ppi[g2.i,1]);  
  }
  if (is.null(g.s)){
    next;
  } else {
    g.s <- unique(g.s);
  }
  g.p.i <- which(pathway_gene[,2] %in% g.s);
  if (length(g.p.i)==0){
    next;
  }
  tmp.dat   <- data.frame(seg1=rep(pathway_gene[i,1], length(g.p.i)),
                          po1=rep(pathway_gene[i,3], length(g.p.i)),
                          gene1=rep(pathway_gene[i,2], length(g.p.i)),
                          seg2=pathway_gene[g.p.i,1],
                          po2=pathway_gene[g.p.i,3],
                          gene2=pathway_gene[g.p.i,2],
                          note=rep("HALLMARK", length(g.p.i))
                          )
  link_data <- rbind(link_data, tmp.dat)
}

cir.col <- rainbow(pathway_num, alpha=0.6);

pdf("do_pathway_link2_0.pdf", 8, 8);
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="", asp=1)
circos(R=345, cir=cir, mapping=cir, type="chr.noB", W=10, col=cir.col);
circos(R=330,  cir=cir, W=100, mapping=link_data, type="link", lwd=1);

dev.off();

link500 <- link_data[sample(1:nrow(link_data), 500),];
p.k     <-  gsub("seg", "", link500[,1]);
col500  <- cir.col[as.numeric(p.k)];

pdf("do_pathway_link2_1.pdf", 8, 8);
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="", asp=1)
circos(R=345, cir=cir, mapping=cir, type="chr.noB", W=10, col=cir.col);
circos(R=330,  cir=cir, W=50, mapping=link500, type="link", lwd=1, col=col500);

dev.off();

y <- 8;
pdf("do_pathway_link2_2.pdf", 8, 8);
plot(c(1,10), c(1,10), type="n", axes=FALSE, xlab="", ylab="", main="", asp=1)
for (i in 1:nrow(pathway_table)){
  points(2, y, pch=15, col=cir.col[i], cex=2)
  text(3, y, pathway_table[i,1], adj=0)
  y <- y - 0.4;
}
dev.off();


