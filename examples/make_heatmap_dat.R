

make.heatmap.dat <- function(cnum=10, rnum=10, clu.num=3, clu.mean=c(0,3,6)) {

  ## matrix
  mnum  <- cnum*rnum;
  
  out <- list();
  for (i in 1:clu.num){
    dat  <- matrix(rnorm(mnum, mean=clu.mean[i]), ncol=cnum);
    out[[i]] <- dat;
  }

  dat.s   <- c();
  for (i in 1:3){
    dat.s  <- rbind(dat.s, cbind(out[[i]], out[[i]], out[[i]]));
  }

  library(combinat);
  perm.i <- permn(c(1:clu.num));
  for (i in 1:length(perm.i)){
    perm.j <- perm.i[[i]]
    dat.s  <- rbind(dat.s, cbind(out[[perm.j[1]]], out[[perm.j[2]]], out[[perm.j[3]]]));
  }
  dat.s
}

