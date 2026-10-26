suppressPackageStartupMessages(library(coda))
dir.create('figures/diagnostics',recursive=TRUE,showWarnings=FALSE)
for(f in list.files('outputs',pattern='_draws.rds$',full.names=TRUE)){
 d<-readRDS(f);vars<-varnames(d);stem<-sub('_draws.rds','',basename(f),fixed=TRUE)
 for(page in seq_len(ceiling(length(vars)/7))){v<-vars[seq((page-1)*7+1,min(page*7,length(vars)))];png(paste0('figures/diagnostics/',stem,'_',page,'.png'),width=2000,height=300*length(v),res=170);par(mfrow=c(length(v),2),mar=c(2.2,3.1,1.3,.5),cex=.8)
 for(p in v){m<-sapply(d,function(z)as.matrix(z)[,p]);ix<-unique(round(seq(1,nrow(m),length.out=min(1000,nrow(m)))));matplot(ix,m[ix,,drop=FALSE],type='l',lty=1,col=c('#225EA8','#D95F02','#238B45','#756BB1'),ylab=p,xlab='Retained draw',main=p,cex.main=.9);ds<-lapply(1:4,function(j)density(m[,j]));plot(ds[[1]],main='',xlim=range(m),ylim=c(0,max(sapply(ds,function(z)max(z$y)))),xlab='Parameter value',col='#225EA8');for(j in 2:4)lines(ds[[j]],col=c('#225EA8','#D95F02','#238B45','#756BB1')[j])};dev.off()}
}

