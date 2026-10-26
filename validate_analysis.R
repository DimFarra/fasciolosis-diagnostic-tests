source('run_analysis.R')
q2<-function(m,cov){c((1-m[1])*(1-m[2])+cov,m[1]*(1-m[2])-cov,(1-m[1])*m[2]-cov,prod(m)+cov)}
q3<-function(m,cv){e12<-m[1]*m[2]+cv[1];e13<-m[1]*m[3]+cv[2];e23<-m[2]*m[3]+cv[3];e123<-prod(m)+cv[1]*m[3]+cv[2]*m[2]+cv[3]*m[1];c(1-sum(m)+e12+e13+e23-e123,m[1]-e12-e13+e123,m[2]-e12-e23+e123,e12-e123,m[3]-e13-e23+e123,e13-e123,e23-e123,e123)}
bits2<-as.matrix(expand.grid(0:1,0:1));bits3<-as.matrix(expand.grid(0:1,0:1,0:1))
set.seed(441)
for(i in 1:1000){m<-runif(3,.1,.9); cv<-runif(3,-.01,.01);q<-q3(m,cv);stopifnot(abs(sum(q)-1)<1e-12,max(abs(crossprod(q,bits3)-m))<1e-12);actual<-c(sum(q*bits3[,1]*bits3[,2]),sum(q*bits3[,1]*bits3[,3]),sum(q*bits3[,2]*bits3[,3]))-c(m[1]*m[2],m[1]*m[3],m[2]*m[3]);stopifnot(max(abs(actual-cv))<1e-12)
 ci<-apply(bits3,1,function(b)prod(m^b*(1-m)^(1-b)));stopifnot(max(abs(q3(m,c(0,0,0))-ci))<1e-12)
}
# Numerical rank is a local diagnostic, not a general identifiability proof.
herdprob<-function(theta){d<-q2(theta[4:5],theta[8]);nd<-q2(1-theta[6:7],theta[9]);unlist(lapply(theta[1:3],function(p)(p*d+(1-p)*nd)[1:3]))}
z<-c(.1,.3,.7,.9,.6,.95,.85,.005,.002);jac<-sapply(seq_along(z),function(j){zp<-zm<-z;zp[j]<-zp[j]+1e-6;zm[j]<-zm[j]-1e-6;(herdprob(zp)-herdprob(zm))/2e-6});sv<-svd(jac)$d;write.csv(data.frame(singular_value=sv),'outputs/herd_jacobian_singular_values.csv',row.names=FALSE)
# Validate probability constraints on every 20th retained posterior draw and PPC.
pp<-list();checks<-list();set.seed(774)
for(f in list.files('outputs',pattern='_draws.rds$',full.names=TRUE)){
 stem<-sub('_draws.rds$','',basename(f));model<-sub('_.*','',stem);scenario<-substring(stem,nchar(model)+2);m<-as.matrix(readRDS(f));idx<-unique(round(seq(1,nrow(m),length.out=min(4000,nrow(m)))));m<-m[idx,,drop=FALSE];dep<-!grepl('independent',scenario);inp<-readRDS(sub('_draws.rds','_inputs.rds',f,fixed=TRUE));y<-inp$data$y;ns<-if(model=='Herd')inp$data$n else inp$data$n
 pvals<-matrix(NA,nrow(m),if(model=='Herd')3 else 1);minq<-1;maxerr<-0
 expcells<-array(NA,c(nrow(m),length(y)));repstats<-obsstats<-numeric(nrow(m))
 for(i in seq_len(nrow(m))){v<-m[i,];get<-function(n)if(n%in%names(v))v[n] else 0
 if(model=='Herd'){d<-q2(v[c('SeS','SeC')],get('covSe'));nd<-q2(1-v[c('SpS','SpC')],get('covSp'));prob<-t(sapply(v[paste0('p[',1:3,']')],function(p)p*d+(1-p)*nd));expected<-prob*ns;yr<-t(sapply(1:3,function(j)rmultinom(1,ns[j],prob[j,])[,1]))}
 else{d<-q3(v[c('SeL','SeC','SeS')],c(get('covSeLC'),get('covSeLS'),get('covSeCS')));nd<-q3(1-v[c('SpL','SpC','SpS')],c(get('covSpLC'),get('covSpLS'),get('covSpCS')));prob<-v['p']*d+(1-v['p'])*nd;expected<-prob*ns;yr<-rmultinom(1,ns,prob)[,1]}
 minq<-min(minq,d,nd);maxerr<-max(maxerr,abs(sum(d)-1),abs(sum(nd)-1));stopifnot(min(d,nd)>-1e-9,v['SeS']+v['SpS']>1)
 expcells[i,]<-as.vector(expected);obsstats[i]<-sum((y-expected)^2/(expected+1e-8));repstats[i]<-sum((yr-expected)^2/(expected+1e-8))
 }
 checks[[stem]]<-data.frame(model,scenario,checked_draws=nrow(m),min_conditional_probability=minq,max_normalization_error=maxerr)
 pp[[stem]]<-data.frame(model,scenario,ppp=mean(repstats>=obsstats),median_observed_discrepancy=median(obsstats))
 qs<-t(apply(expcells,2,quantile,c(.025,.5,.975)));write.csv(data.frame(cell=seq_along(y),observed=as.vector(y),lower=qs[,1],median=qs[,2],upper=qs[,3]),paste0('outputs/',stem,'_expected_cells.csv'),row.names=FALSE)
}
write.csv(do.call(rbind,checks),'outputs/probability_validation.csv',row.names=FALSE);write.csv(do.call(rbind,pp),'outputs/posterior_predictive_checks.csv',row.names=FALSE)
writeLines(c('1000 independent probability identity tests passed.','Every checked posterior draw has valid conditional cell probabilities and the stated class orientation.',paste('Herd Jacobian numerical rank at generic interior point:',sum(sv>1e-7),'of 9.')), 'outputs/validation_report.txt')
# Nominal Beta factors are jointly restricted by orientation and cell validity.
# Prior-only simulation displays this induced marginal distribution explicitly.
for(model in c('Herd','Slaughterhouse'))for(sc in c('baseline','beta22_tests','uniform_tests','uniform_all','baseline_independent','beta22_tests_independent','uniform_tests_independent','uniform_all_independent','uniform_prevalence_independent')){
 if(!sc %in% reported_scenarios(model)) next
 if(file.exists(paste0('outputs/',model,'_',sc,'_effective_prior.csv')))next
 dep<-!grepl('independent',sc)
 f<-paste0('outputs/',model,'_',sc,'_inputs.rds');inp<-readRDS(f);dat<-inp$data;dat$y[]<-0L;dat$n[]<-0L;txt<-readLines(paste0('models/',model,'_',sc,'.jags'))
 jm<-jags.model(textConnection(txt),data=dat,inits=lapply(1:4,function(j)init_for(model,dep,55330+j,j)),n.chains=4,n.adapt=5000,quiet=TRUE);update(jm,5000,progress.bar='none');mon<-if(model=='Herd')ferik_monitor(dep) else slaughter_monitor(dep);dr<-coda.samples(jm,mon,n.iter=50000,thin=5,progress.bar='none');su<-summarize(dr,model,sc);write.csv(su,paste0('outputs/',model,'_',sc,'_effective_prior.csv'),row.names=FALSE)
}


