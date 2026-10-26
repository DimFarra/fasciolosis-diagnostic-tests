source('run_analysis.R')
config <- blcm_config()
files <- unlist(lapply(config$models, function(model)
  file.path('outputs', paste0(model, '_', reported_scenarios(model, config), '_summary.csv'))))
stopifnot(all(file.exists(files)))
all <- do.call(rbind, lapply(files, read.csv))
write.csv(all, 'outputs/comparison_all.csv', row.names=FALSE)
required <- rbind(
  data.frame(model='Herd', scenario=blcm_config()$primary_scenarios),
  data.frame(model='Slaughterhouse', scenario=c(blcm_config()$primary_scenarios, 'liver_sp90_independent'))
)
observed <- unique(all[c('model','scenario')])
stopifnot(all(apply(required, 1, function(z) any(observed$model==z[['model']] & observed$scenario==z[['scenario']]))), all(all$passes))

diag <- aggregate(cbind(rhat,ess,mc_pct_sd,iterations_per_chain)~model+scenario, all,
                  function(x)c(min=min(x),max=max(x)))
write.csv(diag, 'outputs/diagnostics_all.csv', row.names=FALSE)

# Record the actual sequentially transferred test priors used by the slaughterhouse baseline.
tp <- herd_posterior_test_prior('outputs')
transfer_map <- c(se_sero='SeS', sp_sero='SpS', se_copro='SeC', sp_copro='SpC')
transfer <- do.call(rbind, lapply(names(tp), function(n) {
  a <- tp[[n]][1]; b <- tp[[n]][2]
  data.frame(
    slaughter_parameter=n,
    herd_parameter=transfer_map[[n]],
    beta_alpha=a, beta_beta=b,
    mean=a/(a+b),
    sd=sqrt(a*b/((a+b)^2*(a+b+1)))
  )
}))
write.csv(transfer, 'outputs/transferred_test_priors.csv', row.names=FALSE)

pri <- list()
for(mod in c('Herd','Slaughterhouse')) for(sc in reported_scenarios(mod)) {
  if(mod=='Herd' && grepl('liver',sc)) next
  pr <- make_prior(mod, sc, 'outputs')
  for(n in names(pr)) {
    a <- pr[[n]][1]; b <- pr[[n]][2]
    pri[[length(pri)+1]] <- data.frame(
      model=mod, scenario=sc, parameter=n, a, b,
      mean=a/(a+b),
      sd=sqrt(a*b/((a+b)^2*(a+b+1))),
      median=qbeta(.5,a,b), lower=qbeta(.025,a,b), upper=qbeta(.975,a,b)
    )
  }
}
write.csv(do.call(rbind,pri), 'outputs/nominal_priors_all.csv', row.names=FALSE)

write.csv(data.frame(group=rownames(ferik_counts),n=ferik_n,ferik_counts),'outputs/herd_counts.csv',row.names=FALSE)
write.csv(data.frame(pattern=names(slaughter_counts),count=slaughter_counts),'outputs/slaughter_counts.csv',row.names=FALSE)

dir.create('figures',showWarnings=FALSE)
labels <- c(SeS='ELISA sensitivity',SpS='ELISA specificity',SeC='Coprology sensitivity',SpC='Coprology specificity',SeL='Liver sensitivity',SpL='Liver specificity')
plot_sensitivity_original_style <- function(mod) {
  pars <- if(mod=='Herd') names(labels)[1:4] else names(labels)
  scens <- c('baseline_independent','widened_independent','beta22_tests_independent','uniform_tests_independent','uniform_all_independent','uniform_prevalence_independent')
  scLabs <- c('Baseline','Wider priors','Beta(2,2) tests','Uniform tests','Uniform all','Uniform prevalence')
  available <- scens %in% unique(all$scenario[all$model == mod])
  scens <- scens[available]; scLabs <- scLabs[available]
  cols <- c('#1B9E77','#D95F02','#377EB8','#984EA3','#FF7F00','#A65628')
  pchs <- c(21,22,23,24,25,8); offsets <- seq(-.30,.30,length.out=length(scens))
  d <- all[all$model==mod & all$parameter%in%pars,]
  png(paste0('figures/',mod,'_sensitivity.png'),width=2400,height=1550,res=240)
  par(mar=c(8.2,5.4,4.8,1.2),family='sans',cex=.92)
  plot(NA,xlim=c(.45,length(pars)+.55),ylim=c(0,100),xaxt='n',xlab='',ylab='Posterior median (95% CrI), %',bty='l',
       main=if(mod=='Herd') 'Herd model: sensitivity analysis' else 'Slaughterhouse model: sensitivity analysis')
  abline(h=seq(0,100,25),col='#E2E2E2',lty=1)
  abline(v=seq(1.5,length(pars)-.5,1),col='#F0F0F0',lty=1)
  for(i in seq_along(scens)) for(j in seq_along(pars)) {
    v <- d[d$scenario==scens[i] & d$parameter==pars[j],]
    stopifnot(nrow(v)==1)
    x <- j+offsets[i]
    segments(x,100*v$lower,x,100*v$upper,lwd=1.8,col='#303030')
    points(x,100*v$median,pch=pchs[i],bg=cols[i],col='black',cex=1.15)
  }
  text(seq_along(pars),-7,labels[pars],srt=38,adj=1,xpd=NA,cex=.9)
  legend('top',legend=scLabs,col='black',pt.bg=cols,pch=pchs,horiz=TRUE,bty='n',xpd=NA,inset=c(0,-.02),cex=.78)
  mtext('All displayed prior-sensitivity scenarios use conditional independence',side=1,line=6.8,cex=.78)
  dev.off()
}
for(mod in c('Herd','Slaughterhouse')) plot_sensitivity_original_style(mod)
for(mod in c('Herd','Slaughterhouse')) {
  draw_file <- paste0('outputs/',mod,'_baseline_independent_draws.rds')
  if (file.exists(draw_file)) {
    d <- readRDS(draw_file)
    png(paste0('figures/',mod,'_main_diagnostics.png'),width=2400,height=2500,res=200)
    plot(d,auto.layout=TRUE); dev.off()
  }
}


