# Check Beta prior shapes against their stated matching targets.
suppressPackageStartupMessages(library(PriorGen))
x<-data.frame(parameter=c('Arab','Goran','Fulani','SeS','SpS','SeC','SpC','Slaughter'),target_type=c('median','median',rep('mode',6)),target=c(.25,.01,.47,.97,.93,.60,.80,.23),lower_bound=c(.1,.001,.3,.9,.85,.3,.6,.1),prob_above=.95)
shapes<-t(sapply(seq_len(nrow(x)),function(i){args<-list(percentile=x$prob_above[i],lower.v=FALSE,percentile.value=x$lower_bound[i]);args[[paste0('the',x$target_type[i])]]<-x$target[i];do.call(findbeta,args)$parameters}))
write.csv(cbind(x,shapes),'outputs/reconstructed_elicitation.csv',row.names=FALSE)

