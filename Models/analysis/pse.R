setwd(paste0(Sys.getenv('CS_HOME'),'/NetworkGrowth/Models/analysis'))

library(ggplot2)
library(dplyr)
library(GGally)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

resdir = paste0(Sys.getenv('CS_HOME'),'/NetworkGrowth/Results/PSE/');dir.create(resdir,recursive = T)

resdirs = list(
  'random' = '../netlogo/pse/20211026_1258_PSE_RANDOM_LOCAL/',
  'connection' = '../netlogo/pse/20211103_2144_PSE_CONNECTION_LOCAL/',
  'gravity' = '../netlogo/pse/20211028_1101_PSE_GRAVITY_LOCAL/',
  'breakdown' = '../netlogo/pse/20211028_0804_PSE_BREAKDOWN_LOCAL/',
  'cost' = '../netlogo/pse/20211030_1731_PSE_COST_LOCAL',
  'biological' = '../netlogo/pse/20211030_1727_PSE_BIOLOGICAL_LOCAL/'
)

gen<-function(filename){as.integer(sapply(strsplit(sapply(strsplit(filename,"population"),function(s){s[2]}),".csv"),function(s){s[1]}))}

latestgen <- function(dir){
  if(is.null(dir)){return(NULL)}
  else{return(max(sapply(list.files(dir,pattern=".csv"),gen)))}
}

indics = c("meanBwCentrality","meanClosenessCentrality","meanPathLength","meanRelativeSpeed","nwDiameter")#,"nwLength")

res = data.frame()
for(model in names(resdirs)){
  currentdir = resdirs[[model]]
  currentfiles = paste0(currentdir,'/',list.files(currentdir,pattern = '.csv'))
  for(currentfile in currentfiles){
    currentres = as.tbl(read.csv(currentfile));currentgen=gen(currentfile)
    res = rbind(res,cbind(currentres[,indics],model=rep(model,nrow(currentres)),generation=rep(currentgen,nrow(currentres))))
  }
}

# some cleaning
res=res[!is.na(res$meanBwCentrality),]
res=res[res$meanClosenessCentrality<0.2,]



ggsave(plot = ggpairs(res,aes(color=model),columns = indics,
                      lower = list(continuous = wrap("points", alpha = 0.6,size=0.05,shape='.')),
                      diag = list(continuous = wrap("densityDiag", alpha = 0.4))
)+stdtheme,filename = paste0(resdir,'scatter_pse_allmodels.png'),width = 40,height=30,units='cm')


###
# with real data
real = read.csv('data/real.csv')
# renorm - these were not in the old computation (for mmost distance based indics, no norm by world diag -> must be rescaled here)
# ! need to recompute with correct formulas - still try to look what the plot looks like
real$meanPathLength = ((real$meanPathLength - mean(real$meanPathLength))/sd(real$meanPathLength))*sd(res$meanPathLength) + mean(res$meanPathLength)
real = real[real$meanPathLength<max(res$meanPathLength),]
real = real[real$meanRelativeSpeed<quantile(real$meanRelativeSpeed,c(0.8)),]
real$meanRelativeSpeed = ((real$meanRelativeSpeed - mean(real$meanRelativeSpeed))/sd(real$meanRelativeSpeed))*sd(res$meanRelativeSpeed) + mean(res$meanRelativeSpeed)
real$nwDiameter = ((real$nwDiameter - min(real$nwDiameter))/sd(real$nwDiameter))*sd(res$nwDiameter) + min(res$nwDiameter)
real=real[real$nwDiameter<2.7,]
real = real[real$meanClosenessCentrality<0.2,]


res = rbind(res,cbind(real[,indics],generation=rep(0,nrow(real)),model=rep('real',nrow(real))))

ggsave(plot = ggpairs(res,aes(color=model),columns = indics,
                      lower = list(continuous = wrap("points", alpha = 0.6,size=0.05,shape='.')),
                      diag = list(continuous = wrap("densityDiag", alpha = 0.4))
)+stdtheme,filename = paste0(resdir,'scatter_pse_withreal.png'),width = 40,height=30,units='cm')
# this is crap - do not show this



######
## hypervolumes intersections

library(hypervolume)
library(reshape2)

models = names(resdirs)

set.seed(42)

hvs = list()
for(model in models){
  show(model)
  start = as.numeric(Sys.time())
  sample = res[res$model==model,indics]
  sample = sample[sample.int(nrow(sample),size=min(1000,nrow(sample)),replace = F),]
  hvs[[model]] = hypervolume_gaussian(sample)
  show(as.numeric(Sys.time()) - start)
}
sapply(hvs,function(hv){hv@Volume})

overlaps = matrix(data=rep(1,length(models)*length(models)),nrow=length(models))
for(model1 in 1:length(models)){
  for (model2 in 1:length(models)){
    if(model1!=model2){
      show(paste0(models[model1],' / ',models[model2]))
      hvset = hypervolume_set(hvs[[models[model1]]],hvs[[models[model2]]],check.memory=FALSE,num.points.max=100000)
      overlaps[model1,model2] = hvset@HVList$Intersection@Volume / hvs[[models[model2]]]@Volume
    }
  }
}
rownames(overlaps)<-models
colnames(overlaps)<-models
diag(overlaps)=NA
save(overlaps,hvs,file=paste0(resdir,'hvs.RData'))

g=ggplot(melt(overlaps),aes(x=Var1,y=Var2,fill=value))
ggsave(g+geom_raster()+xlab('')+ylab('')+scale_fill_continuous(name='Overlap')+stdtheme,
       file=paste0(resdir,'pointclouds-overlap.png'),width=30,height=26,units='cm')





