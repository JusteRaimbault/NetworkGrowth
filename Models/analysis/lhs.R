
setwd(paste0(Sys.getenv('CS_HOME'),'/NetworkGrowth/Models/analysis'))

library(ggplot2)
library(dplyr)
library(GGally)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

resdir = paste0(Sys.getenv('CS_HOME'),'/NetworkGrowth/Results/LHS/');dir.create(resdir,recursive = T)

res <- as.tbl(read.csv(paste0(Sys.getenv('CN_HOME'),'/Models/MesoCoevol/MesoCoevol/exploration/2017_08_18_08_05_31_NETWORK_LHS_GRID.csv')))

res$heuristic = floor(res$nwHeuristic)
res$meanBwCentrality<-as.numeric(as.character(res$meanBwCentrality))
res$meanPathLength<-as.numeric(as.character(res$meanPathLength))
res$nwDiameter<-as.numeric(as.character(res$nwDiameter))
res$meanRelativeSpeed<-as.numeric(as.character(res$meanRelativeSpeed))
res$meanClosenessCentrality<-as.numeric(as.character(res$meanClosenessCentrality))
res$nwLength<-as.numeric(as.character(res$nwLength))
res=res[!is.na(res$meanBwCentrality),]
res=res[res$meanClosenessCentrality<0.2,]

heuristics = c("random","connexion","det-brkdn","rnd-brkdn","cost","biological")
res$heuristic = heuristics[res$heuristic+1]

indics = c("meanBwCentrality","meanClosenessCentrality","meanPathLength","meanRelativeSpeed","nwDiameter","nwLength")


####
# Scatterplot

ggsave(plot = ggpairs(res,aes(color=heuristic),columns = indics,
                      lower = list(continuous = wrap("points", alpha = 0.6,size=0.05,shape='.')),
                      diag = list(continuous = wrap("densityDiag", alpha = 0.4))
)+stdtheme,filename = paste0(resdir,'scatter_models.png'),width = 40,height=30,units='cm')


######
## hypervolumes intersections

library(hypervolume)
library(reshape2)

models = heuristics

set.seed(42)

hvs = list()
for(model in models){
  show(model)
  start = as.numeric(Sys.time())
  sample = res[res$heuristic==model,indics]
  sample = sample[sample.int(nrow(sample),size=1000,replace = F),]
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
g+geom_raster()+xlab('')+ylab('')+scale_fill_continuous(name='Overlap')+stdtheme
ggsave(file=paste0(resdir,'pointclouds-overlap.png'),width=30,height=26,units='cm')







