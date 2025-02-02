library(VIM)
library(GGally) 
library(ggplot2)

data=read.csv(paste(SPATH, "/Data/source/DMSdata.csv", sep = ""))
colnames(data)=c("Ocean", "Water_Density", "DMHg", "Hg0", "MMHg", "THg")
data$Ocean=factor(data$Ocean)
p=ncol(data)
Y.knn = kNN(data, k = 5)[,1:p]

yy=function(data,mapping,...){
  ggplot(data=Y.knn,mapping=mapping)+ 
    #geom_histogram(position = "identity",alpha=0.5,aes(y = ..density..,fill = factor(V8))) +
    stat_density(geom = "line",position = "identity", aes(colour = factor(Ocean)))
  #geom_density(colour="black")
}
zz=function(data,mapping,...){
  ggplot(data=Y.knn,mapping=mapping)+
    geom_point(aes(shape = factor(Ocean), color = factor(Ocean)))+
    scale_shape_manual(values=c("Current"=16,"Undercurrent"=2))  
  #geom_smooth(method = lm)+
  #geom_smooth(aes(group = 1), method="lm", size = 0.5, se = F,colour="black")
}

cor(data[,-1])
Undercurrent=Y.knn[which(Y.knn$Ocean=="Undercurrent"),]
Undercurrent.cor=cor(Undercurrent[,-1])
Current=Y.knn[which(Y.knn$Ocean=="Current"),]
Current.cor=cor(Current[,-1])

Fig.S1 <- ggpairs(Y.knn, columns = 2:6,
                  aes(colour=factor(Ocean), shape=factor(Ocean)),  #加shape參數
                  upper = list(continuous = wrap("cor", size = 3.2)),
                  lower = list(continuous =zz),diag = list(continuous =yy),title="", 
                  axisLabels = "show", legend = c(1, 1)) +
  theme(text = element_text(size=14, lineheight = 30),
        plot.margin = unit(c(0,0,0,0), "cm"),
        strip.text = element_text(size = 15, margin = margin(r = 10, l = 10, b = 10, t = 10))) +
  # theme_light(10) + 
  theme(legend.position="none")
ggsave(paste(SPATH, '/Results/Fig.S1.eps', sep=''), plot = Fig.S1, width=8, height=8)
