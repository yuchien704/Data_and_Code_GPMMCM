data=read.csv(paste(SPATH, "/Data/source/DMSdata.csv", sep = ""))
true.clus = as.numeric(as.factor(data[,1]))
Y = data[,-1]
g = list()
mean = sd = matrix(NA, 5,2)
for(i in 1:5){
  for(j in 1:2){
    g[[j]] = which(true.clus==j)
    mean[i,j] = mean(Y[g[[j]],i],na.rm = T)
    sd[i,j] = sd(na.omit(Y[g[[j]],i]))
  }
}

n = nrow(data)
cen.rate = round(rbind(0,sum(data[,3]<=2)/n, sum(data[,4]<=4)/n,sum(data[,5]<=11.3)/n,sum(data[,6]<=0.22)/n),4)*100
na.rate = round(colSums(is.na(data[,-1]))/n,4)*100

Table3 = cbind(cen.rate, na.rate, 
               round(mean,2)[,1], round(sd,2)[,1], 
               round(mean,2)[,2], round(sd,2)[,2])
colnames(Table3) = c("Censoring proportion", "Missing proportion", 
                     "California Current-mean","California Current-sd",
                     "California Undercurrent-mean","California Undercurrent-sd")
write.csv(Table3, paste(SPATH,'/Results/Table3.csv',sep=""), row.names = TRUE)
