if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mfuzz")
library("Mfuzz")
df3=read.csv("/Users/huangminjing/Desktop/serum-4group-for-mfuzz.csv",header = T,row.names = 1)
df3a<-as.matrix(df3)
df3Ex<- ExpressionSet(assayData = df3a)
df3F <- filter.NA(df3Ex,thres = 0.25)
df3F <- fill.NA(df3F,mode = 'mean')
df3F <- standardise(df3F)
set.seed(2021)
cl <- mfuzz(df3F,c=8,m=1.25)
pdf("mfuzz.pdf")
mfuzz.plot2(df3F, cl=cl,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()
dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:8){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}
data=read.csv("/Users/huangminjing/Desktop/serum-BC-0.00001-for-depfc2&0.5p0.05-mfuzz.csv",header = T,row.names = 1)
data <- as.matrix(data)
mfuzz_class <- new('ExpressionSet',exprs = data)
gene.s <- standardise(data)
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- standardise(mfuzz_class)
set.seed(123)
cluster_num <- 6
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
mfuzz_cluster$size
mfuzz.plot2(mfuzz_class,cl=mfuzz_cluster,mfrow=c(2,3),ylim.set=c(0,0),centre=T,x11=F,centre.lwd=0.2,Xwidth=1,Xheight=1,new.window= FALSE)
dir.create(path="/Users/Desktop/mfuzz",recursive = TRUE)
for(i in 1:5){
  potname<-names(mfuzz_cluster$cluster[unname(mfuzz_cluster$cluster)==i])
  write.csv(mfuzz_cluster[[4]][potname,i],paste0("/Users/Desktop/mfuzz","/mfuzz_",i,".csv"))
}

mfuzz.plot2(mfuzz_class,cl=mfuzz_cluster,centre=T,x11=F,colo = "fancy", centre.lwd=,single=2,new.window= FALSE)
protein_cluster <- mfuzz_cluster$cluster
protein_cluster <- cbind(data[names(protein_cluster), ], protein_cluster)
write.csv(protein_cluster, '/Users/Desktop/res.csv')
protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)
write.csv(protein_standard_cluster, '/Users/Desktop/res-standerlize.csv')