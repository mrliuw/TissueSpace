setwd("/tissue-space")
options(stringsAsFactors = FALSE)
library(lsa)
data=read.table("GTEx_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",header=T,skip=2,sep = "\t")
data0=data[,-c(1,2)]
aa=rowMeans(data0)
data=data[-which(aa==0),]
data0=data[,-c(1,2)]
aa=apply(data0,1,max)
data=data[-which(aa<10),]
data0=data[,-c(1,2)]
aa=apply(data0,1,sum)
data=data[-which(aa<30),]
#transform to rank,small value small rank
geneid=data[,1]
#geneid=noquote(unlist(strsplit(geneid,"\\."))[seq(1,2*length(geneid),2)])
geneid=unlist(strsplit(geneid,"\\."))[seq(1,2*length(geneid),2)]
genename=data[,2]
data=as.matrix(apply(data[,-c(1,2)],2,rank, ties.method= "max"))
rownames(data)=geneid
dim(data)
##chosing best k
space=lsa(data,dims=dimcalc_share())
D <- diag(space$sk)
#svd analysis,calculate different k and cor,MAE
MAE=0;correlation=0
for (i in c(1:6)){
  pre_rank=space$tk[,1:(50*i)] %*% D[1:(50*i),1:(50*i)] %*% t(space$dk[,1:(50*i)]) 
  MAE[i]<-sum(abs(data-pre_rank))/dim(data)[1]/dim(data)[2]
  correl=lapply(1:dim(data)[2],function(x) cor(data[,x],pre_rank[,x],method="spearman"))
  correlation[i]=paste(mean(unlist(correl)),sd(unlist(correl)))
  print(c(50*i,MAE[i],correlation[i]))}
plot(data.frame(x=seq(50,300,50),y=MAE),ylab="Mean absolute error",xlab="Indices of singular value",
     main="Fold in mean absolute error at different k")

space=lsa(data,dim=100)
D <- diag(space$sk)
rownames(space$tk)=geneid
genevec=space$tk %*% D;write.csv(genevec,"genevec.csv")
tissuevec=D %*% t(space$dk);write.csv(tissuevec,"tissuevec.csv")
saveRDS(space, "space_file")
#space=readRDS("space_file");geneid=rownames(space$tk)
#feature annotation 
library(gProfileR)

for (i in 1:100){
  genes=geneid[order(genevec[,i],decreasing = T)]
  go=gprofiler(as.vector(genes[1:30]), 
               organism = "hsapiens",numeric_ns="")
  write.table(go,"genevec_decoding.csv",append =T,row.names=rep(i,nrow(go)),sep=",")}
for (i in 1:100){
  write.table(which.max(tissuevec[i,]),"tissuevec_max_feature.csv",append =T,
              col.names=F,sep=",")}
for (i in 1:17382){
  write.table(which.max(tissuevec[,i]),"max_tissuevec_feature.csv",append =T,
              col.names=F,row.names=TRUE,sep=",")}
#GTEx fold in and classification
gtex=read.csv("GSE45878.csv",header=T,row.names=1)
#aa=union(gtex[,1],setdiff(geneid,gtex[,1]))
aa=gtex[rownames(gtex) %in% intersect(rownames(data),rownames(gtex)),]
bb=data.frame(cbind(id=rownames(aa),aa))
cc=data.frame(cbind(id=geneid))
new_gtex=dplyr::full_join(bb,cc,by="id") #merge is slow!gtex first
gtex=new_gtex[,-1]
gtex=as.matrix(apply(gtex,2,rank, ties.method= "max",na.last="keep"))
rownames(gtex)=new_gtex[,1]
gtex[is.na(gtex)]=0
gtex=gtex[rownames(data),] #reorder as data space
index=0
for (i in 1:837){
semantic=gtex[,i] %*% space$tk
cos=apply(tissuevec,2,cosine,as.vector(semantic))
index[i]=names(which.max(cos))}
write.csv(index,"validation_index.csv",row.names=colnames(gtex))
#foldin=(data[,2]) %*% space$tk %*% t(space$tk)
#intersect project in correlate compare
expo=read.csv("GSE2109.csv",header=T,row.names=1)
aa=expo[rownames(expo) %in% intersect(rownames(data),rownames(expo)),]
bb=data.frame(cbind(id=rownames(aa),aa))
cc=data.frame(cbind(id=geneid))
new_expo=dplyr::full_join(bb,cc,by="id") #merge is slow!gtex first
expo=new_expo[,-1]
expo=as.matrix(apply(expo,2,rank, ties.method= "max",na.last="keep"))
rownames(expo)=new_expo[,1]
expo[is.na(expo)]=0
expo=expo[rownames(data),] #reorder as data space
index=0
for (i in 1:2158){
  semantic=expo[,i] %*% space$tk
  cos=apply(tissuevec,2,cosine,as.vector(semantic))
  index[i]=names(which.max(cos))}
write.csv(index,"validation_index_expo.csv",row.names=colnames(expo))
###GSE36376
expo=read.csv("GSE36376.csv",header=T,row.names=1)
aa=expo[rownames(expo) %in% intersect(rownames(data),rownames(expo)),]
bb=data.frame(cbind(id=rownames(aa),aa))
cc=data.frame(cbind(id=geneid))
new_expo=dplyr::full_join(bb,cc,by="id") #merge is slow!gtex first
expo=new_expo[,-1]
expo=as.matrix(apply(expo,2,rank, ties.method= "max",na.last="keep"))
rownames(expo)=new_expo[,1]
expo[is.na(expo)]=0
expo=expo[rownames(data),] #reorder as data space
index=0
for (i in 1:dim(expo)[2]){
  semantic=expo[,i] %*% space$tk
  cos=apply(tissuevec,2,cosine,as.vector(semantic))
  index[i]=names(which.max(cos))
  write.table(semantic,"GSE36376vec.csv",row.names=colnames(expo)[i],col.names = F,append=T)
  }
write.csv(index,"validation_index_GSE36376.csv",row.names=colnames(expo))

###GSE130970
expo=read.csv("GSE130970.csv",header=T,row.names=1)
aa=expo[rownames(expo) %in% intersect(rownames(data),rownames(expo)),]
bb=data.frame(cbind(id=rownames(aa),aa))
cc=data.frame(cbind(id=geneid))
new_expo=dplyr::full_join(bb,cc,by="id") #merge is slow!gtex first
expo=new_expo[,-1]
expo=as.matrix(apply(expo,2,rank, ties.method= "max",na.last="keep"))
rownames(expo)=new_expo[,1]
expo[is.na(expo)]=0
expo=expo[rownames(data),] #reorder as data space
index=0
for (i in 1:dim(expo)[2]){
  semantic=expo[,i] %*% space$tk
  cos=apply(tissuevec,2,cosine,as.vector(semantic))
  index[i]=names(which.max(cos))
  write.table(semantic,"GSE130970vec.csv",row.names=colnames(expo)[i],col.names = F,append=T)}
write.csv(index,"validation_index_GSE130970.csv",row.names=colnames(expo))
###correlate samplevec with traits
gse36376trait=read.csv("gse36376-trait.csv",header=T)
gse130970trait=read.csv("gse130970-trait.csv",header=T)
GSE130970vec=read.csv("GSE130970vec.csv",header=F,row.names = 1)
GSE36376vec=read.csv("GSE36376vec.csv",header=F,row.names = 1)
cor=0;shap=0;wilco=0
for (i in 1:100){correl=cor.test(GSE36376vec[,i],gse36376trait[,3],method="spearman")
cor[i]=paste(correl$estimate,correl$p.value,sep=",")
shap[i]=shapiro.test(GSE36376vec[,i])$p.value
wilco[i]=wilcox.test(GSE36376vec[1:193,i],GSE36376vec[194:433,i])$p.value
}
wilco;cor
cor=0;shap=0;wilco=0
for (i in 1:100){correl=cor.test(GSE130970vec[,i],gse130970trait[,9],method="spearman");
cor[i]=paste(correl$estimate,correl$p.value,sep=",")
shap[i]=shapiro.test(GSE130970vec[,i])$p.value
wilco[i]=wilcox.test(GSE130970vec[which(gse130970trait$Sex==1),i],GSE130970vec[which(gse130970trait$Sex==0),i])$p.value
}
wilco;cor

library(survival)
library("survminer")
data("lung")
head(lung)
fit <- survfit(Surv(time, status) ~ sex, data = lung)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
