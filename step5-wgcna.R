## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2019-03-09 19:56:32
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2019-03-09  First version
###
### ---------------

# 这个代码运行结果有点诡异，暂时放在这里。

rm(list = ls())  ## 魔幻操作，一键清空~

library(WGCNA)
## step1:input
if(T){
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL)#toTable提取包中的probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵
  head(ids)
  dat=dat[ids$probe_id,] #过滤，筛出不存在包中的探针
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median)
  #对dat这个矩阵按行操作，取每一行的中位数，将结果添加到ids矩阵median列
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  #对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]
  #去除ids矩阵symbol重复项
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  #把ids矩阵的symbol名作为dat矩阵的行名
  dat[1:4,1:4]  
  dim(dat)
  dat=dat[,group_list=='TNBC']
  
  load(file='trait.Rdata') 
  trait=trait[group_list=='TNBC',]
  head(trait)
  
  WGCNA_matrix = t(dat[order(apply(dat,1,mad), decreasing = T)[1:5000],])
  #mad绝对中位差
  #对dat矩阵每行求绝对中位差，order排序，降序排序，取前5000
  # 这个地方有问题。
  datExpr0 <- WGCNA_matrix  ## top 10000 mad genes
  datExpr <- datExpr0 
  
  datTraits=trait
  datTraits$gsm=rownames(datTraits)#添加gsm的列
  
  ## 下面主要是为了防止临床表型与样本名字对不上
  sampleNames = rownames(datExpr);
  traitRows = match(sampleNames, datTraits$gsm) 
  #match函数，返回sampleNames向量中每个元素在datTraits$gsm向量中的位置
  rownames(datTraits) = datTraits[traitRows, 'gsm'] 
  # pheatmap::pheatmap(datExpr)
}
## 查看是否有离群样品
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
## step2:beta value
if(T){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
sft$powerEstimate
## step3: Net
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate, 
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.1,
    # 调整 mergeCutHeight 参数可以控制最终模块数量。
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = FALSE, 
    verbose = 3
  )
  table(net$colors) 
  # 可以看到过多的基因直接的相似性太高，以至于它们根本就没什么特别的模块，哪怕是调整mergeCutHeight参数。
}
## step4: check sample clinical information
if(T){
  #明确样本数和基因数
  nGenes = ncol(datExpr) #基因名
  nSamples = nrow(datExpr) #样本
  #首先针对样本做个系统聚类树
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  #dist()函数计算变量间距离
  #hclust()函数进行聚类
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
  #针对前面构造的样品矩阵添加对应颜色
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$T)), 
                                  colors = rainbow(5),signed = FALSE)
  ## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
  #  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
  ## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码
  # 当然，这样给的颜色不明显，意义不大。
  # 构造所有样品的系统聚类树及性状热图
  par(mar = c(1,4,3,1),cex=0.8)
  plotDendroAndColors(datExpr_tree, sample_colors,
                      groupLabels = colnames(sample),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  
}
## step5:
if(T){
  datTraits$T=as.factor(datTraits$T)
  design=model.matrix(~0+ datTraits$T)
  colnames(design)=levels(datTraits$T)
  moduleColors <- labels2colors(net$colors)
  #将数值标签的矢量或阵列转换为与标签对应的矢量或颜色阵列
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = levels(datTraits$T),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}

save(net,datExpr,file='TNBC_net_top10000genes.Rdata')

load(file='TNBC_net_top10000genes.Rdata')
table(net$colors) 
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleGenes <- data.frame(module=mergedColors,
                          name=colnames(datExpr)
)
head(moduleGenes) 
library(clusterProfiler)
library(org.Hs.eg.db)
g_diff=unique(moduleGenes[,2])
#unique()去除重复，返回一个把重复元素或行删除的向量
gene.df <- bitr(g_diff, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db) 
head(gene.df)
moduleGenes=merge(moduleGenes,gene.df,by.x='name',by.y='SYMBOL')
head(moduleGenes)
mydf <- moduleGenes[,c(2,4)]  
formula_res <- compareCluster(ENTREZID~module, 
                              data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))
dotplot(formula_res)
library(ggplot2)
ggsave('modules_enrich_top10000.png')











