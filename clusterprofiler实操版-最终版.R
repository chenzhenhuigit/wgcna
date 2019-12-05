getwd()
workingDir = "G:\\RCZH\\gse83514" 
setwd(workingDir)
getwd()
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(org.Bt.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

######数据第一步处理######
keytypes(org.Bt.eg.db) 
data <- read.table("green module.txt",header=FALSE) #输入文件名称
data$V1 <- as.character(data$V1) 
test1 = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Bt.eg.db")
head(test1,2)

#go分析
#group go分析
#ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Bt.eg.db, ont = "CC",level = 3,readable = TRUE)
#enrichGO 富集分析
ego_ALL <- enrichGO(gene = test1$ENTREZID, 
                    OrgDb = org.Bt.eg.db, 
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "none", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种,优先BH
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,#0.05也可
                    readable = TRUE) 
head(ego_ALL,2)

#表格形式输出输出结果
write.table(ego_ALL, file="ALL-enrich.xls",row.names=F, col.names=T,quote=FALSE,sep="\t")
#绘制图形
##可视化--点图#点图，按富集的数从大到小的
dotplot(ego_ALL,title="Enrichment GO_dot")
##可视化--条形图
barplot(ego_ALL, showCategory=20,title="Enrichment GO_bar")#条状图，按p从小到大排，绘制前20个Term
##可视化--
#goplot(ego_ALL,title="Enrichment GO_dot")
##绘制有向无环图###
ego_MF <- enrichGO(gene = test1$ENTREZID,OrgDb = org.Bt.eg.db,ont = "MF", pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = FALSE)
plotGOgraph(ego_MF)
ego_BP <- enrichGO(gene = test1$ENTREZID,OrgDb = org.Bt.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = FALSE)
plotGOgraph(ego_BP)
ego_CC <- enrichGO(gene = test1$ENTREZID,OrgDb = org.Bt.eg.db,ont = "CC", pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = FALSE)
plotGOgraph(ego_CC)
##GO terms之间的基因重叠关系


###基因和go terms之间的对应关系
#cnetplot(ego_ALL,showCategory = 5)

#4.1 候选基因进行通路分析
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'bta', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05)
head(kk,2)

#结果输出
write.table(kk, file="KEGG-enrich.xls",row.names=F, col.names=T,quote=FALSE,sep="\t")

#气泡图
dotplot(kk,title="Enrichment KEGG_dot")
#emapplot(kk, showCategory = 30, color = "p.adjust", layout = "kk")
#查看特定通路
hsa04750 <- pathview(gene.data = geneList,
                     pathway.id = "hsa04750", #上述结果中的hsa04750通路
                     species = "hsa",
                     limit = list(gene=max(abs(geneList)), cpd=1))
