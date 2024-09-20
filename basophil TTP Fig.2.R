library(TCC)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggnewscale)
library(Seurat)

set.seed(1234)

exp_mat = read.table(file="GeneExpressionMatrix.tsv", sep="\t", header=T )
head(exp_mat)
colnames(exp_mat)[1] ="Gene_ID"

GTF = read.table("Mus_musculus.GRCm39.112.chr.gtf.genes.tab_new", header=T, sep="\t") #from Ensembl
GTF = GTF[,c(1:2)]
head(GTF)
View(GTF)
GTF[GTF$Gene_Symbol == "",]$Gene_Symbol <- NA

exp_mat = merge(GTF, exp_mat, by="Gene_ID")
exp_mat <- na.omit(exp_mat)
exp_mat = exp_mat[!duplicated(exp_mat[,'Gene_Symbol']),]
rownames(exp_mat)=exp_mat$Gene_Symbol
exp_mat = exp_mat %>% arrange(Gene_Symbol)
View(exp_mat)

colnames(exp_mat) = 
  c("Gene_ID","Gene_Symbol","WT_unstim_1","WT_unstim_2","WT_unstim3","KO_unstim1","KO_unstim2","KO_unstim3",
         "WT_TNP_2h_1","WT_TNP_2h_2","WT_TNP_2h_3","KO_TNP_2h_1","KO_TNP_2h_2","KO_TNP_2h_3",
         "WT_TNP_4h_1","WT_TNP_4h_2","WT_TNP_4h_3","KO_TNP_4h_1","KO_TNP_4h_2","KO_TNP_4h_3")
View(exp_mat)
exp_mat = exp_mat[,-2]
dim(exp_mat)
View(exp_mat)

#================TNP 4h (Fig.2D)==================
#construct TCC class object (4h)
data_4h = exp_mat[, c(14:19)]
group_4h = c(rep(5,"WT_TNP4h", 3),rep(6,"KO_TNP4h", 3))
tcc_4h = new("TCC", data_4h, group_4h)
head(tcc_4h)

#filter low-count genes
dim(tcc_4h$count)
tcc_4h <- filterLowCountGenes(tcc_4h)
dim(tcc_4h$count)

#normalization
tcc_4h <- calcNormFactors(tcc_4h, norm.method="tmm", test.method="edger",
                          iteration=3, FDR=0.1, floorPDEG=0.05)
tcc_4h$norm.factors
normalized.count_4h <- getNormalizedData(tcc_4h)
head(normalized.count_4h)
dim(normalized.count_4h)

#DE analysis
tcc_4h <- estimateDE(tcc_4h, test.method="edger", FDR=0.1)
result_4h <- getResult(tcc_4h, sort=FALSE)
head(result_4h)
rownames(result_4h) = result_4h$gene_id
tmp_4h <- cbind(normalized.count_4h, result_4h)
dim(tmp_4h)
View(tmp_4h)

#DE number 
sum(tcc_4h$stat$q.value < 0.1) #total DEGs
dim(tmp_4h[tmp_4h$m.value>0 & tmp_4h$estimatedDEG==1,]) #DEGs upregulated by TTP-deficiency

#Volcano plot
interested = c("Il4", "Ccl3", "Areg")

tmp_4h$diffexpressed <- "NO"
# if log2Foldchange (m.value) > 0 and q.value < 0.1, set as "UP" 
tmp_4h$diffexpressed[tmp_4h$m.value > 0 & tmp_4h$estimatedDEG==1] <- "UP"
# if log2Foldchange (m.value) < 0 and q.value < 0.1, set as "DOWN"
tmp_4h$diffexpressed[tmp_4h$m.value < 0 & tmp_4h$estimatedDEG==1] <- "DOWN"

tmp_4h$diffexpressed[tmp_4h$gene_id %in% interested] <- "INTEREST"
tmp_4h$diffexpressed = factor(tmp_4h$diffexpressed, levels = c("INTEREST", "NO", "UP", "DOWN"))

p = ggplot(data=tmp_4h, aes(x=m.value, y=-log10(q.value), col=diffexpressed, size=diffexpressed)) +
  geom_point() + 
  scale_color_manual(values=c("red", "gray", "bisque", "bisque"))+
  scale_size_manual(values=c(2,0.5,0.5,0.5))+
  geom_vline(xintercept = 0, colour="gray40", linetype = 5)+
  geom_hline(yintercept = 1, colour="gray40", linetype = 5)+
  xlim(-5,5)+
  theme_classic()+
  theme(legend.position = "none")
p
p = LabelPoints(plot = p, points = interested, repel = TRUE)
p

#================2h (Fig.2C)==================
#construct TCC class object (2h)
data_2h = exp_mat[, c(8:13)]
group_2h = c(rep(3,"WT_TNP2h", 3),rep(4,"KO_TNP2h", 3))
tcc_2h = new("TCC", data_2h, group_2h)
head(tcc_2h)

#filter low-count genes
dim(tcc_2h$count)
tcc_2h <- filterLowCountGenes(tcc_2h)
dim(tcc_2h$count)

#normalization
tcc_2h <- calcNormFactors(tcc_2h, norm.method="tmm", test.method="edger",
                          iteration=3, FDR=0.1, floorPDEG=0.05)
tcc_2h$norm.factors
normalized.count_2h <- getNormalizedData(tcc_2h)
head(normalized.count_2h)
dim(normalized.count_2h)

#DE analysis
tcc_2h <- estimateDE(tcc_2h, test.method="edger", FDR=0.1)
result_2h <- getResult(tcc_2h, sort=FALSE)
head(result_2h)
rownames(result_2h) = result_2h$gene_id
tmp_2h <- cbind(normalized.count_2h, result_2h)
head(tmp_2h)
View(tmp_2h)

#DE number
sum(tcc_2h$stat$q.value < 0.1) #total DEGs
dim(tmp_2h[tmp_2h$m.value>0 & tmp_2h$estimatedDEG==1,]) #DEGs upregulated by TTP-deficiency

#Volcano plot
interested = c("Il13", "Ccl3", "Cxcl2", "Ptgs2", "Areg")

tmp_2h$diffexpressed <- "NO"
# if log2Foldchange (m.value) > 0 and q.value < 0.1, set as "UP" 
tmp_2h$diffexpressed[tmp_2h$m.value > 0 & tmp_2h$estimatedDEG==1] <- "UP"
# if log2Foldchange (m.value) < 0 and q.value < 0.1, set as "DOWN"
tmp_2h$diffexpressed[tmp_2h$m.value < 0 & tmp_2h$estimatedDEG==1] <- "DOWN"

tmp_2h$diffexpressed[tmp_2h$gene_id %in% interested] <- "INTEREST"
tmp_2h$diffexpressed = factor(tmp_2h$diffexpressed, levels = c("INTEREST", "NO", "UP", "DOWN"))

p = ggplot(data=tmp_2h, aes(x=m.value, y=-log10(q.value), col=diffexpressed, size=diffexpressed)) +
  geom_point() + 
  scale_color_manual(values=c("red", "gray", "wheat1", "wheat1"))+
  scale_size_manual(values=c(2,0.5,0.5,0.5))+
  geom_vline(xintercept = 0, colour="gray40", linetype = 5)+
  geom_hline(yintercept = 1, colour="gray40", linetype = 5)+
  xlim(-5,5)+
  theme_classic()+
  theme(legend.position = "none")
p
p = LabelPoints(plot = p, points = interested, repel = TRUE)
p

#================unstim (Supplementary Figure)==================
#construct TCC class object (unstim)
data_unstim = exp_mat[, c(2:7)]
group_unstim = c(rep(1,"WT_unstim", 3),rep(2,"KO_unstim", 3))
tcc_unstim = new("TCC", data_unstim, group_unstim)
head(tcc_unstim)

#filter low-count genes
dim(tcc_unstim$count)
tcc_unstim <- filterLowCountGenes(tcc_unstim)
dim(tcc_unstim$count)

#normalization
tcc_unstim <- calcNormFactors(tcc_unstim, norm.method="tmm", test.method="edger",
                              iteration=3, FDR=0.1, floorPDEG=0.05)
tcc_unstim$norm.factors
normalized.count_unstim <- getNormalizedData(tcc_unstim)
head(normalized.count_unstim)
dim(normalized.count_unstim)

#DE analysis
tcc_unstim <- estimateDE(tcc_unstim, test.method="edger", FDR=0.1)
result_unstim <- getResult(tcc_unstim, sort=FALSE)
head(result_unstim)
rownames(result_unstim) = result_unstim$gene_id
tmp_unstim <- cbind(normalized.count_unstim, result_unstim)
head(tmp_unstim)
View(tmp_unstim)

#DE number
sum(tcc_unstim$stat$q.value < 0.1) #total DEGs
dim(tmp_unstim[tmp_unstim$m.value>0 & tmp_unstim$estimatedDEG==1,]) #DEGs upregulated by TTP-deficiency

#Volcano plot
tmp_unstim$diffexpressed <- "NO"
# if log2Foldchange (m.value) > 0 and q.value < 0.1, set as "UP" 
tmp_unstim$diffexpressed[tmp_unstim$m.value > 0 & tmp_unstim$estimatedDEG==1] <- "UP"
# if log2Foldchange (m.value) < 0 and q.value < 0.1, set as "DOWN"
tmp_unstim$diffexpressed[tmp_unstim$m.value < 0 & tmp_unstim$estimatedDEG==1] <- "DOWN"
tmp_unstim$diffexpressed = factor(tmp_unstim$diffexpressed, levels = c("NO", "UP", "DOWN"))

p = ggplot(data=tmp_unstim, aes(x=m.value, y=-log10(q.value), col=diffexpressed, size=diffexpressed)) +
  geom_point() + 
  scale_color_manual(values=c("gray", "wheat1", "wheat1"))+
  scale_size_manual(values=c(0.5,0.5,0.5))+
  geom_vline(xintercept = 0, colour="gray40", linetype = 5)+
  geom_hline(yintercept = 1, colour="gray40", linetype = 5)+
  xlim(-5,5)+
  theme_classic()+
  theme(legend.position = "none")
p
p = LabelPoints(plot = p, points = interested, repel = TRUE)
p

#================enrichment analysis (Fig.2E)==================
#extract upregulated genes
upDEGs_unstim = rownames(tmp_unstim %>% filter(m.value> 0) %>% filter(estimatedDEG==1))
upDEGs_2h = rownames(tmp_2h %>% filter(m.value> 0) %>% filter(estimatedDEG==1))
upDEGs_4h = rownames(tmp_4h %>% filter(m.value> 0) %>% filter(estimatedDEG==1))


#convert to ENTREZID
#unstim
upDEGs_unstim <- bitr(upDEGs_unstim, 
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)
#2hr
upDEGs_2h <- bitr(upDEGs_2h, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
#4hr
upDEGs_4h <- bitr(upDEGs_4h, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)

#show dotplot
a = as.list(upDEGs_unstim$ENTREZID)
b = as.list(upDEGs_2h$ENTREZID)
c = as.list(upDEGs_4h$ENTREZID)
gene_list = list("unstim" = a, "2h" = b, "4h" = c)

ck <- compareCluster(geneCluster = gene_list, fun = enrichGO, OrgDb = org.Mm.eg.db, ont = "BP")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 

GO_list = c("regulation of inflammatory response", "production of molecular mediator of immune response", "cytokine production involved in immune response", "cytoplasmic translation", "ribosome biogenesis")
enrichplot::dotplot(ck, showCategory = GO_list,  font.size=8) + geom_count() + scale_size_area()