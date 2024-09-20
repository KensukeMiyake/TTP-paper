library(TCC)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(nichenetr)

set.seed(1234)

exp_mat = read.table(file="GeneExpressionMatrix.tsv", sep="\t", header=T) #from GSE222701
head(exp_mat)
colnames(exp_mat)[1] ="Gene_ID"

GTF = read.table("Mus_musculus.GRCm39.112.chr.gtf.genes.tab_new", header=T, sep="\t") #from Ensembl
GTF = GTF[,c(1:2)]
head(GTF)
View(GTF)
GTF[GTF$Gene_Symbol == "",]$Gene_Symbol <- NA

exp_mat = merge(GTF, exp_mat, by="Gene_ID")
exp_mat = na.omit(exp_mat)
exp_mat = exp_mat[!duplicated(exp_mat[,'Gene_Symbol']),]
rownames(exp_mat)=exp_mat$Gene_Symbol
exp_mat = exp_mat %>% arrange(Gene_Symbol)
View(exp_mat)

colnames(exp_mat) = 
  c("Gene_ID","Gene_Symbol","OVA_Prebaso_1","OVA_Prebaso_2","OVA_Prebaso_3","OVA_Mature_1","OVA_Mature_2","OVA_Mature_3",
         "TNP_Prebaso_1","TNP_Prebaso_2","TNP_Prebaso_3","TNP_Mature_1","TNP_Mature_2","TNP_Mature_3",
         "IL-3_Prebaso_1","IL-3_Prebaso_2","IL-3_Prebaso_3","IL-3_Mature_1","IL-3_Mature_2","IL-3_Mature_3")
View(exp_mat)
exp_mat = exp_mat[,-2]
View(exp_mat)

#================normalization and DEG analyses==================
#construct TCC class object
data = exp_mat[, c(5,6,7,11,12,13)]
group = c(rep(1,"OVA_Mature", 3),rep(2,"TNP_Mature", 3))
tcc = new("TCC", data, group)
head(tcc)

#filter low-count genes
dim(tcc$count)
tcc <- filterLowCountGenes(tcc)
dim(tcc$count)

#normalization of count data
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
tcc$norm.factors
normalized.count <- getNormalizedData(tcc)
head(normalized.count)
dim(normalized.count)

#DE analysis
tcc <- estimateDE(tcc, test.method="edger", FDR=0.1)
result <- getResult(tcc, sort=FALSE)
head(result)
sum(tcc$stat$q.value < 0.1) 
rownames(result) = result$gene_id

#================gene list of CCCH zinc finger proteins==================
#CCCH zinc finger protein list
zinc_list =c("CNOT4","CPSF4","CPSF4L","DHX57","DUS3L","HELZ","HELZ2","LENG9","MBNL1","MBNL2","MBNL3","MKRN1",
             "MKRN2","MKRN3","NUPL2","PAN3","PPP1R10","PRR3","RC3H1","RC3H2","RNF113A","RNF113B","TIPARP","TOE1",
             "TRMT1","U2AF1","U2AF1L4","ZC3HAV1L","PARP12","ZC3HAV1","ZC3H3","ZC3H4","UNK","UNKL","ZC3H6","ZC3H7A",
             "ZC3H7B","ZC3H8","ZGPAT","ZC3H10","ZC3H11A","ZC3H12A","ZC3H12B","ZC3H12C","ZC3H12D","ZC3H13","ZC3H14",
             "ZC3H15","RBM22","RBM26","ZC3H18","ZMAT5","RBM27","ZRSR2","ZFP36","ZFP36L1","ZFP36L2")
length(zinc_list)
#REF: https://www.nature.com/articles/nri.2016.129

#ZC3H1 = PARP12
#ZC3H2 = ZC3HAV1
#ZC3H5 = UNK
#ZC3H5L = UNKL
#ZC3H9 = ZGPAT
#ZC3H16 = RBM22
#ZC3H17 = RBM26
#ZC3H19 = ZMAT5
#ZC3H20 = RBM27
#ZC3H22 = ZRSR2
#RNF113B : not found in NCBI Gene as mouse gene (2024.9.3)


#convert human genes to mouse genes
zinc_list = nichenetr::convert_human_to_mouse_symbols(zinc_list, version = 1)
length(zinc_list) 
zinc_list #RNF113B & U2AF1: could not be converted to mouse genes
zinc_list = zinc_list[is.na(zinc_list)==F]
length(zinc_list)
zinc_list = unname(zinc_list)

#Zc3h21 and Zfp36l3: rodent-specific genes (REF: https://www.nature.com/articles/nri.2016.129)
#Mus musculus Zc3h21 gene was not found in NCBI Gene and in our data, and ZRSR2 pseudogene 1(ZC3H21) was found as Homo sapiens gene (2024.9.13)
#Mus musculus Zrsr2-ps1 was found in NCBI Gene and in our data (2024.9.13)
#U2AF1 could not be converted to murine genes by nichenetr package (2024.9.3), but U2af1 gene was found in our murine data 
zinc_list = c(zinc_list, c("Zrsr2-ps1", "Zfp36l3", "U2af1", "Rnf113b")) 
DEGs_zinc = filter(result, gene_id %in% zinc_list)
dim(DEGs_zinc) #57 genes
View(DEGs_zinc) #not found in our data: RNF113B, Zfp36l3


#================volcano plot of CCCH zinc finger proteins==================
interested = c("Zfp36")

DEGs_zinc$diffexpressed <- "NO"
# if log2Foldchange(m.value) > 0.5 and q.value < 0.1, set as "UP" 
DEGs_zinc$diffexpressed[DEGs_zinc$m.value > 0.5 & DEGs_zinc$q.value < 0.1] <- "UP"
# if log2Foldchange(m.value) < -0.5 and q.value < 0.1, set as "DOWN"
DEGs_zinc$diffexpressed[DEGs_zinc$m.value < -0.5 & DEGs_zinc$q.value < 0.1] <- "DOWN"

DEGs_zinc$diffexpressed[DEGs_zinc$gene_id %in% interested] <- "INTEREST"
DEGs_zinc$diffexpressed = factor(DEGs_zinc$diffexpressed, levels = c("INTEREST", "NO", "UP", "DOWN"))


DEGs_zinc$delabel <- NA
DEGs_zinc$delabel[DEGs_zinc$diffexpressed != "NO"] <- rownames(DEGs_zinc)[DEGs_zinc$diffexpressed != "NO"]


p = ggplot(data=DEGs_zinc, aes(x=m.value, y=-log10(q.value), col=diffexpressed, size=diffexpressed, label = delabel)) +
  geom_point() + 
  geom_text_repel()+
  scale_color_manual(values=c("red", "gray", "orange2", "turquoise4"))+
  scale_size_manual(values=c(3,1,2,2))+
  geom_vline(xintercept = -0.5, colour="gray40", linetype = 5)+
  geom_vline(xintercept = 0.5, colour="gray40", linetype = 5)+
  geom_hline(yintercept = 1, colour="gray40", linetype = 5)+
  xlim(-1.8, 1.8)+
  theme_classic()+
  theme(legend.position = "none")
p
